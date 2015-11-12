#!/usr/bin/env python2.7
#Jeltje van Baren 
"""
Batch CNV script using Toil

Runs Varscan on input bam files.

If the option is specified (s3_dir), the output table will be placed
in S3.  ~/.boto config file and S3AM: https://github.com/BD2KGenomics/s3am
are required for this step.

Dependencies:
Docker  -   apt-get install docker.io
Toil    -   pip install toil
S3AM*   -   pip install --pre S3AM  (optional)
Curl    -   apt-get install curl
"""
import argparse
import base64
from collections import OrderedDict
import hashlib
import os
import subprocess
import multiprocessing
import shutil
import sys
from toil.job import Job


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', default=None, help='configuration file with ID and URLs to bam inputs (control, tumor): uuid,url,url,...')
    parser.add_argument('-s', '--ssec', default=None, help='Path to Key File for SSE-C Encryption')
    parser.add_argument('-r', '--ref', required=True, help='Reference fasta file')
    parser.add_argument('-f', '--fai', required=True, help='Reference fasta file (fai)')
    parser.add_argument('-b', '--cent', required=True, help='centromere locations (bed format)')
    parser.add_argument('-w', '--white', required=True, help='exome whitelist (bed format)')
    parser.add_argument('-o', '--out', default=None, help='full path where final results will be output')
    parser.add_argument('-3', '--s3_dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                             'cgl-driver-projects/ckcc/rna-seq-samples/')
    return parser


# Convenience Functions
def generate_unique_key(master_key_path, url):
    """
    Input1: Path to the BD2K Master Key (for S3 Encryption)
    Input2: S3 URL (e.g. https://s3-us-west-2.amazonaws.com/cgl-driver-projects-encrypted/wcdt/exome_bams/DTB-111-N.bam)

    Returns: 32-byte unique key generated for that URL
    """
    with open(master_key_path, 'r') as f:
        master_key = f.read()
    assert len(master_key) == 32, 'Invalid Key! Must be 32 characters. ' \
                                  'Key: {}, Length: {}'.format(master_key, len(master_key))
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is invalid and is not 32 characters: {}'.format(new_key)
    return new_key


def download_encrypted_file(work_dir, url, key_path, name):
    """
    Downloads encrypted file from S3

    Input1: Working directory
    Input2: S3 URL to be downloaded
    Input3: Path to key necessary for decryption
    Input4: name of file to be downloaded
    """
    file_path = os.path.join(work_dir, name)
    key = generate_unique_key(key_path, url)
    encoded_key = base64.b64encode(key)
    encoded_key_md5 = base64.b64encode(hashlib.md5(key).digest())
    h1 = 'x-amz-server-side-encryption-customer-algorithm:AES256'
    h2 = 'x-amz-server-side-encryption-customer-key:{}'.format(encoded_key)
    h3 = 'x-amz-server-side-encryption-customer-key-md5:{}'.format(encoded_key_md5)
    try:
        subprocess.check_call(['curl', '-fs', '--retry', '5', '-H', h1, '-H', h2, '-H', h3, url, '-o', file_path])
    except OSError:
        raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)

def download_S3_file(work_dir, url, name):
    """
    Downloads file from S3

    Input1: Working directory
    Input2: S3 URL to be downloaded
    Input3: name of file to be downloaded
    """
    file_path = os.path.join(work_dir, name)
    # curl fs means fail quietly but show error when failing (not clear how that works)
    try:
        subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
    except OSError:
        raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)


def download_from_url(job, input_args, ids, name):
    """
    Downloads a file from a URL and places it in the jobStore

    Input1: Toil job instance
    Input2: Input arguments
    Input3: jobstore id dictionary
    Input4: Name of key used to access url in input_args
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, name)
    url = input_args[name]
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    job.fileStore.updateGlobalFile(ids[name], file_path)


def return_input_paths(job, work_dir, ids, *args):
    """
    Returns the paths of files from the FileStore

    Input1: Toil job instance
    Input2: Working directory
    Input3: jobstore id dictionary
    Input4: names of files to be returned from the jobstore

    Returns: path(s) to the file(s) requested -- unpack these!
    """
    paths = OrderedDict()
    for name in args:
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = job.fileStore.readGlobalFile(ids[name], os.path.join(work_dir, name))
        else:
            file_path = os.path.join(work_dir, name)
        paths[name] = file_path
        if len(args) == 1:
            return file_path

    return paths.values()


def move_to_output_dir(work_dir, output_dir, uuid=None, files=list()):
    """
    Moves files from work_dir to output_dir

    Input1: Working directory
    Input2: Output directory
    Input3: UUID to be preprended onto file name
    Input4: list of file names to be moved from working dir to output dir
    """
    for fname in files:
        if uuid is None:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


# Start of Job Functions
def batch_start(job, input_args):
    """
    Downloads and places shared files that are used by all samples for alignment
    """
    input_args['cpu_count'] = multiprocessing.cpu_count()
    shared_files = ['ref.fa', 'ref.fa.fai', 'cent.bed', 'white.bed']
    shared_ids = {x: job.fileStore.getEmptyFileStoreID() for x in shared_files}
    for fname in shared_files:
        job.addChildJobFn(download_from_url, input_args, shared_ids, fname)
    job.addFollowOnJobFn(spawn_batch_jobs, shared_ids, input_args)


def spawn_batch_jobs(job, shared_ids, input_args):
    """
    Spawns a varscan job for every sample in the input configuration file
    """
    samples = []
    config = input_args['config']
    cores = input_args['cpu_count']
    with open(config, 'r') as f_in:
        for line in f_in:
            uuid, c_url, t_url = line.strip().split(',')
            samples.append((uuid, c_url, t_url))
    for sample in samples:
        job.addChildJobFn(varscan, shared_ids, input_args, sample, cores=cores)


def varscan(job, ids, input_args, sample):
    """
    Runs varscan on the input bams for this sample

    Input1: Toil Job instance
    Input2: jobstore id dictionary
    Input3: Input arguments dictionary
    Input4: Sample UUID and urls
    """
    uuid, c_url, t_url = sample
    # TODO How do I do this?
#    output_files = ['copyCalled.recenter', 'output.copynumber']
#    output_ids = {x: job.fileStore.getEmptyFileStoreID() for x in output_files}
    ids['cnv'] = job.fileStore.getEmptyFileStoreID()
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    key_path = input_args['ssec']
    cores = input_args['cpu_count']

    # I/O
    return_input_paths(job, work_dir, ids, 'ref.fa', 'ref.fa.fai', 'cent.bed', 'white.bed') 

    # Get bams associated with this sample
    download_encrypted_file(work_dir, c_url, key_path, uuid + ".control.bam")
    download_encrypted_file(work_dir, t_url, key_path, uuid + ".tumor.bam")
#    for url in urls:
#        download_S3_file(work_dir, url, os.path.basename(url))
    #sam_path=input_args['insam']
    #shutil.copy(sam_path, os.path.join(work_dir, 'input.sam'))

    # Setup docker base and varscan command
    docker_cmd = ['docker', 'run', '--rm', '-v', '{}:/data'.format(work_dir)]


#docker run -v /data/tmp/vtest:/data jeltje/varscan2:0.0.1 -t /data/tumor.10x.bam -c /data/control.10x.bam -q <sample ID> -s /data -i /data/Homo_sapiens_assembly19.fasta -b /data/GRCh37.centromeres -w /data/broadTargets.bed 2> data/err.log > varscan.out
    varscan_command = ["jeltje/varscan", 
		"-t", "/data/{}".format(uuid + '.tumor.bam'),
		"-c", "/data/{}".format(uuid + '.control.bam'),
		"-q", "{}".format(uuid),
		"-s", "/data",
		"-i", "/data/ref.fa",
		"-b", "/data/cent.bed",
		"-w", "/data/white.bed"]

    # Piping the output to a file handle
    # check_call blocks progress until finished
    outfile = uuid + '.cnv'
    with open(os.path.join(work_dir, outfile), 'w') as f_out:
        subprocess.check_call(docker_cmd + varscan_command, stdout=f_out)

    # Save in JobStore
    job.fileStore.updateGlobalFile(ids['cnv'], os.path.join(work_dir, outfile))
    # Move file in output_dir
    if input_args['output_dir']:
        move_to_output_dir(work_dir, output_dir, uuid=None, files=[outfile])
    # Copy file to S3
    if input_args['s3_dir']:
        job.addChildJobFn(upload_file_to_s3, ids, input_args, sample[0], cores=cores)

def upload_file_to_s3(job, ids, input_args, uuid):
    """
    Uploads output file from sample to S3

    Input1: Toil Job instance
    Input2: jobstore id dictionary
    Input3: Input arguments dictionary
    Input4: Sample uuid
    """
    work_dir = job.fileStore.getLocalTempDir()
    key_path = input_args['ssec']
    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.split('/')[0]
    bucket_dir = '/'.join(s3_dir.split('/')[1:])
    base_url = 'https://s3-us-west-2.amazonaws.com/'
    outfile = uuid + '.cnv'
    url = os.path.join(base_url, bucket_name, bucket_dir, outfile)
    #I/O
    job.fileStore.readGlobalFile(ids['cnv'], os.path.join(work_dir, outfile))
    # Generate keyfile for upload
    with open(os.path.join(work_dir, uuid + '.key'), 'wb') as f_out:
        f_out.write(generate_unique_key(key_path, url))
    # Commands to upload to S3 via S3AM
    s3am_command = ['s3am',
                    'upload',
                    '--sse-key-file', os.path.join(work_dir, uuid + '.key'),
                    'file://{}'.format(os.path.join(work_dir, outfile)),
                    bucket_name,
                    os.path.join(bucket_dir, outfile)]

    subprocess.check_call(s3am_command)


if __name__ == "__main__":
    # Define Parser object and add to toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    # Store input_URLs for downloading
    inputs = {'config': args.config,
              'ref.fa': args.ref,
              'ref.fa.fai': args.fai,
              'cent.bed': args.cent,
              'white.bed': args.white,
              'ssec':args.ssec,
              'output_dir': args.out,
              's3_dir': args.s3_dir,
              'cpu_count': None}

    # Launch jobs
    Job.Runner.startToil(Job.wrapJobFn(batch_start, inputs), args)
