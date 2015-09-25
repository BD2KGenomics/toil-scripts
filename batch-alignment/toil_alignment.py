#!/usr/bin/env python2.7
# John Vivian
"""
Batch alignment script using Toil

"""
import argparse
import base64
from collections import OrderedDict
import hashlib
import os
import subprocess
import multiprocessing
import shutil
from toil.job import Job


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', required=True, help='configuration file: uuid,url,url,...')
    parser.add_argument('-r', '--ref', required=True, help='Reference fasta file')
    parser.add_argument('-m', '--amb', required=True, help='Reference fasta file (amb)')
    parser.add_argument('-n', '--ann', required=True, help='Reference fasta file (ann)')
    parser.add_argument('-b', '--bwt', required=True, help='Reference fasta file (bwt)')
    parser.add_argument('-p', '--pac', required=True, help='Reference fasta file (pac)')
    parser.add_argument('-a', '--sa', required=True, help='Reference fasta file (sa)')
    parser.add_argument('-f', '--fai', required=True, help='Reference fasta file (fai)')
    parser.add_argument('-s', '--ssec', help='Path to Key File for SSE-C Encryption')
    parser.add_argument('-o', '--out', required=True, help='full path where final results will be output')
    parser.add_argument('-3', '--s3_dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                             'cgl-driver-projects/ckcc/rna-seq-samples/')
    return parser


# Convenience Functions
def generate_unique_key(master_key_path, url):
    with open(master_key_path, 'r') as f:
        master_key = f.read()
    assert len(master_key) == 32, 'Invalid Key! Must be 32 characters. ' \
                                  'Key: {}, Length: {}'.format(master_key, len(master_key))
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is invalid and is not 32 characters: {}'.format(new_key)
    return new_key


def download_encrypted_file(work_dir, url, key_path, name):
    """
    Downloads encrypted files
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


def download_from_url(job, input_args, ids, name):
    """
    Downloads a URL that was supplied as an argument to running this script in LocalTempDir.
    After downloading the file, it is stored in the FileStore.
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
    Returns the paths of files from the FileStore if they are not present.
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
    A list of files to move from work_dir to output_dir.
    """
    for fname in files:
        if uuid is None:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


# Start of Job Functions
def batch_start(job, input_args):
    shared_files = ['ref.fa', 'ref.fa.amb', 'ref.fa.ann', 'ref.fa.bwt', 'ref.fa.pac', 'ref.fa.sa', 'ref.fa.fai']
    shared_ids = {x: job.fileStore.getEmptyFileStoreID() for x in shared_files}
    # Download shared files used by all samples in the pipeline
    for fname in shared_files:
        job.addChildJobFn(download_from_url, input_args, shared_ids, fname)
    job.addFollowOnJobFn(spawn_batch_jobs, shared_ids, input_args)


def spawn_batch_jobs(job, shared_ids, input_args):
    samples = []
    config = input_args['config']
    with open(config, 'r') as f_in:
        for line in f_in:
            line = line.strip().split(',')
            uuid = line[0]
            urls = line[1:]
            samples.append((uuid, urls))
    for sample in samples:
        job.addChildJobFn(bwa, shared_ids, input_args, sample, cores=32, memory='20 G', disk='100 G')


def bwa(job, ids, input_args, sample):
    uuid, urls = sample
    ids['bam'] = job.fileStore.getEmptyFileStoreID()
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    key_path = input_args['ssec']
    cores = str(multiprocessing.cpu_count())

    # I/O
    return_input_paths(job, work_dir, ids, 'ref.fa', 'ref.fa.amb', 'ref.fa.ann',
                                                     'ref.fa.bwt', 'ref.fa.pac', 'ref.fa.sa', 'ref.fa.fai')
    # Get fastqs associated with this sample
    for url in urls:
        download_encrypted_file(work_dir, url, key_path, os.path.basename(url))

    # Parameters for BWA and Bamsort
    docker_cmd = ['docker', 'run', '--rm', '-v', '{}:/data'.format(work_dir)]

    bwa_command = ["jvivian/bwa",
                   "mem",
                   "-R", "@RG\tID:{0}\tPL:Illumina\tSM:{0}\tLB:KapaHyper".format(uuid),
                   "-T", str(0),
                   "-t", cores,
                   "/data/ref.fa"] + [os.path.join('/data/',  os.path.basename(x)) for x in urls]

    bamsort_command = ["jeltje/biobambam",
                       "/usr/local/bin/bamsort",
                       "inputformat=sam",
                       "level=1",
                       "inputthreads={}".format(cores),
                       "outputthreads={}".format(cores),
                       "calmdnm=1",
                       "calmdnmrecompindetonly=1",
                       "calmdnmreference=/data/ref.fa",
                       "I=/data/{}".format(uuid + '.sam')]
    # Piping the output to a file handle
    with open(os.path.join(work_dir, uuid + '.sam'), 'w') as f_out:
        subprocess.check_call(docker_cmd + bwa_command, stdout=f_out)

    with open(os.path.join(work_dir, uuid + '.bam'), 'w') as f_out:
        subprocess.check_call(docker_cmd + bamsort_command, stdout=f_out)

    # Save in JobStore
    job.fileStore.updateGlobalFile(ids['bam'], os.path.join(work_dir, uuid + '.bam'))
    # Place file in output_dir
    # move_to_output_dir(work_dir, output_dir, uuid=None, files=[uuid + '.bam'])
    # Move file to S3
    if input_args['s3_dir']:
        job.addChildJobFn(upload_bam_to_s3, ids, input_args, sample, cores=4, memory='20 G', disk='30 G')



def upload_bam_to_s3(job, ids, input_args, sample):
    uuid, urls = sample
    key_path = input_args['ssec']
    work_dir = job.fileStore.getLocalTempDir()
    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.split('/')[0]
    bucket_dir = '/'.join(s3_dir.split('/')[1:])
    base_url = 'https://s3-us-west-2.amazonaws.com/'
    url = os.path.join(base_url, bucket_name, bucket_dir, uuid + '.bam')
    #I/O
    job.fileStore.readGlobalFile(ids['bam'], os.path.join(work_dir, uuid + '.bam'))
    # Generate keyfile for upload
    with open(os.path.join(work_dir, uuid + '.key'), 'wb') as f_out:
        f_out.write(generate_unique_key(key_path, url))
    # Commands to upload to S3 via S3AM
    s3am_command = ['s3am',
                    'upload',
                    '--sse-key-file', os.path.join(work_dir, uuid + '.key'),
                    'file://{}'.format(os.path.join(work_dir, uuid + '.bam')),
                    bucket_name,
                    os.path.join(bucket_dir, uuid + '.bam')]

    subprocess.check_call(s3am_command)


if __name__ == "__main__":
    # Define Parser object and add to toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    # Store input_URLs for downloading
    inputs = {'config': args.config,
              'ref.fa': args.ref,
              'ref.fa.amb': args.amb,
              'ref.fa.ann': args.ann,
              'ref.fa.bwt': args.bwt,
              'ref.fa.pac': args.pac,
              'ref.fa.sa': args.sa,
              'ref.fa.fai': args.fai,
              'ssec':args.ssec,
              'output_dir': args.out,
              's3_dir': args.s3_dir,
              'cpu_count': None}

    # Launch jobs
    Job.Runner.startToil(Job.wrapJobFn(batch_start, inputs), args)