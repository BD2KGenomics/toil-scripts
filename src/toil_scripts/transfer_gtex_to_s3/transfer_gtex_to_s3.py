#!/usr/bin/env python2.7
"""
Toil script to move TCGA data into an S3 bucket.

Dependencies
Curl:       apt-get install curl
Docker:     wget -qO- https://get.docker.com/ | sh
Toil:       pip install toil
S3AM:       pip install --pre s3am
"""
import argparse
import glob
import hashlib
import os
import shutil
import subprocess
import tarfile
from toil.job import Job


def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--sra', default=None, required=True,
                        help='Path to a file with one analysis ID per line for data hosted on CGHub.')
    parser.add_argument('-k', '--dbgap_key', default=None, required=True,
                        help='Path to a CGHub key that has access to the TCGA data being requested. An exception will'
                             'be thrown if "-g" is set but not this argument.')
    parser.add_argument('--s3_dir', default=None, required=True, help='S3 Bucket. e.g. tcga-data')
    parser.add_argument('--ssec', default=None, required=True, help='Path to Key File for SSE-C Encryption')
    parser.add_argument('--single_end', default=None, action='store_true', help='Set this flag if data is single-end')
    parser.add_argument('--sudo', dest='sudo', default=None, action='store_true',
                        help='Docker usually needs sudo to execute locally, but not when running Mesos or when '
                             'the user is a member of a Docker group.')
    return parser


# Convenience Functions
def generate_unique_key(master_key_path, url):
    """
    master_key_path: str    Path to the BD2K Master Key (for S3 Encryption)
    url: str                S3 URL (e.g. https://s3-us-west-2.amazonaws.com/bucket/file.txt)

    Returns: str            32-byte unique key generated for that URL
    """
    with open(master_key_path, 'r') as f:
        master_key = f.read()
    assert len(master_key) == 32, 'Invalid Key! Must be 32 characters. ' \
                                  'Key: {}, Length: {}'.format(master_key, len(master_key))
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is invalid and is not 32 characters: {}'.format(new_key)
    return new_key


def docker_call(work_dir, tool_parameters, tool, java_opts=None, sudo=False, outfile=None):
    """
    Makes subprocess call of a command to a docker container.


    tool_parameters: list   An array of the parameters to be passed to the tool
    tool: str               Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools)
    java_opts: str          Optional commands to pass to a java jar execution. (e.g. '-Xmx15G')
    outfile: file           Filehandle that stderr will be passed to
    sudo: bool              If the user wants the docker command executed as sudo
    """
    base_docker_call = 'docker run --log-driver=none --rm -v {}:/data'.format(work_dir).split()
    if sudo:
        base_docker_call = ['sudo'] + base_docker_call
    if java_opts:
        base_docker_call = base_docker_call + ['-e', 'JAVA_OPTS={}'.format(java_opts)]
    try:
        if outfile:
            subprocess.check_call(base_docker_call + [tool] + tool_parameters, stdout=outfile)
        else:
            subprocess.check_call(base_docker_call + [tool] + tool_parameters)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status: {}'
                           ''.format(base_docker_call + [tool] + tool_parameters))
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


def parse_sra(path_to_config):
    """
    Parses genetorrent config file.  Returns list of samples: [ [id1, id1 ], [id2, id2], ... ]
    Returns duplicate of ids to follow UUID/URL standard.
    """
    samples = []
    with open(path_to_config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                samples.append(line.strip())
    return samples


def tarball_files(work_dir, tar_name, uuid=None, files=None):
    """
    Tars a group of files together into a tarball

    work_dir: str       Current Working Directory
    tar_name: str       Name of tarball
    uuid: str           UUID to stamp files with
    files: str(s)       List of filenames to place in the tarball from working directory
    """
    with tarfile.open(os.path.join(work_dir, tar_name), 'w:gz') as f_out:
        for fname in files:
            if uuid:
                f_out.add(os.path.join(work_dir, fname), arcname=uuid + '.' + fname)
            else:
                f_out.add(os.path.join(work_dir, fname), arcname=fname)


# Job Functions
def start_batch(job, input_args):
    """
    This function will administer 5 jobs at a time then recursively call itself until subset is empty
    """
    samples = parse_sra(input_args['sra'])
    # for analysis_id in samples:
    job.addChildJobFn(download_and_transfer_sample, input_args, samples, cores=1, disk='30')


def download_and_transfer_sample(job, input_args, samples):
    """
    Downloads a sample from dbGaP via SRAToolKit, then uses S3AM to transfer it to S3

    input_args: dict        Dictionary of input arguments
    analysis_id: str        An analysis ID for a sample in CGHub
    """
    if len(samples) > 1:
        a = samples[len(samples)/2:]
        b = samples[:len(samples)/2]
        job.addChildJobFn(download_and_transfer_sample, input_args, a, disk='30G')
        job.addChildJobFn(download_and_transfer_sample, input_args, b, disk='30G')
    else:
        analysis_id = samples[0]
        work_dir = job.fileStore.getLocalTempDir()
        sudo = input_args['sudo']
        # Acquire dbgap_key
        shutil.copy(input_args['dbgap_key'], os.path.join(work_dir, 'dbgap.ngc'))
        # Call to fastq-dump to pull down SRA files and convert to fastq
        if input_args['single_end']:
            parameters = [analysis_id]
        else:
            parameters = ['--split-files', analysis_id]
        docker_call(tool='quay.io/ucsc_cgl/fastq-dump:2.5.7--4577a6c1a3c94adaa0c25dd6c03518ee610433d1',
                    work_dir=work_dir, tool_parameters=parameters, sudo=sudo)
        # Collect files and encapsulate into a tarball
        shutil.rmtree(os.path.join(work_dir, 'sra'))
        sample_name = analysis_id + '.tar.gz'
        if input_args['single_end']:
            r = [os.path.basename(x) for x in glob.glob(os.path.join(work_dir, '*.f*'))]
            tarball_files(work_dir, tar_name=sample_name, files=r)
        else:
            r1 = [os.path.basename(x) for x in glob.glob(os.path.join(work_dir, '*_1*'))]
            r2 = [os.path.basename(x) for x in glob.glob(os.path.join(work_dir, '*_2*'))]
            tarball_files(work_dir, tar_name=sample_name, files=r1 + r2)
        # Parse s3_dir to get bucket and s3 path
        key_path = input_args['ssec']
        s3_dir = input_args['s3_dir']
        bucket_name = s3_dir.lstrip('/').split('/')[0]
        base_url = 'https://s3-us-west-2.amazonaws.com/'
        url = os.path.join(base_url, bucket_name, sample_name)
        # Generate keyfile for upload
        with open(os.path.join(work_dir, 'temp.key'), 'wb') as f_out:
            f_out.write(generate_unique_key(key_path, url))
        # Upload to S3 via S3AM
        s3am_command = ['s3am',
                        'upload',
                        '--sse-key-file', os.path.join(work_dir, 'temp.key'),
                        'file://{}'.format(os.path.join(work_dir, sample_name)),
                        's3://' + bucket_name + '/']
        subprocess.check_call(s3am_command)


def main():
    """
    Transfer gTEX data from dbGaP (NCBI) to S3
    """
    # Define Parser object and add to toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # Store inputs from argparse
    inputs = {'sra': args.sra,
              'dbgap_key': args.dbgap_key,
              'ssec': args.ssec,
              's3_dir': args.s3_dir,
              'single_end': args.single_end,
              'sudo': args.sudo}
    # Sanity checks
    if args.ssec:
        assert os.path.isfile(args.ssec)
    if args.sra:
        assert os.path.isfile(args.sra)
    if args.dbgap_key:
        assert os.path.isfile(args.dbgap_key)
    # Start Pipeline
    Job.Runner.startToil(Job.wrapJobFn(start_batch, inputs), args)

if __name__ == '__main__':
    main()