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
from toil.job import Job
from toil_lib.jobs import map_job
from toil_lib.programs import docker_call


def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--genetorrent', default=None, required=True,
                        help='Path to a file with one analysis ID per line for data hosted on CGHub.')
    parser.add_argument('-k', '--genetorrent_key', default=None, required=True,
                        help='Path to a CGHub key that has access to the TCGA data being requested. An exception will'
                             'be thrown if "-g" is set but not this argument.')
    parser.add_argument('--s3_dir', default=None, required=True, help='S3 Bucket. e.g. tcga-data')
    parser.add_argument('--ssec', default=None, required=True, help='Path to Key File for SSE-C Encryption')
    return parser


# Convenience Functions
def generate_unique_key(master_key_path, url):
    """
    master_key_path: str    Path to the 32-byte Master Key (for S3 Encryption)
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


def parse_genetorrent(path_to_config):
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


# Job Functions
def download_and_transfer_sample(job, sample, inputs):

    """
    Downloads a sample from CGHub via GeneTorrent, then uses S3AM to transfer it to S3

    input_args: dict        Dictionary of input arguments
    analysis_id: str        An analysis ID for a sample in CGHub
    """
    analysis_id = sample[0]
    work_dir = job.fileStore.getLocalTempDir()
    folder_path = os.path.join(work_dir, os.path.basename(analysis_id))
    # Acquire genetorrent key and download sample
    shutil.copy(inputs['genetorrent_key'], os.path.join(work_dir, 'cghub.key'))
    parameters = ['-vv', '-c', 'cghub.key', '-d', analysis_id]
    docker_call(tool='quay.io/ucsc_cgl/genetorrent:3.8.7--9911761265b6f08bc3ef09f53af05f56848d805b',
                work_dir=work_dir, parameters=parameters)
    try:
        sample = glob.glob(os.path.join(folder_path, '*tar*'))[0]
    except KeyError as e:
        print 'No tarfile found inside of folder: '.format(e)
        raise
    # Upload sample to S3AM
    key_path = inputs['ssec']
    if sample.endswith('gz'):
        sample_name = analysis_id + '.tar.gz'
        shutil.move(sample, os.path.join(work_dir, sample_name))
    else:
        sample_name = analysis_id + '.tar'
        shutil.move(sample, os.path.join(work_dir, sample_name))
    # Parse s3_dir to get bucket and s3 path
    s3_dir = inputs['s3_dir']
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
    This is a Toil pipeline to transfer TCGA data into an S3 Bucket

    Data is pulled down with Genetorrent and transferred to S3 via S3AM.
    """
    # Define Parser object and add to toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # Store inputs from argparse
    inputs = {'genetorrent': args.genetorrent,
              'genetorrent_key': args.genetorrent_key,
              'ssec': args.ssec,
              's3_dir': args.s3_dir}
    # Sanity checks
    if args.ssec:
        assert os.path.isfile(args.ssec)
    if args.genetorrent:
        assert os.path.isfile(args.genetorrent)
    if args.genetorrent_key:
        assert os.path.isfile(args.genetorrent_key)
    samples = parse_genetorrent(args.genetorrent)
    # Start pipeline
    # map_job accepts a function, an iterable, and *args. The function is launched as a child
    # process with one element from the iterable and *args, which in turn spawns a tree of child jobs.
    Job.Runner.startToil(Job.wrapJobFn(map_job, download_and_transfer_sample, samples, inputs), args)


if __name__ == '__main__':
    main()
