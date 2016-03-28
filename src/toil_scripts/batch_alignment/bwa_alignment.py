#!/usr/bin/env python2.7
# John Vivian
# Date: 11/5/2015

"""
Alignment Pipeline with sorting, header fix, and read group addition.
This pipeline was created from cannabalized portions of Arjun Rao's Precision Immuno Pipeline:
https://github.com/arkal/docker_scripts/blob/master/precision_immuno.py
Thanks bro!

        0 --> 1 --> 2 --> 3 --> 4

0   bwa alignment to a reference
1   samtools sam to bam conversion (and sort)
2   Fix header
3   Add read groups
4   Upload to S3

===================================================================
:Dependencies:
curl            - apt-get install curl
Toil            - pip install --pre toil
Docker          - http://docs.docker.com/engine/installation/

Optional:
S3AM            - pip install --s3am (requires ~/.boto config file)
"""
from __future__ import print_function
import argparse
import base64
from collections import OrderedDict
import hashlib
import multiprocessing
import os
import errno
import subprocess
import shutil
import tarfile
import logging
from toil.job import Job

_log = logging.getLogger(__name__)


def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, add_help=True)
    parser.add_argument('-c', '--config', required=True, help='configuration file. One sample per line: uuid,url,url')
    parser.add_argument('-l', '--lb', required=True, help='the LB (library) entry to go in the BAM header')
    parser.add_argument('-r', '--ref', required=True, help='Reference fasta file')
    parser.add_argument('-m', '--amb', required=True, help='Reference fasta file (amb)')
    parser.add_argument('-n', '--ann', required=True, help='Reference fasta file (ann)')
    parser.add_argument('-b', '--bwt', required=True, help='Reference fasta file (bwt)')
    parser.add_argument('-p', '--pac', required=True, help='Reference fasta file (pac)')
    parser.add_argument('-a', '--sa', required=True, help='Reference fasta file (sa)')
    parser.add_argument('-f', '--fai', required=True, help='Reference fasta file (fai)')
    parser.add_argument('-t', '--alt', required=False, help='Alternate file for reference build (alt). Necessary for alt aware alignment.')
    parser.add_argument('-s', '--ssec', help='Path to Key File for SSE-C Encryption')
    parser.add_argument('-o', '--output_dir', default=None, help='full path where final results will be output')
    parser.add_argument('-u', '--sudo', dest='sudo', action='store_true', help='Docker usually needs sudo to execute '
                                                                               'locally, but not''when running Mesos '
                                                                               'or when a member of a Docker group.')
    parser.add_argument('-3', '--s3_dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                             'cgl-driver-projects/ckcc/rna-seq-samples/')
    parser.add_argument('-k', '--use_bwakit', action='store_true', help='Use bwakit instead of the binary build of bwa')
    parser.add_argument('-se', '--file_size', default='100G', help='Approximate input file size. Should be given as %d[TGMK], e.g., for a 100 gigabyte file, use --file_size 100G')
    parser.add_argument('--trim', action='store_true', help='Trim adapters. Ony works with --use_bwakit.')
    parser.add_argument('--skip_sort', action='store_true', help='Skip sorting. Only works with --use_bwakit.')
    parser.set_defaults(sudo=False)
    return parser


# Convenience functions
def mkdir_p(path):
    """
    It is Easier to Ask for Forgiveness than Permission
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def generate_unique_key(master_key_path, url):
    """
    Given a master key and a url, this function returns a 32-byte unique key

    master_key_path: str    Absolute path to the Master Key (for S3 Encryption)
    url: str                S3 URL (format: https://s3-us-west-2.amazonaws.com/<bucket>/<key.name>)
    """
    with open(master_key_path, 'r') as f:
        master_key = f.read()
    assert len(master_key) == 32, 'Invalid Key! Must be 32 characters. ' \
                                  'Key: {}, Length: {}'.format(master_key, len(master_key))
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is invalid and is not 32 characters: {}'.format(new_key)
    return new_key


def download_encrypted_file(job, url, key_path):
    """
    Downloads encrypted files from S3 via header injection

    url: str        URL to be downloaded
    key_path: str   Path to the master key needed to derive unique encryption keys per file
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, os.path.basename(url))

    with open(key_path, 'r') as f:
        key = f.read()
    if len(key) != 32:
        raise RuntimeError('Invalid Key! Must be 32 bytes: {}'.format(key))

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
    return job.fileStore.writeGlobalFile(file_path)


def download_from_url(job, url):
    """
    Downloads a URL that was supplied as an argument to running this script in LocalTempDir.
    After downloading the file, it is stored in the FileStore.

    url: str        URL to be downloaded. filename is derived from URL
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, os.path.basename(url))
    if not os.path.exists(file_path):
        if url.startswith('s3:'):
            download_from_s3_url(file_path, url)
        else:
            try:
                download_cmd = ['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path]
                _log.info("Downloading file using command %s." % " ".join(download_cmd))
                subprocess.check_call(download_cmd)
            except OSError:
                raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    return job.fileStore.writeGlobalFile(file_path)


def return_input_paths(job, work_dir, ids, *args):
    """
    Returns the paths of files from the FileStore if they are not present.
    This function should be unpacked for every item being returned, unless none
    of the paths are needed in which case it should be unassigned.

    work_dir: str       Path to the current working directory
    ids: dict           Dictionary of fileStore IDs, accessed by filename (str)
    *args: str(s)       Files to be retrieved and placed in the current working directory
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


def docker_call(work_dir,
                tool_parameters,
                tool,
                java_opts=None,
                outfile=None,
                sudo=False,
                docker_parameters=None,
                check_output=False,
                no_rm=False):
    """
    Makes subprocess call of a command to a docker container.

    tool_parameters: list   An array of the parameters to be passed to the tool
    tool: str               Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools)
    java_opts: str          Optional commands to pass to a java jar execution. (e.g. '-Xmx15G')
    outfile: file           Filehandle that stderr will be passed to
    sudo: bool              If the user wants the docker command executed as sudo
    """
    rm = '--rm'
    if no_rm:
        rm = ''

    base_docker_call = ('docker run %s --log-driver=none -v %s:/data' % (rm, work_dir)).split()

    if sudo:
        base_docker_call = ['sudo'] + base_docker_call
    if java_opts:
        base_docker_call = base_docker_call + ['-e', 'JAVA_OPTS={}'.format(java_opts)]
    if docker_parameters:
        base_docker_call = base_docker_call + docker_parameters

    _log.warn("Calling docker with %s." % " ".join(base_docker_call + [tool] + tool_parameters))

    try:
        if outfile:
            subprocess.check_call(base_docker_call + [tool] + tool_parameters, stdout=outfile)
        else:
            if check_output:
                return subprocess.check_output(base_docker_call + [tool] + tool_parameters)
            else:
                subprocess.check_call(base_docker_call + [tool] + tool_parameters)

    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


def copy_to_output_dir(work_dir, output_dir, uuid=None, files=()):
    """
    A list of files to move from work_dir to output_dir.

    work_dir: str       Current working directory
    output_dir: str     Desired output directory
    uuid: str           UUID to be prepended to output files. Optional.
    files: list         List of files (by filename) to be copied to the desired output directory
    """
    for fname in files:
        if uuid is None:
            shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


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


def download_from_s3_url(file_path, url):
    from urlparse import urlparse
    from boto.s3.connection import S3Connection
    s3 = S3Connection()
    try:
        parsed_url = urlparse(url)
        if not parsed_url.netloc or not parsed_url.path.startswith('/'):
            raise RuntimeError("An S3 URL must be of the form s3:/BUCKET/ or "
                               "s3://BUCKET/KEY. '%s' is not." % url)
        bucket = s3.get_bucket(parsed_url.netloc)
        key = bucket.get_key(parsed_url.path[1:])
        _log.info("From URL %s, parsed --> %s, bucket %s, key %s", url, parsed_url, bucket, key)
        key.get_contents_to_filename(file_path)
    finally:
        s3.close()


# Job Functions
def download_shared_files(job, input_args):
    """
    Downloads shared files that are used by all samples for alignment and places them in the jobstore.

    input_args: dict        Input arguments (passed from main())
    """
    shared_files = ['ref.fa', 'ref.fa.amb', 'ref.fa.ann', 'ref.fa.bwt', 'ref.fa.pac', 'ref.fa.sa', 'ref.fa.fai']

    if input_args['ref.fa.alt']:
        shared_files.append('ref.fa.alt')

    shared_ids = {}
    for fname in shared_files:
        url = input_args[fname]
        shared_ids[fname] = job.addChildJobFn(download_from_url, url).rv()
    job.addFollowOnJobFn(parse_config, shared_ids, input_args)


def parse_config(job, shared_ids, input_args):
    """
    Stores the UUID and urls associated with the input files to be retrieved.
    Configuration file has one sample per line, with the following format:  UUID,1st_url,2nd_url

    shared_ids: dict        Dictionary of fileStore IDs for the shared files downloaded in the previous step
    input_args: dict        Input argumentts
    """
    samples = []
    config = input_args['config']

    # does the config file exist locally? if not, try to read from job store
    if not os.path.exists(config):

        work_dir = job.fileStore.getLocalTempDir()
        config_path = os.path.join(work_dir, 'config.txt')
        job.fileStore.readGlobalFile(config, config_path)
        config = config_path

    with open(config, 'r') as f_in:
        for line in f_in:
            if line.strip():
                line = line.strip().split(',')
                assert len(line) == 3, 'Improper formatting'
                uuid = line[0]
                urls = line[1:]
                samples.append((uuid, urls))
    input_args['cpu_count'] = multiprocessing.cpu_count()
    job_vars = (input_args, shared_ids)
    for sample in samples:
        job.addChildJobFn(download_inputs, job_vars, sample, cores=input_args['cpu_count'], memory='20 G', disk=input_args['file_size'])


def download_inputs(job, job_vars, sample):
    """
    Downloads the sample inputs (R1.fq.gz and R2.fq.gz)

    job_vars: tuple         Contains the dictionaries: input_args and ids
    sample: tuple           Contains the uuid (str) and urls (list of strings)
    """
    input_args, ids = job_vars
    uuid, urls = sample
    input_args['uuid'] = uuid
    for i in xrange(2):
        if input_args['ssec']:
            key_path = input_args['ssec']
            ids['r{}.fq.gz'.format(i+1)] = job.addChildJobFn(download_encrypted_file, urls[i], key_path, disk=input_args['file_size']).rv()
        else:
            ids['r{}.fq.gz'.format(i+1)] = job.addChildJobFn(download_from_url, urls[i], disk=input_args['file_size']).rv()
    job.addFollowOnJobFn(static_dag_declaration, job_vars)


def static_dag_declaration(job, job_vars):
    """
    Unnecessary static declaration, but done to test out the promise mechanism

    job_vars: tuple         Contains the dictionaries: input_args and ids
    """
    # Define
    input_args, ids = job_vars
    cores = input_args['cpu_count']

    # if we run with bwakit, we skip conversion and reheadering
    if input_args['use_bwakit']:

        bwa = job.wrapJobFn(run_bwa, job_vars, cores=cores, disk=input_args['file_size'])
        
        # Link
        job.addChild(bwa)
    
    else:

        bwa = job.wrapJobFn(run_bwa, job_vars, cores=cores, disk=input_args['file_size'])
        conversion = job.wrapJobFn(bam_conversion, job_vars, bwa.rv(), disk=input_args['file_size'])
        header = job.wrapJobFn(fix_bam_header, job_vars, conversion.rv(), disk=input_args['file_size'])
        rg = job.wrapJobFn(add_readgroups, job_vars, header.rv(), memory='15G', disk=input_args['file_size'])
        
        # Link
        job.addChild(bwa)
        bwa.addChild(conversion)
        conversion.addChild(header)
        header.addChild(rg)


def copy_or_upload(work_dir, input_args, output_file):
    """
    Looks at what input args are defined and uploads the file to it's proper location.
    """

    # Either write file to local output directory or upload to S3 cloud storage
    if input_args['output_dir']:
        copy_to_output_dir(work_dir=work_dir,
                           output_dir=input_args['output_dir'],
                           files=[output_file])

    if input_args['s3_dir']:
        upload_to_s3(work_dir,
                     input_args,
                     output_file)


def run_bwa(job, job_vars):
    """
    This module aligns two fastqs into a SAMFILE using Burrows-Wheeler Alignment

    job_vars: tuple         Contains the dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    cores = input_args['cpu_count']
    # Retrieve samples
    return_input_paths(job, work_dir, ids, 'r1.fq.gz', 'r2.fq.gz')
    # Retrieve input files
    return_input_paths(job, work_dir, ids, 'ref.fa', 'ref.fa.amb', 'ref.fa.ann',
                       'ref.fa.bwt', 'ref.fa.pac', 'ref.fa.sa', 'ref.fa.fai')

    if input_args['ref.fa.alt']:
        return_input_paths(job, work_dir, ids, 'ref.fa.alt')

    if input_args['use_bwakit']:

        # add a read group line
        # this is an illumina only pipeline, therefore, hardcoding illumina is OK
        # program unit is required by GATK, but we don't know it
        # so, fake a bad value
        uuid = input_args['uuid']
        library = input_args['lb']
        rg = "@RG\\tID:%s\\tLB:%s\\tPL:ILLUMINA\\tPU:12345\\tSM:%s" % (uuid,
                                                                  library,
                                                                  uuid)
        
        # do we want to use bwakit to sort?
        opt_args = []
        if input_args['sort']:
            opt_args.append('-s')

        # do we want to use bwakit to trim adapters?
        if input_args['trim']:
            opt_args.append('-a')

        # Call: bwakit
        parameters = (['-t%d' % cores,
                       '-R%s' % rg] +
                      opt_args +
                      ['-o', '/data/aligned',
                       '/data/ref.fa',
                       '/data/r1.fq.gz',
                       '/data/r2.fq.gz'])

        docker_call(tool='quay.io/ucsc_cgl/bwakit:0.7.12--528bb9bf73099a31e74a7f5e6e3f2e0a41da486e',
                    tool_parameters=parameters, work_dir=work_dir, sudo=sudo)

        # bwa insists on adding an `*.aln.sam` suffix, so rename the output file
        output_file = '{}.bam'.format(uuid)
        os.rename(os.path.join(work_dir, 'aligned.aln.bam'),
                  os.path.join(work_dir, output_file))
        
        # Either write file to local output directory or upload to S3 cloud storage
        copy_or_upload(work_dir, input_args, output_file)
        
    else:
        # Call: BWA
        parameters = ["mem",
                      "-t", str(cores),
                      "/data/ref.fa"] + ['/data/r1.fq.gz', '/data/r2.fq.gz']

        with open(os.path.join(work_dir, 'aligned.sam'), 'w') as samfile:
            docker_call(tool='quay.io/ucsc_cgl/bwa:0.7.12--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                        tool_parameters=parameters, work_dir=work_dir, outfile=samfile, sudo=sudo)

        # write aligned sam file to global file store
        return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'aligned.sam'))


def bam_conversion(job, job_vars, sam_id):
    """
    This module converts SAMFILE to BAMFILE using samtools

    job_vars: tuple         Contains the dictionaries: input_args and ids
    sam_id: str             FileStore ID for the samfile
    """
    # Unpack variables
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # Retrieve input file
    job.fileStore.readGlobalFile(sam_id, os.path.join(work_dir, 'aligned.sam'))
    # Call: Samtools
    parameters = ['view',
                  '-bS',
                  '/data/aligned.sam']

    with open(os.path.join(work_dir, 'aligned.bam'), 'w') as bamfile:
        docker_call(tool='quay.io/ucsc_cgl/samtools:1.2--35ac87df5b21a8e8e8d159f26864ac1e1db8cf86',
                    tool_parameters=parameters, work_dir=work_dir, outfile=bamfile, sudo=sudo)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'aligned.bam'))


def fix_bam_header(job, job_vars, bam_id):
    """
    This module modified the header in BAMFILE using samtools

    job_vars: tuple         Contains the dictionaries: input_args and ids
    bam_id: str             FileStore ID for the BAMFILE
    """
    # Unpack variables
    work_dir = job.fileStore.getLocalTempDir()
    input_args, ids = job_vars
    sudo = input_args['sudo']
    # Retrieve input file
    try:
        job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'aligned.bam'))
    except:
        # remove and retry?
        try:
            os.remove(os.path.join(work_dir, 'aligned.bam'))
        except:
            pass
        
        job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'aligned.bam'))
    # Call: Samtools
    parameters = ['view',
                  '-H',
                  '/data/aligned.bam']

    with open(os.path.join(work_dir, 'aligned_bam.header'), 'w') as headerfile:
        docker_call(tool='quay.io/ucsc_cgl/samtools:1.2--35ac87df5b21a8e8e8d159f26864ac1e1db8cf86',
                    tool_parameters=parameters, work_dir=work_dir, outfile=headerfile, sudo=sudo)
    with open(headerfile.name, 'r') as headerfile, open('/'.join([work_dir, 'output_bam.header']), 'w') as outheaderfile:
            for line in headerfile:
                if line.strip().startswith('@PG'):
                    line = '\t'.join([x for x in line.strip().split('\t') if not x.startswith('CL')])
                print(line.strip(), file=outheaderfile)
    # Call: Samtools
    parameters = ['reheader',
                  '/data/output_bam.header',
                  '/data/aligned.bam']
    # Memory leak, crashes
    with open(os.path.join(work_dir, 'aligned_fixPG.bam'), 'w') as fixpg_bamfile:
        docker_call(tool='quay.io/ucsc_cgl/samtools:1.2--35ac87df5b21a8e8e8d159f26864ac1e1db8cf86',
                    tool_parameters=parameters, work_dir=work_dir, outfile=fixpg_bamfile, sudo=sudo)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'aligned_fixPG.bam'))


def add_readgroups(job, job_vars, bam_id):
    """
    This module adds the appropriate read groups to the BAMFILE using picard-tools

    job_vars: tuple         Contains the dictionaries: input_args and ids
    ban_id: str             FileStore ID for the BAMFILE from the previous step
    """
    # Unpack variables
    work_dir = job.fileStore.getLocalTempDir()
    input_args, ids = job_vars
    sudo = input_args['sudo']
    uuid = input_args['uuid']
    output_file = '{}.bam'.format(uuid)
    # Retrieve input file
    job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'aligned_fixpg.bam'))
    _log.warn("Read global file to %s, adding read groups and saving to %s." % (os.path.join(work_dir, 'aligned_fixpg.bam'), os.path.join(work_dir, output_file)))

    # Call: Samtools
    parameters = ['AddOrReplaceReadGroups',
                  'I=/data/aligned_fixpg.bam',
                  'O=/data/{}'.format(output_file),
                  'SO=coordinate',
                  'ID={}'.format(uuid),
                  'LB={}'.format(input_args['lb']),
                  'PL=ILLUMINA',
                  'PU=12345',
                  'SM={}'.format(uuid)]
    docker_call(tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                tool_parameters=parameters, work_dir=work_dir, java_opts='-Xmx15G', sudo=sudo)

    # Either write file to local output directory or upload to S3 cloud storage
    copy_or_upload(work_dir, input_args, output_file)
        

def upload_to_s3(work_dir, input_args, output_file):
    """
    Uploads a file to S3 via S3AM with encryption.  This pipeline assumes all
    BAMS contain patient data and must be encrypted upon upload.

    job_vars: tuple         Contains the dictionaries: input_args and ids
    """

    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.lstrip('/').split('/')[0]
    bucket_dir = '/'.join(s3_dir.lstrip('/').split('/')[1:])

    # what is the path we are uploading to?
    s3_path = os.path.join('s3://', bucket_name, bucket_dir, output_file)

    # how many cores should we use?
    cores = input_args['cpu_count']

    # does this need to be uploaded with encryption?
    if input_args['ssec']:
        key_path = input_args['ssec']

        base_url = 'https://s3-us-west-2.amazonaws.com'
        url = os.path.join(base_url, bucket_name, bucket_dir, output_file)

        # Generate keyfile for upload
        uuid = input_args['uuid']
        with open(os.path.join(work_dir, uuid + '.key'), 'wb') as f_out:
            f_out.write(generate_unique_key(key_path, url))

        # Upload to S3 via S3AM
        s3am_command = ['s3am',
                        'upload',
                        '--sse-key-file', os.path.join(work_dir, uuid + '.key'),
                        '--part-size=64M',
                        '--upload-slots={}'.format(cores),
                        '--download-slots={}'.format(cores),
                        'file://{}'.format(os.path.join(work_dir, output_file)),
                        s3_path]
    else:
        # Upload to S3 via S3AM
        s3am_command = ['s3am',
                        'upload',
                        '--part-size=64M',
                        '--upload-slots={}'.format(cores),
                        '--download-slots={}'.format(cores),
                        'file://{}'.format(os.path.join(work_dir, output_file)),
                        s3_path]

    # run upload
    _log.info("Calling s3am with %s" % " ".join(s3am_command))
    try:
        subprocess.check_call(s3am_command)
    except:
        _log.error("Upload to %s failed. Cancelling..." % s3_path)
        s3am_cancel = ['s3am', 'cancel', s3_path]
        
        try:
            subprocess.check_call(s3am_cancel)
        except:
            _log.error("Cancelling upload with '%s' failed." % " ".join(s3am_cancel))

        raise


def main():
    """
    This is a Toil pipeline used to perform alignment of fastqs.
    It uses BWA to align, samtools to sort and fix the header,
    and picard to add the read group (RG) information to the BAM
    that is in a format compatible with GATK.

    Config: "bwa_config.csv" is a required file that stores sample information to be run.
    Script: "launch_bwa_hg19.sh" provides all of the inputs needed to run this script locally.
    Script: "bwa_launch_mesos.sh" provides all of the inputs needed to run this script on a mesos cluster in AWS
    """
    # Define Parser object and add to Toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    # validate bwakit specific args
    arg_errors = ["Argument validation failed:"]
    if not args.use_bwakit:
        if args.alt:
            arg_errors.append('Alternate haplotype build supplied, but not using bwakit.')
        if args.skip_sort:
            arg_errors.append('--skip_sort can only be used with --use_bwakit.')
        if args.trim:
            arg_errors.append('--trim can only be used with --use_bwakit.')
    if len(arg_errors) > 1:
        raise ValueError("\n".join(arg_errors))

    # Store input parameters in a dictionary
    inputs = {'config': args.config,
              'ref.fa': args.ref,
              'ref.fa.amb': args.amb,
              'ref.fa.ann': args.ann,
              'ref.fa.bwt': args.bwt,
              'ref.fa.pac': args.pac,
              'ref.fa.sa': args.sa,
              'ref.fa.fai': args.fai,
              'ref.fa.alt': args.alt,
              'ssec': args.ssec,
              'output_dir': args.output_dir,
              'sudo': args.sudo,
              's3_dir': args.s3_dir,
              'lb': args.lb,
              'uuid': None,
              'cpu_count': None,
              'file_size': args.file_size,
              'use_bwakit': args.use_bwakit,
              'sort': not args.skip_sort,
              'trim': args.trim}

    # Launch Pipeline
    Job.Runner.startToil(Job.wrapJobFn(download_shared_files, inputs), args)


if __name__ == '__main__':
    main()
