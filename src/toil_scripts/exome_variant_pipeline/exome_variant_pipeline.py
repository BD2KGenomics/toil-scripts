#!/usr/bin/env python2.7
# John Vivian
# Fall-2015

"""
Toil pipeline for exome variant analysis

Tree Structure of variant pipeline (per sample)

                        /-------> 14 --\
   0------------> 5 ------------> 15 ------> 17 ---> 18
/ / \ \          / \    \-------> 16 --/
1 2 3 4         6   7
                |   |
                8   9
                |   |
                10  11
                |   |
                12  13

0 = Start node
1 = reference index
2 = reference dict
3 = normal bam index
4 = tumor bam index
5 = pre-processing node / DAG declaration
6,7 = RealignerTargetCreator
8,9 = IndelRealigner
10,11 = BaseRecalibration
12,13 = PrintReads
14 = MuTect
15 = Pindel
16 = MuSe
17 = Consolidate Output
18 = Upload to S3
===================================================================
:Dependencies:
curl            - apt-get install curl
docker          - https://docs.docker.com/linux/step_one/
Toil            - pip install toil

Optional:
S3AM            - pip install --pre s3am (requires ~/.boto)
"""
import argparse
import base64
from collections import OrderedDict
from contextlib import closing
import glob
import hashlib
import os
import shutil
import subprocess
import multiprocessing
import errno
import tarfile
from toil.job import Job


def build_parser():
    """
    Contains argparse arguments
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r', '--reference', required=True, help="Reference Genome URL")
    parser.add_argument('-f', '--config', required=True, help="Each line contains (CSV): UUID,Normal_URL,Tumor_URL")
    parser.add_argument('-p', '--phase', required=True, help='1000G_phase1.indels.hg19.sites.fixed.vcf URL')
    parser.add_argument('-m', '--mills', required=True, help='Mills_and_1000G_gold_standard.indels.hg19.sites.vcf URL')
    parser.add_argument('-d', '--dbsnp', required=True, help='dbsnp_132_b37.leftAligned.vcf URL')
    parser.add_argument('-c', '--cosmic', required=True, help='b37_cosmic_v54_120711.vcf URL')
    parser.add_argument('-o', '--output_dir', default=None, help='Full path to final output dir')
    parser.add_argument('-s', '--ssec', help='A key that can be used to fetch encrypted data')
    parser.add_argument('-3', '--s3_dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                             'cgl-driver-projects/ckcc/rna-seq-samples/')
    parser.add_argument('-u', '--sudo', dest='sudo', action='store_true', help='Docker usually needs sudo to execute '
                                                                               'locally, but not''when running Mesos '
                                                                               'or when a member of a Docker group.')
    return parser


# Convenience functions used in the pipeline
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


def flatten(x):
    """
    Flattens a nested array into a single list

    x: list/tuple       The nested list/tuple to be flattened.
    """
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


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


def download_encrypted_file(job, input_args, name):
    """
    Downloads encrypted files from S3 via header injection

    input_args: dict    Input dictionary defined in main()
    name: str           Symbolic name associated with file
    """
    work_dir = job.fileStore.getLocalTempDir()
    key_path = input_args['ssec']
    file_path = os.path.join(work_dir, name)
    url = input_args[name]
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


def copy_to_output_dir(work_dir, output_dir, uuid=None, files=None):
    """
    A list of files to move from work_dir to output_dir.

    work_dir: str       Current working directory
    output_dir: str     Output directory for files to go
    uuid: str           UUID to "stamp" onto output files
    files: list         List of files to iterate through
    """
    for fname in files:
        if uuid is None:
            shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


def download_from_url(job, url, name):
    """
    Simple curl request made for a given url

    url: str    URL to download
    name: str   Name to give downloaded file
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, name)
    if not os.path.exists(file_path):
        if url.startswith('s3:'):
            download_from_s3_url(file_path, url)
        else:
            try:
                subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
            except subprocess.CalledProcessError:
                raise RuntimeError(
                    '\nNecessary file could not be acquired: {}. Check input URL'.format(url))
            except OSError:
                raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    return job.fileStore.writeGlobalFile(file_path)


def return_input_paths(job, work_dir, ids, *args):
    """
    Given one or more strings representing file_names, return the paths to those files. Each item must be unpacked!

    work_dir: str       Current working directory
    ids: dict           Dictionary of fileStore IDs
    *args: str(s)       for every file in *args, place file in work_dir via FileStore
    """
    paths = OrderedDict()
    for name in args:
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = job.fileStore.readGlobalFile(ids[name], os.path.join(work_dir, name))
        else:
            file_path = name
        paths[name] = file_path
        if len(args) == 1:
            paths = file_path

    return paths


def docker_path(file_path):
    """
    Returns the path internal to the docker container (for standard reasons, this is always /data)
    """
    return os.path.join('/data', os.path.basename(file_path))


def docker_call(work_dir, tool_parameters, tool, java_opts=None, outfile=None, sudo=False):
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
    except subprocess.CalledProcessError, e:
        raise RuntimeError('docker command returned a non-zero exit status. {}'.format(e))
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


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


def get_mean_insert_size(work_dir, bam_path):
    cmd = "sudo docker run --log-driver=none --rm -v {}:/data quay.io/ucsc_cgl/samtools " \
          "view -f66 {}".format(work_dir, bam_path)
    process = subprocess.Popen(args=cmd, shell=True, stdout=subprocess.PIPE)
    b_sum = 0L
    b_count = 0L
    while True:
        line = process.stdout.readline()
        if not line:
            break
        tmp = line.split("\t")
        if abs(long(tmp[8])) < 10000:
            b_sum += abs(long(tmp[8]))
            b_count += 1
    process.wait()
    mean = b_sum / b_count
    print "Using insert size: %d" % mean
    return int(mean)


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
        key.get_contents_to_filename(file_path)
    finally:
        s3.close()


# Start of Job Functions
def download_shared_files(job, input_args):
    """
    Downloads files shared by all samples in the pipeline

    input_args: dict        Dictionary of input arguments (from main())
    """
    shared_ids = {}
    for fname in ['ref.fasta', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf', 'cosmic.vcf']:
        shared_ids[fname] = job.addChildJobFn(download_from_url, url=input_args[fname], name=fname).rv()
    job.addFollowOnJobFn(reference_preprocessing, input_args, shared_ids)


def reference_preprocessing(job, input_args, shared_ids):
    """
    Create index and dict file for reference

    input_args: dict        Dictionary of input argumnets
    shared_ids: dict        Dictionary of fileStore IDs
    """
    ref_id = shared_ids['ref.fasta']
    sudo = input_args['sudo']
    shared_ids['ref.fasta.fai'] = job.addChildJobFn(create_reference_index, ref_id, sudo).rv()
    shared_ids['ref.dict'] = job.addChildJobFn(create_reference_dict, ref_id, sudo).rv()
    job.addFollowOnJobFn(spawn_batch_jobs, input_args, shared_ids)


def create_reference_index(job, ref_id, sudo):
    """
    Uses Samtools to create reference index file (.fasta.fai)

    ref_id: str     The fileStore ID of the reference
    sudo: bool      Boolean item to determine whether to invoke sudo with docker
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path to reference
    ref_path = docker_path(job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta')))
    # Call: Samtools
    command = ['faidx', ref_path]
    docker_call(work_dir=work_dir, tool_parameters=command,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e', sudo=sudo)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.fasta.fai'))


def create_reference_dict(job, ref_id, sudo):
    """
    Uses Picardtools to create reference dictionary (.dict) for the sample

    ref_id: str     The fileStore ID of the reference
    sudo: bool      Boolean item to determine whether to invoke sudo with docker
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta'))
    # Call: picardtools
    output = os.path.splitext(docker_path(ref_path))[0]
    command = ['CreateSequenceDictionary', 'R={}'.format(docker_path(ref_path)), 'O={}.dict'.format(output)]
    docker_call(work_dir=work_dir, tool_parameters=command,
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e', sudo=sudo)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.dict'))


def spawn_batch_jobs(job, input_args, shared_ids):
    """
    Spawn a pipeline for each sample in the configuration file

    input_args: dict        Dictionary of input argumnets
    shared_ids: dict        Dictionary of fileStore IDs
    """
    samples = []
    config = input_args['config']
    with open(config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                samples.append(line.strip().split(','))
    for sample in samples:
        job.addChildJobFn(download_samples, shared_ids, input_args, sample)


def download_samples(job, ids, input_args, sample):
    """
    Defines sample variables then downloads the sample.

    ids: dict           Dictionary of fileStore IDs
    input_args: dict    Dictionary of input arguments
    sample: str         Contains uuid, normal url, and tumor url
    """
    uuid, normal_url, tumor_url = sample
    # Create a unique
    sample_input = dict(input_args)
    sample_input['uuid'] = uuid
    sample_input['normal.bam'] = normal_url
    sample_input['tumor.bam'] = tumor_url
    sample_input['cpu_count'] = multiprocessing.cpu_count()
    if sample_input['output_dir']:
        sample_input['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    # Download sample bams and launch pipeline
    if input_args['ssec']:
        ids['normal.bam'] = job.addChildJobFn(download_encrypted_file, sample_input, 'normal.bam').rv()
        ids['tumor.bam'] = job.addChildJobFn(download_encrypted_file, sample_input, 'tumor.bam').rv()
    else:
        ids['normal.bam'] = job.addChildJobFn(download_from_url, url=sample_input['normal.bam'], name='normal.bam').rv()
        ids['tumor.bam'] = job.addChildJobFn(download_from_url, url=sample_input['tumor.bam'], name='tumor.bam').rv()
    job_vars = (sample_input, ids)
    job.addFollowOnJobFn(pipeline_launchpoint, job_vars)


def pipeline_launchpoint(job, job_vars):
    """
    Statically link the rest of the workflow

    job_vars: tuple     Contains the input_args and ids dictionaries
    """
    pre_processing = job.wrapJobFn(index_bams, job_vars).encapsulate()
    run_mutect = job.wrapJobFn(mutect, job_vars, pre_processing.rv(), disk='75G')
    # run_pindel = job.wrapJobFn(pindel, job_vars, pre_processing.rv())
    # run_muse = job.wrapJobFn(muse, job_vars, pre_processing.rv())
    consolidate = job.wrapJobFn(consolidate_output, job_vars, run_mutect.rv())#, run_pindel.rv(), run_muse.rv())
    # Wire up DAG
    job.addChild(pre_processing)
    pre_processing.addChild(run_mutect)
    # pre_processing.addChild(run_pindel)
    # pre_processing.addChild(run_muse)
    pre_processing.addFollowOn(consolidate)


def index_bams(job, job_vars):
    """
    Create index (.bai) files for each sample bam

    job_vars: tuple     Contains the input_args and ids dictionaries
    """
    normal_ids = job.addChildJobFn(index, job_vars, 'normal', cores=1, memory='1 G', disk='8 G').rv()
    tumor_ids = job.addChildJobFn(index, job_vars, 'tumor', cores=1, memory='1 G', disk='8 G').rv()
    return normal_ids, tumor_ids


def index(job, job_vars, sample):
    """
    Runs samtools index to create (.bai) files

    job_vars: tuple     Contains the input_args and ids dictionaries
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    cores = int(input_args['cpu_count'])
    # Retrieve file path
    bam = '{}.bam'.format(sample)
    path = return_input_paths(job, work_dir, ids, bam)
    # Call: index the normal.bam
    parameters = ['index', '{}'.format(docker_path(path))]
    docker_call(work_dir=work_dir, tool_parameters=parameters,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e', sudo=sudo)
    # Write to fileStore
    ids[bam + '.bai'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, bam) + '.bai')
    return job.addChildJobFn(realigner_target_creator, job_vars, sample, cores=cores, memory='10 G', disk='15 G').rv()


def realigner_target_creator(job, job_vars, sample):
    """
    Creates <type>.intervals file needed for indel realignment

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.bam'.format(sample), 'ref.fasta.fai', 'ref.dict',
                       '{}.bam.bai'.format(sample), 'phase.vcf', 'mills.vcf')

    # Output file path
    output = os.path.join(work_dir, '{}.intervals'.format(sample))
    # Call: GATK -- RealignerTargetCreator
    parameters = ['-T', 'RealignerTargetCreator',
                  '-nt', str(cores),
                  '-R', '/data/ref.fasta',
                  '-I', '/data/{}.bam'.format(sample),
                  '-known', '/data/phase.vcf',
                  '-known', '/data/mills.vcf',
                  '--downsampling_type', 'NONE',
                  '-o', docker_path(output)]

    docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                work_dir=work_dir, tool_parameters=parameters, java_opts='-Xmx10g', sudo=sudo)
    # Write to fileStore
    ids['{}.intervals'.format(sample)] = job.fileStore.writeGlobalFile(output)
    return job.addChildJobFn(indel_realignment, job_vars, sample, cores=1, memory='10 G', disk='30 G').rv()


def indel_realignment(job, job_vars, sample):
    """
    Creates realigned bams using <sample>.intervals file from previous step

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.bam'.format(sample), 'phase.vcf', 'mills.vcf',
                       '{}.intervals'.format(sample), 'ref.fasta.fai', 'ref.dict', '{}.bam.bai'.format(sample))
    # Output file path
    output = os.path.join(work_dir, '{}.indel.bam'.format(sample))
    # Call: GATK -- IndelRealigner
    parameters = ['-T', 'IndelRealigner',
                  '-R', '/data/ref.fasta',
                  '-I', '/data/{}.bam'.format(sample),
                  '-known', '/data/phase.vcf',
                  '-known', '/data/mills.vcf',
                  '-targetIntervals', '/data/{}.intervals'.format(sample),
                  '--downsampling_type', 'NONE',
                  '-maxReads', str(720000),
                  '-maxInMemory', str(5400000),
                  '-o', docker_path(output)]
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                work_dir=work_dir, tool_parameters=parameters, java_opts='-Xmx10g', sudo=sudo)
    # Write to fileStore
    ids['{}.indel.bam'.format(sample)] = job.fileStore.writeGlobalFile(output)
    ids['{}.indel.bai'.format(sample)] = job.fileStore.writeGlobalFile(os.path.splitext(output)[0] + '.bai')
    return job.addChildJobFn(base_recalibration, job_vars, sample, cores=cores, memory='15 G', disk='15 G').rv()


def base_recalibration(job, job_vars, sample):
    """
    Creates recal table to perform Base Quality Score Recalibration

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.indel.bam'.format(sample), 'dbsnp.vcf', 'ref.fasta.fai',
                       'ref.dict', '{}.indel.bai'.format(sample))
    # Output file path
    output = os.path.join(work_dir, '{}.recal.table'.format(sample))
    # Call: GATK -- IndelRealigner
    parameters = ['-T', 'BaseRecalibrator',
                  '-nct', str(cores),
                  '-R', '/data/ref.fasta',
                  '-I', '/data/{}.indel.bam'.format(sample),
                  '-knownSites', '/data/dbsnp.vcf',
                  '-o', docker_path(output)]
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                work_dir=work_dir, tool_parameters=parameters, java_opts='-Xmx15g', sudo=sudo)
    # Write to fileStore
    ids['{}.recal.table'.format(sample)] = job.fileStore.writeGlobalFile(output)
    return job.addChildJobFn(print_reads, job_vars, sample, cores=cores, memory='15 G', disk='40 G').rv()


def print_reads(job, job_vars, sample):
    """
    Create bam that has undergone Base Quality Score Recalibration (BQSR)

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.indel.bam'.format(sample), 'ref.fasta.fai',
                       'ref.dict', '{}.indel.bai'.format(sample), '{}.recal.table'.format(sample))
    # Output file
    output = os.path.join(work_dir, '{}.bqsr.bam'.format(sample))
    # Call: GATK -- PrintReads
    parameters = ['-T', 'PrintReads',
                  '-nct', str(cores),
                  '-R', '/data/ref.fasta',
                  '--emit_original_quals',
                  '-I', '/data/{}.indel.bam'.format(sample),
                  '-BQSR', '/data/{}.recal.table'.format(sample),
                  '-o', docker_path(output)]
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                work_dir=work_dir, tool_parameters=parameters, java_opts='-Xmx15g', sudo=sudo)
    # Write to fileStore
    bam_id = job.fileStore.writeGlobalFile(output)
    bai_id = job.fileStore.writeGlobalFile(os.path.splitext(output)[0] + '.bai')
    return bam_id, bai_id


def mutect(job, job_vars, bam_ids):
    """
    Calls MuTect to perform variant analysis

    job_vars: tuple     Contains the input_args and ids dictionaries
    bam_ids: tuple      Contains a tuple of normal/tumor fileStore ids for bams and index files (bai)
    """
    # Unpack convenience variables for job
    normal_ids, tumor_ids = bam_ids
    normal_bam_id, normal_bai_id = normal_ids
    tumor_bam_id, tumor_bai_id = tumor_ids
    input_args, ids = job_vars
    sudo = input_args['sudo']
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve input files
    job.fileStore.readGlobalFile(normal_bam_id, os.path.join(work_dir, 'normal.bam'))
    job.fileStore.readGlobalFile(normal_bai_id, os.path.join(work_dir, 'normal.bai'))
    job.fileStore.readGlobalFile(tumor_bam_id, os.path.join(work_dir, 'tumor.bam'))
    job.fileStore.readGlobalFile(tumor_bai_id, os.path.join(work_dir, 'tumor.bai'))
    return_input_paths(job, work_dir, ids, 'ref.fasta', 'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict', 'cosmic.vcf')
    # Output VCF
    uuid = input_args['uuid']
    output_name = uuid + '.vcf'
    mut_out = docker_path(os.path.join(work_dir, uuid + '.out'))
    mut_cov = docker_path(os.path.join(work_dir, uuid + '.cov'))
    # Call: MuTect
    parameters = ['--analysis_type', 'MuTect',
                  '--reference_sequence', 'ref.fasta',
                  '--cosmic', 'cosmic.vcf',
                  '--dbsnp', 'dbsnp.vcf',
                  '--input_file:normal', 'normal.bam',
                  '--input_file:tumor', 'tumor.bam',
                  '--tumor_lod', str(10),
                  '--initial_tumor_lod', str(4.0),
                  '--out', mut_out,
                  '--coverage_file', mut_cov,
                  '--vcf', docker_path(output_name)]
    docker_call(work_dir=work_dir, tool_parameters=parameters,
                tool='quay.io/ucsc_cgl/mutect:1.1.7--e8bf09459cf0aecb9f55ee689c2b2d194754cbd3', sudo=sudo)
    # Tarball files
    tarball_files(work_dir, 'mutect.tar.gz', files=[uuid + '.vcf', uuid + '.cov', uuid + '.out'])
    # Return fileStore ID
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'mutect.tar.gz'))


def pindel(job, job_vars, bam_ids):
    """
    Calls Pindel to compute indels / deletions in

    job_vars: tuple     Contains the input_args and ids dictionaries
    bam_ids: tuple      Contains a tuple of normal/tumor fileStore ids for bams and index files (bai)
    """
    # Unpack convenience variables for job
    normal_ids, tumor_ids = bam_ids
    normal_bam_id, normal_bai_id = normal_ids
    tumor_bam_id, tumor_bai_id = tumor_ids
    input_args, ids = job_vars
    sudo = input_args['sudo']
    uuid = input_args['uuid']
    cores = input_args['cpu_count']
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve input files
    job.fileStore.readGlobalFile(normal_bam_id, os.path.join(work_dir, 'normal.bam'))
    job.fileStore.readGlobalFile(normal_bai_id, os.path.join(work_dir, 'normal.bai'))
    job.fileStore.readGlobalFile(tumor_bam_id, os.path.join(work_dir, 'tumor.bam'))
    job.fileStore.readGlobalFile(tumor_bai_id, os.path.join(work_dir, 'tumor.bai'))
    return_input_paths(job, work_dir, ids, 'ref.fasta')
    # Call: Pindel
    parameters = ['-f', '/data/ref.fasta',
                  '-i', '/data/config.txt',
                  '--number_of_threads', cores,
                  '--minimum_support_for_event', '3',
                  '--report_long_insertions', 'true',
                  '--report_breakpoints', 'true',
                  '-o', uuid]
    docker_call(tool='quay.io/ucsc_cgl/pindel:0.2.5b6--4e8d1b31d4028f464b3409c6558fb9dfcad73f88',
                work_dir=work_dir, tool_parameters=parameters, sudo=sudo)
    output_files = [f for f in glob.glob(os.path.join(work_dir, '*')) if '.bam' not in f and '.bai' not in f]
    tarball_files(work_dir, 'pindel.tar.gz', files=[output_files])
    # Return fileStore ID
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'pindel.tar.gz'))


def muse(job, job_vars, bam_ids):
    """
    Calls MuSe to find variants

    job_vars: tuple     Contains the input_args and ids dictionaries
    bam_ids: tuple      Contains a tuple of normal/tumor fileStore ids for bams and index files (bai)
    """
    # Unpack convenience variables for job
    normal_ids, tumor_ids = bam_ids
    normal_bam_id, normal_bai_id = normal_ids
    tumor_bam_id, tumor_bai_id = tumor_ids
    input_args, ids = job_vars
    sudo = input_args['sudo']
    uuid = input_args['uuid']
    cores = int(input_args['cpu_count'])
    output_name = uuid + '.muse.vcf'
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve input files
    job.fileStore.readGlobalFile(normal_bam_id, os.path.join(work_dir, 'normal.bam'))
    job.fileStore.readGlobalFile(normal_bai_id, os.path.join(work_dir, 'normal.bai'))
    job.fileStore.readGlobalFile(tumor_bam_id, os.path.join(work_dir, 'tumor.bam'))
    job.fileStore.readGlobalFile(tumor_bai_id, os.path.join(work_dir, 'tumor.bai'))
    return_input_paths(job, work_dir, ids, 'ref.fasta', 'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict')
    # Call: MuSE
    parameters = ['--muse', 'MuSEv1.0rc',
                  '--mode', 'wxs',
                  '--dbsnp', 'dbsnp.vcf',
                  '--fafile', 'ref.fasta',
                  '--tumor-bam', 'tumor.bam',
                  '--tumor-bam-index', 'tumor.bam.bai',
                  '--normal-bam', 'normal.bam',
                  '--normal-bam-index', 'normal.bam.bai',
                  '--outfile', output_name,
                  '--cpus', cores]
    docker_call(tool='quay.io/ucsc_cgl/muse:1.0--6add9b0a1662d44fd13bbc1f32eac49326e48562',
                work_dir=work_dir, tool_parameters=parameters, sudo=sudo)
    tarball_files(work_dir, 'pindel.tar.gz', files=[output_name])
    # Return fileStore ID
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'pindel.tar.gz'))


def consolidate_output(job, job_vars, mutect_id):#, pindel_id, muse_id):
    """
    Combine the contents of separate zipped outputs into one via streaming

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    output_ids: tuple   Nested tuple of all the output fileStore IDs
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    # Retrieve output file paths to consolidate
    mutect_tar = job.fileStore.readGlobalFile(mutect_id, os.path.join(work_dir, 'mutect.tar.gz'))
    # pindel_tar = job.fileStore.readGlobalFile(pindel_id, os.path.join(work_dir, 'pindel.tar.gz'))
    # muse_tar = job.fileStore.readGlobalFile(muse_id, os.path.join(work_dir, 'muse.tar.gz'))
    # I/O
    out_tar = os.path.join(work_dir, uuid + '.tar.gz')
    # Consolidate separate tarballs into one as streams (avoids unnecessary untaring)
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        # for tar in [mutect_tar, pindel_tar, muse_tar]:
        for tar in [mutect_tar]:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == mutect_tar:
                            tarinfo.name = os.path.join(uuid, 'mutect', os.path.basename(tarinfo.name))
                        # elif tar == pindel_tar:
                        #     tarinfo.name = os.path.join(uuid, 'Kallisto', os.path.basename(tarinfo.name))
                        # else:
                        #     tarinfo.name = os.path.join(uuid, 'MuSe', os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
    # Move to output directory of selected
    if input_args['output_dir']:
        output_dir = input_args['output_dir']
        mkdir_p(output_dir)
        copy_to_output_dir(work_dir, output_dir, uuid=None, files=[uuid + '.tar.gz'])
    # Write output file to fileStore
    output_tar = job.fileStore.writeGlobalFile(out_tar)
    # If S3 bucket argument specified, upload to S3
    if input_args['s3_dir']:
        job.addChildJobFn(upload_output_to_s3, input_args, output_tar)


def upload_output_to_s3(job, input_args, output_tar):
    """
    Uploads Mutect files to S3 using S3AM

    mutect_id: str      FileStore ID for the mutect.vcf
    input_args: dict    Dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.lstrip('/').split('/')[0]
    bucket_dir = '/'.join(s3_dir.lstrip('/').split('/')[1:])
    # Retrieve VCF file
    job.fileStore.readGlobalFile(output_tar, os.path.join(work_dir, uuid + '.tar.gz'))
    # Upload to S3 via S3AM
    s3am_command = ['s3am',
                    'upload',
                    'file://{}'.format(os.path.join(work_dir, uuid + '.tar.gz')),
                    os.path.join('s3://', bucket_name, bucket_dir, uuid + '.tar.gz')]
    subprocess.check_call(s3am_command)


def main():
    """
    This is a Toil pipeline used to perform variant analysis (usually on exomes) from Tumor/Normal BAMs.
    All samples are co-cleaned (GATK Indel Realignment (IR) and Base Quality Score Recalibration (BQSR))
    before variant analysis is performed by MuTect.  The final output of this pipeline is a tarball
    containing the output of MuTect (.vcf, .cov, .out).

    Please see the associated README.md for an overview and quickstart walkthrough.
    """
    # Define Parser object and add to jobTree
    argparser = build_parser()
    Job.Runner.addToilOptions(argparser)
    pargs = argparser.parse_args()
    # Variables to pass to initial job
    inputs = {'ref.fasta': pargs.reference,
              'config': pargs.config,
              'phase.vcf': pargs.phase,
              'mills.vcf': pargs.mills,
              'dbsnp.vcf': pargs.dbsnp,
              'cosmic.vcf': pargs.cosmic,
              'output_dir': pargs.output_dir,
              'ssec': pargs.ssec,
              's3_dir': pargs.s3_dir,
              'sudo': pargs.sudo,
              'uuid': None,
              'normal.bam': None,
              'tumor.bam': None,
              'cpu_count': None}

    # Launch Pipeline
    Job.Runner.startToil(Job.wrapJobFn(download_shared_files, inputs), pargs)

if __name__ == '__main__':
    main()
