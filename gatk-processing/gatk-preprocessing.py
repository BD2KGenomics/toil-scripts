#!/usr/bin/env python2.7
"""
@author Jacob Pfeil
@data 01/13/2016



        0



0 = Start node
1 = reference index
2 = reference dict
"""

import argparse
import collections
import multiprocessing
import tarfile
import os
import errno
import hashlib
import base64
import subprocess
import shutil
from toil.job import Job


def build_parser():
    """
    Contains argparse arguments
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r', '--reference', required=True, help="Reference Genome URL")
    parser.add_argument('-f', '--config', required=True, help="Each line contains (CSV): UUID, BAM_URL")
    parser.add_argument('-p', '--phase', required=True, help='1000G_phase1.indels.hg19.sites.fixed.vcf URL')
    parser.add_argument('-m', '--mills', required=True, help='Mills_and_1000G_gold_standard.indels.hg19.sites.vcf URL')
    parser.add_argument('-d', '--dbsnp', required=True, help='dbsnp_132_b37.leftAligned.vcf URL')
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


def move_to_output_dir(work_dir, output_dir, *filenames):
    """`
    Moves files from the working directory to the output directory.

    :param work_dir: the working directory
    :param output_dir: the output directory
    :param filenames: remaining arguments are filenames
    """
    for filename in filenames:
        origin = os.path.join(work_dir, filename)
        dest = os.path.join(output_dir, filename)
        shutil.move(origin, dest)


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
    paths = collections.OrderedDict()
    for name in args:
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = job.fileStore.readGlobalFile(ids[name], os.path.join(work_dir, name))
        else:
            file_path = name
        paths[name] = file_path
        if len(args) == 1:
            paths = file_path

    return paths


def read_from_filestore(job, work_dir, ids, *filenames):
    """
    Reads file from fileStore and writes it to working directory.

    :param job: Job instance
    :param work_dir: working directory
    :param ids: shared file promises, dict
    :param filenames: remaining arguments are filenames
    """
    for filename in filenames:
        if not os.path.exists(os.path.join(work_dir, filename)):
            job.fileStore.readGlobalFile(ids[filename], os.path.join(work_dir, filename))


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


def create_reference_index(job, ref_id, sudo):
    """
    Uses Samtools to create reference index file (.fasta.fai)

    ref_id: str     The fileStore ID of the reference
    sudo: bool      Boolean item to determine whether to invoke sudo with docker
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path to reference
    ref_path = docker_path(job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fa')))
    # Call: Samtools
    command = ['faidx', ref_path]
    docker_call(work_dir=work_dir, tool_parameters=command,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e', sudo=sudo)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.fa.fai'))


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


def index(job, shared_ids, input_args, sample):
    """
    Runs samtools index to create (.bai) files

    job_vars: tuple     Contains the input_args and ids dictionaries
    """
    assert sample.endswith('.bam'), "ERROR: index takes a bam file"
    # Unpack convenience variables for job
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # Retrieve file path
    path = return_input_paths(job, work_dir, shared_ids, sample)
    # Call: index the normal.bam
    parameters = ['index', '{}'.format(docker_path(path))]
    docker_call(work_dir=work_dir, tool_parameters=parameters,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e', sudo=sudo)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, sample) + '.bai')


def sort(job, shared_ids, input_args, sample):
    """
    Uses picardtools SortSam to sort a sample bam file
    """
    assert sample.endswith('.bam'), "ERROR: sort takes a bam file"
    filename, ext = os.path.splitext(os.path.basename(sample))
    output = '{}.sorted.bam'.format(filename)
    work_dir = job.fileStore.getLocalTempDir()

    #Retrieve file path
    read_from_filestore(job, work_dir, sample)
    #Call: picardtools
    command = ['SortSam',
               'INPUT={}'.format(sample),
               'OUTPUT={}'.format(output),
               'SORT_ORDER=coordinate']
    sudo = input_args['sudo']
    docker_call(work_dir=work_dir, tool_parameters=command,
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e', sudo=sudo)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, output))


def mark_dups(job, shared_ids, input_args, sample):
    """
    Uses picardtools MarkDuplicates
    """
    assert sample.endswith('.bam'), "ERROR: mark_dups takes a bam file"
    filename, ext = os.path.splitext(os.path.basename(sample))
    output = '{}.mkdups.bam'.format(filename)
    sudo = input_args['sudo']
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    read_from_filestore(job, work_dir, shared_ids, sample)
    # Call: picardtools
    command = ['MarkDuplicates',
               'INPUT={}'.format(sample),
               'OUTPUT={}'.format(output),
               'SORT_ORDER=coordinate',
               'METRICS_FILE=metrics.txt',
               'ASSUME_SORTED=true']
    docker_call(work_dir=work_dir, tool_parameters=command,
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e', sudo=sudo)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, output))


def check_data(job):
    for f in os.listdir('data/'):
        job.fileStore.writeGlobalFile(os.path.join('data', f))

# TODO Start of Pipeline
def download_shared_files(job, input_args):
    """
    Downloads files shared by all samples in the pipeline

    input_args: dict        Dictionary of input arguments (from main())
    """
    job.addChildJobFn(check_data).rv()
    shared_ids = {}
    for fname in ['ref.fa', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf']:
        shared_ids[fname] = job.addChildJobFn(download_from_url, url=input_args[fname], name=fname).rv()
#    job.addFollowOnJobFn(reference_preprocessing, shared_ids, input_args)

#
def reference_preprocessing(job, shared_ids, input_args):
    """
    Create index and dict file for reference

    input_args: dict        Dictionary of input argumnets
    shared_ids: dict        Dictionary of fileStore IDs
    """
    ref_id = shared_ids['ref.fa']
    sudo = input_args['sudo']
    shared_ids['ref.fa.fai'] = job.addChildJobFn(create_reference_index, ref_id, sudo).rv()
    shared_ids['ref.dict'] = job.addChildJobFn(create_reference_dict, ref_id, sudo).rv()
    # job.addFollowOnJobFn(spawn_batch_jobs, shared_ids, input_args)


# def spawn_batch_jobs(job, shared_ids, input_args):
#     """
#     Spawn a pipeline for each sample in the configuration file
#
#     input_args: dict        Dictionary of input argumnets
#     shared_ids: dict        Dictionary of fileStore IDs
#     """
#     samples = []
#     config = input_args['config']
#     with open(config, 'r') as f:
#         for line in f:
#             if not line.isspace():
#                 sample = line.strip().split(',')
#                 job.addChildJobFn(download_samples, shared_ids, input_args, sample)
#
#
# def download_samples(job, shared_ids, input_args, sample):
#     """
#     Defines sample variables then downloads the sample.
#
#     ids: dict           Dictionary of fileStore IDs
#     input_args: dict    Dictionary of input arguments
#     sample: str         Contains uuid, normal url, and tumor url
#     """
#     uuid, url = sample
#     Create a unique
    # input_args['uuid'] = uuid
    # input_args['sample.bam'] = url
    # cores = multiprocessing.cpu_count()
    # if input_args['output_dir']:
    #     input_args['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    # Download sample bams and launch pipeline
    # if input_args['ssec']:
    #     shared_ids['sample.bam'] = job.addChildJobFn(download_encrypted_file, input_args, 'sample.bam').rv()
    # else:
    #     shared_ids['sample.bam'] = job.addChildJobFn(download_from_url, url=input_args['sample.bam'], name='sample.bam').rv()
    # job_vars = (input_args, shared_ids)
    # job.addChildJobFn(index_sample, job_vars, cores=cores, memory='10 G', disk='15 G')

def sample_preprocessing(job, shared_ids, input_args):
    """
    """
    shared_ids['sample.bam.bai'] = job.addChildJobFn(index, shared_ids,
                                                     input_args,
                                                     'sample.bam').rv()

    shared_ids['sample.sorted.bam'] = job.addChildJobFn(sort, shared_ids,
                                                        input_args,
                                                        'sample.bam').rv()

    shared_ids['sample.mkdups.bam'] = job.addChildJobFn(mark_dups,
                                                        shared_ids,
                                                        input_args,
                                                        'sample.sorted.bam').rv()





# def realigner_target_creator(job, job_vars):
#     """
#     Creates <type>.intervals file needed for indel realignment
#
#     job_vars: tuple     Contains the input_args and ids dictionaries
#     sample: str         Either "normal" or "tumor" to track which one is which
#     """
#     Unpack convenience variables for job
    # input_args, ids = job_vars
    # work_dir = job.fileStore.getLocalTempDir()
    # sudo = input_args['sudo']
    # cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    # return_input_paths(job, work_dir, ids, 'ref.fa', 'sample.mkdups.bam', 'ref.fa.fai', 'ref.dict',
    #                    'sample.mkdups.bam.bai', 'phase.vcf', 'mills.vcf')
    #
    # Output file path
    # output = os.path.join(work_dir, 'sample.intervals')
    # Call: GATK -- RealignerTargetCreator
    # parameters = ['-T', 'RealignerTargetCreator',
    #               '-nt', str(cores),
    #               '-R', '/data/ref.fasta',
    #               '-I', '/data/sample.mkdups.bam',
    #               '-known', '/data/phase.vcf',
    #               '-known', '/data/mills.vcf',
    #               '--downsampling_type', 'NONE',
    #               '-o', docker_path(output)]
    #
    # docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
    #             work_dir=work_dir, tool_parameters=parameters, java_opts='-Xmx10g', sudo=sudo)
    # Write to fileStore
    # ids['sample.intervals'] = job.fileStore.writeGlobalFile(output)
    # return job.addChildJobFn(indel_realignment, job_vars, cores=1, memory='10 G', disk='30 G').rv()
#
#
# def indel_realignment(job, job_vars):
#     """
#     Creates realigned bams using <sample>.intervals file from previous step
#
#     job_vars: tuple     Contains the input_args and ids dictionaries
#     sample: str         Either "normal" or "tumor" to track which one is which
#     """
#     Unpack convenience variables for job
    # input_args, ids = job_vars
    # work_dir = job.fileStore.getLocalTempDir()
    # sudo = input_args['sudo']
    # cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    # return_input_paths(job, work_dir, ids, 'ref.fa', 'sample.mkdups.bam', 'phase.vcf', 'mills.vcf',
    #                    'sample.intervals', 'ref.fa.fai', 'ref.dict', 'sample.mkdups.bam.bai')
  #  Output file path
    # output = os.path.join(work_dir, 'sample.indel.bam')
 #   Call: GATK -- IndelRealigner
    # parameters = ['-T', 'IndelRealigner',
    #               '-R', '/data/ref.fasta',
    #               '-I', '/data/sample.mkdups.bam',
    #               '-known', '/data/phase.vcf',
    #               '-known', '/data/mills.vcf',
    #               '-targetIntervals', '/data/sample.intervals',
    #               '--downsampling_type', 'NONE',
    #               '-maxReads', str(720000),
    #               '-maxInMemory', str(5400000),
    #               '-o', docker_path(output)]
    # docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
    #             work_dir=work_dir, tool_parameters=parameters, java_opts='-Xmx10g', sudo=sudo)
#    Write to fileStore
    # ids['sample.indel.bam'] = job.fileStore.writeGlobalFile(output)
    # ids['sample.indel.bai'] = job.fileStore.writeGlobalFile(os.path.splitext(output)[0] + '.bai')
    # return job.addChildJobFn(base_recalibration, job_vars, cores=cores, memory='15 G', disk='15 G').rv()

#
# def base_recalibration(job, job_vars):
#     """
#     Creates recal table to perform Base Quality Score Recalibration
#
#     job_vars: tuple     Contains the input_args and ids dictionaries
#     sample: str         Either "normal" or "tumor" to track which one is which
#     """
#     Unpack convenience variables for job
    # input_args, ids = job_vars
    # work_dir = job.fileStore.getLocalTempDir()
    # sudo = input_args['sudo']
    # cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    # return_input_paths(job, work_dir, ids, 'ref.fa', 'sample.indel.bam', 'dbsnp.vcf', 'ref.fa.fai',
    #                    'ref.dict', 'sample.indel.bai')
    # Output file path
    # output = os.path.join(work_dir, 'sample.recal.table')
    # Call: GATK -- IndelRealigner
    # parameters = ['-T', 'BaseRecalibrator',
    #               '-nct', str(cores),
    #               '-R', '/data/ref.fasta',
    #               '-I', '/data/sample.indel.bam',
    #               '-knownSites', '/data/dbsnp.vcf',
    #               '-o', docker_path(output)]
    # docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
    #             work_dir=work_dir, tool_parameters=parameters, java_opts='-Xmx15g', sudo=sudo)
    # Write to fileStore
    # ids['sample.recal.table'] = job.fileStore.writeGlobalFile(output)
    # return job.addChildJobFn(print_reads, job_vars, cores=cores, memory='15 G', disk='40 G').rv()
#
#
# def print_reads(job, job_vars, sample):
#     """
#     Create bam that has undergone Base Quality Score Recalibration (BQSR)
#
#     job_vars: tuple     Contains the input_args and ids dictionaries
#     sample: str         Either "normal" or "tumor" to track which one is which
#     """
#     Unpack convenience variables for job
    # input_args, ids = job_vars
    # uuid = input_args['uuid']
    # work_dir = job.fileStore.getLocalTempDir()
    # sudo = input_args['sudo']
    # cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    # return_input_paths(job, work_dir, ids, 'ref.fa', 'sample.indel.bam', 'ref.fa.fai',
    #                    'ref.dict', 'sample.indel.bai', 'sample.recal.table')
    # Output file
    # output = os.path.join(work_dir, '{}.bqsr.bam'.format(uuid))
    # Call: GATK -- PrintReads
    # parameters = ['-T', 'PrintReads',
    #               '-nct', str(cores),
    #               '-R', '/data/ref.fa',
    #               '--emit_original_quals',
    #               '-I', '/data/sample.indel.bam',
    #               '-BQSR', '/data/sample.recal.table',
    #               '-o', docker_path(output)]
    # docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
    #             work_dir=work_dir, tool_parameters=parameters, java_opts='-Xmx15g', sudo=sudo)
    # Write to fileStore

    # bam_id = job.fileStore.writeGlobalFile(output)
    # bai_id = job.fileStore.writeGlobalFile(os.path.splitext(output)[0] + '.bai')


def main():
    """
    GATK Pre-processing Script
    """
    # Define Parser object and add to jobTree
    argparser = build_parser()
    Job.Runner.addToilOptions(argparser)
    pargs = argparser.parse_args()
    # Variables to pass to initial job
    inputs = {'ref.fa': pargs.reference,
              'config': pargs.config,
              'phase.vcf': pargs.phase,
              'mills.vcf': pargs.mills,
              'dbsnp.vcf': pargs.dbsnp,
              'output_dir': pargs.output_dir,
              's3_dir': pargs.s3_dir,
              'sudo': pargs.sudo,
              'cpu_count': multiprocessing.cpu_count()}


    # Launch Pipeline
    Job.Runner.startToil(Job.wrapJobFn(download_shared_files, inputs), pargs)

if __name__ == '__main__':
    main()
