#!/usr/bin/env python2.7
"""
@author Jacob Pfeil
@data 01/13/2016

Toil pipeline for processing bam files for GATK halpotype calling

1 = download shared data
2 = reference preprocessing
3 = download sample 
4 = index
5 = sort
6 = mark duplicates
7 = index
8 = realigner target 
9 = indel realignment
10 = index
11 = base recalibration
12 = ouput bqsr file
"""
import argparse
import collections
import multiprocessing
import os
import hashlib
import base64
import subprocess
import shutil
import sys
import logging
import errno
from toil.job import Job

from toil_scripts import download_from_s3_url
from toil_scripts.lib.urls import s3am_upload
from toil_scripts.lib.programs import docker_call

_log = logging.getLogger(__name__)

debug = False 

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
    parser.add_argument('-x', '--suffix', default=".bqsr", help='additional suffix, if any')
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
    """
    Moves files from the working directory to the output directory.

    :param work_dir: the working directory
    :param output_dir: the output directory
    :param filenames: remaining arguments are filenames
    """
    for filename in filenames:
        origin = os.path.join(work_dir, filename)
        dest = os.path.join(output_dir, filename)
        try:
            shutil.move(origin, dest)
        except IOError:
            mkdir_p(output_dir)
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

def upload_or_move(job, work_dir, input_args, output):

    # are we moving this into a local dir, or up to s3?
    if input_args['output_dir']:
        # get output path and
        output_dir = input_args['output_dir']
        # FIXME: undefined function
        make_directory(output_dir)
        move_to_output_dir(work_dir, output_dir, output)

    elif input_args['s3_dir']:
        s3am_upload(fpath=os.path.join(work_dir, output),
                    s3_dir=input_args['s3_dir'],
                    s3_key_path=input_args['ssec'])

    else:
        raise ValueError('No output_directory or s3_dir defined. Cannot determine where to store %s' % output)

def download_from_url_gatk(job, url, name):
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
                if debug:
                    debug_log = open('download_from_url', 'a')
                    debug_log.write(file_path + '\n')
                    debug_log.close()
                    f = open(file_path, 'w')
                    f.write('debug')
                    f.close()
                else:
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


def create_reference_index(job, ref_id, sudo):
    """
    Uses Samtools to create reference index file (.fasta.fai)

    ref_id: str     The fileStore ID of the reference
    sudo: bool      Boolean item to determine whether to invoke sudo with docker
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path to reference
    try:
        job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fa'))  
    except:
        sys.stderr.write("Failed when reading global file %s to %s. Retrying with dict index." % (ref_id,
                                                                                                  os.path.join(work_dir, 'ref.fa')))
        
        try:
            job.fileStore.readGlobalFile(ref_id['ref.fa'], os.path.join(work_dir, 'ref.fa'))  
        except:
            sys.stderr.write("Reading %s on retry failed." % ref_id['ref.fa'])
            raise

    # Call: Samtools
    command = ['faidx', 'ref.fa']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=['ref.fa'],
                outputs={'ref.fa.fai': None},
                sudo=sudo)
    output = os.path.join(work_dir, 'ref.fa.fai')
    assert os.path.exists(output)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(output)


def create_reference_dict(job, ref_id, input_args):
    """
    Uses Picardtools to create reference dictionary (.dict) for the sample

    ref_id: str     The fileStore ID of the reference
    input_args: dict        Dictionary of input arguments (from main())
    """
    sudo = input_args['sudo']
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fa'))
    # Call: picardtools
    command = ['CreateSequenceDictionary', 'R=ref.fa', 'O=ref.dict']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % input_args['memory']},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=['ref.fa'],
                outputs={'ref.dict': None},
                sudo=sudo)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.dict'))


# TODO Start of Pipeline
def download_gatk_files(job, input_args):
    """
    Downloads files shared by all samples in the pipeline

    input_args: dict        Dictionary of input arguments (from main())
    """
    
    shared_ids = {}
    for fname in ['ref.fa', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf']:
        shared_ids[fname] = job.addChildJobFn(download_from_url_gatk, url=input_args[fname], name=fname).rv()
    job.addFollowOnJobFn(reference_preprocessing, shared_ids, input_args)


def reference_preprocessing(job, shared_ids, input_args):
    """
    Create index and dict file for reference

    input_args: dict        Dictionary of input argumnets
    shared_ids: dict        Dictionary of fileStore IDs
    """
    ref_id = shared_ids['ref.fa']
    if isinstance(ref_id, dict):
        sys.stderr.write("shared_ids['ref.fa'] is a dict. %s['ref_id'] = %s." % (shared_ids, ref_id))
        ref_id = ref_id['ref.fa']

    sudo = input_args['sudo']
    shared_ids['ref.fa.fai'] = job.addChildJobFn(create_reference_index, ref_id, sudo).rv()
    shared_ids['ref.dict'] = job.addChildJobFn(create_reference_dict, ref_id, input_args).rv()
    job.addFollowOnJobFn(spawn_batch_preprocessing, shared_ids, input_args)


def spawn_batch_preprocessing(job, shared_ids, input_args):
    """
    Spawn a pipeline for each sample in the configuration file

    input_args: dict        Dictionary of input argumnets
    shared_ids: dict        Dictionary of fileStore IDs
    """
    samples = []
    config = input_args['config']

    # does the config file exist locally? if not, try to read from job store
    if not os.path.exists(config):

        work_dir = job.fileStore.getLocalTempDir()
        config_path = os.path.join(work_dir, 'config.txt')
        job.fileStore.readGlobalFile(config, config_path)
        config = config_path

    with open(config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                samples.append(line.strip().split(','))
    for sample in samples:
        job.addChildJobFn(download_sample, shared_ids, input_args, sample)


def download_sample(job, shared_ids, input_args, sample):
    """
    Defines sample variables then downloads the sample.

    ids: dict           Dictionary of fileStore IDs
    input_args: dict    Dictionary of input arguments
    sample: str         Contains uuid, normal url, and tumor url
    """
    uuid, url = sample
    # Create a unique
    input_args['uuid'] = uuid
    if input_args['output_dir']:
        input_args['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    # Download sample bams and launch pipeline
    if input_args['ssec']:
        shared_ids['sample.bam'] = job.addChildJobFn(download_encrypted_file, input_args, 'sample.bam').rv()
    else:
        shared_ids['sample.bam'] = job.addChildJobFn(download_from_url_gatk, url=url, name='sample.bam').rv()
    job.addFollowOnJobFn(remove_supplementary_alignments, shared_ids, input_args)

def remove_supplementary_alignments(job, shared_ids, input_args):
    work_dir = job.fileStore.getLocalTempDir()

    #Retrieve file path
    read_from_filestore(job, work_dir, shared_ids, 'sample.bam')
    outpath = os.path.join(work_dir, 'sample.nosuppl.bam')

    command = ['view',
               '-b', '-o', '/data/sample.nosuppl.bam',
               '-F', '0x800',
               '/data/sample.bam']
               
    sudo = input_args['sudo']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c',
                inputs=['sample.bam'],
                outputs={'sample.nosuppl.bam': None},
                sudo=sudo)
    shared_ids['sample.nosuppl.bam'] = job.fileStore.writeGlobalFile(outpath)
    job.addChildJobFn(sort_sample, shared_ids, input_args)


def _move_bai(outpath):

    # bam index gets written without bam extension
    os.rename(outpath.replace('bam', 'bai'), outpath + '.bai')


def sort_sample(job, shared_ids, input_args):
    """
    Uses picardtools SortSam to sort a sample bam file
    """
    work_dir = job.fileStore.getLocalTempDir()

    #Retrieve file path
    read_from_filestore(job, work_dir, shared_ids, 'sample.nosuppl.bam')
    outpath = os.path.join(work_dir, 'sample.sorted.bam')
    #Call: picardtools
    command = ['SortSam',
               'INPUT=sample.nosuppl.bam',
               'OUTPUT=sample.sorted.bam',
               'SORT_ORDER=coordinate',
               'CREATE_INDEX=true']
    sudo = input_args['sudo']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % input_args['memory']},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=['sample.nosuppl.bam'],
                outputs={'sample.sorted.bam': None, 'sample.sorted.bai': None},
                sudo=sudo)
    shared_ids['sample.sorted.bam'] = job.fileStore.writeGlobalFile(outpath)
    job.addChildJobFn(mark_dups_sample, shared_ids, input_args)


def mark_dups_sample(job, shared_ids, input_args):
    """
    Uses picardtools MarkDuplicates
    """
    sudo = input_args['sudo']
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    read_from_filestore(job, work_dir, shared_ids, 'sample.sorted.bam')
    outpath = os.path.join(work_dir, 'sample.mkdups.bam')
    # Call: picardtools
    command = ['MarkDuplicates',
               'INPUT=sample.sorted.bam',
               'OUTPUT=sample.mkdups.bam',
               'METRICS_FILE=metrics.txt',
               'ASSUME_SORTED=true',
               'CREATE_INDEX=true']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % input_args['memory']},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=['sample.sorted.bam'],
                outputs={'sample.mkdups.bam': None, 'sample.mkdups.bai': None},
                sudo=sudo)
    shared_ids['sample.mkdups.bam'] = job.fileStore.writeGlobalFile(outpath)

    # picard writes the index for file.bam at file.bai, not file.bam.bai
    _move_bai(outpath)
    shared_ids['sample.mkdups.bam.bai'] = job.fileStore.writeGlobalFile(outpath + ".bai")
    job.addChildJobFn(realigner_target_creator, shared_ids, input_args, cores = input_args['cpu_count'])


def realigner_target_creator(job, shared_ids, input_args):
    """
    Creates <type>.intervals file needed for indel realignment

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # Retrieve input file paths
    read_from_filestore(job, work_dir, shared_ids, 'ref.fa',
                        'sample.mkdups.bam', 'ref.fa.fai', 'ref.dict',
                        'sample.mkdups.bam.bai', 'phase.vcf', 'mills.vcf')

    # Output file path
    output = os.path.join(work_dir, 'sample.intervals')
    # Call: GATK -- RealignerTargetCreator
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'RealignerTargetCreator',
                  '-nt', str(input_args['cpu_count']),
                  '-R', 'ref.fa',
                  '-I', 'sample.mkdups.bam',
                  '-known', 'phase.vcf',
                  '-known', 'mills.vcf',
                  '--downsampling_type', 'NONE',
                  '-o', 'sample.intervals']

    docker_call(work_dir=work_dir, parameters=parameters,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=['ref.fa','sample.mkdups.bam', 'ref.fa.fai', 'ref.dict',
                        'sample.mkdups.bam.bai', 'phase.vcf', 'mills.vcf'],
                outputs={'sample.intervals': None},
                env={'JAVA_OPTS':'-Xmx%sg' % input_args['memory']},
                sudo=sudo)
    shared_ids['sample.intervals'] = job.fileStore.writeGlobalFile(output)
    job.addChildJobFn(indel_realignment, shared_ids, input_args)


def indel_realignment(job, shared_ids, input_args):
    """
    Creates realigned bams using <sample>.intervals file from previous step

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # Retrieve input file paths
    return_input_paths(job, work_dir, shared_ids, 'ref.fa',
                       'sample.mkdups.bam', 'phase.vcf', 'mills.vcf',
                       'sample.intervals', 'ref.fa.fai', 'ref.dict',
                       'sample.mkdups.bam.bai')
    # Output file path
    outpath = os.path.join(work_dir, 'sample.indel.bam')
    # Call: GATK -- IndelRealigner
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'IndelRealigner',
                  '-R', 'ref.fa',
                  '-I', 'sample.mkdups.bam',
                  '-known', 'phase.vcf',
                  '-known', 'mills.vcf',
                  '-targetIntervals', 'sample.intervals',
                  '--downsampling_type', 'NONE',
                  '-maxReads', str(720000),
                  '-maxInMemory', str(5400000),
                  '-o', 'sample.indel.bam']

    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                work_dir=work_dir, parameters=parameters,
                inputs=['ref.fa', 'sample.mkdups.bam', 'phase.vcf', 'mills.vcf',
                        'sample.intervals', 'ref.fa.fai', 'ref.dict',
                        'sample.mkdups.bam.bai'],
                outputs={'sample.indel.bam': None, 'sample.indel.bai': None},
                env={'JAVA_OPTS':'-Xmx10g'}, sudo=sudo)

    # Write to fileStore
    shared_ids['sample.indel.bam'] = job.fileStore.writeGlobalFile(outpath)
    _move_bai(outpath)
    shared_ids['sample.indel.bam.bai'] = job.fileStore.writeGlobalFile(outpath + ".bai")
    job.addChildJobFn(base_recalibration, shared_ids, input_args, cores = input_args['cpu_count'])


def base_recalibration(job, shared_ids, input_args):
    """
    Creates recal table to perform Base Quality Score Recalibration

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # Retrieve input file paths
    return_input_paths(job, work_dir, shared_ids, 'ref.fa', 'sample.indel.bam',
                       'dbsnp.vcf', 'ref.fa.fai',
                       'ref.dict', 'sample.indel.bam.bai')
    # Output file path
    output = os.path.join(work_dir, 'sample.recal.table')
    # Call: GATK -- IndelRealigner
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'BaseRecalibrator',
                  '-nct', str(input_args['cpu_count']),
                  '-R', 'ref.fa',
                  '-I', 'sample.indel.bam',
                  '-knownSites', 'dbsnp.vcf',
                  '-o', 'sample.recal.table']
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                work_dir=work_dir, parameters=parameters,
                inputs=['ref.fa', 'sample.indel.bam', 'dbsnp.vcf', 'ref.fa.fai',
                        'ref.dict', 'sample.indel.bam.bai'],
                outputs={'sample.recal.table': None},
                env={'JAVA_OPTS':'-Xmx%sg' % input_args['memory']}, sudo=sudo)
    # Write to fileStore
    shared_ids['sample.recal.table'] = job.fileStore.writeGlobalFile(output)
    job.addChildJobFn(print_reads, shared_ids, input_args, cores = input_args['cpu_count'])


def print_reads(job, shared_ids, input_args):
    """
    Create bam that has undergone Base Quality Score Recalibration (BQSR)

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    uuid = input_args['uuid']
    suffix = input_args['suffix']
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # Retrieve input file paths
    return_input_paths(job, work_dir, shared_ids, 'ref.fa', 'sample.indel.bam',
                       'ref.fa.fai', 'ref.dict', 'sample.indel.bam.bai', 'sample.recal.table')
    # Output file
    outfile = '{}{}.bam'.format(uuid, suffix)
    gatk_outfile_idx = '{}{}.bai'.format(uuid, suffix)
    outfile_idx = '{}{}.bam.bai'.format(uuid, suffix)
    outpath = os.path.join(work_dir, outfile)
    # Call: GATK -- PrintReads
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'PrintReads',
                  '-nct', str(input_args['cpu_count']),
                  '-R', 'ref.fa',
                  '--emit_original_quals',
                  '-I', 'sample.indel.bam',
                  '-BQSR', 'sample.recal.table',
                  '-o', outfile]
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                work_dir=work_dir, parameters=parameters,
                inputs=['ref.fa', 'sample.indel.bam', 'ref.fa.fai', 'ref.dict', 
                        'sample.indel.bam.bai', 'sample.recal.table'],
                outputs={outfile: None, gatk_outfile_idx: None},
                env={'JAVA_OPTS':'-Xmx%sg' % input_args['memory']}, sudo=sudo)
    
    upload_or_move(job, work_dir, input_args, outfile)
    _move_bai(outpath)
    upload_or_move(job, work_dir, input_args, outfile_idx)

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
              'ssec': pargs.ssec,
              'suffix': pargs.suffix,
              'cpu_count': multiprocessing.cpu_count(), # FIXME: should not be called from toil-leader, see #186
              'memory': '15'}

    # Launch Pipeline
    Job.Runner.startToil(Job.wrapJobFn(download_gatk_files, inputs), pargs)

if __name__ == '__main__':
    main()
