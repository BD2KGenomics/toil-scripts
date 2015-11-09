#!/usr/bin/env python2.7
# John Vivian
# 7-17-15

"""
Tree Structure of GATK Pipeline (per sample)
   0------------> 5 ----------> 14
/ / \ \          / \
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
5 = pre-processing node
6, 7 = RealignerTargetCreator
8,9 = IndelRealigner
10,11 = BaseRecalibration
12, 13 = PrintReads
14 = MuTect

0-4, 6-13 are "children jobs"
5 and 14 are "follow-on jobs"
===================================================================
:Dependencies:
curl            - apt-get install curl
docker          - apt-get install docker (or 'docker.io' for linux)
Toil            - pip install --pre toil
"""
import argparse
import base64
from collections import OrderedDict
import hashlib
import os
import shutil
import subprocess
import multiprocessing
import errno
from toil.job import Job


def build_parser():
    """
    Contains arguments for the all of necessary input files
    """
    parser = argparse.ArgumentParser()
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
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


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


def download_encrypted_file(job, input_args, name):
    """
    Downloads encrypted files from S3 via header injection
    Input1: input arguments defined in main()
    Input2: dictionary of jobStore IDs
    Input3: symbolic name associated with file
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
    encoded_key_md5 = base64.b64encode( hashlib.md5(key).digest() )
    h1 = 'x-amz-server-side-encryption-customer-algorithm:AES256'
    h2 = 'x-amz-server-side-encryption-customer-key:{}'.format(encoded_key)
    h3 = 'x-amz-server-side-encryption-customer-key-md5:{}'.format(encoded_key_md5)
    try:
        subprocess.check_call(['curl', '-fs', '--retry', '5', '-H', h1, '-H', h2, '-H', h3, url, '-o', file_path])
    except OSError:
        raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    # job.fileStore.updateGlobalFile(ids[name], file_path)
    return job.fileStore.writeGlobalFile(file_path)


def copy_to_output_dir(work_dir, output_dir, uuid=None, files=None):
    """
    A list of files to move from work_dir to output_dir.
    """
    for fname in files:
        if uuid is None:
            shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


def download_from_URL(job, url, name):
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
    # job.fileStore.updateGlobalFile(ids[name], file_path)
    return job.fileStore.writeGlobalFile(file_path)


def return_input_paths(job, work_dir, ids, *args):
    """
    Given one or more strings representing file_names, return the paths to those files.
    """
    # for every file in *args, place in work_dir via the FileStore and then return the mounted docker path.
    paths = OrderedDict()
    for name in args:
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = docker_path(job.fileStore.readGlobalFile(ids[name], os.path.join(work_dir, name)))
        else:
            file_path = docker_path(name)
        paths[name] = file_path
        if len(args) == 1:
            paths = file_path

    return paths


def docker_path(file_path):
    """
    Returns file names to the perceived path inside of the tools' Docker containers
    """
    return os.path.join('/data', os.path.basename(file_path))


def docker_call(work_dir, tool_parameters, tool, java_opts=None):
    """
    Makes subprocess call of a command to a docker container.
    :parameter tool: name of the Docker image to be used (e.g. computationalgenomicslab/samtools)
    """
    if java_opts:
        base_docker_call = 'docker run -e JAVA_OPTS={} -v {}:/data'.format(java_opts, work_dir)
    else:
        base_docker_call = 'docker run -v {}:/data'.format(work_dir)
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_parameters)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


# Start of Job Functions
def download_shared_files(job, input_args):
    """
    Downloads files shared by all samples in the pipeline
    """
    shared_ids = {}
    for fname in ['ref.fasta', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf', 'cosmic.vcf']:
        shared_ids[fname] = job.addChildJobFn(download_from_URL, url=input_args[fname], name=fname).rv()
    job.addFollowOnJobFn(reference_preprocessing, input_args, shared_ids)


def reference_preprocessing(job, input_args, shared_ids):
    """
    Create index and dict file for reference
    """
    ref_id = shared_ids['ref.fasta']
    shared_ids['ref.fasta.fai'] = job.addChildJobFn(create_reference_index, ref_id).rv()
    shared_ids['ref.dict'] = job.addChildJobFn(create_reference_dict, ref_id).rv()
    job.addFollowOnJobFn(spawn_batch_jobs, input_args, shared_ids)


def create_reference_index(job, ref_id):
    """
    Uses Samtools to create reference index file (.fasta.fai)
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path to reference
    ref_path = docker_path(job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta')))
    # Call: Samtools
    command = ['faidx', ref_path]
    docker_call(work_dir, command, tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.fasta.fai'))


def create_reference_dict(job, ref_id):
    """
    Uses Picardtools to create reference dictionary (.dict) for the sample
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta'))
    # Call: picardtools
    output = os.path.splitext(docker_path(ref_path))[0]
    command = ['CreateSequenceDictionary', 'R={}'.format(docker_path(ref_path)), 'O={}.dict'.format(output)]
    docker_call(work_dir, command, tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.dict'))


def spawn_batch_jobs(job, input_args, shared_ids):
    """
    Spawn a pipeline for each sample in the configuration file
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
        ids['normal.bam'] = job.addChildJobFn(download_from_URL, url=sample_input['normal.bam'], name='normal.bam').rv()
        ids['tumor.bam'] = job.addChildJobFn(download_from_URL, url=sample_input['tumor.bam'], name='tumor.bam').rv()
    job_vars = (sample_input, ids)
    job.addFollowOnJobFn(pipeline_launchpoint, job_vars)


def pipeline_launchpoint(job, job_vars):
    """
    Statically link the rest of the workflow
    """
    A = job.wrapJobFn(index_bams, job_vars).encapsulate()
    B = job.wrapJobFn(mutect, job_vars, A.rv())
    # Since A encapsulates all of the pre-processing, mutect can be added as a child of A
    job.addChild(A)
    A.addChild(B)


def index_bams(job, job_vars):
    """
    Create index (.bai) files for each sample bam
    """
    normal_ids = job.addChildJobFn(index, job_vars, 'normal', cores=1, memory='1 G', disk='8 G').rv()
    tumor_ids = job.addChildJobFn(index, job_vars, 'tumor', cores=1, memory='1 G', disk='8 G').rv()
    return (normal_ids, tumor_ids)


def index(job, job_vars, sample):
    """
    Runs samtools index to create (.bai) files
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = int(input_args['cpu_count'])
    # Retrieve file path
    bam = '{}.bam'.format(sample)
    path = return_input_paths(job, work_dir, ids, bam)
    # Call: index the normal.bam
    parameters = ['index', '{}'.format(docker_path(path))]
    docker_call(work_dir, tool_parameters=parameters, tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    # Write to fileStore
    ids[bam + '.bai'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, bam) + '.bai')
    return job.addChildJobFn(realigner_target_creator, job_vars, sample, cores=cores, memory='10 G', disk='15 G').rv()


def realigner_target_creator(job, job_vars, sample):
    """
    Creates <type>.intervals file needed for indel realignment
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    ref_fasta, bam, ref_fai, ref_dict, bam_bai,\
    phase_vcf, mills_vcf = return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.bam'.format(sample),
                                              'ref.fasta.fai', 'ref.dict', '{}.bam.bai'.format(sample), 'phase.vcf', 'mills.vcf')

    # Output file path
    output = os.path.join(work_dir, '{}.intervals'.format(sample))
    # Call: GATK -- RealignerTargetCreator
    parameters = ['-T', 'RealignerTargetCreator',
                  '-nt', str(cores),
                  '-R', docker_path(ref_fasta),
                  '-I', docker_path(bam),
                  '-known', docker_path(phase_vcf),
                  '-known', docker_path(mills_vcf),
                  '--downsampling_type', 'NONE',
                  '-o', docker_path(output)]

    docker_call(work_dir, tool_parameters=parameters, tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e', java_opts='-Xmx10g')
    # Write to fileStore
    ids['{}.intervals'.format(sample)] = job.fileStore.writeGlobalFile(output)
    return job.addChildJobFn(indel_realignment, job_vars, sample, cores=1, memory='10 G', disk='30 G').rv()


def indel_realignment(job, job_vars, sample):
    """
    Creates realigned bams using <sample>.intervals file from previous step
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    ref_fasta, bam, phase_vcf, \
    mills_vcf, intervals, ref_fai, \
    ref_dict, bam_bai = return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.bam'.format(sample), 'phase.vcf',
                                           'mills.vcf', '{}.intervals'.format(sample), 'ref.fasta.fai', 'ref.dict',
                                           '{}.bam.bai'.format(sample))
    # Output file path
    output = os.path.join(work_dir, '{}.indel.bam'.format(sample))
    # Call: GATK -- IndelRealigner
    command = ['-T', 'IndelRealigner',
              '-R', ref_fasta,
              '-I', bam,
              '-known', phase_vcf,
              '-known', mills_vcf,
              '-targetIntervals', intervals,
              '--downsampling_type', 'NONE',
              '-maxReads', str(720000),
              '-maxInMemory', str(5400000),
              '-o', docker_path(output)]
    docker_call(work_dir, command, tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e', java_opts='-Xmx10g')
    # Write to fileStore
    ids['{}.indel.bam'.format(sample)] = job.fileStore.writeGlobalFile(output)
    ids['{}.indel.bai'.format(sample)] = job.fileStore.writeGlobalFile(os.path.splitext(output)[0] + '.bai')
    return job.addChildJobFn(base_recalibration, job_vars, sample, cores=cores, memory='15 G', disk='15 G').rv()


def base_recalibration(job, job_vars, sample):
    """
    Creates recal table to perform Base Quality Score Recalibration
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    ref_fasta, indel_bam, dbsnp_vcf, ref_fai, \
    ref_dict, bam_bai = return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.indel.bam'.format(sample),
                                           'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict', '{}.indel.bai'.format(sample))
    # Output file path
    output = os.path.join(work_dir, '{}.recal.table'.format(sample))
    # Call: GATK -- IndelRealigner
    command = ['-T', 'BaseRecalibrator',
              '-nct', str(cores),
              '-R', ref_fasta,
              '-I', indel_bam,
              '-knownSites', dbsnp_vcf,
              '-o', docker_path(output)]
    docker_call(work_dir, command, tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e', java_opts='-Xmx15g')
    # Write to fileStore
    ids['{}.recal.table'.format(sample)] = job.fileStore.writeGlobalFile(output)
    return job.addChildJobFn(print_reads, job_vars, sample, cores=cores, memory='15 G', disk='40 G').rv()


def print_reads(job, job_vars, sample):
    """
    Create bam that has undergone Base Quality Score Recalibration (BQSR)
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = int(input_args['cpu_count'])
    # Retrieve input file paths
    ref_fasta, indel_bam, dbsnp_vcf, ref_fai, \
    ref_dict, bam_bai, recal_table = return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.indel.bam'.format(sample),
                                                        'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict', '{}.indel.bai'.format(sample),
                                                        '{}.recal.table'.format(sample))
    # Output file
    output = os.path.join(work_dir, '{}.bqsr.bam'.format(sample))
    # Call: GATK -- PrintReads
    command = ['-T', 'PrintReads',
              '-nct', str(cores),
              '-R', ref_fasta,
              '--emit_original_quals',
              '-I', indel_bam,
              '-BQSR', recal_table,
              '-o', docker_path(output)]
    docker_call(work_dir, command, tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e', java_opts='-Xmx15g')
    # Write to fileStore
    bam_id = job.fileStore.writeGlobalFile(output)
    bai_id = job.fileStore.writeGlobalFile(os.path.splitext(output)[0] + '.bai')
    return (bam_id, bai_id)


def mutect(job, job_vars, bam_ids):
    """
    Calls MuTect to perform variant analysis
    """
    # Unpack convenience variables for job
    normal_ids, tumor_ids = bam_ids
    normal_bam_id, normal_bai_id = normal_ids
    tumor_bam_id, tumor_bai_id = tumor_ids
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve input files
    normal_bam = job.fileStore.readGlobalFile(normal_bam_id, os.path.join(work_dir, 'normal.bam'))
    job.fileStore.readGlobalFile(normal_bai_id, os.path.join(work_dir, 'normal.bai'))
    tumor_bam = job.fileStore.readGlobalFile(tumor_bam_id, os.path.join(work_dir, 'tumor.bam'))
    job.fileStore.readGlobalFile(tumor_bai_id, os.path.join(work_dir, 'tumor.bai'))
    return_input_paths(job, work_dir, ids, 'ref.fasta', 'dbsnp.vcf', 'ref.fasta.fai',
                                              'ref.dict', 'cosmic.vcf')
    # Output VCF
    uuid = input_args['uuid']
    output_name = uuid + '.vcf'
    mut_out = docker_path(os.path.join(work_dir, 'mutect.out'))
    mut_cov = docker_path(os.path.join(work_dir, 'mutect.cov'))
    # Call: MuTect
    command = ['--analysis_type', 'MuTect',
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
    docker_call(work_dir, command, tool='quay.io/ucsc_cgl/mutect:1.1.7--e8bf09459cf0aecb9f55ee689c2b2d194754cbd3')
    # Write to fileStore
    mutect_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, output_name))
    if input_args['output_dir']:
        output_dir = input_args['output_dir']
        mkdir_p(output_dir)
        copy_to_output_dir(work_dir, output_dir, uuid=None, files=[output_name])
    if input_args['s3_dir']:
        job.addChildJobFn(upload_vcf_to_s3, mutect_id, input_args)


def upload_vcf_to_s3(job, mutect_id, input_args):
    """
    Uploads VCF file to S3 using S3AM
    """
    uuid = input_args['uuid']
    work_dir = job.fileStore.getLocalTempDir()
    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.lstrip('/').split('/')[0]
    bucket_dir = '/'.join(s3_dir.lstrip('/').split('/')[1:])
    # Retrieve VCF file
    job.fileStore.readGlobalFile(mutect_id, os.path.join(work_dir, uuid + '.vcf'))
    # Upload to S3 via S3AM
    s3am_command = ['s3am',
                    'upload',
                    'file://{}'.format(os.path.join(work_dir, uuid + '.vcf')),
                    bucket_name,
                    os.path.join(bucket_dir, uuid + '.bam')]
    subprocess.check_call(s3am_command)


if __name__ == '__main__':
    # Define Parser object and add to jobTree
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # Variables to pass to initial job
    input_args = {'ref.fasta': args.reference,
                  'config': args.config,
                  'phase.vcf': args.phase,
                  'mills.vcf': args.mills,
                  'dbsnp.vcf': args.dbsnp,
                  'cosmic.vcf': args.cosmic,
                  'output_dir': args.output_dir,
                  'ssec': args.ssec,
                  's3_dir': args.s3_dir,
                  'uuid': None,
                  'normal.bam': None,
                  'tumor.bam': None,
                  'cpu_count': None }

    # Launch Pipeline
    Job.Runner.startToil(Job.wrapJobFn(download_shared_files, input_args), args)
