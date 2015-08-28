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
jobTree         - https://github.com/benedictpaten/jobTree
"""
import argparse
from collections import OrderedDict
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
    parser.add_argument('-o', '--output_dir', required=True, help='Full path to final output dir')
    parser.add_argument('-s', '--ssec', help='A key that can be used to fetch encrypted data')
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


def move_to_output_dir(work_dir, output_dir, uuid=None, files=None):
    """
    A list of files to move from work_dir to output_dir.
    """
    for fname in files:
        if uuid is None:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


def download_from_URL(job, input_args, ids, name):
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, name)
    url = input_args[name]
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
        except subprocess.CalledProcessError:
            raise RuntimeError(
                '\nNecessary file could not be acquired: {}. Check input URL'.format(url))
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    job.fileStore.updateGlobalFile(ids[name], file_path)


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


def docker_call(work_dir, tool_parameters, tool):
    """
    Makes subprocess call of a command to a docker container.
    :parameter tool: name of the Docker image to be used (e.g. computationalgenomicslab/samtools)
    """
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_parameters)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


# Start of Target Functions
def batch_start(job, input_args):
    # Create IDs for the shared files in the pipeline
    shared_ids = {x: job.fileStore.getEmptyFileStoreID() for x in ['ref.fasta', 'phase.vcf', 'mills.vcf', 'cosmic.vcf',
                                                                      'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict']}
    for file_name in ['ref.fasta', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf', 'cosmic.vcf']:
        job.addChildJobFn(download_from_URL, input_args, shared_ids, file_name)
    job.addFollowOnJobFn(batch_setup, input_args, shared_ids)


def batch_setup(job, input_args, shared_ids):
    # Docker tools used in the pipeline
    job.addChildJobFn(create_reference_index, (input_args, shared_ids))
    job.addChildJobFn(create_reference_dict, (input_args, shared_ids))
    job.addFollowOnJobFn(spawn_batch_jobs, input_args, shared_ids)


def spawn_batch_jobs(job, input_args, shared_ids):
    # Names for every input file used in the pipeline by each sample
    samples = []
    config = input_args['config']
    with open(config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                samples.append(line.strip().split(','))
    for sample in samples:
        job.addChildJobFn(start, shared_ids, input_args, sample)


def start(job, shared_ids, input_args, sample):
    uuid, normal_url, tumor_url = sample
    symbolic_inputs = ['normal.bam.bai', 'tumor.bam.bai', 'normal.intervals', 'tumor.intervals',
                   'normal.indel.bam', 'tumor.indel.bam', 'normal.indel.bai', 'tumor.indel.bai',
                   'normal.recal.table', 'tumor.recal.table', 'normal.bqsr.bam', 'tumor.bqsr.bam',
                   'normal.bqsr.bai', 'tumor.bqsr.bai', 'mutect.vcf', 'tumor.bam', 'normal.bam']

    ids = shared_ids.copy()
    ids.update( {x: job.fileStore.getEmptyFileStoreID() for x in symbolic_inputs} )
    sample_input = dict(input_args)
    # Update input
    sample_input['uuid'] = uuid
    sample_input['normal.bam'] = normal_url
    sample_input['tumor.bam'] = tumor_url
    sample_input['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    job_vars = (sample_input, ids)

    # Download sample bams and launch pipeline
    if input_args['ssec'] is None:
        job.addChildJobFn(download_from_URL, sample_input, ids, 'normal.bam')
        job.addChildJobFn(download_from_URL, sample_input, ids, 'tumor.bam')
    else:
        # TODO: Add capability to download bams with encryption
        pass
    job.addFollowOnJobFn(start_preprocessing, job_vars)


def create_reference_index(job, job_vars):
    """
    Uses Samtools to create reference index file (.fasta.fai)
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = docker_path(job.fileStore.readGlobalFile(ids['ref.fasta'], os.path.join(work_dir, 'ref.fasta')))
    # Call: Samtools
    command = ['faidx', ref_path]
    docker_call(work_dir, command, tool='computationalgenomicslab/samtools')
    # Update fileStore for output
    job.fileStore.updateGlobalFile(ids['ref.fasta.fai'], os.path.join(work_dir, 'ref.fasta.fai'))


def create_reference_dict(job, job_vars):
    """
    Uses Picardtools to create reference dictionary (.dict) for the sample
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = job.fileStore.readGlobalFile(ids['ref.fasta'], os.path.join(work_dir, 'ref.fasta'))
    # Call: picardtools
    output = os.path.splitext(docker_path(ref_path))[0]
    command = ['CreateSequenceDictionary', 'R={}'.format(docker_path(ref_path)), 'O={}.dict'.format(output)]
    docker_call(work_dir, command, tool='computationalgenomicslab/picardtools')
    # Update fileStore for output
    job.fileStore.updateGlobalFile(ids['ref.dict'], os.path.join(work_dir, 'ref.dict'))


def index(job, job_vars, sample):
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    bam = '{}.bam'.format(sample)
    path = return_input_paths(job, work_dir, ids, bam)
    # Call: index the normal.bam
    parameters = ['index', '{}'.format(docker_path(path))]
    docker_call(work_dir, tool_parameters=parameters, tool='computationalgenomicslab/samtools')
    # Update FileStore and call child
    job.fileStore.updateGlobalFile(ids[bam + '.bai'], os.path.join(work_dir, bam) + '.bai')
    job.addChildJobFn(rtc, job_vars, sample, cores=int(input_args['cpu_count']), memory='10 G', disk='15 G')


def start_preprocessing(job, job_vars):
    # Add children and followOn jobs
    job.addChildJobFn(index, job_vars, 'normal', cores=1, memory='1 G', disk='8 G')
    job.addChildJobFn(index, job_vars, 'tumor', cores=1, memory='1 G', disk='8 G')
    job.addFollowOnJobFn(mutect, job_vars, cores=1, memory='3 G', disk='30 G')


def rtc(job, job_vars, sample):
    """
    Creates <type>.intervals file needed for indel realignment
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve input file paths
    ref_fasta, bam, ref_fai, ref_dict, bam_bai,\
    phase_vcf, mills_vcf = return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.bam'.format(sample),
                                              'ref.fasta.fai', 'ref.dict', '{}.bam.bai'.format(sample), 'phase.vcf', 'mills.vcf')

    # Output file path
    output = os.path.join(work_dir, '{}.intervals'.format(sample))
    # Call: GATK -- RealignerTargetCreator
    parameters = ['-T', 'RealignerTargetCreator',
                  '-nt', input_args['cpu_count'],
                  '-R', docker_path(ref_fasta),
                  '-I', docker_path(bam),
                  '-known', docker_path(phase_vcf),
                  '-known', docker_path(mills_vcf),
                  '--downsampling_type', 'NONE',
                  '-o', docker_path(output)]

    docker_call(work_dir, tool_parameters=parameters, tool='computationalgenomicslab/gatk')
    # Update fileStore and spawn child job
    job.fileStore.updateGlobalFile(ids['{}.intervals'.format(sample)], output)
    job.addChildJobFn(ir, job_vars, sample, cores=1, memory='4 G', disk='12 G')


def ir(job, job_vars, sample):
    """
    Creates realigned bams using <sample>.intervals file from previous step
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
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
    docker_call(work_dir, command, tool='computationalgenomicslab/gatk')
    # Update fileStore and spawn child job
    job.fileStore.updateGlobalFile(ids['{}.indel.bam'.format(sample)], output)
    job.fileStore.updateGlobalFile(ids['{}.indel.bai'.format(sample)], os.path.splitext(output)[0] + '.bai')
    job.addChildJobFn(br, job_vars, sample, cores=int(input_args['cpu_count']), memory='10 G', disk='15 G')


def br(job, job_vars, sample):
    """
    Creates recal table to perform Base Quality Score Recalibration
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve input file paths
    ref_fasta, indel_bam, dbsnp_vcf, ref_fai, \
    ref_dict, bam_bai = return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.indel.bam'.format(sample),
                                           'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict', '{}.indel.bai'.format(sample))
    # Output file path
    output = os.path.join(work_dir, '{}.recal.table'.format(sample))
    # Call: GATK -- IndelRealigner
    command = ['-T', 'BaseRecalibrator',
              '-nct', input_args['cpu_count'],
              '-R', ref_fasta,
              '-I', indel_bam,
              '-knownSites', dbsnp_vcf,
              '-o', docker_path(output)]
    docker_call(work_dir, command, tool='computationalgenomicslab/gatk')
    # Update fileStore and spawn child job
    job.fileStore.updateGlobalFile(ids['{}.recal.table'.format(sample)], output)
    job.addChildJobFn(pr, job_vars, sample, cores=int(input_args['cpu_count']), memory='15 G', disk='30 G')


def pr(job, job_vars, sample):
    """
    Create bqsr bam
    """
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve input file paths
    ref_fasta, indel_bam, dbsnp_vcf, ref_fai, \
    ref_dict, bam_bai, recal_table = return_input_paths(job, work_dir, ids, 'ref.fasta', '{}.indel.bam'.format(sample),
                                                        'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict', '{}.indel.bai'.format(sample),
                                                        '{}.recal.table'.format(sample))
    # Output file
    output = os.path.join(work_dir, '{}.bqsr.bam'.format(sample))
    # Call: GATK -- PrintReads
    command = ['-T', 'PrintReads',
              '-nct', input_args['cpu_count'],
              '-R', ref_fasta,
              '--emit_original_quals',
              '-I', indel_bam,
              '-BQSR', recal_table,
              '-o', docker_path(output)]
    docker_call(work_dir, command, tool='computationalgenomicslab/gatk')
    # Update GlobalFileStore
    job.fileStore.updateGlobalFile(ids['{}.bqsr.bam'.format(sample)], output)
    job.fileStore.updateGlobalFile(ids['{}.bqsr.bai'.format(sample)], os.path.splitext(output)[0] + '.bai')


def mutect(job, job_vars):
    # Unpack convenience variables for job
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)
    # Retrieve input files
    ref_fasta, normal_bqsr_bam, tumor_bqsr_bam, \
    dbsnp_vcf, normal_bqsr_bai, tumor_bqsr_bai, \
    ref_fai, ref_dict, cosmic_vcf = return_input_paths(job, work_dir, ids, 'ref.fasta', 'normal.bqsr.bam', 'tumor.bqsr.bam',
                                           'dbsnp.vcf', 'normal.bqsr.bai', 'tumor.bqsr.bai', 'ref.fasta.fai', 'ref.dict', 'cosmic.vcf')
    # Output VCF
    uuid = input_args['uuid']
    output_name = '{}_mutect.vcf'.format(uuid)
    mut_out = docker_path(os.path.join(work_dir, 'mutect.out'))
    mut_cov = docker_path(os.path.join(work_dir, 'mutect.cov'))
    # Call: MuTect
    command = ['--analysis_type', 'MuTect',
              '--reference_sequence', ref_fasta,
              '--cosmic', cosmic_vcf,
              '--dbsnp', dbsnp_vcf,
              '--input_file:normal', normal_bqsr_bam,
              '--input_file:tumor', tumor_bqsr_bam,
              '--tumor_lod', str(10),
              '--out', mut_out,
              '--coverage_file', mut_cov,
              '--vcf', docker_path(output_name)]
    docker_call(work_dir, command, tool='computationalgenomicslab/mutect')
    # Update FileStoreID
    job.fileStore.updateGlobalFile(ids['mutect.vcf'], os.path.join(work_dir, output_name))
    move_to_output_dir(work_dir, output_dir, uuid=None, files=[os.path.join(work_dir, output_name)] )


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
                  'uuid': None,
                  'normal.bam': None,
                  'tumor.bam': None,
                  'cpu_count': str(multiprocessing.cpu_count())}

    # Launch jobs
    Job.Runner.startToil(Job.wrapJobFn(batch_start, input_args), args)
