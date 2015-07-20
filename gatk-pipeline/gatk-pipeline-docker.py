#!/usr/bin/env python2.7
# John Vivian
# 7-17-15

"""
Tree Structure of GATK Pipeline
   0------------> 5 ----------> 14
/ / \ \          / \
1 2 3 4         6   7
                |   |
                8   9
                |   |
                10  11
                |   |
                12  13

0  = Start node
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

0-4, 6-13 are "children targets"
5 and 14 are "follow-on targets"
=========================================================================
:Directory Structure:
work_dir = <args.work_dir>/<script_name>/<UUID4>
=========================================================================
:Dependencies:
curl            - apt-get install curl
docker          - apt-get install docker (or 'docker.io' for linux)
jobTree         - https://github.com/benedictpaten/jobTree
"""
import argparse
from collections import OrderedDict
import os
import subprocess
import uuid
import multiprocessing
from jobTree.target import Target


def build_parser():
    """
    Contains arguments for the all of necessary input files
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', required=True, help="Reference Genome URL")
    parser.add_argument('-n', '--normal', required=True, help='Normal BAM URL. Format: UUID.normal.bam')
    parser.add_argument('-t', '--tumor', required=True, help='Tumor BAM URL. Format: UUID.tumor.bam')
    parser.add_argument('-p', '--phase', required=True, help='1000G_phase1.indels.hg19.sites.fixed.vcf URL')
    parser.add_argument('-m', '--mills', required=True, help='Mills_and_1000G_gold_standard.indels.hg19.sites.vcf URL')
    parser.add_argument('-d', '--dbsnp', required=True, help='dbsnp_132_b37.leftAligned.vcf URL')
    parser.add_argument('-c', '--cosmic', required=True, help='b37_cosmic_v54_120711.vcf URL')
    parser.add_argument('-g', '--gatk', required=True, help='GenomeAnalysisTK.jar')
    parser.add_argument('-u', '--mutect', required=True, help='Mutect.jar')
    parser.add_argument('-w', '--work_dir', required=True, help='Where you wanna work from? (full path please)')
    return parser


# Convenience functions used in the pipeline
def download_files(target, input_args, ids, *args):
    """
    Given one or more strings representing file names, download those files and return the paths.
    """
    work_dir = input_args['work_dir']
    # for every file in *args, download to the working directory and update fileStore, returning the docker path.
    paths = OrderedDict()
    urls_and_names = [(input_args[name], os.path.join(work_dir, name)) for name in args]
    pool = multiprocessing.Pool(processes=int(input_args['cpu_count']))
    pool.map(download_worker, urls_and_names)
    for file_path in [os.path.join(work_dir, name) for name in args]:
        paths[name] = docker_path(file_path)
        target.fileStore.updateGlobalFile(ids[name], file_path)
        if len(args) == 1:
            paths = docker_path(file_path)

    return paths


def download_worker(url_and_path):
    url, file_path = url_and_path
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--create-dir', url, '-o', file_path])
        except subprocess.CalledProcessError:
            raise RuntimeError(
                '\nNecessary file could not be acquired: {}. Check input URL'.format(url))
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)


def return_input_paths(target, input_args, ids, *args):
    """
    Given one or more strings representing file_names, return the paths to those files.
    """
    work_dir = input_args['work_dir']
    # for every file in *args, place in work_dir via the FileStore and then return the mounted docker path.
    paths = OrderedDict()
    for name in args:
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = docker_path(target.fileStore.readGlobalFile(ids[name], os.path.join(work_dir, name)))
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


def docker_call(input_args, tool_command, tool):
    """
    Makes subprocess call of a command to a docker container.
    :parameter tool: name of the Docker image to be used (e.g. computationalgenomicslab/samtools)
    """
    work_dir = input_args['work_dir']
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_command.split())
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


# Start of Target Functions
def start(target, input_args):
    # Names for every input file used in the pipeline
    symbolic_inputs = input_args.keys() + \
                      ['ref.fai', 'ref.dict', 'normal.bam.bai', 'tumor.bam.bai', 'normal.intervals', 'tumor.intervals',
                       'normal.indel.bam', 'tumor.indel.bam', 'normal.indel.bai', 'tumor.indel.bai',
                       'normal.recal.table', 'tumor.recal.table', 'normal.bqsr.bam', 'tumor.bqsr.bam',
                       'normal.bqsr.bai', 'tumor.bqsr.bai', 'mutect.vcf']
    # FileStore IDs for all files used in the pipeline
    ids = {x: target.fileStore.getEmptyFileStoreID() for x in symbolic_inputs}
    # Docker tools used in the pipeline
    tools = {'samtools': 'computationalgenomicslab/samtools',
             'picard': 'computationalgenomicslab/picardtools',
             'mutect': 'computationalgenomicslab/mutect',
             'gatk': 'computationalgenomicslab/gatk'}
    # Tuple encapsulation of variables used by targets
    target_vars = (input_args, symbolic_inputs, ids, tools)
    # Download file shared by parallel targets to avoid conflicts in fileStore
    download_files(target, input_args, ids, 'ref.fasta')
    # Add children and followOn targets
    target.addChildTargetFn(create_reference_index, target_vars)
    target.addChildTargetFn(create_reference_dict, target_vars)
    target.addChildTargetFn(index, target_vars, 'normal')
    target.addChildTargetFn(index, target_vars, 'tumor')
    target.addFollowOnTargetFn(start_preprocessing, target_vars)


def create_reference_index(target, target_vars):
    """
    Uses Samtools to create reference index file (.fasta.fai)
    """
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # Retrieve file path
    ref_path = return_input_paths(target, input_args, ids, 'ref.fasta')
    # Call: Samtools
    command = 'faidx {}'.format(ref_path)
    docker_call(input_args, command, tool=tools['samtools'])
    # Update fileStore for output
    target.fileStore.updateGlobalFile(ids['ref.fai'], os.path.join(work_dir, 'ref.fasta') + '.fai')


def create_reference_dict(target, target_vars):
    """
    Uses Picardtools to create reference dictionary (.dict) for the sample
    """
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # Retrieve file path
    ref_path = return_input_paths(target, input_args, ids, 'ref.fasta')
    # Call: picardtools
    output = os.path.splitext(docker_path(ref_path))[0]
    command = 'CreateSequenceDictionary R={} O={}.dict'.format(docker_path(ref_path), output)
    docker_call(input_args, command, tool=tools['picard'])
    # Update fileStore for output
    target.fileStore.updateGlobalFile(ids['ref.dict'], os.path.join(work_dir, 'ref') + '.dict')


def index(target, target_vars, sample):
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # Retrieve file path
    bam = '{}.bam'.format(sample)
    path = download_files(target, input_args, ids, bam)
    # Call: index the normal.bam
    command = 'index {}'.format(docker_path(path))
    docker_call(input_args, command, tool=tools['samtools'])
    # Update FileStore for output
    target.fileStore.updateGlobalFile(ids[bam], os.path.join(work_dir, bam) + '.bai')


def start_preprocessing(target, target_vars):
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    # Download files shared by parallel targets to avoid conflicts in fileStore
    download_files(target, input_args, ids, 'phase.vcf', 'mills.vcf', 'dbsnp.vcf')
    # Add children and followOn targets
    target.addChildTargetFn(rtc, target_vars, 'normal')
    target.addChildTargetFn(rtc, target_vars, 'tumor')
    target.addFollowOnTargetFn(mutect, target_vars)


def rtc(target, target_vars, sample):
    """
    Creates <type>.intervals file needed for indel realignment
    """
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # Retrieve input file paths
    ref_fasta, bam, ref_fai, ref_dict, bam_bai,\
    phase_vcf, mills_vcf = return_input_paths(target, input_args, ids, 'ref.fasta', '{}.bam'.format(sample),
                                              'ref.fai', 'ref.dict', '{}.bam.bai'.format(sample), 'phase.vcf', 'mills.vcf')
    # Output file path
    output = os.path.join(work_dir, '{}.intervals'.format(sample))
    # Call: GATK -- RealignerTargetCreator
    command = '-T RealignerTargetCreator ' \
              '-nt {0} ' \
              '-R {1} ' \
              '-I {2} ' \
              '-known {3} ' \
              '-known {4} ' \
              '--downsampling_type NONE ' \
              '-o {5} '.format(input_args['cpu_count'], ref_fasta, bam, phase_vcf, mills_vcf, docker_path(output))
    docker_call(input_args, command, tool=tools['gatk'])
    # Update fileStore and spawn child target
    target.fileStore.updateGlobalFile(ids['{}.intervals'.format(sample)], output)
    target.addChildTargetFn(ir, target_vars, sample)


def ir(target, target_vars, sample):
    """
    Creates realigned bams using <sample>.intervals file from previous step
    """
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # Retrieve input file paths
    ref_fasta, bam, phase_vcf, \
    mills_vcf, intervals, ref_fai, \
    ref_dict, bam_bai = return_input_paths(target, input_args, ids, 'ref.fasta', '{}.bam'.format(sample), 'phase.vcf',
                                           'mills.vcf', '{}.intervals'.format(sample), 'ref.fai', 'ref.dict',
                                           '{}.bam.bai'.format(sample))
    # Output file path
    output = os.path.join(work_dir, '{}.indel.bam'.format(sample))
    # Call: GATK -- IndelRealigner
    command = '-T IndelRealigner ' \
              '-R {0} ' \
              '-I {1} ' \
              '-known {2} ' \
              '-known {3} ' \
              '-targetIntervals {4} ' \
              '--downsampling_type NONE ' \
              '-maxReads {5} ' \
              '-maxInMemory {6} ' \
              '-o {7}'.format(ref_fasta, bam, phase_vcf, mills_vcf, intervals,
                              str(720000), str(5400000), docker_path(output))
    docker_call(input_args, command, tool=tools['gatk'])
    # Update fileStore and spawn child target
    target.fileStore.updateGlobalFile(ids['{}.indel.bam'.format(sample)], output)
    target.fileStore.updateGlobalFile(ids['{}.indel.bai'.format(sample)], os.path.splitext(output)[0] + '.bai')
    target.addChildTargetFn(br, target_vars, sample)


def br(target, target_vars, sample):
    """
    Creates recal table to perform Base Quality Score Recalibration
    """
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # Retrieve input file paths
    ref_fasta, indel_bam, dbsnp_vcf, ref_fai, \
    ref_dict, bam_bai = return_input_paths(target, input_args, ids, 'ref.fasta', '{}.indel.bam'.format(sample),
                                           'dbsnp.vcf', 'ref.fai', 'ref.dict', '{}.indel.bai'.format(sample))
    # Output file path
    output = os.path.join(work_dir, '{}.recal.table'.format(sample))
    # Call: GATK -- IndelRealigner
    command = '-T BaseRecalibrator ' \
              '-nct {0} ' \
              '-R {1} ' \
              '-I {2} ' \
              '-knownSites {3} ' \
              '-o {4}'.format(input_args['cpu_count'], ref_fasta, indel_bam, dbsnp_vcf, docker_path(output))
    docker_call(input_args, command, tool=tools['gatk'])
    # Update fileStore and spawn child target
    target.fileStore.updateGlobalFile(ids['{}.recal.table'.format(sample)], output)
    target.addChildTargetFn(pr, target_vars, sample)


def pr(target, target_vars, sample):
    """
    Create bqsr bam
    """
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # Retrieve input file paths
    ref_fasta, indel_bam, dbsnp_vcf, ref_fai, \
    ref_dict, bam_bai, recal_table = return_input_paths(target, input_args, ids, 'ref.fasta', '{}.indel.bam'.format(sample),
                                                        'dbsnp.vcf', 'ref.fai', 'ref.dict', '{}.indel.bai'.format(sample),
                                                        '{}.recal.table'.format(sample))
    # Output file
    output = os.path.join(work_dir, '{}.bqsr.bam'.format(sample))
    # Call: GATK -- PrintReads
    command = '-T PrintReads ' \
              '-nct {0} ' \
              '-R {1} ' \
              '--emit_original_quals ' \
              '-I {2} ' \
              '-BQSR {3} ' \
              '-o {4} '.format(input_args['cpu_count'], ref_fasta, indel_bam, recal_table, docker_path(output))
    docker_call(input_args, command, tool=tools['gatk'])
    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['{}.bqsr.bam'.format(sample)], output)
    target.fileStore.updateGlobalFile(ids['{}.bqsr.bai'.format(sample)], os.path.splitext(output)[0] + '.bai')


def mutect(target, target_vars):
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # Retrieve input files
    cosmic_vcf = docker_path(download_files(target, input_args, ids, 'cosmic.vcf'))
    ref_fasta, normal_bqsr_bam, tumor_bqsr_bam, \
    dbsnp_vcf, normal_bqsr_bai, tumor_bqsr_bai, \
    ref_fai, ref_dict = return_input_paths(target, input_args, ids, 'ref.fasta', 'normal.bqsr.bam', 'tumor.bqsr.bam',
                                           'dbsnp.vcf', 'normal.bqsr.bai', 'tumor.bqsr.bai', 'ref.fai', 'ref.dict')
    # Output VCF
    normal_uuid = input_args['normal.bam'].split('/')[-1].split('.')[0]
    tumor_uuid = input_args['tumor.bam'].split('/')[-1].split('.')[0]
    output_name = '{}-normal:{}-tumor.vcf'.format(normal_uuid, tumor_uuid)
    output = os.path.join(work_dir, output_name)
    mut_out = docker_path(os.path.join(input_args['work_dir'], 'mutect.out'))
    mut_cov = docker_path(os.path.join(input_args['work_dir'], 'mutect.cov'))
    # Call: MuTect
    command = '--analysis_type MuTect ' \
              '--reference_sequence {0} ' \
              '--cosmic {1} ' \
              '--dbsnp {2} ' \
              '--input_file:normal {3} ' \
              '--input_file:tumor {4} ' \
              '--tumor_lod 10 ' \
              '--out {5} ' \
              '--coverage_file {6} ' \
              '--vcf {7} '.format(ref_fasta, cosmic_vcf, dbsnp_vcf, normal_bqsr_bam,
                                  tumor_bqsr_bam, mut_out, mut_cov, docker_path(output))
    docker_call(input_args, command, tool=tools['mutect'])
    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids['mutect.vcf'], os.path.join(work_dir,  output_name))
    target.addChildTargetFn(teardown, target_vars)


def teardown(target, target_vars):
    # Unpack convenience variables for target
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']
    # For every file in work_dir (that is not final output), remove.
    files = [os.path.join(work_dir, f) for f in os.listdir(work_dir) if 'tumor.vcf' not in f]
    for f in files:
        os.remove(f)


if __name__ == '__main__':
    # Define Parser object and add to jobTree
    parser = build_parser()
    Target.Runner.addJobTreeOptions(parser)
    args = parser.parse_args()
    # Variables to pass to initial target
    work_dir = os.path.join(str(args.work_dir),
                            'bd2k-{}'.format(os.path.basename(__file__).split('.')[0]), str(uuid.uuid4()))
    input_args = {'ref.fasta': args.reference,
                  'normal.bam': args.normal,
                  'tumor.bam': args.tumor,
                  'phase.vcf': args.phase,
                  'mills.vcf': args.mills,
                  'dbsnp.vcf': args.dbsnp,
                  'cosmic.vcf': args.cosmic,
                  'gatk.jar': args.gatk,
                  'mutect.jar': args.mutect,
                  'work_dir': work_dir,
                  'cpu_count': str(multiprocessing.cpu_count())}
    # Ensure BAMs are in the appropriate format ( UUID.<Sample>.bam)
    for bam in [args.normal, args.tumor]:
        if len(bam.split('/')[-1].split('.')) != 3:
            raise RuntimeError('{} BAM is not in the appropriate format: \
            UUID.normal.bam or UUID.tumor.bam'.format(str(bam).split('.')[1]))
    # Launch jobs
    Target.Runner.startJobTree(Target.wrapTargetFn(start, input_args), args)
    Target.Runner.cleanup(args)
