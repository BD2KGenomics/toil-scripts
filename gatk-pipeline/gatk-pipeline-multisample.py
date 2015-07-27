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
import shutil
import subprocess
import multiprocessing
from jobTree.target import Target


def build_parser():
    """
    Contains arguments for the all of necessary input files
    """
    # TODO: Change to read in from S3 bucket for samples
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', required=True, help="Reference Genome URL")
    parser.add_argument('-n6', '--n6', required=True, help='Normal BAM URL. Format: UUID.normal.bam')
    parser.add_argument('-t6', '--t6', required=True, help='Tumor BAM URL. Format: UUID.tumor.bam')
    parser.add_argument('-n7', '--n7', required=True, help='Normal BAM URL. Format: UUID.normal.bam')
    parser.add_argument('-t7', '--t7', required=True, help='Tumor BAM URL. Format: UUID.tumor.bam')
    parser.add_argument('-n8', '--n8', required=True, help='Normal BAM URL. Format: UUID.normal.bam')
    parser.add_argument('-t8', '--t8', required=True, help='Tumor BAM URL. Format: UUID.tumor.bam')
    parser.add_argument('-n9', '--n9', required=True, help='Normal BAM URL. Format: UUID.normal.bam')
    parser.add_argument('-t9', '--t9', required=True, help='Tumor BAM URL. Format: UUID.tumor.bam')
    parser.add_argument('-p', '--phase', required=True, help='1000G_phase1.indels.hg19.sites.fixed.vcf URL')
    parser.add_argument('-m', '--mills', required=True, help='Mills_and_1000G_gold_standard.indels.hg19.sites.vcf URL')
    parser.add_argument('-d', '--dbsnp', required=True, help='dbsnp_132_b37.leftAligned.vcf URL')
    parser.add_argument('-c', '--cosmic', required=True, help='b37_cosmic_v54_120711.vcf URL')
    parser.add_argument('-o', '--output_dir', required=True, help='Full path to final output dir')
    return parser


# Convenience functions used in the pipeline
def download_from_URL(target, input_args, ids, name):
    work_dir = target.fileStore.getLocalTempDir()
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
    target.fileStore.updateGlobalFile(ids[name], file_path)


def return_input_paths(target, work_dir, ids, *args):
    """
    Given one or more strings representing file_names, return the paths to those files.
    """
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


def docker_call(work_dir, tool_command, tool):
    """
    Makes subprocess call of a command to a docker container.
    :parameter tool: name of the Docker image to be used (e.g. computationalgenomicslab/samtools)
    """
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_command.split())
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


# Start of Target Functions
def batch_start(target, input_args):
    # Create IDs for the shared files in the pipeline
    shared_ids = {x: target.fileStore.getEmptyFileStoreID() for x in ['ref.fasta', 'phase.vcf', 'mills.vcf', 'cosmic.vcf',
                                                                      'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict']}
    for file_name in ['ref.fasta', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf', 'cosmic.vcf']:
        target.addChildTargetFn(download_from_URL, input_args, shared_ids, file_name)
    target.addFollowOnTargetFn(batch_setup, input_args, shared_ids)


def batch_setup(target, input_args, shared_ids):
    # Docker tools used in the pipeline
    tools = {'samtools': 'computationalgenomicslab/samtools',
             'picard': 'computationalgenomicslab/picardtools',
             'mutect': 'computationalgenomicslab/mutect',
             'gatk': 'computationalgenomicslab/gatk'}
    target.addChildTargetFn(create_reference_index, (input_args, shared_ids, tools), cpu=1, memory=375000000, disk=375000000)
    target.addChildTargetFn(create_reference_dict, (input_args, shared_ids, tools), cpu=1, memory=375000000, disk=375000000)
    target.addFollowOnTargetFn(spawn_batch_jobs, input_args, shared_ids, tools)


def spawn_batch_jobs(target, input_args, shared_ids, tools):
    # Names for every input file used in the pipeline by each sample
    symbolic_inputs = ['normal.bam.bai', 'tumor.bam.bai', 'normal.intervals', 'tumor.intervals',
                       'normal.indel.bam', 'tumor.indel.bam', 'normal.indel.bai', 'tumor.indel.bai',
                       'normal.recal.table', 'tumor.recal.table', 'normal.bqsr.bam', 'tumor.bqsr.bam',
                       'normal.bqsr.bai', 'tumor.bqsr.bai', 'mutect.vcf', 'tumor.bam', 'normal.bam']
    for i in xrange(6,10):
        # Create unique IDs for each item created by each sample
        ids = dict(shared_ids.items() + {x: target.fileStore.getEmptyFileStoreID() for x in symbolic_inputs}.items())
        new_input = dict(input_args)
        new_input['normal.bam'] = input_args['n{}'.format(i)]
        new_input['tumor.bam'] = input_args['t{}'.format(i)]
        target.addChildTargetFn(download_bam, (new_input, ids, tools), cpu=.2, memory=1e8, disk=6.25e8)



def download_bam(target, target_vars):
    input_args, ids, tools = target_vars
    target.addChildTargetFn(download_from_URL, input_args, ids, 'normal.bam')
    target.addChildTargetFn(download_from_URL, input_args, ids, 'tumor.bam')
    target.addFollowOnTargetFn(start, target_vars)


def start(target, target_vars):
    # Add children and followOn targets
    target.addChildTargetFn(index, target_vars, 'normal', cpu=1, memory=1e9, disk=8e9)
    target.addChildTargetFn(index, target_vars, 'tumor', cpu=1, memory=1e9, disk=8e9)
    target.addFollowOnTargetFn(start_preprocessing, target_vars)


def create_reference_index(target, target_vars):
    """
    Uses Samtools to create reference index file (.fasta.fai)
    """
    # Unpack convenience variables for target
    input_args, ids, tools = target_vars
    work_dir = target.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = docker_path(target.fileStore.readGlobalFile(ids['ref.fasta'], os.path.join(work_dir, 'ref.fasta')))
    # Call: Samtools
    command = 'faidx {}'.format(ref_path)
    docker_call(work_dir, command, tool=tools['samtools'])
    # Update fileStore for output
    target.fileStore.updateGlobalFile(ids['ref.fasta.fai'], os.path.join(work_dir, 'ref.fasta.fai'))


def create_reference_dict(target, target_vars):
    """
    Uses Picardtools to create reference dictionary (.dict) for the sample
    """
    # Unpack convenience variables for target
    input_args, ids, tools = target_vars
    work_dir = target.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = target.fileStore.readGlobalFile(ids['ref.fasta'], os.path.join(work_dir, 'ref.fasta'))
    # Call: picardtools
    output = os.path.splitext(docker_path(ref_path))[0]
    command = 'CreateSequenceDictionary R={} O={}.dict'.format(docker_path(ref_path), output)
    docker_call(work_dir, command, tool=tools['picard'])
    # Update fileStore for output
    target.fileStore.updateGlobalFile(ids['ref.dict'], os.path.join(work_dir, 'ref.dict'))


def index(target, target_vars, sample):
    # Unpack convenience variables for target
    input_args, ids, tools = target_vars
    work_dir = target.fileStore.getLocalTempDir()
    # Retrieve file path
    bam = '{}.bam'.format(sample)
    path = return_input_paths(target, work_dir, ids, bam)
    # Call: index the normal.bam
    command = 'index {}'.format(docker_path(path))
    docker_call(work_dir, command, tool=tools['samtools'])
    # Update FileStore for output
    target.fileStore.updateGlobalFile(ids[bam + '.bai'], os.path.join(work_dir, bam) + '.bai')


def start_preprocessing(target, target_vars):
    # Unpack convenience variables for target
    # input_args, symbolic_inputs, ids, tools = target_vars
    # work_dir = target.fileStore.getLocalTempDir()
    # Download files shared by parallel targets to avoid conflicts in fileStore
    # download_files(target, input_args, work_dir,  ids, 'phase.vcf', 'mills.vcf', 'dbsnp.vcf')
    # Add children and followOn targets
    target.addChildTargetFn(rtc, target_vars, 'normal', cpu=4, memory=10e9, disk=15e9)
    target.addChildTargetFn(rtc, target_vars, 'tumor', cpu=4, memory=10e9, disk=15e9)
    target.addFollowOnTargetFn(mutect, target_vars, cpu=1, memory=3e9, disk=30e9)


def rtc(target, target_vars, sample):
    """
    Creates <type>.intervals file needed for indel realignment
    """
    # Unpack convenience variables for target
    input_args, ids, tools = target_vars
    work_dir = target.fileStore.getLocalTempDir()
    # Retrieve input file paths
    ref_fasta, bam, ref_fai, ref_dict, bam_bai,\
    phase_vcf, mills_vcf = return_input_paths(target, work_dir, ids, 'ref.fasta', '{}.bam'.format(sample),
                                              'ref.fasta.fai', 'ref.dict', '{}.bam.bai'.format(sample), 'phase.vcf', 'mills.vcf')

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
              '-o {5} '.format(input_args['cpu_count'], docker_path(ref_fasta), docker_path(bam), docker_path(phase_vcf),
                               docker_path(mills_vcf), docker_path(output))
    docker_call(work_dir, command, tool=tools['gatk'])
    # Update fileStore and spawn child target
    target.fileStore.updateGlobalFile(ids['{}.intervals'.format(sample)], output)
    target.addChildTargetFn(ir, target_vars, sample, cpu=1, memory=3.5e9, disk=12e9)


def ir(target, target_vars, sample):
    """
    Creates realigned bams using <sample>.intervals file from previous step
    """
    # Unpack convenience variables for target
    input_args, ids, tools = target_vars
    work_dir = target.fileStore.getLocalTempDir()
    # Retrieve input file paths
    ref_fasta, bam, phase_vcf, \
    mills_vcf, intervals, ref_fai, \
    ref_dict, bam_bai = return_input_paths(target, work_dir, ids, 'ref.fasta', '{}.bam'.format(sample), 'phase.vcf',
                                           'mills.vcf', '{}.intervals'.format(sample), 'ref.fasta.fai', 'ref.dict',
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
    docker_call(work_dir, command, tool=tools['gatk'])
    # Update fileStore and spawn child target
    target.fileStore.updateGlobalFile(ids['{}.indel.bam'.format(sample)], output)
    target.fileStore.updateGlobalFile(ids['{}.indel.bai'.format(sample)], os.path.splitext(output)[0] + '.bai')
    target.addChildTargetFn(br, target_vars, sample, cpu=4, memory=10e9, disk=15e9)


def br(target, target_vars, sample):
    """
    Creates recal table to perform Base Quality Score Recalibration
    """
    # Unpack convenience variables for target
    input_args, ids, tools = target_vars
    work_dir = target.fileStore.getLocalTempDir()
    # Retrieve input file paths
    ref_fasta, indel_bam, dbsnp_vcf, ref_fai, \
    ref_dict, bam_bai = return_input_paths(target, work_dir, ids, 'ref.fasta', '{}.indel.bam'.format(sample),
                                           'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict', '{}.indel.bai'.format(sample))
    # Output file path
    output = os.path.join(work_dir, '{}.recal.table'.format(sample))
    # Call: GATK -- IndelRealigner
    command = '-T BaseRecalibrator ' \
              '-nct {0} ' \
              '-R {1} ' \
              '-I {2} ' \
              '-knownSites {3} ' \
              '-o {4}'.format(input_args['cpu_count'], ref_fasta, indel_bam, dbsnp_vcf, docker_path(output))
    docker_call(work_dir, command, tool=tools['gatk'])
    # Update fileStore and spawn child target
    target.fileStore.updateGlobalFile(ids['{}.recal.table'.format(sample)], output)
    target.addChildTargetFn(pr, target_vars, sample, cpu=4, memory=15e9, disk=30e9)


def pr(target, target_vars, sample):
    """
    Create bqsr bam
    """
    # Unpack convenience variables for target
    input_args, ids, tools = target_vars
    work_dir = target.fileStore.getLocalTempDir()
    # Retrieve input file paths
    ref_fasta, indel_bam, dbsnp_vcf, ref_fai, \
    ref_dict, bam_bai, recal_table = return_input_paths(target, work_dir, ids, 'ref.fasta', '{}.indel.bam'.format(sample),
                                                        'dbsnp.vcf', 'ref.fasta.fai', 'ref.dict', '{}.indel.bai'.format(sample),
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
    docker_call(work_dir, command, tool=tools['gatk'])
    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['{}.bqsr.bam'.format(sample)], output)
    target.fileStore.updateGlobalFile(ids['{}.bqsr.bai'.format(sample)], os.path.splitext(output)[0] + '.bai')


def mutect(target, target_vars):
    # Unpack convenience variables for target
    input_args, ids, tools = target_vars
    work_dir = target.fileStore.getLocalTempDir()
    # Retrieve input files
    ref_fasta, normal_bqsr_bam, tumor_bqsr_bam, \
    dbsnp_vcf, normal_bqsr_bai, tumor_bqsr_bai, \
    ref_fai, ref_dict, cosmic_vcf = return_input_paths(target, work_dir, ids, 'ref.fasta', 'normal.bqsr.bam', 'tumor.bqsr.bam',
                                           'dbsnp.vcf', 'normal.bqsr.bai', 'tumor.bqsr.bai', 'ref.fasta.fai', 'ref.dict', 'cosmic.vcf')
    # Output VCF
    normal_uuid = input_args['normal.bam'].split('/')[-1].split('.')[0]
    tumor_uuid = input_args['tumor.bam'].split('/')[-1].split('.')[0]
    output_name = '{}-normal:{}-tumor.vcf'.format(normal_uuid, tumor_uuid)
    mut_out = docker_path(os.path.join(work_dir, 'mutect.out'))
    mut_cov = docker_path(os.path.join(work_dir, 'mutect.cov'))
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
                                  tumor_bqsr_bam, mut_out, mut_cov, docker_path(output_name))
    docker_call(work_dir, command, tool=tools['mutect'])
    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids['mutect.vcf'], os.path.join(work_dir,  output_name))
    shutil.copyfile(os.path.join(work_dir, output_name), os.path.join(input_args['output_dir'], output_name))


if __name__ == '__main__':
    # Define Parser object and add to jobTree
    parser = build_parser()
    Target.Runner.addJobTreeOptions(parser)
    args = parser.parse_args()
    # Variables to pass to initial target
    input_args = {'ref.fasta': args.reference,
                  'n6' : args.n6,
                  't6' : args.t6,
                  'n7' : args.n7,
                  't7' : args.t7,
                  'n8' : args.n8,
                  't8' : args.t8,
                  'n9' : args.n9,
                  't9' : args.t9,
                  'phase.vcf': args.phase,
                  'mills.vcf': args.mills,
                  'dbsnp.vcf': args.dbsnp,
                  'cosmic.vcf': args.cosmic,
                  'output_dir': args.output_dir,
                  'cpu_count': str(multiprocessing.cpu_count())}
    # Ensure BAMs are in the appropriate format ( UUID.<Tumor/Normal>.bam)
    for i in xrange(6,10):
        for bam in [input_args['n{}'.format(i)], input_args['t{}'.format(i)]]:
            if len(bam.split('/')[-1].split('.')) != 3:
                raise RuntimeError('{} BAM is not in the appropriate format: \
                UUID.normal.bam or UUID.tumor.bam'.format(str(bam).split('.')[1]))
    # Launch jobs
    Target.Runner.startJobTree(Target.wrapTargetFn(batch_start, input_args), args)
    Target.Runner.cleanup(args)
