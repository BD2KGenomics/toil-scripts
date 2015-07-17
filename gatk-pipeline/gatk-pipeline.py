#!/usr/bin/env python2.7
# John Vivian
# 7-15-15

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
jobTree         - https://github.com/benedictpaten/jobTree
"""
import argparse
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
def download_input(target, input_args, ids, name):
    """
    Downloads file if not present from supplied URL.
    Updates the FileStoreID to point to a file.
    :name: Key from self.input_urls.
    :returns: Path to file (work_dir path)
    :rtype: str
    """
    # Get path to file
    work_dir = input_args['work_dir']
    file_path = os.path.join(work_dir, name)

    # Check if file exists, download if not present
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--create-dir', input_args[name], '-o', file_path])
        except subprocess.CalledProcessError:
            raise RuntimeError('\nNecessary file could not be acquired: {}. Check input URL')
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')

    assert os.path.exists(file_path)

    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids[name], file_path)

    return file_path


def read_and_rename_global_file(target, input_args, file_store_id, name):
    """
    Given a FileStoreID, returns the filepath linked to it.
    :new_extension: Adds an extension to the file at the file_path.
    :alternate_name: A path to a filename that you want to use. (allows directory control as well as name)
    """
    work_dir = input_args['work_dir']
    if not os.path.exists(os.path.join(work_dir, name)):
        name = target.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))

    return os.path.join(work_dir, name)


def docker_path(filepath):
    return os.path.join('/data', os.path.basename(filepath))


def docker_call(input_args, tool_command, tool):
    """
    Makes subprocess call of a command to a docker container.
    Abstracts away the docker commands needed to run the tool.
    :type tool_command: str
    :type tool: str
    :param tool: the docker image for a given tool
    """
    work_dir = input_args['work_dir']
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_command.split())
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


# Start of Target Functions
def start(target, input_args):
    # Create dictionary to hold empty FileStoreIDs
    symbolic_inputs = input_args.keys() + ['ref.fai', 'ref.dict', 'normal.bai', 'tumor.bai', 'normal.intervals',
                                           'tumor.intervals', 'normal.indel.bam', 'tumor.indel.bam', 'normal.indel.bai',
                                           'tumor.indel.bai', 'normal.recal.table', 'tumor.recal.table',
                                           'normal.bqsr.bam', 'tumor.bqsr.bam', 'normal.bqsr.bai', 'tumor.bqsr.bai',
                                           'mutect.vcf']
    ids = {x: target.fileStore.getEmptyFileStoreID() for x in symbolic_inputs}
    tools = {'samtools': 'jvivian/samtools:1.2',
             'picard': 'jvivian/picardtools:1.113',
             'mutect': 'computationalgenomicslab/mutect'}
    target_vars = (input_args, symbolic_inputs, ids, tools)

    # Download ref.fasta so it isn't updated in parallel
    download_input(target, input_args, ids, 'ref.fasta')

    target.addChildTargetFn(create_reference_index, target_vars)
    target.addChildTargetFn(create_reference_dict, target_vars)
    target.addChildTargetFn(create_normal_index, target_vars)
    target.addChildTargetFn(create_tumor_index, target_vars)

    target.addFollowOnTargetFn(start_preprocessing, target_vars)


def create_reference_index(target, target_vars):
    """
    Uses Samtools to create reference index file (.fasta.fai)
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve reference & store in FileStoreID
    ref_path = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')

    # Tool call
    command = 'samtools faidx {}'.format(docker_path(ref_path))
    docker_call(input_args, command, tool=tools['samtools'])

    # Update FileStoreID of output
    target.fileStore.updateGlobalFile(ids['ref.fai'], ref_path + '.fai')


def create_reference_dict(target, target_vars):
    """
    Uses Picardtools to create reference dictionary (.dict)
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve reference & store in FileStoreID
    ref_path = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')

    # Tool call
    output = os.path.splitext(docker_path(ref_path))[0]
    command = 'picard-tools CreateSequenceDictionary R={} O={}.dict'.format(docker_path(ref_path), output)
    docker_call(input_args, command, tool=tools['picard'])

    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids['ref.dict'], os.path.splitext(ref_path)[0] + '.dict')


def create_normal_index(target, target_vars):
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve normal bam
    normal_path = download_input(target, input_args, ids, 'normal.bam')

    # Tool call
    command = 'samtools index {}'.format(docker_path(normal_path))
    docker_call(input_args, command, tool=tools['samtools'])

    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids['normal.bai'], normal_path + '.bai')


def create_tumor_index(target, target_vars):
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve tumor bam
    tumor_path = download_input(target, input_args, ids, 'tumor.bam')

    # Tool call
    command = 'samtools index {}'.format(docker_path(tumor_path))
    docker_call(input_args, command, tool=tools['samtools'])

    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids['tumor.bai'], tumor_path + '.bai')


def start_preprocessing(target, target_vars):
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    download_input(target, input_args, ids, 'phase.vcf')
    download_input(target, input_args, ids, 'mills.vcf')
    download_input(target, input_args, ids, 'gatk.jar')
    download_input(target, input_args, ids, 'dbsnp.vcf')

    target.addChildTargetFn(normal_rtc, target_vars)
    target.addChildTargetFn(tumor_rtc, target_vars)
    target.addFollowOnTargetFn(mutect, target_vars)


def normal_rtc(target, target_vars):
    """
    Creates normal.intervals file
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    phase_vcf = read_and_rename_global_file(target, input_args, ids['phase.vcf'], name='phase.vcf')
    mills_vcf = read_and_rename_global_file(target, input_args, ids['mills.vcf'], name='mills.vcf')
    gatk_jar = read_and_rename_global_file(target, input_args, ids['gatk.jar'], name='gatk.jar')
    ref_fasta = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')
    normal_bam = read_and_rename_global_file(target, input_args, ids['normal.bam'], name='normal.bam')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')
    read_and_rename_global_file(target, input_args, ids['normal.bai'], name='normal.bai')

    # Output File
    output = os.path.join(input_args['work_dir'], 'normal.intervals')

    # Create interval file
    try:
        subprocess.check_call(['java', '-Xmx15g', '-jar', gatk_jar, '-T', 'RealignerTargetCreator',
                               '-nt', input_args['cpu_count'], '-R', ref_fasta, '-I', normal_bam, '-known', phase_vcf,
                               '-known', mills_vcf, '--downsampling_type', 'NONE', '-o', output])
    except subprocess.CalledProcessError:
        raise RuntimeError('RealignerTargetCreator failed to finish')
    except OSError:
        raise RuntimeError('Failed to find "java" or gatk_jar')

    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['normal.intervals'], output)

    # Spawn Child
    target.addChildTargetFn(normal_ir, target_vars)


def tumor_rtc(target, target_vars):
    """
    Creates tumor.intervals file
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    phase_vcf = read_and_rename_global_file(target, input_args, ids['phase.vcf'], name='phase.vcf')
    mills_vcf = read_and_rename_global_file(target, input_args, ids['mills.vcf'], name='mills.vcf')
    gatk_jar = read_and_rename_global_file(target, input_args, ids['gatk.jar'], name='gatk.jar')
    ref_fasta = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')
    tumor_bam = read_and_rename_global_file(target, input_args, ids['tumor.bam'], name='tumor.bam')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')
    read_and_rename_global_file(target, input_args, ids['tumor.bai'], name='tumor.bai')

    # Output File
    output = os.path.join(input_args['work_dir'], 'tumor.intervals')

    # Create interval file
    try:
        subprocess.check_call(['java', '-Xmx15g', '-jar', gatk_jar, '-T', 'RealignerTargetCreator',
                               '-nt', input_args['cpu_count'], '-R', ref_fasta, '-I', tumor_bam, '-known', phase_vcf,
                               '-known', mills_vcf, '--downsampling_type', 'NONE', '-o', output])
    except subprocess.CalledProcessError:
        raise RuntimeError('RealignerTargetCreator failed to finish')
    except OSError:
        raise RuntimeError('Failed to find "java" or gatk_jar')

    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['tumor.intervals'], output)

    # Spawn Child
    target.addChildTargetFn(tumor_ir, target_vars)


def normal_ir(target, target_vars):
    """
    Creates realigned normal bams
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    gatk_jar = read_and_rename_global_file(target, input_args, ids['gatk.jar'], name='gatk.jar')
    ref_fasta = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')
    normal_bam = read_and_rename_global_file(target, input_args, ids['normal.bam'], name='normal.bam')
    phase_vcf = read_and_rename_global_file(target, input_args, ids['phase.vcf'], name='phase.vcf')
    mills_vcf = read_and_rename_global_file(target, input_args, ids['mills.vcf'], name='mills.vcf')
    normal_intervals = read_and_rename_global_file(target, input_args, ids['normal.intervals'], name='normal.intervals')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')
    read_and_rename_global_file(target, input_args, ids['normal.bai'], name='normal.bai')

    # Output file
    output = os.path.join(input_args['work_dir'], 'normal.indel.bam')

    # Create interval file
    try:
        subprocess.check_call(['java', '-Xmx15g', '-jar', gatk_jar, '-T', 'IndelRealigner',
                               '-R', ref_fasta, '-I', normal_bam, '-known', phase_vcf, '-known', mills_vcf,
                               '-targetIntervals', normal_intervals, '--downsampling_type', 'NONE',
                               '-maxReads', str(720000), '-maxInMemory', str(5400000), '-o', output])
    except subprocess.CalledProcessError:
        raise RuntimeError('IndelRealignment failed to finish')
    except OSError:
        raise RuntimeError('Failed to find "java" or gatk_jar')

    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['normal.indel.bam'], output)
    target.fileStore.updateGlobalFile(ids['normal.indel.bai'], os.path.splitext(output)[0] + '.bai')

    # Spawn Child
    target.addChildTargetFn(normal_br, target_vars)


def tumor_ir(target, target_vars):
    """
    Creates realigned tumor bams
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    gatk_jar = read_and_rename_global_file(target, input_args, ids['gatk.jar'], name='gatk.jar')
    ref_fasta = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')
    tumor_bam = read_and_rename_global_file(target, input_args, ids['tumor.bam'], name='tumor.bam')
    phase_vcf = read_and_rename_global_file(target, input_args, ids['phase.vcf'], name='phase.vcf')
    mills_vcf = read_and_rename_global_file(target, input_args, ids['mills.vcf'], name='mills.vcf')
    tumor_intervals = read_and_rename_global_file(target, input_args, ids['tumor.intervals'], name='tumor.intervals')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')
    read_and_rename_global_file(target, input_args, ids['tumor.bai'], name='tumor.bai')

    # Output file
    output = os.path.join(input_args['work_dir'], 'tumor.indel.bam')

    # Create interval file
    try:
        subprocess.check_call(['java', '-Xmx15g', '-jar', gatk_jar, '-T', 'IndelRealigner',
                               '-R', ref_fasta, '-I', tumor_bam, '-known', phase_vcf, '-known', mills_vcf,
                               '-targetIntervals', tumor_intervals, '--downsampling_type', 'NONE',
                               '-maxReads', str(720000), '-maxInMemory', str(5400000), '-o', output])
    except subprocess.CalledProcessError:
        raise RuntimeError('IndelRealignment failed to finish')
    except OSError:
        raise RuntimeError('Failed to find "java" or gatk_jar')

    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['tumor.indel.bam'], output)
    target.fileStore.updateGlobalFile(ids['tumor.indel.bai'], os.path.splitext(output)[0] + '.bai' )

    # Spawn Child
    target.addChildTargetFn(tumor_br, target_vars)


def normal_br(target, target_vars):
    """
    Creates normal recal table
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    dbsnp_vcf = read_and_rename_global_file(target, input_args, ids['dbsnp.vcf'], name='dbsnp.vcf')
    gatk_jar = read_and_rename_global_file(target, input_args, ids['gatk.jar'], name='gatk.jar')
    ref_fasta = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')
    normal_indel_bam = read_and_rename_global_file(target, input_args, ids['normal.indel.bam'], name='normal.indel.bam')
    read_and_rename_global_file(target, input_args, ids['normal.indel.bai'], name='normal.indel.bai')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')

    # Output file
    output = os.path.join(input_args['work_dir'], 'normal.recal.table')

    # Create interval file
    try:
        subprocess.check_call(['java', '-Xmx7g', '-jar', gatk_jar, '-T', 'BaseRecalibrator',
                               '-nct', input_args['cpu_count'], '-R', ref_fasta, '-I', normal_indel_bam,
                               '-knownSites', dbsnp_vcf, '-o', output])
    except subprocess.CalledProcessError:
        raise RuntimeError('BaseRecalibrator failed to finish')
    except OSError:
        raise RuntimeError('Failed to find "java" or gatk_jar')

    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['normal.recal.table'], output)

    # Spawn Child
    target.addChildTargetFn(normal_pr, target_vars)


def tumor_br(target, target_vars):
    """
    Creates tumor recal table
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    dbsnp_vcf = read_and_rename_global_file(target, input_args, ids['dbsnp.vcf'], name='dbsnp.vcf')
    gatk_jar = read_and_rename_global_file(target, input_args, ids['gatk.jar'], name='gatk.jar')
    ref_fasta = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')
    tumor_indel_bam = read_and_rename_global_file(target, input_args, ids['tumor.indel.bam'], name='tumor.indel.bam')
    read_and_rename_global_file(target, input_args, ids['tumor.indel.bai'], name='tumor.indel.bai')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')

    # Output file
    output = os.path.join(input_args['work_dir'], 'tumor.recal.table')

    # Create interval file
    try:
        subprocess.check_call(['java', '-Xmx7g', '-jar', gatk_jar, '-T', 'BaseRecalibrator',
                               '-nct', input_args['cpu_count'], '-R', ref_fasta, '-I', tumor_indel_bam,
                               '-knownSites', dbsnp_vcf, '-o', output])
    except subprocess.CalledProcessError:
        raise RuntimeError('BaseRecalibrator failed to finish')
    except OSError:
        raise RuntimeError('Failed to find "java" or gatk_jar')

    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['tumor.recal.table'], output)

    # Spawn Child
    target.addChildTargetFn(tumor_pr, target_vars)


def normal_pr(target, target_vars):
    """
    Create normal.bqsr.bam
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    gatk_jar = read_and_rename_global_file(target, input_args, ids['gatk.jar'], name='gatk.jar')
    ref_fasta = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')
    normal_indel_bam = read_and_rename_global_file(target, input_args, ids['normal.indel.bam'], name='normal.indel.bam')
    normal_recal = read_and_rename_global_file(target, input_args, ids['normal.recal.table'], name='normal.recal.table')
    read_and_rename_global_file(target, input_args, ids['normal.indel.bai'], name='normal.indel.bai')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')

    # Output file
    output = os.path.join(input_args['work_dir'], 'normal.bqsr.bam')

    # Create interval file
    try:
        subprocess.check_call(['java', '-Xmx7g', '-jar', gatk_jar, '-T', 'PrintReads',
                               '-nct', input_args['cpu_count'], '-R', ref_fasta, '--emit_original_quals',
                               '-I', normal_indel_bam, '-BQSR', normal_recal, '-o', output])
    except subprocess.CalledProcessError:
        raise RuntimeError('PrintReads failed to finish')
    except OSError:
        raise RuntimeError('Failed to find "java" or gatk_jar')

    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['normal.bqsr.bam'], output)
    target.fileStore.updateGlobalFile(ids['normal.bqsr.bai'], os.path.splitext(output)[0] + '.bai')


def tumor_pr(target, target_vars):
    """
    Create tumor.bqsr.bam
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    gatk_jar = read_and_rename_global_file(target, input_args, ids['gatk.jar'], name='gatk.jar')
    ref_fasta = read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta')
    tumor_indel_bam = read_and_rename_global_file(target, input_args, ids['tumor.indel.bam'], name='tumor.indel.bam')
    tumor_recal = read_and_rename_global_file(target, input_args, ids['tumor.recal.table'], name='tumor.recal.table')
    read_and_rename_global_file(target, input_args, ids['tumor.indel.bai'], name='tumor.indel.bai')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')

    # Output file
    output = os.path.join(input_args['work_dir'], 'tumor.bqsr.bam')

    # Create interval file
    try:
        subprocess.check_call(['java', '-Xmx7g', '-jar', gatk_jar, '-T', 'PrintReads',
                               '-nct', input_args['cpu_count'], '-R', ref_fasta, '--emit_original_quals',
                               '-I', tumor_indel_bam, '-BQSR', tumor_recal, '-o', output])
    except subprocess.CalledProcessError:
        raise RuntimeError('PrintReads failed to finish')
    except OSError:
        raise RuntimeError('Failed to find "java" or gatk_jar')

    # Update GlobalFileStore
    target.fileStore.updateGlobalFile(ids['tumor.bqsr.bam'], output)
    target.fileStore.updateGlobalFile(ids['tumor.bqsr.bai'], os.path.splitext(output)[0] + '.bai')


def mutect(target, target_vars):
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve input files
    mutect_jar = docker_path(download_input(target, input_args, ids, name='mutect.jar'))
    cosmic_vcf = docker_path(download_input(target, input_args, ids, name='cosmic.vcf'))

    dbsnp_vcf = docker_path(read_and_rename_global_file(target, input_args, ids['dbsnp.vcf'], name='dbsnp.vcf'))
    normal_bqsr_bam = docker_path(read_and_rename_global_file(target, input_args, ids['normal.bqsr.bam'], name='normal.bqsr.bam'))
    tumor_bqsr_bam = docker_path(read_and_rename_global_file(target, input_args, ids['tumor.bqsr.bam'], name='tumor.bqsr.bam'))
    ref_fasta = docker_path(read_and_rename_global_file(target, input_args, ids['ref.fasta'], name='ref.fasta'))
    read_and_rename_global_file(target, input_args, ids['normal.bqsr.bai'], name='normal.bqsr.bai')
    read_and_rename_global_file(target, input_args, ids['tumor.bqsr.bai'], name='tumor.bqsr.bai')
    read_and_rename_global_file(target, input_args, ids['ref.fai'], name='ref.fasta.fai')
    read_and_rename_global_file(target, input_args, ids['ref.dict'], name='ref.dict')

    # Output VCF
    normal_uuid = input_args['normal.bam'].split('/')[-1].split('.')[0]
    tumor_uuid = input_args['tumor.bam'].split('/')[-1].split('.')[0]
    output = docker_path(os.path.join(input_args['work_dir'], '{}-normal:{}-tumor.vcf'.format(normal_uuid, tumor_uuid)))
    mut_out = docker_path(os.path.join(input_args['work_dir'], 'mutect.out'))
    mut_cov = docker_path(os.path.join(input_args['work_dir'], 'mutect.cov'))

    # Tool call
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
                                  tumor_bqsr_bam, mut_out, mut_cov, output)

    docker_call(input_args, command, tool=tools['mutect'])

    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids['mutect.vcf'],
                            os.path.join(input_args['work_dir'], '{}-normal:{}-tumor.vcf'.format(normal_uuid, tumor_uuid)))

    target.addChildTargetFn(teardown, target_vars)


def teardown(target, target_vars):
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars
    work_dir = input_args['work_dir']

    files = [os.path.join(work_dir, f) for f in os.listdir(work_dir) if 'tumor.vcf' not in f]
    for f in files:
        os.remove(f)


if __name__ == '__main__':
    parser = build_parser()
    Target.Runner.addJobTreeOptions(parser)
    args = parser.parse_args()

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

    # Ensure user supplied URLs to files and that BAMs are in the appropriate format
    for bam in [args.normal, args.tumor]:
        if len(bam.split('/')[-1].split('.')) != 3:
            raise RuntimeError('{} BAM is not in the appropriate format: \
            UUID.normal.bam or UUID.tumor.bam'.format(str(bam).split('.')[1]))

    Target.Runner.startJobTree(Target.wrapTargetFn(start, input_args), args)
    # Target.Runner.cleanup(args)
