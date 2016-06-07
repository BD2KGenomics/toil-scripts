#!/usr/bin/env python2.7
"""
@author Jacob Pfeil
@data 01/13/2016

Toil pipeline for processing bam files for GATK halpotype calling

1 = download shared data
2 = reference preprocessing
3 = download sample
4 = remove secondary alignments
5 = index sample
6 = sort sample
7 = mark duplicates
8 = index
9 = realigner target
10 = indel realignment
11 = index
12 = base recalibration
13 = output bqsr fil6
"""
import argparse
import multiprocessing
import os
import sys
import logging
import textwrap
import yaml

from copy import deepcopy
from urlparse import urlparse

from toil.job import Job
from toil_scripts.lib import require
from toil_scripts.lib.files import copy_to
from toil_scripts.lib.urls import download_url_job, s3am_upload
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import s3am_upload
from toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline import generate_file

from bd2k.util.processes import which

_log = logging.getLogger(__name__)

# Convenience functions used in the pipeline
def upload_or_move(filename, work_dir=None, output_dir=None, s3_key_path=None):
    if output_dir.startswith('S://'):
        s3am_upload(fpath=os.path.join(work_dir, filename),
                    s3_dir=output_dir,
                    s3_key_path=s3_key_path)

    elif output_dir:
        # We must copy files to the dest because of caching
        copy_to(filename, output_dir, work_dir)

    else:
        raise ValueError('No output_directory or s3_dir defined. Cannot determine where to store %s' % filename)


def get_files_from_filestore(job, workDir, inputDict):
    """
    Retrieves files from filestore and copies them to the work directory

    work_dir: str       Current working directory
    ids: dict           Dictionary of fileStore IDs {name: fileStoreID}
    """
    for name, fileStoreID in inputDict.iteritems():
        if not os.path.exists(os.path.join(workDir, name)):
            file_path = job.fileStore.readGlobalFile(fileStoreID, os.path.join(workDir, name))
        else:
            file_path = name
        inputDict[name] = file_path
    return inputDict


def samtools_faidx(job, fileStoreID, mock=False):
    """
    Uses Samtools faidx to make a fasta index file (.fasta.fai)

    fileStoreID str: fasta fileStore ID
    :returns fileStoreID to index file
    """
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': fileStoreID}
    get_files_from_filestore(job, work_dir, inputs)

    outputs = {'genome.fa.fai': None}
    # Call: Samtools
    command = ['faidx', 'genome.fa']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=inputs.keys(), outputs=outputs,
                mock=mock)
    outpath = os.path.join(work_dir, 'genome.fa.fai')
    # Write to fileStore
    return job.fileStore.writeGlobalFile(outpath)


def picard_create_sequence_dictionary(job, fileStoreID, memory=8, mock=False):
    """
    Uses Picardtools to create reference dictionary (.dict)

    :param fileStoreID str: fasta fileStoreID
    :param memory int: allocated memory resource in GB
    :param mock bool: runs tool in debug mode
    :return fileStoreID str: fasta.dict fileStoreID
    """
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': fileStoreID}
    outputs = {'genome.dict': None}
    get_files_from_filestore(job, work_dir, inputs)

    # Call: picardtools
    command = ['CreateSequenceDictionary', 'R=genome.fa', 'O=genome.dict']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=mock)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'genome.dict'))


def download_shared_files(job, config):
    """
    Downloads files shared by all samples in the pipeline

    :param config dict: Dictionary of input arguments (from main())
    :return dict: updated config dictionary with shared fileStoreIDS
    """
    job.fileStore.logToMaster('Downloading Shared Reference Files')
    reference_names = ['genome.fa', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf']
    for name in reference_names:
        key, _ = name.split('.')
        config[name] = job.addChildJobFn(download_url_job, url=config[key], name=name,
                                         s3_key_path=config['ssec']).rv()
    return config


def reference_preprocessing(job, config, mock=False):
    """
    Create index and dict file for reference

    :param config dict: Dictionary of shared fileStore IDs
    :return dict: updated config dictionary with shared fileStoreIDs
    """
    job.fileStore.logToMaster('Preparing Reference Files')
    fileStoreID = config['genome.fa']
    if 'genome.fa.fai' not in config:
        config['genome.fa.fai'] = job.addChildJobFn(samtools_faidx, fileStoreID, mock=mock).rv()
    if 'genome.dict' not in config:
        config['genome.dict'] = job.addChildJobFn(picard_create_sequence_dictionary, fileStoreID, mock=mock).rv()
    return config


def samtools_view(job, bamFileStoreID, flag='0', mock=False):
    """
    Outputs bam file using a samtools view flag value
    :param job: Job instance
    :param bamFileStoreID str: bam fileStoreID
    :param flag str: samtools defined flags
    :param mock bool: If True, run in mock mode
    :return str: bam fileStoreID

    '0x800'
    """
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'sample.bam': bamFileStoreID}
    outputs = {'sample.output.bam': None},
    get_files_from_filestore(job, work_dir, inputs)
    outpath = os.path.join(work_dir, 'sample.output.bam')

    command = ['view',
               '-b', '-o', '/data/sample.output.bam',
               '-F', flag,
               '/data/sample.bam']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=mock)
    return job.fileStore.writeGlobalFile(outpath)


def picard_sort_sam(job, fileStoreID, config, xmx=8):
    """
    Uses picardtools SortSam to sort a sample bam file

    :param job:
    :param fileStoreID:
    :param config:
    :param xmx:
    :return:
    """
    job.fileStore.logToMaster('Running Picard SortSam: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'sample.bam': fileStoreID}
    outputs={'sample.sorted.bam': None, 'sample.sorted.bai': None},
    get_files_from_filestore(job, work_dir, inputs)
    outpath_bam = os.path.join(work_dir, 'sample.sorted.bam')
    outpath_bai = os.path.join(work_dir, 'sample.sorted.bai')

    #Call: picardtools
    command = ['SortSam',
               'INPUT=sample.bam',
               'OUTPUT=sample.sorted.bam',
               'SORT_ORDER=coordinate',
               'CREATE_INDEX=true']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % xmx},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock'])
    bam_id = job.fileStore.writeGlobalFile(outpath_bam)
    bai_id = job.fileStore.writeGlobalFile(outpath_bai)
    return bam_id, bai_id


def picard_mark_duplicates(job, bamFileStoreID, baiFileStoreID, config, xmx=8):
    """
    Uses picardtools MarkDuplicates

    :param job:
    :param bamFileStoreID:
    :param baiFileStoreID:
    :param config:
    :param xmx:
    :return:
    """
    job.fileStore.logToMaster('Running Picard MarkDuplicates: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'sample.sorted.bam': bamFileStoreID,
              'sample.sorted.bai': baiFileStoreID}
    outputs={'sample.mkdups.bam': None, 'sample.mkdups.bai': None}
    outpath_bam = os.path.join(work_dir, 'sample.mkdups.bam')
    outpath_bai = os.path.join(work_dir, 'sample.mkdups.bai')

    # Retrieve file path
    get_files_from_filestore(job, work_dir, inputs)

    # Call: picardtools
    command = ['MarkDuplicates',
               'INPUT=sample.sorted.bam',
               'OUTPUT=sample.mkdups.bam',
               'METRICS_FILE=metrics.txt',
               'ASSUME_SORTED=true',
               'CREATE_INDEX=true']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % xmx},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=inputs.keys(), mock=config['mock'],
                outputs=outputs)

    bam_id = job.fileStore.writeGlobalFile(outpath_bam)
    bai_id = job.fileStore.writeGlobalFile(outpath_bai)
    return bam_id, bai_id


def gatk_realigner_target_creator(job, bamFileStoreID, baiFileStoreID, config, xmx=8):
    """
    Creates <type>.intervals file needed for indel realignment

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    job.fileStore.logToMaster('Running GATK RealignerTargetCreator: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'sample.sorted.bam': bamFileStoreID,
              'sample.sorted.bai': baiFileStoreID}
    references = ['genome.fa','genome.fa.fai', 'genome.dict', 'phase.vcf', 'mills.vcf']
    inputs.update({key: config[key] for key in references})
    get_files_from_filestore(job, work_dir, inputs)
    # Output file path
    outputs={'sample.intervals': None},
    output = os.path.join(work_dir, 'sample.intervals')
    # Call: GATK -- RealignerTargetCreator
    cores = multiprocessing.cpu_count()
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'RealignerTargetCreator',
                  '-nt', str(cores),
                  '-R', 'ref.fa',
                  '-I', 'sample.mkdups.bam',
                  '-known', 'phase.vcf',
                  '-known', 'mills.vcf',
                  '--downsampling_type', 'NONE',
                  '-o', 'sample.intervals']

    docker_call(work_dir=work_dir, parameters=parameters,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(), outputs=outputs,
                env={'JAVA_OPTS':'-Xmx%sg' % xmx},
                mock=config['mock'])
    return job.fileStore.writeGlobalFile(output)


def gatk_indel_realigner(job, bamFileStoreID, baiFileStoreID, intervalsFileStoreID, config, memory=8):
    """
    Creates realigned bams using <sample>.intervals file from previous step

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    job.fileStore.logToMaster('Running GATK IndelRealigner: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve input file paths
    references = ['genome.fa','genome.fa.fai', 'genome.dict', 'phase.vcf', 'mills.vcf']
    inputs = {'sample.mkdups.bam': bamFileStoreID, 'sample.mkdups.bam.bai': baiFileStoreID,
              'sample.intervals': intervalsFileStoreID}
    inputs.update({key: config[key] for key in references})
    get_files_from_filestore(job, work_dir, inputs)
    # Output file path
    outpath_bam = os.path.join(work_dir, 'sample.indel.bam')
    outpath_bai = os.path.join(work_dir, 'sample.indel.bai')
    outputs={'sample.indel.bam': None, 'sample.indel.bai': None},
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
                inputs=inputs, outputs=outputs,
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                mock=config['mock'])

    # Write to fileStore
    bam_id = job.fileStore.writeGlobalFile(outpath_bam)
    bai_id = job.fileStore.writeGlobalFile(outpath_bai)
    return bam_id, bai_id


def gatk_base_recalibrator(job, bamFileStoreID, baiFileStoreID, config, xmx=8):
    """
    Creates recal table to perform Base Quality Score Recalibration

    :param job:
    :param bamFileStoreID:
    :param baiFileStoreID:
    :param config:
    :param xmx:
    :return:
    """
    job.fileStore.logToMaster('Running GATK BaseRecalibrator: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()
    # Retrieve input file paths
    references = ['genome.fa','genome.fa.fai', 'genome.dict', 'dbsnp.vcf']
    inputs = {'sample.indel.bam': bamFileStoreID, 'sample.indel.bam.bai': baiFileStoreID}
    inputs.update({key: config[key] for key in references})
    get_files_from_filestore(job, work_dir, inputs)
    # Output file path
    output = os.path.join(work_dir, 'sample.recal.table')

    # Call: GATK -- IndelRealigner
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'BaseRecalibrator',
                  '-nct', str(cores),
                  '-R', 'genome.fa',
                  '-I', 'sample.indel.bam',
                  '-knownSites', 'dbsnp.vcf',
                  '-o', 'sample.recal.table']
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                work_dir=work_dir, parameters=parameters,
                inputs=inputs.keys(),
                outputs={'sample.recal.table': None},
                env={'JAVA_OPTS':'-Xmx%sg' % xmx},
                mock=config['mock'])
    # Write to fileStore
    return job.fileStore.writeGlobalFile(output)


def gatk_print_reads(job, bamFileStoreID, baiFileStoreID, recalFileStoreID, config, xmx=8):
    """
    Create bam that has undergone Base Quality Score Recalibration (BQSR)

    :param job:
    :param bamFileStoreID:
    :param baiFileStoreID:
    :param recalFileStoreID:
    :param config:
    :param xmx:
    :return:

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    job.fileStore.logToMaster('Running GATK PrintReads: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()
    # Retrieve input file paths
    inputs = {'sample.indel.bam': bamFileStoreID, 'sample.indel.bam.bai': baiFileStoreID,
              'recal.table': recalFileStoreID}
    inputs.update({key: config[key] for key in ['genome.fa','genome.fa.fai', 'genome.dict']})
    get_files_from_filestore(job, work_dir, inputs)
    # Output file
    # Unpack convenience variables for job
    uuid = config['uuid']
    suffix = config['suffix'] if config['suffix'] is not None else ''
    outfile = '{}{}.bam'.format(uuid, suffix)
    gatk_outfile_idx = '{}{}.bai'.format(uuid, suffix)
    outpath_bam = os.path.join(work_dir, outfile)
    outpath_bai = os.path.join(work_dir, gatk_outfile_idx)
    # Call: GATK -- PrintReads
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'PrintReads',
                  '-nct', str(cores),
                  '-R', 'ref.fa',
                  '--emit_original_quals',
                  '-I', 'sample.indel.bam',
                  '-BQSR', 'sample.recal.table',
                  '-o', outfile]
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                work_dir=work_dir, parameters=parameters,
                inputs=inputs,
                outputs={outfile: None, gatk_outfile_idx: None},
                env={'JAVA_OPTS':'-Xmx%sg' % xmx},
                mock=config['mock'])
    
    upload_or_move(job, work_dir, config['output_dir'], outfile)
    outfile_base, bai = os.path.splitext(outfile)
    outfile_idx = outfile_base + '.bam' + bai
    upload_or_move(job, work_dir, config['output_dir'], outfile_idx)
    # Write to fileStore
    bam_id = job.fileStore.writeGlobalFile(outpath_bam)
    bai_id = job.fileStore.writeGlobalFile(outpath_bai)
    return bam_id, bai_id


def generate_config():
    return textwrap.dedent("""
        # GATK Preprocessing Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        genome:                   # Required: Reference Genome URL
        phase:                    # Required: URL (1000G_phase1.indels.hg19.sites.fixed.vcf)
        mills:                    # Required: URL (Mills_and_1000G_gold_standard.indels.hg19.sites.vcf)
        dbsnp:                    # Required: URL (dbsnp_132_b37.leftAligned.vcf URL)
        ssec:                     # Optional: (string) Path to Key File for SSE-C Encryption
        mock:                     # Optional: (bool) If True, run in mock mode
    """[1:])

def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 2 tab-separated columns: UUID, URL
        #
        #   UUID        This should be a unique identifier for the sample to be processed
        #   URL         A URL (http://, ftp://, file://, s3://) pointing to the input SAM or BAM file
        #
        #   Example below. Lines beginning with # are ignored.
        #
        #   UUID_1    file:///path/to/sample.bam
        #
        #   Place your samples below, one per line.
    """[1:])


def parse_manifest(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if not line.isspace() and not line.startswith('#'):
                sample = line.strip().split('\t')
                require(len(sample) == 2, 'Bad manifest format! '
                                          'Expected 2 tab separated columns, got: {}'.format(sample))
                uuid, url = sample
                require(urlparse(url).scheme and urlparse(url), 'Invalid URL passed for {}'.format(url))
                samples.append(sample)
    return samples


def gatk_preprocessing_pipeline(job, uuid, url, config):
    config = deepcopy(config)
    config['uuid'] = uuid
    sample_name = os.path.basename(url)
    download = job.wrapJobFn(download_url_job, url, name=sample_name, s3_key_path=config['ssec'])
    rm_secondary = job.wrapJobFn(samtools_view, download.rv(), '0x800', mock=config['mock'])
    picard_sort = job.wrapJobFn(picard_sort_sam, rm_secondary.rv(), config)
    mdups = job.wrapJobFn(picard_mark_duplicates, picard_sort.rv(0), picard_sort.rv(1), config)
    target = job.wrapJobFn(gatk_realigner_target_creator, mdups.rv(0), mdups.rv(1), config)
    indel = job.wrapJobFn(gatk_indel_realigner, mdups.rv(0), mdups.rv(1), target.rv(), config)
    base = job.wrapJobFn(gatk_base_recalibrator, indel.rv(0), indel.rv(1), config)
    print_reads = job.wrapJobFn(gatk_print_reads, indel.rv(0), indel.rv(1), base.rv(), config)

    job.addChild(download)
    download.addChild(rm_secondary)
    rm_secondary.addChild(picard_sort)
    picard_sort.addChild(mdups)
    mdups.addChild(target)
    target.addChild(indel)
    indel.addChild(base)
    base.addChild(print_reads)


def main():
    """
    GATK Pre-processing Script
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')


    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the CGL germline pipeline')
    parser_run.add_argument('--config', default='config-toil-germline.yaml', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--manifest', default='manifest-toil-germline.tsv', type=str,
                            help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--bam', default=None, type=str,
                            help='URL for the sample BAM. URLs can take the form: http://, file://, s3://, '
                                 'and gnos://. The UUID for the sample must be given with the "--uuid" flag.')
    parser_run.add_argument('--uuid', default=None, type=str, help='Provide the UUID of a sample when using the'
                                                                   '"--bam" option')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'config-toil-germline.yaml'), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-toil-germline.tsv'), generate_manifest)
    if 'generate' in args.command:
        sys.exit()
    if args.command == 'run':
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program)), program + ' must be installed on every node.'.format(program))

        config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}

        require(set(config) > {'genome', 'mills', 'dbsnp', 'phase'},
                'Missing inputs for preprocessing, check config file')

        if args.bam or args.uuid:
            require(args.bam and args.uuid, '"--bam" and "--uuid" must all be supplied')
            samples = [[args.uuid, args.normal, args.tumor]]
        else:
            samples = parse_manifest(args.manifest)

        root = Job.wrapJobFn(download_shared_files, config)
        ref = Job.wrapJobFn(reference_preprocessing, root.rv(), mock=config['mock'])

        root.addFollowOn(ref)

        for uuid, url in samples:
            ref.addFollowOnJobFn(gatk_preprocessing_pipeline, uuid, url, ref.rv())

        Job.Runner.startToil(root, args)

if __name__ == '__main__':
    main()
