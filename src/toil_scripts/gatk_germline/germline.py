#!/usr/bin/env python2.7

from __future__ import print_function
import argparse
import multiprocessing
import os
import textwrap
import yaml
from copy import deepcopy
from urlparse import urlparse
from bd2k.util.processes import which

from toil.job import Job
from toil_scripts.lib.urls import download_url_job
from toil_scripts.lib import require
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.files import upload_or_move_job
from toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline import generate_file
from toil_scripts.gatk_processing.gatk_preprocessing import samtools_faidx, picard_create_sequence_dictionary


def gatk_germline_pipeline(job, uuid, url, config, bai_url=None):
    """
    Runs GATK Germline Pipeline on a single sample. Writes gvcf and vqsr.vcf to output directory.

    :param job: Toil Job instance
    :param uuid str: Unique identifier for the sample
    :param url str: URL to sample bam
    :param config dict: Configuration options for pipeline
    """
    config = deepcopy(config)
    config['uuid'] = uuid

    if config ['xmx'] is None:
        config['xmx'] = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024 ** 3)
    download_bam = job.wrapJobFn(download_url_job, url, name='toil.bam', s3_key_path=config['ssec'])

    if bai_url:
        bam_index = job.wrapJobFn(download_url_job, bai_url, name='toil.bam.bai', s3_key_path=config['ssec'])
        download_bam.addChild(bam_index)
    else:
        bam_index = job.wrapJobFn(samtools_index, download_bam.rv(), config)
        download_bam.addChild(bam_index)

    haplotype_caller = job.wrapJobFn(gatk_haplotype_caller, download_bam.rv(), bam_index.rv(), config)
    genotype_gvcf = job.wrapJobFn(gatk_genotype_gvcf, haplotype_caller.rv(), config)

    snp_recal = job.wrapJobFn(gatk_variant_recalibrator_snp, genotype_gvcf.rv(), config)
    apply_snp_recal = job.wrapJobFn(gatk_apply_variant_recalibration_snp, genotype_gvcf.rv(),
                                    snp_recal.rv(0), snp_recal.rv(1), config)

    indel_recal = job.wrapJobFn(gatk_variant_recalibrator_indel, genotype_gvcf.rv(), config)
    apply_indel_recal = job.wrapJobFn(gatk_apply_variant_recalibration_indel, apply_snp_recal.rv(),
                                      indel_recal.rv(0), indel_recal.rv(1), config)
    # Output vcfs
    raw_name = '{}.raw{}.gvcf'.format(uuid, config['suffix'])
    vqsr_name = '{}.vqsr{}.vcf'.format(uuid, config['suffix'])
    output_raw = job.wrapJobFn(upload_or_move_job, raw_name, genotype_gvcf.rv(), config['output_dir'])
    output_vqsr = job.wrapJobFn(upload_or_move_job, vqsr_name, apply_indel_recal.rv(), config['output_dir'])

    # Create DAG
    job.addChild(download_bam)
    # Wait for bam index
    download_bam.addFollowOn(haplotype_caller)
    haplotype_caller.addChild(genotype_gvcf)
    genotype_gvcf.addChild(snp_recal)
    genotype_gvcf.addChild(indel_recal)
    genotype_gvcf.addChild(output_raw)
    genotype_gvcf.addFollowOn(apply_snp_recal)
    apply_snp_recal.addChild(apply_indel_recal)
    apply_indel_recal.addChild(output_vqsr)


def get_files_from_filestore(job, work_dir, input_dict):
    """
    Copies files from the file store to a work directory

    :param job: Toil Job instance
    :param work_dir: current working directory
    :param input_dict: {filename: fileStoreID}
    :return dict: {filename: filepath}
    """
    for name, fileStoreID in input_dict.iteritems():
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = job.fileStore.readGlobalFile(fileStoreID, os.path.join(work_dir, name))
        else:
            file_path = name
        input_dict[name] = file_path
    return input_dict


def download_shared_files(job, config):
    """
    Downloads reference files shared by all samples in the pipeline

    :param config dict: pipeline configuration options
    :return dict: updated config dictionary with shared fileStoreIDS
    """
    job.fileStore.logToMaster('Downloading shared reference files')
    reference_names = ['genome.fa', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf', 'hapmap.vcf', 'omni.vcf']
    for name in reference_names:
        key, _ = name.split('.')
        config[name] = job.addChildJobFn(download_url_job, config[key], name=name,
                                         s3_key_path=config['ssec']).rv()
    return config


def reference_preprocessing(job, config, mock=False):
    """
    Creates a genome.fa.fai and genome.dict file for reference genome.fa
    if one is not already present in the config dictionary

    :param config dict: pipeline configuration options and shared files
                        requires a key called 'genome.fa'
    :return dict: updated config dictionary with reference index files
    """
    job.fileStore.logToMaster('Preparing Reference Files')
    genome_id = config['genome.fa']
    if 'genome.fa.fai' not in config:
        config['genome.fa.fai'] = job.addChildJobFn(samtools_faidx, genome_id, mock=mock).rv()
    if 'genome.dict' not in config:
        config['genome.dict'] = job.addChildJobFn(picard_create_sequence_dictionary, genome_id, mock=mock).rv()
    return config


def samtools_index(job, bam_id, config):
    """
    Index sample bam using samtools index

    :param job: Job instance
    :param bam_id str: bam file store ID
    :param config dict: pipeline configuration options and shared files
    :return:
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    inputs = {'toil.bam': bam_id}
    outputs={'toil.bam.bai': None}
    get_files_from_filestore(job, work_dir, inputs)
    # Call: index the normal.bam
    parameters = ['index', 'toil.bam']
    docker_call(work_dir = work_dir,
                parameters = parameters,
                tool = 'quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])
    output_path = os.path.join(work_dir, 'toil.bam.bai')
    return job.fileStore.writeGlobalFile(output_path)


def gatk_haplotype_caller(job, bam_id, bai_id, config):
    """
    Uses GATK HaplotypeCaller to identify SNPs and Indels and writes a gVCF.

    :param job: Job instance
    :param bam_id str: bam file store ID
    :param bai_id str: bai file store ID
    :param config dict: pipeline configuration options and shared files
    :return str: gvcf file store ID
    """
    job.fileStore.logToMaster('Running GATK HaplotypeCaller: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    inputs['toil.bam'] = bam_id
    inputs['toil.bam.bai'] = bai_id
    get_files_from_filestore(job, work_dir, inputs)

    # Call GATK -- HaplotypeCaller
    command = ['-nct', str(cores),
               '-R', 'genome.fa',
               '-T', 'HaplotypeCaller',
               '--genotyping_mode', 'Discovery',
               '--emitRefConfidence', 'GVCF',
               '-I', 'toil.bam',
               '-o', 'toil.gvcf',
               '-variant_index_type', 'LINEAR',
               '-variant_index_parameter', '128000',
               '--annotation', 'QualByDepth',
               '--annotation', 'DepthPerSampleHC',
               '--annotation', 'FisherStrand',
               '--annotation', 'ReadPosRankSumTest']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs={'toil.gvcf': None}
    docker_call(work_dir = work_dir,
                env={'JAVA_OPTS':'-Xmx%sg' % config['xmx']},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])

    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'toil.gvcf'))


def gatk_genotype_gvcf(job, gvcf_id, config):
    """
    Genotypes the gVCF generated by HaplotypeCaller.

    :param job: Job instance
    :param gvcf_id str: file store ID for gvcf
    :param config dict: pipeline configuration options and shared files
    :return str: vcf file store ID
    """
    job.fileStore.logToMaster('Running GATK GenotypeGVCFs: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    inputs['toil.vcf'] = gvcf_id
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-nt', str(cores),
               '-R', 'genome.fa',
               '-T', 'GenotypeGVCFs',
               '--variant', 'toil.vcf',
               '--out', 'toil.vcf',
               '-stand_emit_conf', '10.0',
               '-stand_call_conf', '30.0']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs = {'toil.vcf': None}
    docker_call(work_dir = work_dir,
                env={'JAVA_OPTS':'-Xmx%sg' % config['xmx']},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])

    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'toil.vcf'))


def gatk_variant_recalibrator_snp(job, vcf_id, config):
    """
    Variant quality score recalibration for SNP variants.

    :param job: Job instance
    :param vcf_id str: vcf file store ID
    :param config dict: pipeline configuration options and shared files
    :return tuple: recalibration table, tranches file, plots file
    """
    job.fileStore.logToMaster('Running GATK VariantRecalibrator (SNP Mode): {}'.format(config['uuid']))
    cores = multiprocessing.cpu_count()
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict', 'hapmap.vcf', 'omni.vcf', 'dbsnp.vcf', 'phase.vcf']
    inputs = {key: config[key] for key in references}
    inputs['toil.vcf'] = vcf_id
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-T', 'VariantRecalibrator',
               '-R', 'genome.fa',
               '-input', 'toil.vcf',
               '-nt', str(cores),
               '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', 'hapmap.vcf',
               '-resource:omni,known=false,training=true,truth=false,prior=12.0', 'omni.vcf',
               '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0', 'dbsnp.vcf',
               '-resource:1000G,known=false,training=true,truth=false,prior=10.0', 'phase.vcf',
               '-an', 'QD', '-an', 'DP', '-an', 'FS', '-an', 'ReadPosRankSum',
               '-mode', 'SNP', '-minNumBad', '1000',
               '-recalFile', 'HAPSNP.recal',
               '-tranchesFile', 'HAPSNP.tranches',
               '-rscriptFile', 'HAPSNP.plots']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs = {'HAPSNP.recal': None, 'HAPSNP.tranches': None, 'HAPSNP.plots': None}
    docker_call(work_dir = work_dir,
                env={'JAVA_OPTS':'-Xmx%sg' % config['xmx']},
                parameters = command,
                tool ='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])

    recal_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPSNP.recal'))
    tranches_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPSNP.tranches'))
    plots_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPSNP.plots'))
    return recal_id, tranches_id, plots_id


def gatk_apply_variant_recalibration_snp(job, vcf_id, recal_id, tranches_id, config):
    """
    Apply variant quality score recalibration for SNP variants.

    :param job: Job instance
    :param vcf_id str: vcf file store ID
    :param recal_id: recalibration table file store ID
    :param tranches_id: tranches file store ID
    :param config dict: pipeline configuration options and shared files
    :return str: SNP recalibrated VCF file store ID
    """
    job.fileStore.logToMaster('Running GATK ApplyRecalibration (SNP Mode): {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {'toil.vcf': vcf_id,
              'HAPSNP.tranches': tranches_id,
              'HAPSNP.recal': recal_id}
    inputs.update({key: config[key] for key in references})
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-T', 'ApplyRecalibration',
               '-input', 'toil.vcf',
               '-o', 'toil.vqsr.vcf',
               '-R', 'genome.fa',
               '-nt', '1',
               '-ts_filter_level', '99.0',
               '-tranchesFile', 'HAPSNP.tranches',
               '-recalFile', 'HAPSNP.recal',
               '-mode', 'SNP']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs={'toil.vqsr.vcf': None}
    docker_call(work_dir = work_dir,
                env={'JAVA_OPTS':'-Xmx%sg' % config['xmx']},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'toil.vqsr.vcf'))


def gatk_variant_recalibrator_indel(job, vcf_id, config):
    """
    Variant quality score recalibration for INDEL variants.

    :param job: Job instance
    :param vcf_id str: vcf file store ID
    :param config dict: pipeline configuration options and shared files
    :return tuple: recalibration table, tranches file, plots file
    """
    job.fileStore.logToMaster('Running GATK VariantRecalibrator (INDEL Mode): {}'.format(config['uuid']))
    cores = multiprocessing.cpu_count()
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict', 'mills.vcf']
    inputs = {key: config[key] for key in references}
    inputs['toil.vcf'] = vcf_id
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-T', 'VariantRecalibrator',
               '-R', 'genome.fa',
               '-input', 'toil.vcf',
               '-nt', str(cores),
               '-resource:mills,known=true,training=true,truth=true,prior=12.0', 'mills.vcf',
               '-an', 'DP', '-an', 'FS', '-an', 'ReadPosRankSum',
               '-mode', 'INDEL',
               '-minNumBad', '1000',
               '-recalFile', 'HAPINDEL.recal',
               '-tranchesFile', 'HAPINDEL.tranches',
               '-rscriptFile', 'HAPINDEL.plots',
               '--maxGaussians', '4']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs = {'HAPINDEL.recal': None, 'HAPINDEL.tranches': None, 'HAPINDEL.plots': None}
    docker_call(work_dir = work_dir,
                env={'JAVA_OPTS':'-Xmx%sg' % config['xmx']},
                parameters = command,
                tool ='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])
    recal_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPINDEL.recal'))
    tranches_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPINDEL.tranches'))
    plots_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPINDEL.plots'))
    return recal_id, tranches_id, plots_id


def gatk_apply_variant_recalibration_indel(job, vcf_id, recal_id, tranches_id, config):
    """
    Apply variant quality score recalibration for indel variants.

    :param job: Job instance
    :param vcf_id str: vcf file store ID
    :param recal_id: recalibration table file store ID
    :param tranches_id: tranches file store ID
    :param config dict: pipeline configuration options and shared files
    :return str: INDEL recalibrated VCF file store ID
    """
    job.fileStore.logToMaster('Running GATK ApplyRecalibration (INDEL Mode): {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    uuid = config['uuid']
    suffix = config['suffix']

    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {'toil.vcf': vcf_id,
              'HAPINDEL.recal': recal_id,
              'HAPINDEL.tranches': tranches_id}
    inputs.update({key: config[key] for key in references})
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-T', 'ApplyRecalibration',
               '-input', 'unified.raw.BOTH.gatk.vcf',
               '-o', 'toil.vqsr.vcf',
               '-R', 'genome.fa',
               '-nt', '1',
               '-ts_filter_level', '99.0',
               '-tranchesFile', 'HAPINDEL.tranches',
               '-recalFile', 'HAPINDEL.recal',
               '-mode', 'INDEL']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs={'toil.vqsr.vcf': None}
    docker_call(work_dir = work_dir,
                env={'JAVA_OPTS':'-Xmx%sg' % config['xmx']},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs = inputs.keys(),
                outputs = outputs,
                mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'toil.vqsr.vcf'))


def generate_config():
    return textwrap.dedent("""
        # GATK Germline Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        # Required: Reference Genome URL
        genome:
        # Required: URL (1000G_phase1.indels.hg19.sites.fixed.vcf)
        phase:
        # Required: URL (Mills_and_1000G_gold_standard.indels.hg19.sites.vcf)
        mills:
        # Required: URL (dbsnp_132_b37.leftAligned.vcf URL)
        dbsnp:
        # Required: URL (hapmap_3.3.b37.vcf)
        hapmap:
        # Required: URL (1000G_omni.5.b37.vcf)
        omni:
        # Approximate input file size. Should be given as %d[TGMK], e.g.,
        # for a 100 gigabyte file, use file-size: 100G
        file-size: 100G
        # Memory allocation for Java option Xmx
        xmx: 100G
        # Set to True if the input file is already indexed and a .bam.bai file is
        # present at the same url as the .bam file
        indexed': False
        # Optional: (string) Path to Key File for SSE-C Encryption
        ssec:
        # Optional: (string) Suffix to be added to final output
        suffix:
        # Optional: (string) Path to output directory
        output_dir:
        # Optional: (boolean) Set to True to allow seq dict incompatibility
        unsafe_mode:
        # Optional: (boolean) Set to True to run pipeline in mock mode
        mock_mode:
    """[1:])

def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 2 tab-separated columns: UUID and URL
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
                sample = line.strip().split()
                require(len(sample) == 2, 'Bad manifest format! '
                                          'Expected 2 tab separated columns, got: {}'.format(sample))
                uuid, url = sample
                require(urlparse(url).scheme and urlparse(url), 'Invalid URL passed for {}'.format(url))
                samples.append(sample)
    return samples


def main():
    """
    GATK HaplotypeCaller with variant recalibration.
    Writes gvcf and vqsr.vcf to output directory

                Tree Structure of GATK Pipeline
                0 --> 1 --> 2 --> 3 --> 4 --> 5 ----> 8
                                             / \      |
                                            6   7     9
    0 = Start Node
    1 = Download References
    2 = Prepare References
    4 = Download Sample
    5 = HaplotypeCaller SNP & Indel
    6 = VariantRecalibrator SNPs
    7 = VariantRecalibrator Indels
    8 = ApplyRecalibration SNPs
    9 = ApplyRecalibration Indels

    ===================================================================
    :Dependencies:
    curl            - apt-get install curl
    docker          - apt-get install docker (or 'docker.io' for linux)
    toil            - pip install --pre toil
    """
    # Define Parser object and add to jobTree
    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    # Generate subparsers
    subparsers = parser.add_subparsers(dest='command')
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the GATK germline pipeline')
    parser_run.add_argument('--config', default='config-toil-germline.yaml', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config".')
    parser_run.add_argument('--manifest', default='manifest-toil-germline.tsv', type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s".')
    parser_run.add_argument('--sample', default=None, nargs=2, type=str,
                       help='Space delimited sample UUID and BAM file in the format: uuid url')
    parser_run.add_argument('--output-dir', default=None, help='Full path to directory or filename where '
                                                               'final results will be output')
    parser_run.add_argument('-s', '--suffix', default=None,
                            help='Additional suffix to add to the names of the output files')
    Job.Runner.addToilOptions(parser_run)
    options = parser.parse_args()

    cwd = os.getcwd()
    if options.command == 'generate-config' or options.command == 'generate':
        generate_file(os.path.join(cwd, 'config-toil-germline.yaml'), generate_config)
    if options.command == 'generate-manifest' or options.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-toil-germline.tsv'), generate_manifest)
    # Pipeline execution
    elif options.command == 'run':
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program)), program + ' must be installed on every node.'.format(program))

        require(os.path.exists(options.config), '{} not found. Please run '
                                                '"generate-config"'.format(options.config))
        samples = []
        if options.sample:
            uuid, bam = options.sample.split()
            samples.append((uuid, bam))
        elif options.manifest:
            samples.extend(parse_manifest(options.manifest))

        # Parse config
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(options.config).read()).iteritems()}
        config_fields = set(config)
        required_fields = {'genome', 'mills', 'dbsnp', 'phase', 'hapmap', 'omni',
                           'output_dir', 'mock_mode', 'unsafe_mode', 'xmx', 'suffix'}
        require(config_fields > required_fields,
                'Missing config parameters:\n{}'.format(required_fields - config_fields))

        if config['output_dir'] is  None:
            config['output_dir'] = options.output_dir if options.output_dir else os.getcwd()

        if config['suffix'] is  None:
            config['suffix'] = options.suffix if options.suffix else ''

        root = Job.wrapJobFn(download_shared_files, config)
        ref = Job.wrapJobFn(reference_preprocessing, root.rv(), mock=config['mock_mode'])

        root.addFollowOn(ref)

        for uuid, url in samples:
            ref.addFollowOnJobFn(gatk_germline_pipeline, uuid, url, ref.rv())

        Job.Runner.startToil(root, options)

if __name__ == '__main__':
    main()
