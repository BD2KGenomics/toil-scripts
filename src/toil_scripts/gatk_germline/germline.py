#!/usr/bin/env python2.7
from __future__ import print_function
import argparse
from copy import deepcopy
import os
import re
from urlparse import urlparse

from bd2k.util.humanize import human2bytes
from bd2k.util.processes import which
import synapseclient
from toil.job import Job, PromisedRequirement
import yaml

from toil_scripts.gatk_germline.germline_config import generate_config, generate_manifest
from toil_scripts.gatk_germline.hard_filter import hard_filter_pipeline
from toil_scripts.gatk_germline.vqsr import vqsr_pipeline
from toil_scripts.lib import require
from toil_scripts.lib.files import get_files_from_filestore, upload_or_move_job
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url_job
from toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline import generate_file
from toil_scripts.tools.aligners import run_bwakit
from toil_scripts.tools.indexing import run_samtools_faidx
from toil_scripts.tools.preprocessing import run_gatk_preprocessing, run_samtools_sort, \
    run_samtools_index, run_samtools_view, run_picard_create_sequence_dictionary
from toil_scripts.tools.variant_annotation import run_oncotator
from toil_scripts.tools.variant_filters import split_vcf_by_name


def parse_manifest(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param path_to_manifest str: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    bam_re = r"^(?P<uuid>\S+)\s(?P<url>\S+[bsc][r]?am)"
    fq_re = r"^(?P<uuid>\S+)\s(?P<url>\S+)\s(?P<url2>\S+)?\s?(?P<rg_line>@RG\S+)"
    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith('#'):
                continue
            bam_match = re.match(bam_re, line)
            fastq_match = re.match(fq_re, line)
            if bam_match:
                uuid = bam_match.group('uuid')
                url = bam_match.group('url')
                url2 = None
                rg_line = None
                require('.bam' in url.lower() or '.fastq' in url.lower(),
                        'Expected .bam extension:\n{}:\t{}'.format(uuid, url))
            elif fastq_match:
                uuid = fastq_match.group('uuid')
                url = fastq_match.group('url')
                url2 = fastq_match.group('url2')
                rg_line = fastq_match.group('rg_line')
                require('.fq' in url.lower() or '.fastq' in url.lower(),
                        'Expected .fq extension:\n{}:\t{}'.format(uuid, url))
            else:
                msg = 'Could not parse entry in manifest: %s\n%s' % (f.path, line)
                raise ValueError(msg)
            require(urlparse(url).scheme or url.startswith('syn'),
                    'Invalid URL passed for {}'.format(url))
            samples.append((uuid, url, url2, rg_line))
    return samples


def convert_space_resource(attribute, default):
    """
    Converts resource requirement from human readable bytes to raw integer value.

    :param str|int|None attribute: Space resource value
    :param str default: Human readable parameter size
    :return: bytes
    :rtype: int
    """
    if attribute is None:
        return human2bytes(default)
    elif isinstance(attribute, str) and re.match('\d+\s*[MKG]', attribute.upper()):
        return human2bytes(attribute)
    elif isinstance(attribute, int):
        return attribute
    else:
        raise ValueError('Could not convert resource requirement: %s' % attribute)


def download_shared_files(job, config):
    """
    Downloads reference files shared by all samples in the Toil Germline pipeline

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Pipeline configuration options
    :return: Updated config with shared fileStoreIDS
    :rtype: Namespace
    """
    job.fileStore.logToMaster('Downloading shared reference files')
    shared_files = {'genome_fasta', 'genome_fai', 'genome_dict'}
    nonessential_files = {'genome_fai', 'genome_dict'}

    # Corrects naming convention in other Toil scripts
    if getattr(config, 'ref', None): config.genome_fasta = config.ref
    if getattr(config, 'fai', None): config.genome_fai = config.fai
    if getattr(config, 'dict', None): config.genome_dict = config.dict

    # Download necessary files for pipeline configuration
    if config.run_bwa:
        shared_files |= {'amb', 'ann', 'bwt', 'pac', 'sa', 'alt'}
        nonessential_files.add('alt')
    if config.preprocess:
        shared_files |= {'phase', 'mills', 'dbsnp'}
    if config.run_vqsr:
        shared_files |= {'phase', 'mills', 'dbsnp', 'hapmap', 'omni'}
    if config.run_oncotator:
        shared_files.add('oncotator_db')
    for name in shared_files:
        try:
            url = getattr(config, name, None)
            if url is None:
                continue
            setattr(config, name, job.addChildJobFn(download_url_job,
                                                    url,
                                                    name=name,
                                                    synapse_login=config.synapse_login,
                                                    s3_key_path=config.ssec,
                                                    disk='5G').rv())
        finally:
            if getattr(config, name, None) is None and name not in nonessential_files:
                raise ValueError("Necessary configuration parameter is missing:\n{}".format(name))
    return job.addFollowOnJobFn(reference_preprocessing, config).rv()


def reference_preprocessing(job, config):
    """
    Creates a genome.fa.fai and genome.dict file for reference genome.fa
    if one is not already present in the config dictionary

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Pipeline configuration options and shared files.
                             Requires genome_fasta attribute
    :return: Updated config with reference index files
    :rtype: Namespace
    """
    job.fileStore.logToMaster('Preparing Reference Files')
    genome_id = config.genome_fasta
    if getattr(config, 'genome_fai', None) is None:
        config.genome_fai = job.addChildJobFn(run_samtools_faidx,
                                              genome_id,
                                              cores=config.cores).rv()
    if getattr(config, 'genome_dict', None) is None:
        config.genome_dict = job.addChildJobFn(run_picard_create_sequence_dictionary,
                                               genome_id,
                                               cores=config.cores).rv()
    return config


def prepare_bam(job, uuid, url, url2, config, rg_line=None):
    """
    Prepares BAM file for GATK Germline pipeline.

    0: Align FASTQ or Download BAM
    1: Remove secondary alignments
    2: Sort
    3: Index

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: Unique identifier for the sample
    :param str url: URL to BAM file or FASTQs
    :param Namespace config: Configuration options for pipeline
    :param str rg_line: RG line for BWA alignment. Default is None.
    :return: BAM and BAI FileStoreIDs
    :rtype: tuple
    """
    # 0: Align FASQ or realign BAM
    if config.run_bwa:
        get_bam = job.wrapJobFn(setup_and_run_bwa_kit, url, url2, rg_line, config).encapsulate()

    # 0: Download BAM
    elif '.bam' in url.lower():
        job.fileStore.logToMaster("Downloading BAM: %s" % uuid)
        get_bam = job.wrapJobFn(download_url_job,
                                url,
                                name='toil.bam',
                                synapse_login=config.synapse_login,
                                s3_key_path=config.ssec,
                                disk=config.file_size).encapsulate()
    else:
        raise ValueError('Could not generate BAM file for %s\n'
                         'Either provide a FASTQ URL and set run-bwa or provide a BAM URL.' % uuid)

    # 1: Remove secondary alignments
    rm_secondary = job.wrapJobFn(run_samtools_view,
                                 get_bam.rv(),
                                 flag='0x800',
                                 ncores=config.cores,
                                 disk=2*config.file_size, cores=config.cores)

    # 2: Sort BAM file if necessary
    if config.sorted:
        sorted_bam = rm_secondary

    else:
        sorted_bam_disk = PromisedRequirement(lambda x: 3 * x.size + human2bytes('10G'),
                                              rm_secondary.rv())
        sorted_bam = job.wrapJobFn(run_samtools_sort,
                                   rm_secondary.rv(),
                                   ncores=config.cores,
                                   cores=config.cores,
                                   disk=sorted_bam_disk)

    # 3: Index BAM
    index_bam_disk = PromisedRequirement(lambda x: x.size, sorted_bam.rv())
    index_bam = job.wrapJobFn(run_samtools_index, sorted_bam.rv(), disk=index_bam_disk)

    # Wiring
    job.addChild(get_bam)
    get_bam.addChild(rm_secondary)
    rm_secondary.addChild(sorted_bam)
    sorted_bam.addChild(index_bam)
    return sorted_bam.rv(), index_bam.rv()


def setup_and_run_bwa_kit(job, url, url2, rg_line, config):
    """
    Downloads and aligns FASTQs or realigns BAM using BWA.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str url: FASTQ or BAM sample URL. Must have .fq or .bam extension.
    :param Namespace config: Input parameters and shared FileStoreIDs
    :return: BAM FileStoreID
    :rtype: str
    """
    bwa_config = deepcopy(config)
    bwa_config.rg_line = rg_line

    # bwa_alignment uses different naming convention,
    bwa_config.ref = config.genome_fasta
    bwa_config.fai = config.genome_fai

    basename, ext = os.path.splitext(url)
    ext = ext.lower()
    if ext == '.gz':
        _, ext = os.path.splitext(basename)
        ext = ext.lower()

    require(ext in ['.fq', '.fastq', '.bam'],
            'Please use .fq or .bam file extensions')

    # Download fastq files
    samples = []
    input1 = job.addChildJobFn(download_url_job,
                               url,
                               name='file1',
                               synapse_login=config.synapse_login,
                               s3_key_path=config.ssec,
                               disk=config.file_size)
    samples.append(input1.rv())
    if '.bam' in url.lower():
        bwa_config.bam = input1.rv()
    else:
        bwa_config.r1 = input1.rv()

    # Uses first url to deduce paired FASTQ URL, if url looks like a paired file.
    if url2 is None and '1.fq' in url:
        url2 = url.replace('1.fq', '2.fq')
        job.fileStore.logToMaster(
            'Trying to find paired reads using FASTQ URL:\n%s\n%s' % (url, url2))

    if url2:
        input2 = job.addChildJobFn(download_url_job,
                                   url2,
                                   name='file2',
                                   synapse_login=config.synapse_login,
                                   s3_key_path=config.ssec,
                                   disk=config.file_size)
        samples.append(input2.rv())
        bwa_config.r2 = input2.rv()

    # Disk requirement depends on whether there are paired reads
    bwakit_disk = PromisedRequirement(lambda lst: 3 * sum(x.size for x in lst) + human2bytes('8G'),
                                      samples)
    return job.addFollowOnJobFn(run_bwakit,
                                bwa_config,
                                config.cores,
                                sort=False,         # BAM files are sorted later in the pipeline
                                trim=config.trim,
                                cores=config.cores,
                                disk=bwakit_disk).rv()


def gatk_haplotype_caller(job, bam_id, bai_id, config, emit_threshold=10.0, call_threshold=30.0):
    """
    Uses GATK HaplotypeCaller to identify SNPs and INDELs. Outputs variants in a Genomic VCF file.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str bam_id: BAM FileStoreID
    :param str bai_id: BAI FileStoreID
    :param Namespace config: pipeline configuration options and shared files
    :param float emit_threshold: Minimum phred-scale confidence threshold for a variant to be emitted
                               Default: 10
    :param float call_threshold: Minimum phred-scale confidence threshold for a variant to be called
                               Default: 30
    :return: GVCF FileStoreID
    :rtype: str
    """
    job.fileStore.logToMaster('Running GATK HaplotypeCaller: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': config.genome_fasta,
              'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict,
              'input.bam': bam_id, 'input.bam.bai': bai_id}
    get_files_from_filestore(job, work_dir, inputs)

    # Call GATK -- HaplotypeCaller with best practices parameters:
    # https://software.broadinstitute.org/gatk/documentation/article?id=2803
    command = ['-T', 'HaplotypeCaller',
               '-nct', str(config.cores),
               '-R', 'genome.fa',
               '-I', 'input.bam',
               '-o', 'output.gvcf',
               '-stand_call_conf', str(call_threshold),
               '-stand_emit_conf', str(emit_threshold),
               '-variant_index_type', 'LINEAR',
               '-variant_index_parameter', '128000',
               '--genotyping_mode', 'Discovery',
               '--emitRefConfidence', 'GVCF']

    if config.unsafe_mode:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    for annotation in config.annotations:
        if annotation is not None:
            command.extend(['-A', annotation])

    outputs = {'output.gvcf': None}
    docker_call(work_dir=work_dir,
                env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(config.xmx)},
                parameters=command,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.gvcf'))


def annotate_vcfs(job, vcfs, config):
    """
    Annotates VCF files using Oncotator.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param dict vcfs: Dictionary of UUIDs and VCF FileStoreIDs
    :param Namespace config: Input parameters and shared FileStoreIDs
    :return: None
    """
    job.fileStore.logToMaster('Annotating VCFs')
    for uuid, vcf_id in vcfs.iteritems():
        annotated_vcf = job.addChildJobFn(run_oncotator, vcf_id, config)

        output_dir = os.path.join(config.output_dir, uuid)
        filename = '{}.annotated{}.vcf'.format(uuid, config.suffix)
        annotated_vcf.addChildJobFn(upload_or_move_job,
                                    filename,
                                    annotated_vcf.rv(),
                                    output_dir)


def download_and_run_sample(job, uuid, url, config, url2=None, rg_line=None):
    """
    Downloads reference files and runs the GATK best practices germline pipeline for one sample.
    See `toil_scripts.gatk_germline.germline_config.py` file for necessary parameters.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: Unique identifier for the sample
    :param str url: URL to BAM or FASTQ file
    :param Namespace config: Configuration options for pipeline
    :param str|None url2: URL to paired FASTQ file
    :param str|None rg_line: RG line for BWA alignment. Default is None.
    :return: None
    """
    get_shared_files = job.wrapJobFn(download_shared_files, config).encapsulate()
    run_pipeline = Job.wrapJobFn(run_gatk_germline_pipeline,
                                 [(uuid, url, url2, rg_line)],
                                 get_shared_files.rv())
    get_shared_files.addChild(run_pipeline)


def run_gatk_germline_pipeline(job, samples, config):
    """
    Configures and runs the GATK germline pipeline and filters the output. VQSR is done if the
    run_vqsr config parameter is set to True or there are at least 30 samples. If joint genotyping
    is done, the VCF is split per sample after filtering. Otherwise, each sample VCF is hard
    filtered using recommended GATK parameters.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list samples: List of tuples (uuid, url, url2, rg_line) for joint genotyping.
    :param Namespace config: Input parameters and reference FileStoreIDs
    :return: Dictionary of filtered VCF FileStoreIDs
    :rtype: dict
    """
    if len(samples) == 0:
        raise ValueError('No samples were provided!')

    # Generate per sample gvcfs {uuid: gvcf_id}
    gvcfs = {}
    for uuid, url, url2, rg_line in samples:
        gvcfs[uuid] = job.addChildJobFn(gatk_germline_pipeline,
                                        uuid,
                                        url,
                                        config,
                                        url2=url2,
                                        rg_line=rg_line).rv()

    # Stop after preprocessing
    if config.preprocess_only:
        return

    # VQSR requires many variants in order to train a decent model. GATK recommends a minimum of
    # 30 exomes or one large WGS sample:
    # https://software.broadinstitute.org/gatk/documentation/article?id=3225
    filtered_vcfs = {}
    if config.run_vqsr or len(gvcfs) > 30:
        if config.joint:
            joint_vcf = job.addFollowOnJobFn(vqsr_pipeline, gvcfs, config)
            filtered_vcfs = joint_vcf.addChildJobFn(split_vcf_by_name, gvcfs.keys(),
                                                    joint_vcf.rv(), config).rv()
        else:
            for uuid, gvcf in gvcfs.iteritems():
                filtered_vcfs[uuid] = job.addFollowOnJobFn(vqsr_pipeline,
                                                           dict(uuid=gvcf),
                                                           config).rv()
    else:
        for uuid, gvcf in gvcfs.iteritems():
            filtered_vcfs[uuid] = job.addFollowOnJobFn(hard_filter_pipeline,
                                                       uuid, gvcf,
                                                       config).rv()
    return filtered_vcfs


def gatk_germline_pipeline(job, uuid, url, config, url2=None, rg_line=None):
    """
    Runs GATK Germline Pipeline on a single sample. Writes gvcf to output directory
    defined in the config dictionary. Removes secondary alignments.

    0: Align Input or Download BAM file
    1: Generate BAM Index
    2: GATK Preprocessing (Optional)
    3: GATK HaplotypeCaller
    4: Write GVCF to output directory

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: Unique identifier for the sample
    :param str url: URL to sample bam or fq file
    :param Namespace config: Configuration options for pipeline
    :param str url2: Paired FASTQ URL. Default is None.
    :param str rg_line: RG line for BWA alignment. Default is None.
    :return: GVCF FileStoreID
    :rtype: str
    """
    config = deepcopy(config)
    config.uuid = uuid
    config.url = url

    get_bam = job.wrapJobFn(prepare_bam, uuid, url, url2, config, rg_line=rg_line).encapsulate()
    job.addChild(get_bam)
    input_bam = get_bam.rv(0)
    input_bai = get_bam.rv(1)

    if config.preprocess:
        preprocess = job.wrapJobFn(run_gatk_preprocessing,
                                   config.cores,
                                   get_bam.rv(0),
                                   get_bam.rv(1),
                                   config.genome_fasta,
                                   config.genome_dict,
                                   config.genome_fai,
                                   config.phase,
                                   config.mills,
                                   config.dbsnp,
                                   xmx=config.xmx).encapsulate()
        get_bam.addChild(preprocess)

        if config.preprocess_only:
            output_dir = os.path.join(config.output_dir, uuid)
            filename = '{}.{}.bam'.format(uuid, config.suffix)
            output_bam = job.wrapJobFn(upload_or_move_job,
                                       filename,
                                       preprocess.rv(0),
                                       output_dir)
            preprocess.addChild(output_bam)
            return preprocess.rv(0), preprocess.rv(1)
        else:
            input_bam = preprocess.rv(0)
            input_bai = preprocess.rv(1)

    hc_disk = PromisedRequirement(lambda x: 2 * x.size + human2bytes('5G'), input_bam)
    haplotype_caller = job.wrapJobFn(gatk_haplotype_caller,
                                     input_bam,
                                     input_bai,
                                     config,
                                     cores=config.cores,
                                     disk=hc_disk)

    get_bam.addFollowOn(haplotype_caller)
    output_dir = os.path.join(config.output_dir, uuid)
    vqsr_name = '{}.raw{}.gvcf'.format(uuid, config.suffix)
    output_gvcf = job.wrapJobFn(upload_or_move_job,
                                vqsr_name,
                                haplotype_caller.rv(),
                                output_dir)
    haplotype_caller.addChild(output_gvcf)
    return haplotype_caller.rv()


def main():
    """
    GATK Germline Pipeline with variant recalibration and annotation.

    Steps in Pipeline:
        1: Download Shared Files
        2: Prepare References
        3: Download Sample
        4: Align FASTQ (Optional)
        5: Call SNPs & INDELs
        6: Filter Variants
        7: Genotype and Annotate Variants

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
    subparsers.add_parser('generate-inputs',
                          help='Generates an editable inputs in the current working directory.')
    subparsers.add_parser('generate-manifest',
                          help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate',
                          help='Generates a inputs and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run',
                                       help='Runs the GATK germline pipeline')
    parser_run.add_argument('--config',
                            required=True,
                            type=str,
                            help='Path to the (filled in) config file, generated with '
                                 '"generate-config".')
    parser_run.add_argument('--manifest',
                            type=str,
                            help='Path to the (filled in) manifest file, generated with '
                                 '"generate-manifest". '
                                 '\nDefault value: "%(default)s".')
    parser_run.add_argument('--sample',
                            default=None,
                            nargs=2,
                            type=str,
                            help='Space delimited sample UUID and BAM file in the format: uuid url')
    parser_run.add_argument('--output-dir',
                            default=None,
                            help='Path/URL to output directory (Default: cwd).')
    parser_run.add_argument('-s', '--suffix',
                            default=None,
                            help='Additional suffix to add to the names of the output files')
    parser_run.add_argument('--preprocess-only',
                            action='store_true',
                            help='Only runs preprocessing steps')
    Job.Runner.addToilOptions(parser_run)
    options = parser.parse_args()

    cwd = os.getcwd()
    if options.command == 'generate-inputs' or options.command == 'generate':
        generate_file(os.path.join(cwd, 'inputs-toil-germline.yaml'), generate_config)
    if options.command == 'generate-manifest' or options.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-toil-germline.tsv'), generate_manifest)
    # Pipeline execution
    elif options.command == 'run':
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program)),
                    program + ' must be installed on every node.'.format(program))

        require(os.path.exists(options.config), '{} not found. Please run '
                                                '"generate-inputs"'.format(options.config))

        require(options.manifest or options.sample, 'Must provide path to manifest or '
                                                    'sample information at the command line')

        # Read sample manifest
        samples = []
        if options.manifest:
            samples.extend(parse_manifest(options.manifest))

        # Add BAM sample from command line
        if options.sample:
            uuid, url = options.sample
            samples.append((uuid, url, None, None))

        if len(samples) == 0:
            raise ValueError('No samples were detected in the manifest or on the command line')

        # Parse inputs
        inputs = {x.replace('-', '_'): y for x, y in
                  yaml.load(open(options.config).read()).iteritems()}

        required_fields = {'genome_fasta',
                           'assembly',
                           'output_dir',
                           'run_bwa',
                           'sorted',
                           'preprocess',
                           'preprocess_only',
                           'run_vqsr',
                           'joint',
                           'run_oncotator',
                           'cores',
                           'file_size',
                           'xmx',
                           'suffix'}

        input_fields = set(inputs.keys())
        require(input_fields > required_fields,
                'Missing config parameters:\n{}'.format(', '.join(required_fields - input_fields)))

        if inputs['output_dir'] is None:
            inputs['output_dir'] = options.output_dir \
                if options.output_dir else os.path.join(os.getcwd(), 'output')

        if inputs['suffix'] is None:
            inputs['suffix'] = options.suffix if options.suffix else ''

        if inputs['preprocess_only'] is None:
            inputs['preprocess_only'] = options.preprocess_only

        if inputs['run_vqsr']:
            vqsr_fields = {'phase', 'mills', 'dbsnp', 'hapmap', 'omni'}
            require(input_fields > vqsr_fields,
                    'Missing parameters for VQSR:\n{}'
                    .format(', '.join(vqsr_fields - input_fields)))

        # Set resource parameters
        inputs['xmx'] = convert_space_resource(inputs['xmx'], '10G')
        inputs['file_size'] = convert_space_resource(inputs['file_size'], '10G')
        inputs['cores'] = inputs.get('cores', 8)

        # GATK recommended annotations:
        # https://software.broadinstitute.org/gatk/documentation/article?id=2805
        inputs['annotations'] = ['QualByDepth',
                                 'FisherStrand',
                                 'StrandOddsRatio',
                                 'ReadPosRankSumTest',
                                 'MappingQualityRankSumTest',
                                 'RMSMappingQuality',
                                 'InbreedingCoeff']

        # Login to Synapse using the Python API
        if inputs.get('synapse_credentials', None) is not None:
            synapse_name, synapse_pwd = open(inputs['synapse_credentials'],
                                             'r').readline().strip().split(',')
            inputs['synapse_login'] = synapseclient.login(synapse_name,
                                                          synapse_pwd,
                                                          rememberMe=False,
                                                          silent=True)
        else:
            inputs['synapse_login'] = None

        # It is a toil-scripts convention to store input parameters in a Namespace object
        config = argparse.Namespace(**inputs)

        # Download shared files and start germline pipeline
        shared_files = Job.wrapJobFn(download_shared_files, config).encapsulate()
        run_pipeline = Job.wrapJobFn(run_gatk_germline_pipeline,
                                     samples,
                                     shared_files.rv()).encapsulate()
        shared_files.addChild(run_pipeline)

        if inputs['run_oncotator']:
            require(inputs['assembly'].lower() == 'hg19',
                    'Check your reference assembly. We recommend using hg19 when using Oncotator')
            annotate = Job.wrapJobFn(annotate_vcfs, run_pipeline.rv(), shared_files.rv())
            run_pipeline.addChild(annotate)

        Job.Runner.startToil(shared_files, options)


if __name__ == '__main__':
    main()
