#!/usr/bin/env python2.7
from __future__ import print_function
import argparse
from copy import deepcopy
import multiprocessing
import os
from urlparse import urlparse

from bd2k.util.humanize import human2bytes
from bd2k.util.processes import which
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
from toil_scripts.tools.preprocessing import run_germline_preprocessing, run_samtools_sort, \
    run_samtools_index, run_samtools_view, run_picard_create_sequence_dictionary


def download_and_run(job, uuid, url, config, rg_line=None):
    """
    Downloads reference files and runs the GATK best practices germline pipeline for a single
    samples.

    :param JobFunctionWrappingJob job: Toil Job instance
    :param str uuid: Unique identifier for the sample
    :param str url: URL to BAM file
    :param Namespace config: Configuration options for pipeline
    :param str rg_line: RG line for BWA alignment. Default is None.
    :return: None
    """
    get_shared_files = job.wrapJobFn(download_shared_files, config).encapsulate()
    samples = [(uuid, url, rg_line)]
    run_pipeline = Job.wrapJobFn(run_gatk_germline_pipeline, samples, get_shared_files.rv())
    get_shared_files.addChild(run_pipeline)


def run_gatk_germline_pipeline(job, samples, config):
    """
    Configures and runs the GATK germline pipeline with filtering.

    VQSR is done if the run_vqsr parameter is set to True or there are at least
    30 samples. Otherwise, each sample is hard filtered using recommended GATK parameters.

    :param list samples: List of tuples (uuid, url, rg_line) for joint variant calling.
                         (i.e. samples are grouped together for VQSR). Otherwise,
                         samples are filtered using GATK recommended hard filters.
    :param Namespace config: Input parameters and reference files
    :param str suffix: Adds suffix to end of filename.
    :return: None
    """
    if len(samples) == 0:
        raise ValueError('No samples were provided!')

    cores = multiprocessing.cpu_count()

    # Generate per sample gvcfs {uuid: gvcf_id}
    gvcfs = {}
    for uuid, url, rg_line in samples:
        gvcfs[uuid] = job.addChildJobFn(gatk_germline_pipeline,
                                        uuid,
                                        url,
                                        config,
                                        rg_line=rg_line).rv()

    # Stop after preprocessing
    if config.preprocess_only:
        return

    # VQSR requires many variants in order to train a decent model. GATK recommends a minimum of
    # 30 exomes or one large WGS sample:
    # https://software.broadinstitute.org/gatk/documentation/article?id=3225
    if config.run_vqsr or len(gvcfs) > 30:
        # TODO split joint VCF by sample
        if config.joint:
            joint_vcf = job.addFollowOnJobFn(vqsr_pipeline, gvcfs, config)

        else:
            for uuid, gvcf in gvcfs.iteritems():
                job.addFollowOnJobFn(vqsr_pipeline, dict(uuid=gvcf), config)
    else:
        for uuid, gvcf in gvcfs.iteritems():
            job.addFollowOnJobFn(hard_filter_pipeline,
                                 uuid, gvcf,
                                 config)


def gatk_germline_pipeline(job, uuid, url, config, rg_line=None):
    """
    Runs GATK Germline Pipeline on a single sample. Writes gvcf to output directory
    defined in the config dictionary. Removes secondary alignments.

    0: Align FASTQ or Download Bam
    1: Generate BAI
    2: GATK Preprocessing (Optional)
    3: GATK HaplotypeCaller
    4: Write GVCF to output directory

    :param JobFunctionWrappingJob job: Toil Job instance
    :param str uuid: Unique identifier for the sample
    :param str url: URL to sample bam or fq file
    :param Namespace config: Configuration options for pipeline
    :param str rg_line: RG line for BWA alignment. Default is None.
    :return: GVCF FileStoreID
    :rtype: str
    """
    config = deepcopy(config)
    config.uuid = uuid
    config.url = url
    config.rg_line = rg_line

    cores = multiprocessing.cpu_count()

    # Determine how much RAM is available
    if config.xmx is None:
        config.xmx = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')

    get_bam = job.wrapJobFn(prepare_bam, uuid, url, config, rg_line=rg_line).encapsulate()

    job.addChild(get_bam)
    input_bam = get_bam.rv(0)
    input_bai = get_bam.rv(1)

    if config.preprocess:
        preprocess = job.wrapJobFn(run_germline_preprocessing,
                                   cores,
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

    hc_disk = PromisedRequirement(lambda x: 2*x.size + human2bytes('5G'), input_bam)
    haplotype_caller = job.wrapJobFn(gatk_haplotype_caller,
                                     input_bam,
                                     input_bai,
                                     config,
                                     memory=config.xmx,
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


def download_shared_files(job, config):
    """
    Downloads reference files shared by all samples in the pipeline

    :param Namespace config: Pipeline configuration options
    :return: Updated config with shared fileStoreIDS
    :rtype: Namespace
    """
    job.fileStore.logToMaster('Downloading shared reference files')
    references = {'genome_fasta'}
    unrequired_references = {'genome_fai', 'genome_dict'}
    references |= unrequired_references
    if hasattr(config, 'ref'): config.genome_fasta = config.ref
    if hasattr(config, 'fai'): config.genome_fai = config.fai
    if hasattr(config, 'dict'): config.genome_dict = config.dict
    if config.run_bwa:
        bwa_references = {'amb', 'ann', 'bwt', 'pac', 'sa', 'alt'}
        unrequired_references.add('alt')
        references |= bwa_references
    if config.preprocess:
        preprocess_references = {'mills', 'phase', 'dbsnp'}
        references |= preprocess_references
    if config.run_vqsr:
        vqsr_references = {'phase', 'mills', 'dbsnp', 'hapmap', 'omni'}
        references |= vqsr_references
    for name in references:
        try:
            url = getattr(config, name)
            if url is None:
                continue
            setattr(config, name, job.addChildJobFn(download_url_job,
                                                    url,
                                                    name=name,
                                                    s3_key_path=config.ssec).rv())
        except AttributeError:
            setattr(config, name, None)
        finally:
            if getattr(config, name) is None and name not in unrequired_references:
                raise ValueError("Necessary parameter is missing:\n{}".format(name))
    return job.addFollowOnJobFn(reference_preprocessing, config).rv()


def reference_preprocessing(job, config):
    """
    Creates a genome.fa.fai and genome.dict file for reference genome.fa
    if one is not already present in the config dictionary

    :param JobFunctionWrapping Job: passed by Toil
    :param Namespace config: Pipeline configuration options and shared files.
                             Requires genome_fasta attribute
    :return: Updated config with reference index files
    :rtype: Namespace
    """
    job.fileStore.logToMaster('Preparing Reference Files')
    genome_id = config.genome_fasta
    if not hasattr(config, 'genome_fai') or config.genome_fai is None:
        config.genome_fai = job.addChildJobFn(run_samtools_faidx, genome_id).rv()
    if not hasattr(config, 'genome_dict') or config.genome_dict is None:
        config.genome_dict = job.addChildJobFn(run_picard_create_sequence_dictionary,
                                               genome_id).rv()
    return config


def prepare_bam(job, uuid, url, config, rg_line=None):
    """
    Prepares BAM file for GATK Germline pipeline.

    0: Align FASTQ or Download BAM
    1: Remove secondary alignments
    2: Sort BAM
    3: Index BAM

    :param JobFunctionWrappingJob job: Toil Job instance
    :param str uuid: Unique identifier for the sample
    :param str url: URL to BAM file
    :param Namespace config: Configuration options for pipeline
    :param str rg_line: RG line for BWA alignment. Default is None.
    :return: BAM and BAI FileStoreIDs
    :rtype: tuple
    """

    cores = multiprocessing.cpu_count()

    # Determine extension on input file
    base, ext1 = os.path.splitext(url)
    _, ext2 = os.path.splitext(base)

    # Produce or download BAM file
    if {ext1, ext2} & {'.fq', '.fastq'}:
        if not config.run_bwa:
            raise ValueError('Not configured to run BWA. Please set run-bwa: True')
        elif not rg_line:
            raise ValueError('FASTQ file requires an RG line for BWA alignment')
        else:
            get_bam = job.wrapJobFn(setup_and_run_bwa_kit, url, config).encapsulate()

    elif {ext1, ext2} & {'.bam'}:
        job.fileStore.logToMaster("Downloading BAM: %s" % uuid)
        get_bam = job.wrapJobFn(download_url_job,
                                url,
                                name='toil.bam',
                                s3_key_path=config.ssec,
                                disk=config.file_size).encapsulate()
    else:
        raise ValueError('Accepted formats: FASTQ or BAM\n'
                         'Could not determine {} file type from URL:\n{}'.format(uuid, url))

    rm_secondary = job.wrapJobFn(run_samtools_view,
                                 get_bam.rv(),
                                 flag='0x800',
                                 ncores=cores,
                                 disk=config.file_size, cores=cores)

    sort_bam_disk = PromisedRequirement(lambda x: 3*x.size + human2bytes('10G'), rm_secondary.rv())
    sort_bam = job.wrapJobFn(run_samtools_sort,
                             rm_secondary.rv(),
                             ncores=cores,
                             disk=sort_bam_disk)

    index_bam_disk = PromisedRequirement(lambda x: int(1.5 * x.size), sort_bam.rv())
    index_bam = job.wrapJobFn(run_samtools_index, sort_bam.rv(), disk=index_bam_disk)

    job.addChild(get_bam)
    get_bam.addChild(rm_secondary)
    rm_secondary.addChild(sort_bam)
    sort_bam.addChild(index_bam)
    return sort_bam.rv(), index_bam.rv()


def setup_and_run_bwa_kit(job, url, config):
    """
    Downloads and aligns paired FASTQs to reference genome using BWA.

    :param JobFunctionWrapping Job: passed by Toil
    :param str url: FASTQ sample URL
    :param Namespace config: Input parameters and shared FileStoreIDs
    :return: BAM FileStoreID
    :rtype: str
    """
    job.fileStore.logToMaster("Aligning sample: {}".format(config.uuid))
    cores = multiprocessing.cpu_count()
    # Download fastq files
    fq1 = job.addChildJobFn(download_url_job,
                            url,
                            name='toil.1.fq',
                            s3_key_path=config.ssec,
                            disk='50G')
    # Assumes second fastq url is identical
    fq2_url = url.replace('1.fq', '2.fq')
    fq2 = job.addChildJobFn(download_url_job,
                            fq2_url,
                            name='toil.2.fq',
                            s3_key_path=config.ssec,
                            disk='50G')

    bwa_config = deepcopy(config)
    bwa_config.r1 = fq1.rv()
    bwa_config.r2 = fq2.rv()
    # bwa_alignment uses a different naming convention
    bwa_config.ref = config.genome_fasta
    bwa_config.fai = config.genome_fai
    bwakit_disk = PromisedRequirement(lambda x, y: 3*(x.size + y.size) + human2bytes('8G'),
                                      fq1.rv(), fq2.rv())
    return job.addFollowOnJobFn(run_bwakit,
                                bwa_config,
                                cores,
                                sort=False,
                                trim=config.trim,
                                cores=cores,
                                disk=bwakit_disk).rv()


def gatk_haplotype_caller(job, bam_id, bai_id, config, emit_threshold=10.0, call_threshold=30.0):
    """
    Uses GATK HaplotypeCaller to identify SNPs and Indels and writes a gVCF.

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
    inputs = {'genome.fa': config.genome_fasta, 'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict, 'input.bam': bam_id, 'input.bam.bai': bai_id}
    get_files_from_filestore(job, work_dir, inputs)

    cores = multiprocessing.cpu_count()
    # Call GATK -- HaplotypeCaller with best practices parameters:
    # https://software.broadinstitute.org/gatk/documentation/article?id=2803
    command = ['-nct', str(cores),
               '-R', 'genome.fa',
               '-T', 'HaplotypeCaller',
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

    outputs={'output.gvcf': None}
    docker_call(work_dir = work_dir,
                env={'JAVA_OPTS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config.xmx)},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.gvcf'))


def parse_manifest(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param path_to_manifest str: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if not line.isspace() and not line.startswith('#'):
                sample = line.strip().split()
                if len(sample) == 2:
                    uuid, url = sample
                    rg_line = None
                    if not url.endswith('bam'):
                        raise ValueError("Expected a BAM file:\n{}:\t{}".format(uuid, url))
                elif len(sample ) == 3:
                    uuid, url, rg_line = sample
                    if not url.endswith('fq') or url.endswith('fastq'):
                        raise ValueError("Expected a FASTA file:\n{}:\t{}".format(uuid, url))
                else:
                    msg = 'Bad manifest format!\n{}'.format(sample)
                    raise ValueError(msg)
                require(urlparse(url).scheme and urlparse(url),
                        'Invalid URL passed for {}'.format(url))
                samples.append((uuid, url, rg_line))
    return samples


def main():
    """
    GATK HaplotypeCaller with variant recalibration.

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
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')

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
                            help='Path/URL to directory or filename where '
                                 'final results will be output')
    parser_run.add_argument('-s', '--suffix',
                            default=None,
                            help='Additional suffix to add to the names of the output files')
    parser_run.add_argument('--preprocess-only',
                            action='store_true',
                            help='Only runs preprocessing steps')
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

        require(options.manifest or options.sample, 'Must provide path to manifest or '
                                                    'sample information at the command line')

        # Obtain list of samples
        samples = []
        if options.manifest:
            samples.extend(parse_manifest(options.manifest))

        if options.sample:
            uuid, url = options.sample
            samples.append((uuid, url, None))

        if len(samples) == 0:
            raise ValueError('No samples were detected in the manifest or on the command line')

        # Parse config
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(options.config).read()).iteritems()}
        config_fields = set(config)
        required_fields = {'genome_fasta', 'output_dir', 'file_size', 'run_bwa', 'preprocess',
                           'run_vqsr', 'joint', 'xmx'}
        require(config_fields > required_fields,
                'Missing config parameters:\n{}'.format(', '.join(required_fields - config_fields)))

        config['file_size'] = human2bytes(config['file_size'])

        if config['output_dir'] is None:
            config['output_dir'] = options.output_dir \
                if options.output_dir else os.path.join(os.getcwd(), 'output')

        if config['suffix'] is None:
            config['suffix'] = options.suffix if options.suffix else ''

        if config['xmx'] is None:
            pass
        elif not config['xmx'][-1].isdigit():
            config['xmx'] = human2bytes(config['xmx'])
        else:
            config['xmx'] = int(config['xmx'])

        if 'preprocess_only' not in config or config['preprocess_only'] is None:
           config['preprocess_only'] = options.preprocess_only

        if config['run_vqsr']:
            vqsr_fields = {'phase', 'mills', 'dbsnp', 'hapmap', 'omni'}
            require(config_fields > vqsr_fields,
                    'Missing parameters for VQSR:\n{}'
                    .format(', '.join(vqsr_fields - config_fields)))

        # GATK recommended annotations:
        # https://software.broadinstitute.org/gatk/documentation/article?id=2805
        config['annotations'] = ['QualByDepth', 'FisherStrand', 'StrandOddsRatio',
                                 'ReadPosRankSumTest', 'MappingQualityRankSumTest',
                                 'RMSMappingQuality', 'InbreedingCoeff']

        # It is a convention to store configuration attributes in a Namespace object
        config = argparse.Namespace(**config)

        shared_files = Job.wrapJobFn(download_shared_files, config).encapsulate()
        run_pipeline = Job.wrapJobFn(run_gatk_germline_pipeline,samples, shared_files.rv())
        shared_files.addChild(run_pipeline)

        Job.Runner.startToil(shared_files, options)

if __name__ == '__main__':
    main()
