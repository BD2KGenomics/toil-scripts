#!/usr/bin/env python2.7
"""
Pipeline takes a FASTQ or BAM file and generates a filtered VCF file using GATK 3.5

BWA Alignment
0: Download FASTQ(s) or BAM
1: Align to reference
2: Sort BAM
3: Index Bam

GATK Preprocessing
4: Mark duplicates
5: Indel Realignment
6: Base Quality Score Recalibration
7: Apply Recalibration

GATK Variant Discovery
8: Call variants
9: Genotype variants

GATK Hard Filtering
10: Isolate SNPs
11: Isolate Indels
12: Apply SNP filter
13: Apply Indel filter
14: Merge SNP and Indel VCFs

GATK VQSR and Filtering (Optional)
15: Recalibrate SNPs
16: Recalibrate Indels
17: Apply SNP recalibration
18: Apply Indel recalibration
"""
import argparse
from collections import namedtuple
from copy import deepcopy
import logging
import math
from multiprocessing import cpu_count
import os
import re
from urlparse import urlparse

from bd2k.util.humanize import human2bytes
from bd2k.util.processes import which
from toil.job import Job, PromisedRequirement
import yaml

from toil_scripts.gatk_germline.common import output_file_job
from toil_scripts.gatk_germline.germline_config import generate_config, generate_manifest
from toil_scripts.gatk_germline.hard_filter import hard_filter_pipeline
from toil_scripts.gatk_germline.vqsr import vqsr_pipeline
from toil_scripts.lib import require
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url_job
from toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline import generate_file
from toil_scripts.tools.aligners import run_bwakit
from toil_scripts.tools.indexing import run_samtools_faidx
from toil_scripts.tools.preprocessing import run_gatk_preprocessing, \
    run_picard_create_sequence_dictionary, run_samtools_index, run_samtools_sort
from toil_scripts.tools.variant_annotation import gatk_genotype_gvcfs, run_oncotator
from toil_scripts.tools.variant_filters import gatk_combine_gvcfs


logging.basicConfig(level=logging.INFO)


class GermlineSample(namedtuple('GermlineSample', 'uuid url paired_url rg_line')):
    """
    Namedtuple subclass for Toil Germline samples. Stores essential sample information.

    Attributes
    uuid: unique sample identifier that matches sample name in read group
    url: URL/PATH to FASTQ or BAM file
    paired_url: URL/PATH to paired FASTQ file
    rg_line: Read group information (i.e. @RG\tID:foo\tSM:bar)
    """


def run_gatk_germline_pipeline(job, samples, config):
    """
    Downloads shared files and runs the GATK Best Practices Germline pipeline for a list of samples.
    See `toil_scripts.gatk_germline.germline_config.py` file for necessary parameters.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list[GermlineSample] samples: List of GermlineSample namedtuples
    :param Namespace config: Configuration options for pipeline
    """
    # Download shared files and start germline pipeline
    shared_files = Job.wrapJobFn(download_shared_files, config).encapsulate()
    job.addChild(shared_files)

    if config.preprocess_only:
        for germline_sample in samples:
            shared_files.addChildJobFn(prepare_bam,
                                       germline_sample.uuid,
                                       germline_sample.url,
                                       shared_files.rv(),
                                       paired_url=germline_sample.paired_url,
                                       rg_line=germline_sample.rg_line)
    else:
        run_pipeline = Job.wrapJobFn(gatk_germline_pipeline,
                                     samples,
                                     shared_files.rv()).encapsulate()
        shared_files.addChild(run_pipeline)

        if config.run_oncotator:
            annotate = Job.wrapJobFn(annotate_vcfs, run_pipeline.rv(), shared_files.rv())
            run_pipeline.addChild(annotate)


def gatk_germline_pipeline(job, samples, config):
    """
    Configures and runs the GATK Best Practices Germline pipeline.

    Steps in Pipeline
    0: Generate and preprocess BAM
        - Uploads processed BAM to output directory
    1: Call Variants using HaplotypeCaller
        - Uploads GVCF to output directory
    2: Genotype VCF or Joint Genotype
        - Uploads genotyped VCF
    3: Filter Variants using either hard filter or VQSR
        - Uploads filtered VCF

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list samples: List of tuples (uuid, url, paired_url, rg_line)
    :param Namespace config: Input parameters and reference FileStoreIDs
    :return: Dictionary of filtered VCF FileStoreIDs
    :rtype: dict
    """
    require(len(samples) > 0, 'No samples were provided!')

    # 0: Generate processed BAM and BAI files for each sample
    # group preprocessing and variant calling steps in empty Job instance
    group_bam_jobs = Job()
    gvcfs = {}
    for uuid, url, paired_url, rg_line in samples:
        # 0: Generate processed BAM and BAI files for each sample
        get_bam = group_bam_jobs.addChildJobFn(prepare_bam,
                                               uuid,
                                               url,
                                               config,
                                               paired_url=paired_url,
                                               rg_line=rg_line)

        # 1: Generate per sample gvcfs {uuid: gvcf_id}
        # Estimate disk resource as twice the input bam plus 5Gs of reference data
        hc_disk = PromisedRequirement(lambda bam: 2 * bam.size + human2bytes('5G'), get_bam.rv(0))
        get_gvcf = get_bam.addFollowOnJobFn(gatk_haplotype_caller,
                                            get_bam.rv(0),
                                            get_bam.rv(1),
                                            config,
                                            cores=config.cores,
                                            disk=hc_disk,
                                            memory=config.xmx)
        gvcfs[uuid] = get_gvcf.rv()

        # Upload gvcf before genotyping
        vqsr_name = '{}{}.g.vcf'.format(uuid, config.suffix)
        get_gvcf.addChildJobFn(output_file_job,
                               vqsr_name,
                               gvcfs[uuid],
                               os.path.join(config.output_dir, uuid),
                               s3_key_path=config.ssec)

    # VQSR requires many variants in order to train a decent model. GATK recommends a minimum of
    # 30 exomes or one large WGS sample:
    # https://software.broadinstitute.org/gatk/documentation/article?id=3225

    filtered_vcfs = {}
    if config.joint:
        filtered_vcfs = job.addChildJobFn(joint_genotype_and_filter,
                                          gvcfs.items(),
                                          config,
                                          cores=cpu_count()  # Reserve entire node
                                          ).rv()

    else:
        for uuid, gvcf_id in gvcfs.iteritems():
            filtered_vcfs[uuid] = job.addChildJobFn(genotype_and_filter,
                                                    uuid,
                                                    gvcf_id,
                                                    config).rv()

    job.addChild(group_bam_jobs)
    return filtered_vcfs


def joint_genotype_and_filter(job, iterable, config, disk_fraction=0.5):
    """
    Batches GVCFs into one GVCF file and runs the genotype and filter pipeline.
    Request all available cores when calling this function to get a better estimate of the
    available disk space on the worker node.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list[tuples] iterable: Paired UUID and FileStoreIDs
    :param float disk_fraction: Fraction of available disk to reserve for VCFs (Default: 0.5)
    :returns: Yields a batch of paired UUIDs and FileStoreIDs
    :rtype: dict
    """
    work_dir = job.fileStore.getLocalTempDir()
    st = os.statvfs(work_dir)
    free_disk = st.f_bavail * st.f_frsize
    usable_disk = disk_fraction * free_disk

    total_size = sum(elem[1].size for elem in iterable)

    # Number of times the usable disk space goes into the total file size
    times_in = total_size / usable_disk

    # Determine how many batches there will be
    batches = [[] for _ in range(int(math.ceil(times_in)))]
    batch_sums = {index: 0 for index in xrange(len(batches))}
    # Sort iterable by the file size
    sorted_iterable = sorted(iterable, key=lambda elem: elem[1].size, reverse=True)
    # Iterate over sorted list of FileStoreIDs in descending order. Then, assign a GVCF
    # to the smallest GVCF batch. Finally, merge all merged batches into a single GVCF
    for elem in sorted_iterable:
        # Add VCF to batch with the smallest disk footprint
        min_batch_sum = min(batch_sums, key=batch_sums.get)
        batches[min_batch_sum].append(elem)
        # Unpack file store ID
        _, FileStoreID = elem
        batch_sums[min_batch_sum] += FileStoreID.size

    merged_batches = {}
    for index, batch in enumerate(batches):
        batch_name = 'joint-batch-%d' % index
        batch_dict = dict(batch)
        merged_batch_disk = PromisedRequirement(lambda lst: int(1.5 * sum(x.size for x in lst)),
                                               batch_dict.values())
        merged_batches[batch_name] = job.addChildJobFn(gatk_combine_gvcfs,
                                                       batch_dict,
                                                       config,
                                                       disk=merged_batch_disk,
                                                       memory=config.xmx).rv()

    total_merge_disk = PromisedRequirement(lambda lst: int(1.5 * sum(x.size for x in lst)),
                                           merged_batches.values())
    total_merge_job = job.addFollowOnJobFn(gatk_combine_gvcfs,
                                           merged_batches,
                                           config,
                                           disk=total_merge_disk,
                                           memory=config.xmx)
    return total_merge_job.addChildJobFn(genotype_and_filter,
                                         'joint',
                                         total_merge_job.rv(),
                                         config).rv()


def genotype_and_filter(job, uuid, gvcf, config):
    """
    Genotypes and filters a single VCF file using VQSR or hard filters

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: sample unique identifier
    :param str gvcf: GVCF FileStoreIDs
    :param Namespace config: Input parameters and shared FileStoreIDs
    :return: Genotyped and filtered VCF FileStoreID
    :rtype: str
    """
    # Estimate disk requirement as 50% more than the sum of the VCF file sizes plus 5G for
    # reference data.
    genotype_gvcf_disk = PromisedRequirement(lambda x: int(1.5 * x.size + human2bytes('5G')),
                                                           gvcf)
    genotype_gvcf = job.addChildJobFn(gatk_genotype_gvcfs,
                                      {uuid: gvcf},
                                      config,
                                      cores=config.cores,
                                      disk=genotype_gvcf_disk,
                                      memory=config.xmx)

    genotyped_filename = '%s.genotyped%s.vcf' % (uuid, config.suffix)
    genotyped_gvcf_disk = PromisedRequirement(lambda x: x.size, genotype_gvcf.rv())
    genotype_gvcf.addChildJobFn(output_file_job,
                                genotyped_filename,
                                genotype_gvcf.rv(),
                                config.output_dir,
                                s3_key_path=config.ssec,
                                disk=genotyped_gvcf_disk)

    if config.run_vqsr:
        if not config.joint:
            job.fileStore.logToMaster('WARNING: Running VQSR without joint genotyping!')
        joint_vcf = genotype_gvcf.addFollowOnJobFn(vqsr_pipeline,
                                                   uuid,
                                                   genotype_gvcf.rv(),
                                                   config)
    else:
        joint_vcf = genotype_gvcf.addFollowOnJobFn(hard_filter_pipeline,
                                                   uuid,
                                                   genotype_gvcf.rv(),
                                                   config).rv()
    return joint_vcf.rv()


def annotate_vcfs(job, vcfs, config):
    """
    Runs Oncotator for a group of VCF files. Each sample is annotated individually.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param dict vcfs: Dictionary of UUIDs and VCF FileStoreIDs
    :param Namespace config: Input parameters and shared FileStoreIDs
    """
    job.fileStore.logToMaster('Annotating VCFs')
    for uuid, vcf_id in vcfs.iteritems():
        # Estimate disk space as twice the input VCF size plus the database size
        onco_disk = PromisedRequirement(lambda vcf, db: 2 * vcf.size + db.size,
                                        vcf_id,
                                        config.oncotator_db)
        annotated_vcf = job.addChildJobFn(run_oncotator, vcf_id, config,
                                          disk=onco_disk, cores=config.cores, memory=config.xmx)

        output_dir = os.path.join(config.output_dir, uuid)
        filename = '{}.annotated{}.vcf'.format(uuid, config.suffix)
        annotated_vcf.addChildJobFn(output_file_job,
                                    filename,
                                    annotated_vcf.rv(),
                                    output_dir,
                                    s3_key_path=config.ssec)


#####################################
#  Pipeline convenience functions   #
#####################################


def parse_manifest(path_to_manifest):
    """
    Parses manifest file for Toil Germline Pipeline

    :param str path_to_manifest: Path to sample manifest file
    :return: List of GermlineSample namedtuples
    :rtype: list[GermlineSample]
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
                require('.bam' in url.lower(),
                        'Expected .bam extension:\n{}:\t{}'.format(uuid, url))
            elif fastq_match:
                uuid = fastq_match.group('uuid')
                url = fastq_match.group('url')
                url2 = fastq_match.group('url2')
                rg_line = fastq_match.group('rg_line')
                require('.fq' in url.lower() or '.fastq' in url.lower(),
                        'Expected .fq extension:\n{}:\t{}'.format(uuid, url))
            else:
                msg = 'Could not parse entry in manifest: %s\n%s' % (f.name, line)
                raise ValueError(msg)
            require(urlparse(url).scheme or url.startswith('syn'),
                    'Invalid URL passed for {}'.format(url))
            samples.append(GermlineSample(uuid, url, url2, rg_line))
    return samples


def convert_space_resource(value, default):
    """
    Converts resource requirement from human readable bytes to raw integer value.

    :param str|int|None value: Space resource value
    :param str default: Human readable parameter size
    :return: bytes
    :rtype: int
    """
    if value is None:
        return human2bytes(default)
    elif isinstance(value, str) and re.match('\d+\s*[MKG]', value.upper()):
        return human2bytes(value)
    elif isinstance(value, int):
        return value
    else:
        raise ValueError('Could not convert resource requirement: %s' % value)


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
                                                    s3_key_path=config.ssec,
                                                    disk='15G'   # Estimated reference file size
                                                    ).rv())
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
                                               cores=config.cores,
                                               memory=config.xmx).rv()
    return config


def prepare_bam(job, uuid, url, config, paired_url=None, rg_line=None):
    """
    Prepares BAM file for GATK Germline pipeline.

    Steps in pipeline
    0: Align FASTQ or Download BAM
    1: Sort BAM
    2: Index BAM
    3: GATK Preprocessing (Optional)
        - Uploads preprocessed BAM to output directory

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: Unique identifier for the sample
    :param str url: URL to BAM file or FASTQs
    :param Namespace config: Configuration options for pipeline
    :param str|None paired_url: URL to paired FASTQ file
    :param str|None rg_line: RG line for BWA alignment (i.e. @RG\tID:foo\tSM:bar). Default is None.
    :return: BAM and BAI FileStoreIDs
    :rtype: tuple
    """
    # 0: Align FASQ or realign BAM
    if config.run_bwa:
        get_bam = job.wrapJobFn(setup_and_run_bwa_kit,
                                uuid,
                                url,
                                rg_line,
                                config,
                                paired_url=paired_url).encapsulate()

    # 0: Download BAM
    elif '.bam' in url.lower():
        job.fileStore.logToMaster("Downloading BAM: %s" % uuid)
        get_bam = job.wrapJobFn(download_url_job,
                                url,
                                name='toil.bam',
                                s3_key_path=config.ssec,
                                disk=config.file_size).encapsulate()
    else:
        raise ValueError('Could not generate BAM file for %s\n'
                         'Provide a FASTQ URL and set run-bwa or '
                         'provide a BAM URL that includes .bam extension.' % uuid)

    # 1: Sort BAM file if necessary
    # Realigning BAM file shuffles read order
    if config.sorted and not config.run_bwa:
        sorted_bam = get_bam

    else:
        sorted_bam_disk = PromisedRequirement(lambda bam: 3 * bam.size, get_bam.rv())
        sorted_bam = get_bam.addChildJobFn(run_samtools_sort,
                                           get_bam.rv(),
                                           cores=config.cores,
                                           disk=sorted_bam_disk)

    # 2: Index BAM
    index_bam_disk = PromisedRequirement(lambda bam: bam.size, sorted_bam.rv())
    index_bam = job.wrapJobFn(run_samtools_index, sorted_bam.rv(), disk=index_bam_disk)

    job.addChild(get_bam)
    sorted_bam.addChild(index_bam)

    if config.preprocess:
        preprocess = job.wrapJobFn(run_gatk_preprocessing,
                                   sorted_bam.rv(),
                                   index_bam.rv(),
                                   config.genome_fasta,
                                   config.genome_dict,
                                   config.genome_fai,
                                   config.phase,
                                   config.mills,
                                   config.dbsnp,
                                   memory=config.xmx, cores=config.cores).encapsulate()
        sorted_bam.addChild(preprocess)
        index_bam.addChild(preprocess)

        # Update output BAM promises
        output_bam_promise = preprocess.rv(0)
        output_bai_promise = preprocess.rv(1)

        # Save processed BAM
        output_dir = os.path.join(config.output_dir, uuid)
        filename = '{}.processed{}.bam'.format(uuid, config.suffix)
        output_bam = job.wrapJobFn(output_file_job,
                                   filename,
                                   preprocess.rv(0),
                                   output_dir,
                                   s3_key_path=config.ssec)
        preprocess.addChild(output_bam)

    else:
        output_bam_promise = sorted_bam.rv()
        output_bai_promise = index_bam.rv()

    return output_bam_promise, output_bai_promise


def setup_and_run_bwa_kit(job, uuid, url, rg_line, config, paired_url=None):
    """
    Downloads and aligns FASTQs or realigns BAM using BWA.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str url: FASTQ or BAM sample URL. BAM alignment URL must have .bam extension.
    :param Namespace config: Input parameters and shared FileStoreIDs
    :param str|None paired_url: URL to paired FASTQ
    :param str|None rg_line: Read group line (i.e. @RG\tID:foo\tSM:bar)
    :return: BAM FileStoreID
    :rtype: str
    """
    bwa_config = deepcopy(config)
    bwa_config.uuid = uuid
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
            'Please use .fq or .bam file extensions:\n%s' % url)

    # Download fastq files
    samples = []
    input1 = job.addChildJobFn(download_url_job,
                               url,
                               name='file1',
                               s3_key_path=config.ssec,
                               disk=config.file_size)
    samples.append(input1.rv())
    if '.bam' in url.lower():
        bwa_config.bam = input1.rv()
    else:
        bwa_config.r1 = input1.rv()

    # Uses first url to deduce paired FASTQ URL, if url looks like a paired file.
    if paired_url is None and '1.fq' in url:
        paired_url = url.replace('1.fq', '2.fq')
        job.fileStore.logToMaster(
            'Trying to find paired reads using FASTQ URL:\n%s\n%s' % (url, paired_url))

    if paired_url:
        input2 = job.addChildJobFn(download_url_job,
                                   paired_url,
                                   name='file2',
                                   s3_key_path=config.ssec,
                                   disk=config.file_size)
        samples.append(input2.rv())
        bwa_config.r2 = input2.rv()

    # Disk requirement depends on whether there are paired reads, so take the sum of all
    # elements in the samples list and scale it by a factor of 3. Then, add 8G for index data
    bwakit_disk = PromisedRequirement(lambda lst: 3 * sum(x.size for x in lst) + human2bytes('8G'),
                                      samples)
    return job.addFollowOnJobFn(run_bwakit,
                                bwa_config,
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
    job.fileStore.logToMaster('Running GATK HaplotypeCaller')
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': config.genome_fasta,
              'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict,
              'input.bam': bam_id, 'input.bam.bai': bai_id}
    for name, file_store_id in inputs.iteritems():
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))

    # Call GATK -- HaplotypeCaller with parameters to produce a genomic VCF file:
    # https://software.broadinstitute.org/gatk/documentation/article?id=2803
    command = ['-T', 'HaplotypeCaller',
               '-nct', str(job.cores),
               '-R', 'genome.fa',
               '-I', 'input.bam',
               '-o', 'output.g.vcf',
               '-stand_call_conf', str(call_threshold),
               '-stand_emit_conf', str(emit_threshold),
               '-variant_index_type', 'LINEAR',
               '-variant_index_parameter', '128000',
               '--genotyping_mode', 'Discovery',
               '--emitRefConfidence', 'GVCF']

    # Workaround for Adam/GATK pipeline
    if config.unsafe_mode:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    for annotation in config.annotations:
        if annotation is not None:
            command.extend(['-A', annotation])

    outputs = {'output.g.vcf': None}
    docker_call(work_dir=work_dir,
                env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters=command,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.g.vcf'))


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
    Python 2.7
    curl            - apt-get install curl
    virtualenv      - apt-get install python-virtualenv
    pip             - apt-get install python-pip
    toil            - pip install toil
    docker          - http://docs.docker.com/engine/installation/
    """
    # Define Parser object and add to jobTree
    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    # Generate subparsers
    subparsers = parser.add_subparsers(dest='command')
    subparsers.add_parser('generate-config',
                          help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest',
                          help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate',
                          help='Generates a config and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the GATK germline pipeline')
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
                            help='Path/URL to output directory')
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
            require(next(which(program)),
                    program + ' must be installed on every node.'.format(program))

        require(os.path.exists(options.config), '{} not found. Please run '
                                                '"generate-config"'.format(options.config))

        # Read sample manifest
        samples = []
        if options.manifest:
            samples.extend(parse_manifest(options.manifest))

        # Add BAM sample from command line
        if options.sample:
            uuid, url = options.sample
            # samples tuple: (uuid, url, paired_url, rg_line)
            # BAM samples should not have as paired URL or read group line
            samples.append(GermlineSample(uuid, url, None, None))

        require(len(samples) > 0,
                'No samples were detected in the manifest or on the command line')

        # Parse inputs
        inputs = {x.replace('-', '_'): y for x, y in
                  yaml.load(open(options.config).read()).iteritems()}

        required_fields = {'genome_fasta',
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
            inputs['output_dir'] = options.output_dir

        require(inputs['output_dir'] is not None,
                'Missing output directory PATH/URL')

        if inputs['suffix'] is None:
            inputs['suffix'] = options.suffix if options.suffix else ''

        if inputs['preprocess_only'] is None:
            inputs['preprocess_only'] = options.preprocess_only

        if inputs['run_vqsr']:
            vqsr_fields = {'phase', 'mills', 'dbsnp', 'hapmap', 'omni'}
            require(input_fields > vqsr_fields,
                    'Missing parameters for VQSR:\n{}'
                    .format(', '.join(vqsr_fields - input_fields)))

        logging.info('Configuring resource parameters')
        # Set resource parameters
        inputs['xmx'] = convert_space_resource(inputs['xmx'], '10G')
        logging.info('Setting Java memory allocation to %d bytes.' % inputs['xmx'])

        inputs['file_size'] = convert_space_resource(inputs['file_size'], '10G')
        logging.info('Setting input file size to %d bytes' % inputs['file_size'])

        inputs['cores'] = inputs.get('cores', 8)
        logging.info('Allocating %d cores per job.' % inputs['cores'])

        # GATK recommended annotations:
        # https://software.broadinstitute.org/gatk/documentation/article?id=2805
        inputs['annotations'] = ['QualByDepth',
                                 'FisherStrand',
                                 'StrandOddsRatio',
                                 'ReadPosRankSumTest',
                                 'MappingQualityRankSumTest',
                                 'RMSMappingQuality',
                                 'InbreedingCoeff']

        # It is a toil-scripts convention to store input parameters in a Namespace object
        config = argparse.Namespace(**inputs)

        root = Job.wrapJobFn(run_gatk_germline_pipeline, samples, config)
        Job.Runner.startToil(root, options)


if __name__ == '__main__':
    main()
