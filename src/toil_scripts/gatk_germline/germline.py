#!/usr/bin/env python2.7
"""
Runs GATK best practices pipeline for germline SNP and INDEL discovery.

BWA Alignment
0: Download FASTQ(s) or BAM
1: Align to reference
2: Sort BAM
3: Index Bam

GATK Preprocessing
4: Mark duplicates
5: Indel realignment
6: Base quality score recalibration
7: Apply recalibration

GATK Variant Discovery
8:  Call variants
9:  Merge variant calls
10: Genotype variants

GATK Hard Filtering or VQSR
11: Isolate SNPs
12: Isolate Indels
13: Apply SNP filter
14: Apply Indel filter
15: Merge SNP and Indel VCFs

11: Recalibrate SNPs
12: Recalibrate Indels
13: Apply SNP recalibration
14: Apply Indel recalibration

===================================================================

Dependencies
Python 2.7
curl            - apt-get install curl
virtualenv      - apt-get install python-virtualenv
pip             - apt-get install python-pip
toil            - pip install toil
docker          - http://docs.docker.com/engine/installation/
"""
from __future__ import print_function
import argparse
from collections import namedtuple
from copy import deepcopy
import logging
import os
import re
from urlparse import urlparse

from bd2k.util.humanize import human2bytes
from bd2k.util.processes import which
from toil.job import Job, PromisedRequirement
from toil_lib import require
from toil_lib.files import generate_file
from toil_lib.programs import docker_call
from toil_lib.tools.aligners import run_bwakit
from toil_lib.tools.indexing import run_samtools_faidx
from toil_lib.tools.preprocessing import run_gatk_preprocessing, \
    run_picard_create_sequence_dictionary, run_samtools_index, run_samtools_sort
from toil_lib.tools.variant_annotation import gatk_genotype_gvcfs, run_oncotator
from toil_lib.urls import download_url_job
import yaml

from toil_scripts.gatk_germline.common import output_file_job
from toil_scripts.gatk_germline.germline_config_manifest import generate_config, generate_manifest
from toil_scripts.gatk_germline.hard_filter import hard_filter_pipeline
from toil_scripts.gatk_germline.vqsr import vqsr_pipeline


logging.basicConfig(level=logging.INFO)


class GermlineSample(namedtuple('GermlineSample', 'uuid url paired_url rg_line')):
    """
    Namedtuple subclass for Toil Germline samples.

    Attributes
    uuid: unique sample identifier
    url: URL/PATH to FASTQ or BAM file
    paired_url: URL/PATH to paired FASTQ file, or None if BAM file
    rg_line: Read group information (i.e. @RG\tID:foo\tSM:bar), or None if BAM file
    """


def run_gatk_germline_pipeline(job, samples, config):
    """
    Downloads shared files and calls the GATK best practices germline pipeline for a cohort of samples

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list[GermlineSample] samples: List of GermlineSample namedtuples
    :param Namespace config: Configuration options for pipeline
        Requires the following config attributes:
        config.preprocess_only      If True, then stops pipeline after preprocessing steps
        config.joint_genotype       If True, then joint genotypes cohort
        config.run_oncotator        If True, then adds Oncotator to pipeline
        Additional parameters are needed for downstream steps. Refer to pipeline README for more information.
    """
    # Determine the available disk space on a worker node before any jobs have been run.
    work_dir = job.fileStore.getLocalTempDir()
    st = os.statvfs(work_dir)
    config.available_disk = st.f_bavail * st.f_frsize

    # Check that there is a reasonable number of samples for joint genotyping
    num_samples = len(samples)
    if config.joint_genotype and not 30 < num_samples < 200:
        job.fileStore.logToMaster('WARNING: GATK recommends batches of '
                                  '30 to 200 samples for joint genotyping. '
                                  'The current cohort has %d samples.' % num_samples)

    shared_files = Job.wrapJobFn(download_shared_files, config).encapsulate()
    job.addChild(shared_files)

    if config.preprocess_only:
        for sample in samples:
            shared_files.addChildJobFn(prepare_bam,
                                       sample.uuid,
                                       sample.url,
                                       shared_files.rv(),
                                       paired_url=sample.paired_url,
                                       rg_line=sample.rg_line)
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
    Runs the GATK best practices pipeline for germline SNP and INDEL discovery.

    Steps in Pipeline
    0: Generate and preprocess BAM
        - Uploads processed BAM to output directory
    1: Call Variants using HaplotypeCaller
        - Uploads GVCF
    2: Genotype VCF
        - Uploads VCF
    3: Filter Variants using either "hard filters" or VQSR
        - Uploads filtered VCF

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list[GermlineSample] samples: List of GermlineSample namedtuples
    :param Namespace config: Input parameters and reference FileStoreIDs
        Requires the following config attributes:
        config.genome_fasta         FilesStoreID for reference genome fasta file
        config.genome_fai           FilesStoreID for reference genome fasta index file
        config.genome_dict          FilesStoreID for reference genome sequence dictionary file
        config.cores                Number of cores for each job
        config.xmx                  Java heap size in bytes
        config.suffix               Suffix added to output filename
        config.output_dir           URL or local path to output directory
        config.ssec                 Path to key file for SSE-C encryption
        config.joint_genotype       If True, then joint genotype and filter cohort
        config.hc_output            URL or local path to HaplotypeCaller output for testing
    :return: Dictionary of filtered VCF FileStoreIDs
    :rtype: dict
    """
    require(len(samples) > 0, 'No samples were provided!')

    # Get total size of genome reference files. This is used for configuring disk size.
    genome_ref_size = config.genome_fasta.size + config.genome_fai.size + config.genome_dict.size

    # 0: Generate processed BAM and BAI files for each sample
    # group preprocessing and variant calling steps in empty Job instance
    group_bam_jobs = Job()
    gvcfs = {}
    for sample in samples:
        # 0: Generate processed BAM and BAI files for each sample
        get_bam = group_bam_jobs.addChildJobFn(prepare_bam,
                                               sample.uuid,
                                               sample.url,
                                               config,
                                               paired_url=sample.paired_url,
                                               rg_line=sample.rg_line)

        # 1: Generate per sample gvcfs {uuid: gvcf_id}
        # The HaplotypeCaller disk requirement depends on the input bam, bai, the genome reference
        # files, and the output GVCF file. The output GVCF is smaller than the input BAM file.
        hc_disk = PromisedRequirement(lambda bam, bai, ref_size:
                                      2 * bam.size + bai.size + ref_size,
                                      get_bam.rv(0),
                                      get_bam.rv(1),
                                      genome_ref_size)

        get_gvcf = get_bam.addFollowOnJobFn(gatk_haplotype_caller,
                                            get_bam.rv(0),
                                            get_bam.rv(1),
                                            config.genome_fasta, config.genome_fai, config.genome_dict,
                                            annotations=config.annotations,
                                            cores=config.cores,
                                            disk=hc_disk,
                                            memory=config.xmx,
                                            hc_output=config.hc_output)
        # Store cohort GVCFs in dictionary
        gvcfs[sample.uuid] = get_gvcf.rv()

        # Upload individual sample GVCF before genotyping to a sample specific output directory
        vqsr_name = '{}{}.g.vcf'.format(sample.uuid, config.suffix)
        get_gvcf.addChildJobFn(output_file_job,
                               vqsr_name,
                               get_gvcf.rv(),
                               os.path.join(config.output_dir, sample.uuid),
                               s3_key_path=config.ssec,
                               disk=PromisedRequirement(lambda x: x.size, get_gvcf.rv()))

    # VQSR requires many variants in order to train a decent model. GATK recommends a minimum of
    # 30 exomes or one large WGS sample:
    # https://software.broadinstitute.org/gatk/documentation/article?id=3225

    filtered_vcfs = {}
    if config.joint_genotype:
        # Need to configure joint genotype in a separate function to resolve promises
        filtered_vcfs = group_bam_jobs.addFollowOnJobFn(joint_genotype_and_filter,
                                                        gvcfs,
                                                        config).rv()

    # If not joint genotyping, then iterate over cohort and genotype and filter individually.
    else:
        for uuid, gvcf_id in gvcfs.iteritems():
            filtered_vcfs[uuid] = group_bam_jobs.addFollowOnJobFn(genotype_and_filter,
                                                                  {uuid: gvcf_id},
                                                                  config).rv()

    job.addChild(group_bam_jobs)
    return filtered_vcfs


def joint_genotype_and_filter(job, gvcfs, config):
    """
    Checks for enough disk space for joint genotyping, then calls the genotype and filter pipeline function.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param dict gvcfs: Dictionary of GVCFs {Sample ID: FileStoreID}
    :param Namespace config: Input parameters and reference FileStoreIDs
        Requires the following config attributes:
        config.genome_fasta         FilesStoreID for reference genome fasta file
        config.genome_fai           FilesStoreID for reference genome fasta index file
        config.genome_dict          FilesStoreID for reference genome sequence dictionary file
        config.available_disk       Total available disk space
    :returns: FileStoreID for the joint genotyped and filtered VCF file
    :rtype: str
    """
    # Get the total size of genome reference files
    genome_ref_size = config.genome_fasta.size + config.genome_fai.size + config.genome_dict.size

    # Require at least 2.5x the sum of the individual GVCF files
    cohort_size = sum(gvcf.size for gvcf in gvcfs.values())
    require(int(2.5 * cohort_size + genome_ref_size) < config.available_disk,
            'There is not enough disk space to joint '
            'genotype samples:\n{}'.format('\n'.join(gvcfs.keys())))

    job.fileStore.logToMaster('Merging cohort into a single GVCF file')

    return job.addChildJobFn(genotype_and_filter, gvcfs, config).rv()


def genotype_and_filter(job, gvcfs, config):
    """
    Genotypes one or more GVCF files and runs either the VQSR or hard filtering pipeline. Uploads the genotyped VCF file
    to the config output directory.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param dict gvcfs: Dictionary of GVCFs {Sample ID: FileStoreID}
    :param Namespace config: Input parameters and shared FileStoreIDs
        Requires the following config attributes:
        config.genome_fasta         FilesStoreID for reference genome fasta file
        config.genome_fai           FilesStoreID for reference genome fasta index file
        config.genome_dict          FilesStoreID for reference genome sequence dictionary file
        config.suffix               Suffix added to output filename
        config.output_dir           URL or local path to output directory
        config.ssec                 Path to key file for SSE-C encryption
        config.cores                Number of cores for each job
        config.xmx                  Java heap size in bytes
        config.unsafe_mode          If True, then run GATK tools in UNSAFE mode
    :return: FileStoreID for genotyped and filtered VCF file
    :rtype: str
    """
    # Get the total size of the genome reference
    genome_ref_size = config.genome_fasta.size + config.genome_fai.size + config.genome_dict.size

    # GenotypeGVCF disk requirement depends on the input GVCF, the genome reference files, and
    # the output VCF file. The output VCF is smaller than the input GVCF.
    genotype_gvcf_disk = PromisedRequirement(lambda gvcf_ids, ref_size:
                                             2 * sum(gvcf_.size for gvcf_ in gvcf_ids) + ref_size,
                                             gvcfs.values(),
                                             genome_ref_size)

    genotype_gvcf = job.addChildJobFn(gatk_genotype_gvcfs,
                                      gvcfs,
                                      config.genome_fasta,
                                      config.genome_fai,
                                      config.genome_dict,
                                      annotations=config.annotations,
                                      unsafe_mode=config.unsafe_mode,
                                      cores=config.cores,
                                      disk=genotype_gvcf_disk,
                                      memory=config.xmx)

    # Determine if output GVCF has multiple samples
    if len(gvcfs) == 1:
        uuid = gvcfs.keys()[0]
    else:
        uuid = 'joint_genotyped'

    genotyped_filename = '%s.genotyped%s.vcf' % (uuid, config.suffix)
    genotype_gvcf.addChildJobFn(output_file_job,
                                genotyped_filename,
                                genotype_gvcf.rv(),
                                os.path.join(config.output_dir, uuid),
                                s3_key_path=config.ssec,
                                disk=PromisedRequirement(lambda x: x.size, genotype_gvcf.rv()))

    if config.run_vqsr:
        if not config.joint_genotype:
            job.fileStore.logToMaster('WARNING: Running VQSR without joint genotyping.')
        joint_genotype_vcf = genotype_gvcf.addFollowOnJobFn(vqsr_pipeline,
                                                            uuid,
                                                            genotype_gvcf.rv(),
                                                            config)
    else:
        joint_genotype_vcf = genotype_gvcf.addFollowOnJobFn(hard_filter_pipeline,
                                                            uuid,
                                                            genotype_gvcf.rv(),
                                                            config)
    return joint_genotype_vcf.rv()


def annotate_vcfs(job, vcfs, config):
    """
    Runs Oncotator for a group of VCF files. Each sample is annotated individually.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param dict vcfs: Dictionary of VCF FileStoreIDs {Sample identifier: FileStoreID}
    :param Namespace config: Input parameters and shared FileStoreIDs
        Requires the following config attributes:
        config.oncotator_db         FileStoreID to Oncotator database
        config.suffix               Suffix added to output filename
        config.output_dir           URL or local path to output directory
        config.ssec                 Path to key file for SSE-C encryption
        config.cores                Number of cores for each job
        config.xmx                  Java heap size in bytes
    """
    job.fileStore.logToMaster('Running Oncotator on the following samples:\n%s' % '\n'.join(vcfs.keys()))
    for uuid, vcf_id in vcfs.iteritems():
        # The Oncotator disk requirement depends on the input VCF, the Oncotator database
        # and the output VCF. The annotated VCF will be significantly larger than the input VCF.
        onco_disk = PromisedRequirement(lambda vcf, db: 3 * vcf.size + db.size,
                                        vcf_id,
                                        config.oncotator_db)

        annotated_vcf = job.addChildJobFn(run_oncotator,
                                          vcf_id,
                                          config.oncotator_db,
                                          disk=onco_disk,
                                          cores=config.cores,
                                          memory=config.xmx)

        output_dir = os.path.join(config.output_dir, uuid)
        filename = '{}.oncotator{}.vcf'.format(uuid, config.suffix)
        annotated_vcf.addChildJobFn(output_file_job,
                                    filename,
                                    annotated_vcf.rv(),
                                    output_dir,
                                    s3_key_path=config.ssec,
                                    disk=PromisedRequirement(lambda x: x.size, annotated_vcf.rv()))


# Pipeline convenience functions


def parse_manifest(path_to_manifest):
    """
    Parses manifest file for Toil Germline Pipeline

    :param str path_to_manifest: Path to sample manifest file
    :return: List of GermlineSample namedtuples
    :rtype: list[GermlineSample]
    """
    bam_re = r"^(?P<uuid>\S+)\s(?P<url>\S+[bsc][r]?am)"
    fq_re = r"^(?P<uuid>\S+)\s(?P<url>\S+)\s(?P<paired_url>\S+)?\s?(?P<rg_line>@RG\S+)"
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
                paired_url = None
                rg_line = None
                require('.bam' in url.lower(),
                        'Expected .bam extension:\n{}:\t{}'.format(uuid, url))
            elif fastq_match:
                uuid = fastq_match.group('uuid')
                url = fastq_match.group('url')
                paired_url = fastq_match.group('paired_url')
                rg_line = fastq_match.group('rg_line')
                require('.fq' in url.lower() or '.fastq' in url.lower(),
                        'Expected .fq extension:\n{}:\t{}'.format(uuid, url))
            else:
                raise ValueError('Could not parse entry in manifest: %s\n%s' % (f.name, line))
            # Checks that URL has a scheme
            require(urlparse(url).scheme, 'Invalid URL passed for {}'.format(url))
            samples.append(GermlineSample(uuid, url, paired_url, rg_line))
    return samples


def download_shared_files(job, config):
    """
    Downloads shared reference files for Toil Germline pipeline

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
        shared_files |= {'g1k_indel', 'mills', 'dbsnp'}
    if config.run_vqsr:
        shared_files |= {'g1k_snp', 'mills', 'dbsnp', 'hapmap', 'omni'}
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
    Creates a genome fasta index and sequence dictionary file if not already present in the pipeline config.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Pipeline configuration options and shared files.
                             Requires FileStoreID for genome fasta file as config.genome_fasta
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
    Prepares BAM file for Toil germline pipeline.

    Steps in pipeline
    0: Download and align BAM or FASTQ sample
    1: Sort BAM
    2: Index BAM
    3: Run GATK preprocessing pipeline (Optional)
        - Uploads preprocessed BAM to output directory

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: Unique identifier for the sample
    :param str url: URL or local path to BAM file or FASTQs
    :param Namespace config: Configuration options for pipeline
        Requires the following config attributes:
        config.genome_fasta         FilesStoreID for reference genome fasta file
        config.genome_fai           FilesStoreID for reference genome fasta index file
        config.genome_dict          FilesStoreID for reference genome sequence dictionary file
        config.g1k_indel            FileStoreID for 1000G INDEL resource file
        config.mills                FileStoreID for Mills resource file
        config.dbsnp                FileStoreID for dbSNP resource file
        config.suffix               Suffix added to output filename
        config.output_dir           URL or local path to output directory
        config.ssec                 Path to key file for SSE-C encryption
        config.cores                Number of cores for each job
        config.xmx                  Java heap size in bytes
    :param str|None paired_url: URL or local path to paired FASTQ file, default is None
    :param str|None rg_line: RG line for BWA alignment (i.e. @RG\tID:foo\tSM:bar), default is None
    :return: BAM and BAI FileStoreIDs
    :rtype: tuple
    """
    # 0: Align FASTQ or realign BAM
    if config.run_bwa:
        get_bam = job.wrapJobFn(setup_and_run_bwakit,
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
        # The samtools sort disk requirement depends on the input bam, the tmp files, and the
        # sorted output bam.
        sorted_bam_disk = PromisedRequirement(lambda bam: 3 * bam.size, get_bam.rv())
        sorted_bam = get_bam.addChildJobFn(run_samtools_sort,
                                           get_bam.rv(),
                                           cores=config.cores,
                                           disk=sorted_bam_disk)

    # 2: Index BAM
    # The samtools index disk requirement depends on the input bam and the output bam index
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
                                   config.g1k_indel,
                                   config.mills,
                                   config.dbsnp,
                                   memory=config.xmx,
                                   cores=config.cores).encapsulate()
        sorted_bam.addChild(preprocess)
        index_bam.addChild(preprocess)

        # Update output BAM promises
        output_bam_promise = preprocess.rv(0)
        output_bai_promise = preprocess.rv(1)

        # Save processed BAM
        output_dir = os.path.join(config.output_dir, uuid)
        filename = '{}.preprocessed{}.bam'.format(uuid, config.suffix)
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


def setup_and_run_bwakit(job, uuid, url, rg_line, config, paired_url=None):
    """
    Downloads and runs bwakit for BAM or FASTQ files

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: Unique sample identifier
    :param str url: FASTQ or BAM file URL. BAM alignment URL must have .bam extension.
    :param Namespace config: Input parameters and shared FileStoreIDs
        Requires the following config attributes:
        config.genome_fasta         FilesStoreID for reference genome fasta file
        config.genome_fai           FilesStoreID for reference genome fasta index file
        config.cores                Number of cores for each job
        config.trim                 If True, trim adapters using bwakit
        config.amb                  FileStoreID for BWA index file prefix.amb
        config.ann                  FileStoreID for BWA index file prefix.ann
        config.bwt                  FileStoreID for BWA index file prefix.bwt
        config.pac                  FileStoreID for BWA index file prefix.pac
        config.sa                   FileStoreID for BWA index file prefix.sa
        config.alt                  FileStoreID for alternate contigs file or None
    :param str|None paired_url: URL to paired FASTQ
    :param str|None rg_line: Read group line (i.e. @RG\tID:foo\tSM:bar)
    :return: BAM FileStoreID
    :rtype: str
    """
    bwa_config = deepcopy(config)
    bwa_config.uuid = uuid
    bwa_config.rg_line = rg_line

    # bwa_alignment uses a different naming convention
    bwa_config.ref = config.genome_fasta
    bwa_config.fai = config.genome_fai

    # Determine if sample is a FASTQ or BAM file using the file extension
    basename, ext = os.path.splitext(url)
    ext = ext.lower()
    if ext == '.gz':
        _, ext = os.path.splitext(basename)
        ext = ext.lower()

    # The pipeline currently supports FASTQ and BAM files
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

    # If the extension is for a BAM  file, then configure bwakit to realign the BAM file.
    if ext == '.bam':
        bwa_config.bam = input1.rv()
    else:
        bwa_config.r1 = input1.rv()

    # Download the paired FASTQ URL
    if paired_url:
        input2 = job.addChildJobFn(download_url_job,
                                   paired_url,
                                   name='file2',
                                   s3_key_path=config.ssec,
                                   disk=config.file_size)
        samples.append(input2.rv())
        bwa_config.r2 = input2.rv()

    # The bwakit disk requirement depends on the size of the input files and the index
    # Take the sum of the input files and scale it by a factor of 4
    bwa_index_size = sum([getattr(config, index_file).size
                          for index_file in ['amb', 'ann', 'bwt', 'pac', 'sa', 'alt']
                          if getattr(config, index_file, None) is not None])

    bwakit_disk = PromisedRequirement(lambda lst, index_size:
                                      int(4 * sum(x.size for x in lst) + index_size),
                                      samples,
                                      bwa_index_size)

    return job.addFollowOnJobFn(run_bwakit,
                                bwa_config,
                                sort=False,         # BAM files are sorted later in the pipeline
                                trim=config.trim,
                                cores=config.cores,
                                disk=bwakit_disk).rv()


def gatk_haplotype_caller(job,
                          bam, bai,
                          ref, fai, ref_dict,
                          annotations=None,
                          emit_threshold=10.0, call_threshold=30.0,
                          unsafe_mode=False,
                          hc_output=None):
    """
    Uses GATK HaplotypeCaller to identify SNPs and INDELs. Outputs variants in a Genomic VCF file.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str bam: FileStoreID for BAM file
    :param str bai: FileStoreID for BAM index file
    :param str ref: FileStoreID for reference genome fasta file
    :param str ref_dict: FileStoreID for reference sequence dictionary file
    :param str fai: FileStoreID for reference fasta index file
    :param list[str] annotations: List of GATK variant annotations, default is None
    :param float emit_threshold: Minimum phred-scale confidence threshold for a variant to be emitted, default is 10.0
    :param float call_threshold: Minimum phred-scale confidence threshold for a variant to be called, default is 30.0
    :param bool unsafe_mode: If True, runs gatk UNSAFE mode: "-U ALLOW_SEQ_DICT_INCOMPATIBILITY"
    :param str hc_output: URL or local path to pre-cooked VCF file, default is None
    :return: FileStoreID for GVCF file
    :rtype: str
    """
    job.fileStore.logToMaster('Running GATK HaplotypeCaller')

    inputs = {'genome.fa': ref,
              'genome.fa.fai': fai,
              'genome.dict': ref_dict,
              'input.bam': bam,
              'input.bam.bai': bai}

    work_dir = job.fileStore.getLocalTempDir()
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

    if unsafe_mode:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    if annotations:
        for annotation in annotations:
            command.extend(['-A', annotation])

    # Uses docker_call mock mode to replace output with hc_output file
    outputs = {'output.g.vcf': hc_output}
    docker_call(work_dir=work_dir,
                env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters=command,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=True if outputs['output.g.vcf'] else False)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.g.vcf'))


def main():
    """
    GATK germline pipeline with variant filtering and annotation.
    """
    # Define Parser object and add to jobTree
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

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
                                 '"generate-manifest".\nDefault value: "%(default)s".')
    parser_run.add_argument('--sample',
                            default=None,
                            nargs=2,
                            type=str,
                            help='Input sample identifier and BAM file URL or local path')
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
    elif options.command == 'run':
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program)),
                    program + ' must be installed on every node.'.format(program))

        require(os.path.exists(options.config), '{} not found. Please run "generate-config"'.format(options.config))

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
                           'snp_filter_annotations',
                           'indel_filter_annotations',
                           'preprocess',
                           'preprocess_only',
                           'run_vqsr',
                           'joint_genotype',
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
            # Check that essential VQSR parameters are present
            vqsr_fields = {'g1k_snp', 'mills', 'dbsnp', 'hapmap', 'omni'}
            require(input_fields > vqsr_fields,
                    'Missing parameters for VQSR:\n{}'.format(', '.join(vqsr_fields - input_fields)))

        # Check that hard filtering parameters are present. If only running preprocessing steps, then we do
        # not need filtering information.
        elif not inputs['preprocess_only']:
            hard_filter_fields = {'snp_filter_name', 'snp_filter_expression',
                                  'indel_filter_name', 'indel_filter_expression'}
            require(input_fields > hard_filter_fields,
                    'Missing parameters for hard filtering:\n{}'.format(', '.join(hard_filter_fields - input_fields)))

            # Check for falsey hard filtering parameters
            for hard_filter_field in hard_filter_fields:
                require(inputs[hard_filter_field], 'Missing %s value for hard filtering, '
                                                   'got %s.' % (hard_filter_field, inputs[hard_filter_field]))

        # Set resource parameters
        inputs['xmx'] = human2bytes(inputs['xmx'])
        inputs['file_size'] = human2bytes(inputs['file_size'])
        inputs['cores'] = int(inputs['cores'])

        inputs['annotations'] = set(inputs['snp_filter_annotations'] + inputs['indel_filter_annotations'])

        # HaplotypeCaller test data for testing
        inputs['hc_output'] = inputs.get('hc_output', None)

        # It is a toil-scripts convention to store input parameters in a Namespace object
        config = argparse.Namespace(**inputs)

        root = Job.wrapJobFn(run_gatk_germline_pipeline, samples, config)
        Job.Runner.startToil(root, options)


if __name__ == '__main__':
    main()
