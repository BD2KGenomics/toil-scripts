#!/usr/bin/env python2.7
import argparse
from glob import glob
import multiprocessing
import os
import subprocess
import sys
import tarfile
import textwrap
from contextlib import closing
from urlparse import urlparse

import yaml
from toil.job import Job

from toil_scripts.lib import require
from toil_scripts.lib.files import tarball_files, mkdir_p, move_files
from toil_scripts.lib.jobs import map_job
from toil_scripts.lib.programs import which, docker_call
from toil_scripts.lib.urls import download_url_job, s3am_upload


# Start of Job Functions
def download_shared_files(job, samples, config):
    """
    Downloads files shared by all samples in the pipeline

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param list[list] samples: A nested list of samples containing sample information
    """
    job.fileStore.logToMaster('Downloaded shared files')
    file_names = ['reference', 'phase', 'mills', 'dbsnp', 'cosmic']
    urls = [config.reference, config.phase, config.mills, config.dbsnp, config.cosmic]
    for name, url in zip(file_names, urls):
        if url:
            vars(config)[name] = job.addChildJobFn(download_url_job, url=url, s3_key_path=config.ssec).rv()
    job.addFollowOnJobFn(reference_preprocessing, samples, config)


def reference_preprocessing(job, samples, config):
    """
    Spawn the jobs that create index and dict file for reference

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param list[list] samples: A nested list of samples containing sample information
    """
    job.fileStore.logToMaster('Processed reference files')
    config.fai = job.addChildJobFn(samtools_reference_index, config.reference).rv()
    config.dict = job.addChildJobFn(picard_reference_dict, config.reference).rv()
    job.addFollowOnJobFn(map_job, download_sample, samples, config)


def samtools_reference_index(job, ref_id):
    """
    Use Samtools to create reference index file

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str ref_id: FileStoreID for the reference genome
    :return: FileStoreID for reference index
    :rtype: str
    """
    job.fileStore.logToMaster('Created reference index')
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta'))
    command = ['faidx', 'ref.fasta']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.fasta.fai'))


def picard_reference_dict(job, ref_id):
    """
    Use Picard-tools to create reference dictionary

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str ref_id: FileStoreID for the reference genome
    :return: FileStoreID for reference dictionary
    :rtype: str
    """
    job.fileStore.logToMaster('Created reference dictionary')
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta'))
    command = ['CreateSequenceDictionary', 'R=ref.fasta', 'O=ref.dict']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.dict'))


def download_sample(job, sample, config):
    """
    Download sample and store sample specific attributes

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list sample: Contains uuid, normal URL, and tumor URL
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param Namespace shared_ids: Argparse Namespace object containing fileStoreIDs of globally shared inputs
    """
    # Create copy of config that is sample specific
    config = argparse.Namespace(**vars(config))
    uuid, normal_url, tumor_url = sample
    job.fileStore.logToMaster('Downloaded sample: {}'.format(uuid))
    config.uuid = uuid
    config.normal = normal_url
    config.tumor = tumor_url
    config.cores = int(multiprocessing.cpu_count())
    disk = '1G' if config.ci_test else '20G'
    # Download sample bams and launch pipeline
    config.normal_bam = job.addChildJobFn(download_url_job, url=config.normal, s3_key_path=config.ssec,
                                          cghub_key_path=config.gtkey, disk=disk).rv()
    config.tumor_bam = job.addChildJobFn(download_url_job, url=config.tumor, s3_key_path=config.ssec,
                                         cghub_key_path=config.gtkey, disk=disk).rv()
    job.addFollowOnJobFn(index_bams, config)


def index_bams(job, config):
    """
    Convenience job for handling bam indexing to make the static workflow declaration cleaner

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    """
    job.fileStore.logToMaster('Indexed sample BAMS')
    disk = '1G' if config.ci_test else '20G'
    config.normal_bai = job.addChildJobFn(index_bam, config.normal_bam, cores=1, disk=disk).rv()
    config.tumor_bai = job.addChildJobFn(index_bam, config.tumor_bam, cores=1, disk=disk).rv()
    job.addFollowOnJobFn(pre_processing_declaration, config)


def index_bam(job, bam_id):
    """
    Runs samtools index to create (.bai) files

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str bam_id: FileStoreID of the bam file
    :return: BAM index FileStoreID
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'sample.bam'))
    # Call: index the bam
    parameters = ['index', '/data/sample.bam']
    docker_call(work_dir=work_dir, parameters=parameters,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.bam.bai'))


def pre_processing_declaration(job, config):
    if config.preprocessing:
        disk = '1G' if config.ci_test else '20G'
        mem = '2G' if config.ci_test else '10G'
        processed_normal = job.wrapJobFn(start_bam_preprocessing, config.cores, config.normal_bam, config.normal_bai,
                                         config.reference, config.dict, config.fai, config.phase, config.mills,
                                         config.dbsnp, mem, cores=1, memory=mem, disk=disk).encapsulate()
        processed_tumor = job.wrapJobFn(start_bam_preprocessing, config.cores, config.tumor_bam, config.tumor_bai,
                                        config.reference, config.dict, config.fai, config.phase, config.mills,
                                        config.dbsnp, mem, cores=1, memory=mem, disk=disk).encapsulate()
        job.addChild(processed_normal)
        job.addChild(processed_tumor)
        normal_bam = processed_normal.rv()
        tumor_bam = processed_tumor.rv()
        normal_bai = processed_normal.addChildJobFn(index_bam, normal_bam, cores=1, disk=disk).rv()
        tumor_bai = processed_tumor.addChildJobFn(index_bam, tumor_bam, cores=1, disk=disk).rv()
    else:
        normal_bam, normal_bai = config.normal_bam, config.normal_bai
        tumor_bam, tumor_bai = config.tumor_bam, config.tumor_bai
    job.addFollowOnJobFn(static_workflow_declaration, config, normal_bam, normal_bai, tumor_bam, tumor_bai)


def static_workflow_declaration(job, config, normal_bam, normal_bai, tumor_bam, tumor_bai):
    """
    Statically declare workflow so sections can be modularly repurposed

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str normal_bam: Normal BAM FileStoreID
    :param str normal_bai: Normal BAM index FileStoreID
    :param str tumor_bam: Tumor BAM FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    """
    # Mutation and indel tool wiring
    memory = '1G' if config.ci_test else '10G'
    disk = '1G' if config.ci_test else '75G'
    mutect_results, pindel_results, muse_results = None, None, None
    mutect = job.wrapJobFn(run_mutect, normal_bam, normal_bai, tumor_bam, tumor_bai, config.reference, config.dict,
                           config.fai, config.cosmic, config.dbsnp, memory=memory, disk=disk)
    pindel = job.wrapJobFn(run_pindel, config.cores, normal_bam, normal_bai, tumor_bam, tumor_bai, config.reference,
                           config.fai, cores=config.cores,  memory=memory, disk=disk)
    muse = job.wrapJobFn(run_muse, config.cores, normal_bam, normal_bai, tumor_bam, tumor_bai, config.reference,
                         config.dict, config.fai, config.dbsnp, cores=config.cores, memory=memory, disk=disk)
    if config.run_mutect:
        job.addChild(mutect)
        mutect_results = mutect.rv()
    if config.run_pindel:
        job.addChild(pindel)
        pindel_results = pindel.rv()
    if config.run_muse:
        job.addChild(muse)
        muse_results = muse.rv()
    # Pass tool results (whether None or a promised return value) to consolidation step
    consolidation = job.wrapJobFn(consolidate_output, config, mutect_results, pindel_results, muse_results)
    job.addFollowOn(consolidation)


def start_bam_preprocessing(job, cores, bam, bai, ref, ref_dict, fai, phase, mills, dbsnp, mem='10G'):
    """
    Creates intervals file needed for indel realignment

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on a worker node
    :param str bam: Sample BAM FileStoreID
    :param str bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str phase: Phase VCF FileStoreID
    :param str mills: Mills VCF FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :return: FileStoreID for the processed bam
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [ref, fai, ref_dict, bam, bai, phase, mills]
    file_names = ['ref.fasta', 'ref.fasta.fai', 'ref.dict', 'sample.bam', 'sample.bam.bai', 'phase.vcf', 'mills.vcf']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: GATK -- RealignerTargetCreator
    parameters = ['-T', 'RealignerTargetCreator',
                  '-nt', str(cores),
                  '-R', '/data/ref.fasta',
                  '-I', '/data/sample.bam',
                  '-known', '/data/phase.vcf',
                  '-known', '/data/mills.vcf',
                  '--downsampling_type', 'NONE',
                  '-o', '/data/sample.intervals']
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                work_dir=work_dir, parameters=parameters, env=dict(JAVA_OPTS='-Xmx{}'.format(mem)))
    # Write to fileStore
    intervals = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.intervals'))
    return job.addChildJobFn(indel_realignment, cores, intervals, bam, bai, ref, ref_dict, fai, phase,
                             mills, dbsnp, mem, cores=1, memory=mem, disk='25G').rv()


def indel_realignment(job, cores, intervals, bam, bai, ref, ref_dict, fai, phase, mills, dbsnp, mem):
    """
    Creates realigned bams using the intervals file from previous step

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on a worker node
    :param str intervals: Indel interval FileStoreID
    :param str bam: Sample BAM FileStoreID
    :param str bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str phase: Phase VCF FileStoreID
    :param str mills: Mills VCF FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :return: FileStoreID for the processed bam
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [ref, fai, ref_dict, intervals, bam, bai, phase, mills]
    file_names = ['ref.fasta', 'ref.fasta.fai', 'ref.dict', 'sample.intervals',
                  'sample.bam', 'sample.bam.bai', 'phase.vcf', 'mills.vcf']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: GATK -- IndelRealigner
    parameters = ['-T', 'IndelRealigner',
                  '-R', '/data/ref.fasta',
                  '-I', '/data/sample.bam',
                  '-known', '/data/phase.vcf',
                  '-known', '/data/mills.vcf',
                  '-targetIntervals', '/data/sample.intervals',
                  '--downsampling_type', 'NONE',
                  '-maxReads', str(720000),
                  '-maxInMemory', str(5400000),
                  '-o', '/data/sample.indel.bam']
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                work_dir=work_dir, parameters=parameters, env=dict(JAVA_OPTS='-Xmx{}'.format(mem)))
    # Write to fileStore
    indel_bam = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.indel.bam'))
    indel_bai = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.indel.bai'))
    return job.addChildJobFn(base_recalibration, cores, indel_bam, indel_bai, ref, ref_dict, fai, dbsnp, mem,
                             cores=cores, memory=mem, disk='25G').rv()


def base_recalibration(job, cores, indel_bam, indel_bai, ref, ref_dict, fai, dbsnp, mem):
    """
    Creates recal table used in Base Quality Score Recalibration

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on a worker node
    :param str indel_bam: Indel interval FileStoreID
    :param str indel_bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :return: FileStoreID for the processed bam
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [ref, fai, ref_dict, indel_bam, indel_bai, dbsnp]
    file_names = ['ref.fasta', 'ref.fasta.fai', 'ref.dict', 'sample.indel.bam', 'sample.indel.bai', 'dbsnp.vcf']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: GATK -- IndelRealigner
    parameters = ['-T', 'BaseRecalibrator',
                  '-nct', str(cores),
                  '-R', '/data/ref.fasta',
                  '-I', '/data/sample.indel.bam',
                  '-knownSites', '/data/dbsnp.vcf',
                  '-o', '/data/sample.recal.table']
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                work_dir=work_dir, parameters=parameters, env=dict(JAVA_OPTS='-Xmx{}'.format(mem)))
    # Write output to file store
    table = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.recal.table'))
    return job.addChildJobFn(print_reads, cores, table, indel_bam, indel_bai, ref, ref_dict, fai, mem,
                             cores=cores, memory=mem, disk='25G').rv()


def print_reads(job, cores, table, indel_bam, indel_bai, ref, ref_dict, fai, mem):
    """
    Creates BAM that has had the base quality scores recalibrated

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on host node
    :param str table: Recalibration table FileStoreID
    :param str indel_bam: Indel interval FileStoreID
    :param str indel_bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :return: FileStoreID for the processed bam
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [ref, fai, ref_dict, table, indel_bam, indel_bai]
    file_names = ['ref.fasta', 'ref.fasta.fai', 'ref.dict', 'sample.recal.table',
                  'sample.indel.bam', 'sample.indel.bai']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: GATK -- PrintReads
    parameters = ['-T', 'PrintReads',
                  '-nct', str(cores),
                  '-R', '/data/ref.fasta',
                  '--emit_original_quals',
                  '-I', '/data/sample.indel.bam',
                  '-BQSR', '/data/sample.recal.table',
                  '-o', '/data/sample.bqsr.bam']
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.4--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                work_dir=work_dir, parameters=parameters, env=dict(JAVA_OPTS='-Xmx{}'.format(mem)))
    # Write ouptut to file store
    bam_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.bqsr.bam'))
    bai_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.bqsr.bai'))
    return bam_id


def run_mutect(job, normal_bam, normal_bai, tumor_bam, tumor_bai, ref, ref_dict, fai, cosmic, dbsnp):
    """
    Calls MuTect to perform variant analysis

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str normal_bam: Normal BAM FileStoreID
    :param str normal_bai: Normal BAM index FileStoreID
    :param str tumor_bam: Tumor BAM FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str cosmic: Cosmic VCF FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :return: MuTect output (tarball) FileStoreID
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [normal_bam, normal_bai, tumor_bam, tumor_bai, ref, fai, ref_dict, cosmic, dbsnp]
    file_names = ['normal.bam', 'normal.bai', 'tumor.bam', 'tumor.bai', 'ref.fasta',
                  'ref.fasta.fai', 'ref.dict', 'cosmic.vcf', 'dbsnp.vcf']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: MuTect
    parameters = ['--analysis_type', 'MuTect',
                  '--reference_sequence', 'ref.fasta',
                  '--cosmic', '/data/cosmic.vcf',
                  '--dbsnp', '/data/dbsnp.vcf',
                  '--input_file:normal', '/data/normal.bam',
                  '--input_file:tumor', '/data/tumor.bam',
                  '--tumor_lod', str(10),
                  '--initial_tumor_lod', str(4.0),
                  '--out', 'mutect.out',
                  '--coverage_file', 'mutect.cov',
                  '--vcf', 'mutect.vcf']
    docker_call(work_dir=work_dir, parameters=parameters,
                tool='quay.io/ucsc_cgl/mutect:1.1.7--e8bf09459cf0aecb9f55ee689c2b2d194754cbd3')
    # Write output to file store
    output_file_names = ['mutect.vcf', 'mutect.cov', 'mutect.out']
    output_file_paths = [os.path.join(work_dir, x) for x in output_file_names]
    tarball_files('mutect.tar.gz', file_paths=output_file_paths, output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'mutect.tar.gz'))


def run_pindel(job, cores, normal_bam, normal_bai, tumor_bam, tumor_bai, ref, fai):
    """
    Calls Pindel to compute indels / deletions

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on host node
    :param str normal_bam: Normal BAM FileStoreID
    :param str normal_bai: Normal BAM index FileStoreID
    :param str tumor_bam: Tumor BAM FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str fai: Reference index FileStoreID
    :return: Pindel output (tarball) FileStoreID
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [normal_bam, normal_bai, tumor_bam, tumor_bai, ref, fai]
    file_names = ['normal.bam', 'normal.bai', 'tumor.bam', 'tumor.bai', 'ref.fasta', 'ref.fasta.fai']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Create Pindel config
    with open(os.path.join(work_dir, 'pindel-config.txt'), 'w') as f:
        for bam in ['normal.bam', 'tumor.bam']:
            f.write('/data/{} {} {}\n'.format(bam, get_mean_insert_size(work_dir, bam), os.path.splitext(bam)[0]))
    # Call: Pindel
    parameters = ['-f', '/data/ref.fasta',
                  '-i', '/data/pindel-config.txt',
                  '--number_of_threads', str(cores),
                  '--minimum_support_for_event', '3',
                  '--report_long_insertions', 'true',
                  '--report_breakpoints', 'true',
                  '-o', 'pindel']
    docker_call(tool='quay.io/ucsc_cgl/pindel:0.2.5b6--4e8d1b31d4028f464b3409c6558fb9dfcad73f88',
                work_dir=work_dir, parameters=parameters)
    # Collect output files and write to file store
    output_files = glob(os.path.join(work_dir, 'pindel*'))
    tarball_files('pindel.tar.gz', file_paths=output_files, output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'pindel.tar.gz'))


def run_muse(job, cores, normal_bam, normal_bai, tumor_bam, tumor_bai, ref, ref_dict, fai, dbsnp):
    """
    Calls MuSe to find variants

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on host node
    :param str normal_bam: Normal BAM FileStoreID
    :param str normal_bai: Normal BAM index FileStoreID
    :param str tumor_bam: Tumor BAM FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference genome dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :return: MuSe output (tarball) FileStoreID
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [normal_bam, normal_bai, tumor_bam, tumor_bai, ref, ref_dict, fai, dbsnp]
    file_names = ['normal.bam', 'normal.bai', 'tumor.bam', 'tumor.bai',
                  'ref.fasta', 'ref.dict', 'ref.fasta.fai', 'dbsnp.vcf']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: MuSE
    parameters = ['--mode', 'wxs',
                  '--dbsnp', '/data/dbsnp.vcf',
                  '--fafile', '/data/ref.fasta',
                  '--tumor-bam', '/data/tumor.bam',
                  '--tumor-bam-index', '/data/tumor.bai',
                  '--normal-bam', '/data/normal.bam',
                  '--normal-bam-index', '/data/normal.bai',
                  '--outfile', '/data/muse.vcf',
                  '--cpus', str(cores)]
    docker_call(tool='quay.io/ucsc_cgl/muse:1.0--6add9b0a1662d44fd13bbc1f32eac49326e48562',
                work_dir=work_dir, parameters=parameters)
    # Return fileStore ID
    tarball_files('muse.tar.gz', file_paths=[os.path.join(work_dir, 'muse.vcf')], output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'muse.tar.gz'))


def consolidate_output(job, config, mutect, pindel, muse):
    """
    Combine the contents of separate tarball outputs into one via streaming

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str mutect: MuTect tarball FileStoreID
    :param str pindel: Pindel tarball FileStoreID
    :param str muse: MuSe tarball FileStoreID
    """
    work_dir = job.fileStore.getLocalTempDir()
    mutect_tar, pindel_tar, muse_tar = None, None, None
    if mutect:
        mutect_tar = job.fileStore.readGlobalFile(mutect, os.path.join(work_dir, 'mutect.tar.gz'))
    if pindel:
        pindel_tar = job.fileStore.readGlobalFile(pindel, os.path.join(work_dir, 'pindel.tar.gz'))
    if muse:
        muse_tar = job.fileStore.readGlobalFile(muse, os.path.join(work_dir, 'muse.tar.gz'))
    out_tar = os.path.join(work_dir, config.uuid + '.tar.gz')
    # Consolidate separate tarballs into one as streams (avoids unnecessary untaring)
    tar_list = [x for x in [mutect_tar, pindel_tar, muse_tar] if x is not None]
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar in tar_list:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar is mutect_tar:
                            tarinfo.name = os.path.join(config.uuid, 'mutect', os.path.basename(tarinfo.name))
                        elif tar is pindel_tar:
                            tarinfo.name = os.path.join(config.uuid, 'pindel', os.path.basename(tarinfo.name))
                        else:
                            tarinfo.name = os.path.join(config.uuid, 'muse', os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
    # Move to output directory of selected
    if config.output_dir:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(config.uuid, config.output_dir))
        mkdir_p(config.output_dir)
        move_files(file_paths=[out_tar], output_dir=config.output_dir)
    if config.s3_dir:
        job.fileStore.logToMaster('Uploading {} to S3: {}'.format(config.uuid, config.s3_dir))
        s3am_upload(fpath=out_tar, s3_dir=config.s3_dir, num_cores=config.cores)


def upload_output_to_s3(job, input_args, output_tar):
    """
    Uploads Mutect files to S3 using S3AM

    mutect_id: str      FileStore ID for the mutect.vcf
    input_args: dict    Dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.lstrip('/').split('/')[0]
    bucket_dir = '/'.join(s3_dir.lstrip('/').split('/')[1:])
    # Retrieve VCF file
    job.fileStore.readGlobalFile(output_tar, os.path.join(work_dir, uuid + '.tar.gz'))
    # Upload to S3 via S3AM
    s3am_command = ['s3am',
                    'upload',
                    'file://{}'.format(os.path.join(work_dir, uuid + '.tar.gz')),
                    os.path.join('s3://', bucket_name, bucket_dir, uuid + '.tar.gz')]
    subprocess.check_call(s3am_command)


# Pipeline specific functions
def get_mean_insert_size(work_dir, bam_name):
    cmd = "docker run --log-driver=none --rm -v {}:/data quay.io/ucsc_cgl/samtools " \
          "view -f66 {}".format(work_dir, os.path.join(work_dir, bam_name))
    process = subprocess.Popen(args=cmd, shell=True, stdout=subprocess.PIPE)
    b_sum = 0L
    b_count = 0L
    while True:
        line = process.stdout.readline()
        if not line:
            break
        tmp = line.split("\t")
        if abs(long(tmp[8])) < 10000:
            b_sum += abs(long(tmp[8]))
            b_count += 1
    process.wait()
    try:
        mean = b_sum / b_count
    except ZeroDivisionError:
        mean = 150
    print "Using insert size: %d" % mean
    return int(mean)


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
                require(len(sample) == 3, 'Bad manifest format! '
                                          'Expected 3 tab separated columns, got: {}'.format(sample))
                uuid, normal, tumor = sample
                for url in [normal, tumor]:
                    require(urlparse(url).scheme and urlparse(url), 'Invalid URL passed for {}'.format(url))
                samples.append(sample)
    return samples


def generate_config():
    return textwrap.dedent("""
    # CGL Exome Pipeline configuration file
    # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
    # Edit the values in this configuration file and then rerun the pipeline: "toil-variant run"
    # URLs can take the form: http://, file://, s3://, gnos://.
    # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank
    ####################################################################################################################
    reference:              # Required: URL to reference genome
                            # Example: s3://cgl-pipeline-inputs/variant_hg19/hg19.fa\n
    phase:                  # Required: URL to phase indels VCF
                            # Example: s3://cgl-pipeline-inputs/variant_hg19/1000G_phase1.indels.hg19.sites.vcf\n
    mills:                  # Required: URL to Mills indel VCF
                            # Example: s3://cgl-pipeline-inputs/variant_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\n
    dbsnp:                  # Required: URL to dbsnp VCF
                            # Example: s3://cgl-pipeline-inputs/variant_hg19/dbsnp_138.hg19.vcf\n
    cosmic:                 # Required: URL to cosmic VCF
                            # Example: s3://cgl-pipeline-inputs/variant_hg19/cosmic.hg19.vcf\n\n
    run-mutect: true        # Optional: If true, will run MuTect to do mutation calls
    run-pindel: true        # Optional: Iff true, will run pindel to analyze indels
    run-muse: true          # Optional: If true, will run MuSe to do mutation calls
    preprocessing: true     # Optional: If true, will perform indel realignment and base quality score recalibration\n
    output-dir:             # Optional: Provide a full path to where results will appear\n
    s3-dir:                 # Optional: Provide an s3 path (s3://bucket/dir) where results will appear\n
    ssec:                   # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon\n
    gtkey:                  # Optional: Provide a full path to a CGHub Key used to access GNOS hosted data\n
    ci-test:                # Optional: If true, uses resource requirements appropriate for continuous integration\n
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample pair to be run.
        #   There are 3 tab-separated columns: UUID, Normal BAM URL, Tumor BAM URL
        #
        #   UUID            This should be a unique identifier for the sample to be processed
        #   Normal URL      A URL (http://, file://, s3://, gnos://) pointing to the normal bam
        #   Tumor URL       A URL (http://, file://, s3://, gnos://) pointing to the tumor bam
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  file:///path/to/normal.bam  file:///path/to/tumor.bam
        #   UUID_2  http://sample-depot.com/normal.bam  http://sample-depot.com/tumor.bam
        #   UUID_3  s3://my-bucket-name/directory/normal.bam    file:///path/to/tumor.bam
        #
        #   Place your samples below, one per line.
        """[1:])


def generate_file(file_path, generate_func):
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))


def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    Toil exome variant pipeline

    Perform variant analysis given a pair of tumor/normal BAM files.
    Samples are optionally preprocessed (indel realignment and base quality score recalibration)
    before variant analysis is performed using MuTect.  The final output of this pipeline is a tarball
    containing the output of output of MuTect.

    General usage:
    1. Type "toil-variant generate" to create an editable manifest and config in the current working directory.
    2. Parameterize the pipeline by editing the config.
    3. Fill in the manifest with information pertaining to your samples.
    4. Type "toil-variant run [jobStore]" to execute the pipeline.

    Please read the README.md located in the source directory or at:
    https://github.com/BD2KGenomics/toil-scripts/tree/master/src/toil_scripts/exome_variant_pipeline

    Structure of variant pipeline (per sample)

           1 2 3 4          14 -------
           | | | |          |        |
        0 --------- 5 ----- 15 -------- 17
                    |       |        |
                   ---      16 -------
                   | |
                   6 7
                   | |
                   8 9
                   | |
                  10 11
                   | |
                  12 13

    0 = Start node
    1 = reference index
    2 = reference dict
    3 = normal bam index
    4 = tumor bam index
    5 = pre-processing node / DAG declaration
    6,7 = RealignerTargetCreator
    8,9 = IndelRealigner
    10,11 = BaseRecalibration
    12,13 = PrintReads
    14 = MuTect
    15 = Pindel
    16 = MuSe
    17 = Consolidate Output and move/upload results
    ==================================================
    Dependencies
    Curl:       apt-get install curl
    Docker:     wget -qO- https://get.docker.com/ | sh
    Toil:       pip install toil
    Boto:       pip install boto (OPTIONAL)
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the CGL exome pipeline')
    parser_run.add_argument('--config', default='toil-exome.config', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--manifest', default='toil-exome-manifest.tsv', type=str,
                            help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')
    parser_run.add_argument('--normal', default=None, type=str,
                            help='URL for the normal BAM. URLs can take the form: http://, file://, s3://, '
                            'and gnos://. The UUID for the sample must be given with the "--uuid" flag.')
    parser_run.add_argument('--tumor', default=None, type=str,
                            help='URL for the tumor BAM. URLs can take the form: http://, file://, s3://, '
                                 'and gnos://. The UUID for the sample must be given with the "--uuid" flag.')
    parser_run.add_argument('--uuid', default=None, type=str, help='Provide the UUID of a sample when using the'
                                                                   '"--tumor" and "--normal" option')
    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # Add Toil options
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()
    # Parse subparsers related to generation of config and manifest
    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'toil-rnaseq.config'), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'toil-rnaseq-manifest.tsv'), generate_manifest)
    # Pipeline execution
    elif args.command == 'run':
        require(os.path.exists(args.config), '{} not found. Please run '
                                             '"toil-rnaseq generate-config"'.format(args.config))
        if args.normal or args.tumor or args.uuid:
            require(args.normal and args.tumor and args.uuid, '"--tumor", "--normal" and "--uuid" must all be supplied')
            samples = [[args.uuid, args.normal, args.tumor]]
        else:
            samples = parse_manifest(args.manifest)
        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        # Exome pipeline sanity checks
        if config.preprocessing:
            require(config.reference and config.phase and config.mills and config.dbsnp,
                    'Missing inputs for preprocessing, check config file.')
        if config.run_mutect:
            require(config.reference and config.dbsnp and config.cosmic,
                    'Missing inputs for MuTect, check config file.')
        if config.run_pindel:
            require(config.reference, 'Missing input (reference) for Pindel.')
        if config.run_muse:
            require(config.reference and config.dbsnp,
                    'Missing inputs for MuSe, check config file.')
        # Program checks
        for program in ['curl', 'docker']:
            require(which(program), program + ' must be installed on every node.'.format(program))

        # Launch Pipeline
        Job.Runner.startToil(Job.wrapJobFn(download_shared_files, samples, config), args)

if __name__ == '__main__':
    main()
