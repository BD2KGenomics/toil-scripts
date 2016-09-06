#!/usr/bin/env python2.7
import argparse
import multiprocessing
import os
import sys
import tarfile
import textwrap
from contextlib import closing
from urlparse import urlparse

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.job import Job
from toil_lib import require
from toil_lib.files import copy_files
from toil_lib.jobs import map_job
from toil_lib.tools.mutation_callers import run_muse
from toil_lib.tools.mutation_callers import run_mutect
from toil_lib.tools.mutation_callers import run_pindel
from toil_lib.tools.preprocessing import run_gatk_preprocessing
from toil_lib.tools.preprocessing import run_picard_create_sequence_dictionary
from toil_lib.tools.preprocessing import run_samtools_faidx
from toil_lib.tools.preprocessing import run_samtools_index
from toil_lib.urls import download_url_job, s3am_upload


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
            vars(config)[name] = job.addChildJobFn(download_url_job, url=url).rv()
    job.addFollowOnJobFn(reference_preprocessing, samples, config)


def reference_preprocessing(job, samples, config):
    """
    Spawn the jobs that create index and dict file for reference

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param list[list] samples: A nested list of samples containing sample information
    """
    job.fileStore.logToMaster('Processed reference files')
    config.fai = job.addChildJobFn(run_samtools_faidx, config.reference).rv()
    config.dict = job.addChildJobFn(run_picard_create_sequence_dictionary, config.reference).rv()
    job.addFollowOnJobFn(map_job, download_sample, samples, config)


def download_sample(job, sample, config):
    """
    Download sample and store sample specific attributes

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list sample: Contains uuid, normal URL, and tumor URL
    :param Namespace config: Argparse Namespace object containing argument inputs
    """
    # Create copy of config that is sample specific
    config = argparse.Namespace(**vars(config))
    uuid, normal_url, tumor_url = sample
    job.fileStore.logToMaster('Downloaded sample: ' + uuid)
    config.uuid = uuid
    config.normal = normal_url
    config.tumor = tumor_url
    config.cores = min(config.maxCores, int(multiprocessing.cpu_count()))
    disk = '1G' if config.ci_test else '20G'
    # Download sample bams and launch pipeline
    config.normal_bam = job.addChildJobFn(download_url_job, url=config.normal, s3_key_path=config.ssec,
                                          cghub_key_path=config.gtkey, disk=disk).rv()
    config.tumor_bam = job.addChildJobFn(download_url_job, url=config.tumor, s3_key_path=config.ssec,
                                         cghub_key_path=config.gtkey, disk=disk).rv()
    job.addFollowOnJobFn(index_bams, config)


def index_bams(job, config):
    """
    Convenience job for handling bam indexing to make the workflow declaration cleaner

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    """
    job.fileStore.logToMaster('Indexed sample BAMS: ' + config.uuid)
    disk = '1G' if config.ci_test else '20G'
    config.normal_bai = job.addChildJobFn(run_samtools_index, config.normal_bam, cores=1, disk=disk).rv()
    config.tumor_bai = job.addChildJobFn(run_samtools_index, config.tumor_bam, cores=1, disk=disk).rv()
    job.addFollowOnJobFn(preprocessing_declaration, config)


def preprocessing_declaration(job, config):
    """
    Declare jobs related to preprocessing

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    """
    if config.preprocessing:
        job.fileStore.logToMaster('Ran preprocessing: ' + config.uuid)
        disk = '1G' if config.ci_test else '20G'
        mem = '2G' if config.ci_test else '10G'
        processed_normal = job.wrapJobFn(run_gatk_preprocessing, config.normal_bam, config.normal_bai,
                                         config.reference, config.dict, config.fai, config.phase, config.mills,
                                         config.dbsnp, mem, cores=1, memory=mem, disk=disk)
        processed_tumor = job.wrapJobFn(run_gatk_preprocessing, config.tumor_bam, config.tumor_bai,
                                        config.reference, config.dict, config.fai, config.phase, config.mills,
                                        config.dbsnp, mem, cores=1, memory=mem, disk=disk)
        static_workflow = job.wrapJobFn(static_workflow_declaration, config, processed_normal.rv(0),
                                        processed_normal.rv(1), processed_tumor.rv(0), processed_tumor.rv(1))
        job.addChild(processed_normal)
        job.addChild(processed_tumor)
        job.addFollowOn(static_workflow)
    else:
        job.addFollowOnJobFn(static_workflow_declaration, config, config.normal_bam, config.normal_bai,
                             config.tumor_bam, config.tumor_bai)


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
    if config.run_mutect:
        mutect_results = job.addChildJobFn(run_mutect, normal_bam, normal_bai, tumor_bam, tumor_bai, config.reference,
                                           config.dict, config.fai, config.cosmic, config.dbsnp,
                                           cores=1, memory=memory, disk=disk).rv()
    if config.run_pindel:
        pindel_results = job.addChildJobFn(run_pindel, normal_bam, normal_bai, tumor_bam, tumor_bai,
                                           config.reference, config.fai,
                                           cores=config.cores,  memory=memory, disk=disk).rv()
    if config.run_muse:
        muse_results = job.addChildJobFn(run_muse, normal_bam, normal_bai, tumor_bam, tumor_bai,
                                         config.reference, config.dict, config.fai, config.dbsnp,
                                         cores=config.cores, memory=memory, disk=disk).rv()
    # Pass tool results (whether None or a promised return value) to consolidation step
    consolidation = job.wrapJobFn(consolidate_output, config, mutect_results, pindel_results, muse_results)
    job.addFollowOn(consolidation)


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
    # Move to output location
    if urlparse(config.output_dir).scheme == 's3':
        job.fileStore.logToMaster('Uploading {} to S3: {}'.format(config.uuid, config.output_dir))
        s3am_upload(fpath=out_tar, s3_dir=config.output_dir, num_cores=config.cores)
    else:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(config.uuid, config.output_dir))
        mkdir_p(config.output_dir)
        copy_files(file_paths=[out_tar], output_dir=config.output_dir)


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
    #
    # URLs can take the form: http://, ftp://, file://, s3://, gnos://
    # Local inputs follow the URL convention: file:///full/path/to/input
    # S3 URLs follow the convention: s3://bucket/directory/file.txt
    #
    # Comments (beginning with #) do not need to be removed. Optional parameters left blank are treated as false.
    ####################################################################################################################
    # Required: URL to reference genome
    reference: s3://cgl-pipeline-inputs/variant_hg19/hg19.fa

    # Required: URL to phase indels VCF
    phase: s3://cgl-pipeline-inputs/variant_hg19/1000G_phase1.indels.hg19.sites.vcf

    # Required: URL to Mills indel VCF
    mills: s3://cgl-pipeline-inputs/variant_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf

    # Required: URL to dbsnp VCF
    dbsnp: s3://cgl-pipeline-inputs/variant_hg19/dbsnp_138.hg19.vcf

    # Required: URL to cosmic VCF
    cosmic: s3://cgl-pipeline-inputs/variant_hg19/cosmic.hg19.vcf

    # Required: Output location of sample. Can be full path to a directory or an s3:// URL
    # Warning: S3 buckets must exist prior to upload or it will fail.
    output-dir:

    # Optional: If true, will run MuTect to do mutation calls
    run-mutect: true

    # Optional: If true, will run pindel to analyze indel
    run-pindel: true

    # Optional: If true, will run MuSe to do mutation calls
    run-muse: true

    # Optional: If true, will perform indel realignment and base quality score recalibration
    preprocessing: true

    # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
    ssec:

    # Optional: Provide a full path to a CGHub Key used to access GNOS hosted data
    gtkey:

    # Optional: If true, uses resource requirements appropriate for continuous integration
    ci-test: 
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample pair to be run.
        #   There are 3 tab-separated columns: UUID, Normal BAM URL, Tumor BAM URL
        #
        #   UUID            This should be a unique identifier for the sample to be processed
        #   Normal URL      A URL (http://, ftp://, file://, s3://, gnos://) pointing to the normal bam
        #   Tumor URL       A URL (http://, ftp://, file://, s3://, gnos://) pointing to the tumor bam
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
    Toil exome pipeline

    Perform variant / indel analysis given a pair of tumor/normal BAM files.
    Samples are optionally preprocessed (indel realignment and base quality score recalibration)
    The output of this pipeline is a tarball containing results from MuTect, MuSe, and Pindel.

    General usage:
    1. Type "toil-exome generate" to create an editable manifest and config in the current working directory.
    2. Parameterize the pipeline by editing the config.
    3. Fill in the manifest with information pertaining to your samples.
    4. Type "toil-exome run [jobStore]" to execute the pipeline.

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
    parser_run.add_argument('--config', default='config-toil-exome.yaml', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--manifest', default='manifest-toil-exome.tsv', type=str,
                            help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--normal', default=None, type=str,
                            help='URL for the normal BAM. URLs can take the form: http://, ftp://, file://, s3://, '
                                 'and gnos://. The UUID for the sample must be given with the "--uuid" flag.')
    parser_run.add_argument('--tumor', default=None, type=str,
                            help='URL for the tumor BAM. URLs can take the form: http://, ftp://, file://, s3://, '
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
        generate_file(os.path.join(cwd, 'config-toil-exome.yaml'), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-toil-exome.tsv'), generate_manifest)
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
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxint
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
        require(config.output_dir, 'No output location specified: {}'.format(config.output_dir))
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Launch Pipeline
        Job.Runner.startToil(Job.wrapJobFn(download_shared_files, samples, config), args)


if __name__ == '__main__':
    main()
