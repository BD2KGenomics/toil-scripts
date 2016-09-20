#!/usr/bin/env python2.7
import argparse
import multiprocessing
import os
import sys
import textwrap
from urlparse import urlparse

import yaml
from bd2k.util.files import mkdir_p
from toil.job import Job
from toil_lib import require, required_length
from toil_lib.files import copy_file_job
from toil_lib.files import generate_file
from toil_lib.jobs import map_job
from toil_lib.tools.aligners import run_bwakit
from toil_lib.tools.indexing import run_samtools_faidx, run_bwa_index
from toil_lib.urls import download_url_job, s3am_upload_job


def download_reference_files(job, inputs, samples):
    """
    Downloads shared files that are used by all samples for alignment, or generates them if they were not provided.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Input arguments (see main)
    :param list[list[str, list[str, str]]] samples: Samples in the format [UUID, [URL1, URL2]]
    """
    # Create dictionary to store FileStoreIDs of shared input files
    shared_ids = {}
    urls = [('amb', inputs.amb), ('ann', inputs.ann), ('bwt', inputs.bwt),
            ('pac', inputs.pac), ('sa', inputs.sa)]
    # Alt file is optional and can only be provided, not generated
    if inputs.alt:
        urls.append(('alt', inputs.alt))
    # Download reference
    download_ref = job.wrapJobFn(download_url_job, inputs.ref, disk='3G')  # Human genomes are typically ~3G
    job.addChild(download_ref)
    shared_ids['ref'] = download_ref.rv()
    # If FAI is provided, download it. Otherwise, generate it
    if inputs.fai:
        shared_ids['fai'] = job.addChildJobFn(download_url_job, inputs.fai).rv()
    else:
        faidx = job.wrapJobFn(run_samtools_faidx, download_ref.rv())
        shared_ids['fai'] = download_ref.addChild(faidx).rv()
    # If all BWA index files are provided, download them. Otherwise, generate them
    if all(x[1] for x in urls):
        for name, url in urls:
            shared_ids[name] = job.addChildJobFn(download_url_job, url).rv()
    else:
        job.fileStore.logToMaster('BWA index files not provided, creating now')
        bwa_index = job.wrapJobFn(run_bwa_index, download_ref.rv())
        download_ref.addChild(bwa_index)
        for x, name in enumerate(['amb', 'ann', 'bwt', 'pac', 'sa']):
            shared_ids[name] = bwa_index.rv(x)

    # Map_job distributes one sample in samples to the downlaod_sample_and_align function
    job.addFollowOnJobFn(map_job, download_sample_and_align, samples, inputs, shared_ids)


def download_sample_and_align(job, sample, inputs, ids):
    """
    Downloads the sample and runs BWA-kit

    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param tuple(str, list) sample: UUID and URLS for sample
    :param Namespace inputs: Contains input arguments
    :param dict ids: FileStore IDs for shared inputs
    """
    uuid, urls = sample
    r1_url, r2_url = urls if len(urls) == 2 else (urls[0], None)
    job.fileStore.logToMaster('Downloaded sample: {0}. R1 {1}\nR2 {2}\nStarting BWA Run'.format(uuid, r1_url, r2_url))
    # Read fastq samples from file store
    ids['r1'] = job.addChildJobFn(download_url_job, r1_url, s3_key_path=inputs.ssec, disk=inputs.file_size).rv()
    if r2_url:
        ids['r2'] = job.addChildJobFn(download_url_job, r2_url, s3_key_path=inputs.ssec, disk=inputs.file_size).rv()
    else:
        ids['r2'] = None
    # Create config for bwakit
    inputs.cores = min(inputs.maxCores, multiprocessing.cpu_count())
    inputs.uuid = uuid
    config = dict(**vars(inputs))  # Create config as a copy of inputs since it has values we want
    config.update(ids)  # Overwrite attributes with the FileStoreIDs from ids
    config = argparse.Namespace(**config)
    # Define and wire job functions
    bam_id = job.wrapJobFn(run_bwakit, config, sort=inputs.sort, trim=inputs.trim,
                           disk=inputs.file_size, cores=inputs.cores)
    job.addFollowOn(bam_id)
    output_name = uuid + '.bam' + str(inputs.suffix) if inputs.suffix else uuid + '.bam'
    if urlparse(inputs.output_dir).scheme == 's3':
        bam_id.addChildJobFn(s3am_upload_job, file_id=bam_id.rv(), file_name=output_name, s3_dir=inputs.output_dir,
                             s3_key_path=inputs.ssec, cores=inputs.cores, disk=inputs.file_size)
    else:
        mkdir_p(inputs.ouput_dir)
        bam_id.addChildJobFn(copy_file_job, name=output_name, file_id=bam_id.rv(), output_dir=inputs.output_dir,
                                    disk=inputs.file_size)


def generate_config():
    return textwrap.dedent("""
        # BWA Alignment Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline: "toil-bwa run"
        #
        # URLs can take the form: http://, ftp://, file://, s3://, gnos://
        # Local inputs follow the URL convention: file:///full/path/to/input
        # S3 URLs follow the convention: s3://bucket/directory/file.txt
        #
        # Comments (beginning with #) do not need to be removed. Optional parameters left blank are treated as false.
        ##############################################################################################################
        # Required: Reference fasta file
        ref: s3://cgl-pipeline-inputs/alignment/hg19.fa

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        output-dir:

        # Required: The library entry to go in the BAM read group.
        library: Illumina

        # Required: Platform to put in the read group
        platform: Illumina

        # Required: Program Unit for BAM header. Required for use with GATK.
        program_unit: 12345

        # Required: Approximate input file size. Provided as a number followed by (base-10) [TGMK]. E.g. 10M, 150G
        file-size: 50G

        # Optional: If true, sorts bam
        sort: True

        # Optional. If true, trims adapters
        trim: false

        # Optional: Reference fasta file (amb) -- if not present will be generated
        amb: s3://cgl-pipeline-inputs/alignment/hg19.fa.amb

        # Optional: Reference fasta file (ann) -- If not present will be generated
        ann: s3://cgl-pipeline-inputs/alignment/hg19.fa.ann

        # Optional: Reference fasta file (bwt) -- If not present will be generated
        bwt: s3://cgl-pipeline-inputs/alignment/hg19.fa.bwt

        # Optional: Reference fasta file (pac) -- If not present will be generated
        pac: s3://cgl-pipeline-inputs/alignment/hg19.fa.pac

        # Optional: Reference fasta file (sa) -- If not present will be generated
        sa: s3://cgl-pipeline-inputs/alignment/hg19.fa.sa

        # Optional: Reference fasta file (fai) -- If not present will be generated
        fai: s3://cgl-pipeline-inputs/alignment/hg19.fa.fai

        # Optional: (string) Path to Key File for SSE-C Encryption
        ssec:

        # Optional: Use instead of library, program_unit, and platform.
        rg-line:

        # Optional: Alternate file for reference build (alt). Necessary for alt aware alignment
        alt:

        # Optional: If true, runs the pipeline in mock mode, generating a fake output bam
        mock-mode:

        # Optional: Optional suffix to add to sample output
        suffix:
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information for each sample to be run.
        #   Lines beginning with # are ignored.
        #
        #   There are 2 or 3 tab-separated columns: UUID, 1st FASTQ URL, 2nd FASTQ URL (if paired)
        #   If a sample is paired end: UUID    URL1    URL2
        #   If a sample is single-ended: UUID    URL
        #
        #   UUID            This should be a unique identifier for the sample.
        #   URL1/URL2       A URL (http://, ftp://, file://, s3://, gnos://) pointing to the sample fastq files.
        #
        #   Examples below:
        #
        #   Paired_UUID    file:///path/to/R1.fq.gz    file:///path/to/R2.fq.gz
        #   Unpaired_UUID    file:///path/to/unpaired.fq.gz
        #
        #   Place your samples below, one sample per line.
        """[1:])


def parse_manifest(manifest_path):
    """
    Parse manifest file

    :param str manifest_path: Path to manifest file
    :return: samples
    :rtype: list[str, list]
    """
    samples = []
    with open(manifest_path, 'r') as f:
        for line in f:
            if not line.isspace() and not line.startswith('#'):
                sample = line.strip().split('\t')
                require(2 <= len(sample) <= 3, 'Bad manifest format! '
                                               'Expected UUID\tURL1\t[URL2] (tab separated), got: {}'.format(sample))
                uuid = sample[0]
                urls = sample[1:]
                for url in urls:
                    require(urlparse(url).scheme and urlparse(url), 'Invalid URL passed for {}'.format(url))
                samples.append([uuid, urls])
    return samples


def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    Toil BWA pipeline

    Alignment of fastq reads via BWA-kit

    General usage:
    1. Type "toil-bwa generate" to create an editable manifest and config in the current working directory.
    2. Parameterize the pipeline by editing the config.
    3. Fill in the manifest with information pertaining to your samples.
    4. Type "toil-bwa run [jobStore]" to execute the pipeline.

    Please read the README.md located in the source directory or at:
    https://github.com/BD2KGenomics/toil-scripts/tree/master/src/toil_scripts/bwa_alignment

    Structure of the BWA pipeline (per sample)

        0 --> 1

    0 = Download sample
    1 = Run BWA-kit
    ===================================================================
    :Dependencies:
    cURL:       apt-get install curl
    Toil:       pip install toil
    Docker:     wget -qO- https://get.docker.com/ | sh

    Optional:
    S3AM:       pip install --s3am (requires ~/.boto config file)
    Boto:       pip install boto
    """
    # Define Parser object and add to Toil
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the BWA alignment pipeline')
    group = parser_run.add_mutually_exclusive_group()
    parser_run.add_argument('--config', default='config-toil-bwa.yaml', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config".')
    group.add_argument('--manifest', default='manifest-toil-bwa.tsv', type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s".')
    group.add_argument('--sample', nargs='+', action=required_length(2, 3),
                       help='Space delimited sample UUID and fastq files in the format: uuid url1 [url2].')
    # Print docstring help if no arguments provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()
    # Parse subparsers related to generation of config and manifest
    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'config-toil-bwa.yaml'), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-toil-bwa.tsv'), generate_manifest)
    # Pipeline execution
    elif args.command == 'run':
        require(os.path.exists(args.config), '{} not found. Please run generate-config'.format(args.config))
        if not args.sample:
            args.sample = None
            require(os.path.exists(args.manifest), '{} not found and no sample provided. '
                                                   'Please run "generate-manifest"'.format(args.manifest))
        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxint
        samples = [args.sample[0], args.sample[1:]] if args.sample else parse_manifest(args.manifest)
        # Sanity checks
        require(config.ref, 'Missing URL for reference file: {}'.format(config.ref))
        require(config.output_dir, 'No output location specified: {}'.format(config.output_dir))
        # Launch Pipeline
        Job.Runner.startToil(Job.wrapJobFn(download_reference_files, config, samples), args)


if __name__ == "__main__":
    main()
