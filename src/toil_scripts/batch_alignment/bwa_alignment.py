"""
UC Santa Cruz Computational Genomics Lab
Updated: 4/21/16
Authors: John Vivian, Frank Nothaft

Alignment of fastq reads via BWAkit

    0 --> 1 -->(batching)--> 2 --> 3

0   Download shared files
1   Parse config
2   Download sample
3   Run BWAkit

===================================================================
:Dependencies:
curl            - apt-get install curl
Toil            - pip install toil
Docker          - http://docs.docker.com/engine/installation/

Optional:
S3AM            - pip install --s3am (requires ~/.boto config file)
"""
import argparse
import logging
import multiprocessing
import os
import sys
import textwrap

import yaml
from toil.job import Job
from toil_scripts.lib import require
from toil_scripts.lib.files import move_files
from toil_scripts.lib.jobs import map_job
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url_job, s3am_upload
from toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline import generate_file

_log = logging.getLogger(__name__)

def download_shared_files(job, inputs, sample, output_dir, rg_line=None):
    """
    Downloads shared files that are used by all samples for alignment

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param Namespace inputs: Input arguments (see main)
    """
    inputs.sample = sample
    inputs.rg_line = rg_line
    inputs.output_dir = output_dir
    job.fileStore.logToMaster('Downloading shared files for alignment.')

    shared_files = [inputs.ref, inputs.amb, inputs.ann, inputs.bwt, inputs.pac, inputs.sa, inputs.fai]
    if getattr(inputs, 'alt', None):
        shared_files.append(inputs.alt)

    shared_ids = []
    for url in shared_files:
        shared_ids.append(job.addChildJobFn(download_url_job, url).rv())
    job.addFollowOnJobFn(parse_config, shared_ids, inputs)


def parse_config(job, shared_ids, inputs):
    """
    Stores the UUID and urls associated with the input files to be retrieved.
    Configuration file has one sample per line, with the following format:  UUID\t1st_url\t2nd_url

    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param dict shared_ids: stores FileStoreIDs
    :param Namespace inputs: Input arguments
    """
    if inputs.sample:
        uuid, url1, url2 = inputs.sample
        mock_bam = '.'.join(url1.split('.')[:-2])[:-2] + ".bam"
        samples = [[uuid, [url1, url2], mock_bam]]
    else:
        samples = []

        with open(inputs.manifest, 'r') as f_in:
            for line in f_in:
                if not line.isspace() and not line.startswith('#'):
                    line = line.strip().split('\t')
                    assert len(line) == 3, 'Improper formatting. Expected UUID\tURl1\tURL2. Received: {}'.format(line)
                    uuid = line[0]
                    urls = line[1:]
                    mock_bam = '.'.join(line[1].split('.')[:-2])[:-2] + ".bam"
                    samples.append((uuid, urls, mock_bam))
    inputs.maxCores = int(inputs.maxCores) if hasattr(inputs, "maxCores") else sys.maxint
    inputs.cores = min(inputs.maxCores, multiprocessing.cpu_count())
    job.fileStore.logToMaster('Parsed configuration file.')
    job.addChildJobFn(map_job, download_sample, samples, inputs, shared_ids, cores=1, disk=inputs.file_size)


def download_sample(job, sample, inputs, ids):
    """
    Downloads the sample

    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param tuple[str] sample: uuid and URLS for sample
    :param Namespace inputs: Contains input arguments
    :param list ids: list of FileStore IDs for shared inputs
    """
    uuid, urls, mock_bam = sample
    inputs.uuid = uuid
    inputs.mock_bam = mock_bam
    job.fileStore.logToMaster('Downloading sample: {0}\nStarting BWA Run: {0}'.format(uuid))
    ids.insert(0, job.addChildJobFn(download_url_job, urls[0], s3_key_path=inputs.ssec, disk=inputs.file_size).rv())
    ids.insert(1, job.addChildJobFn(download_url_job, urls[1], s3_key_path=inputs.ssec, disk=inputs.file_size).rv())
    job.addFollowOnJobFn(run_bwa, inputs, ids, cores=inputs.cores, disk=inputs.file_size)


def run_bwa(job, inputs, ids):
    """
    Aligns two fastqs into a BAMFILE via BWA

    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param Namespace inputs: Input arguments (see main)
    :param list ids: list of FileStore IDs (R1, R2, reference inputs)
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_names = ['r1.fq.gz', 'r2.fq.gz', 'ref.fa', 'ref.fa.amb', 'ref.fa.ann',
                  'ref.fa.bwt', 'ref.fa.pac', 'ref.fa.sa', 'ref.fa.fai']
    if inputs.alt:
        file_names.append('ref.fa.alt')

    for fileStoreID, name in zip(ids, file_names):
        job.fileStore.readGlobalFile(fileStoreID, os.path.join(work_dir, name))

    rg = ''
    if inputs.rg_line:
        rg = inputs.rg_line
    else:
        rg = "@RG\\tID:{0}".format(inputs.uuid) # '\' character is escaped so bwakit gets individual '\' and 't' characters 
        for tag, info in zip(['LB', 'PL', 'PU', 'SM'], ['library', 'platform', 'program_unit', 'uuid']):
            if hasattr(inputs, info):
                rg = rg + '\\t{0}:{1}'.format(tag, getattr(inputs, info))

    # BWA Options
    opt_args = []
    if inputs.sort:
        opt_args.append('-s')
    if inputs.trim_adapters:
        opt_args.append('-a')
    # Call: bwakit
    parameters = (['-t', str(inputs.cores),
                   '-R', rg] +
                  opt_args +
                  ['-o', '/data/aligned',
                   '/data/ref.fa',
                   '/data/r1.fq.gz',
                   '/data/r2.fq.gz'])
    outputs = {'aligned.aln.bam': inputs.mock_bam}

    docker_call(tool='quay.io/ucsc_cgl/bwakit:0.7.12--528bb9bf73099a31e74a7f5e6e3f2e0a41da486e',
                parameters=parameters, inputs=file_names, outputs=outputs, work_dir=work_dir)

    # BWA insists on adding an `*.aln.sam` suffix, so rename the output file
    output_file = os.path.join(work_dir, '{}.bam'.format(inputs.uuid))
    os.rename(os.path.join(work_dir, 'aligned.aln.bam'),
              output_file)

    # Either write file to local output directory or upload to S3 cloud storage
    job.fileStore.logToMaster('Aligned sample: {}'.format(inputs.uuid))

    if inputs.output_dir.startswith('s3://'):
        s3am_upload(output_file, inputs.output_dir, s3_key_path=inputs.ssec)
    else:
        move_files([output_file], inputs.output_dir)

def generate_config():
    return textwrap.dedent("""
        # BWA Alignment Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        ssec:                     # Optional: (string) Path to Key File for SSE-C Encryption
        library:                  # Optional: the LB (library) entry to go in the BAM header. Exactly one of library
                                  # and rg_line must be provided
        program_unit: 12345       #  Program Unit for BAM header
        platform: ILLUMINA        #  Platform to put in the read group
        rg_line:                  # Optional: Use instead of library, program_unit, and platform. Exactly one of library
                                  # and rg_line must be provided.
        ref:                      # Required: Reference fasta file
        amb:                      # Required: Reference fasta file (amb)
        ann:                      # Required: Reference fasta file (ann)
        bwt:                      # Required: Reference fasta file (bwt)
        pac:                      # Required: Reference fasta file (pac)
        sa:                       # Required: Reference fasta file (sa)
        fai:                      # Required: Reference fasta file (fai)
        alt:                      # Optional: Alternate file for reference build (alt). Necessary for alt aware alignment.
        file_size: 100G           # Approximate input file size. Should be given as %d[TGMK], e.g.,
                                  # for a 100 gigabyte file, use --file_size 100G
        sort: True                # Whether to sort
        trim_adapters: False      # Trim adapters.
    """[1:])

def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 3 tab-separated columns: UUID, URL1, URL2
        #
        #   UUID        This should be a unique identifier for the sample to be processed
        #   URL1/URL2         A URL (http://, ftp://, file://, s3://, gnos://) pointing to the paired sample fastq files
        #
        #   Example below. Lines beginning with # are ignored.
        #
        #   UUID_1    file:///path/to/R1.fq    file:///path/to/R2.fq
        #
        #   Place your samples below, one per line.
        """[1:])

def main():
    """
    This is a Toil pipeline used to perform alignment of fastqs.
    """
    # Define Parser object and add to Toil
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the BWA alignment pipeline')
    group = parser_run.add_mutually_exclusive_group(required=True)
    parser_run.add_argument('--config', default='bwa-alignment.config', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config".')
    group.add_argument('--manifest', default='bwa-alignment-manifest.csv', type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s".')
    group.add_argument('--sample', default=None, nargs='3', type=str,
                       help='Space delimited sample UUID and paired fastq files in the format: uuid url1 url2')
    parser_run.add_argument('--output-dir', default=None, help='Full path to directory or filename where '
                                                               'final results will be output')    
    parser_run.add_argument('-s', '--suffix', default='.bqsr', help='Additional suffix to add to the names of the output files.')
    Job.Runner.addToilOptions(parser_run)
    options = parser.parse_args()

    cwd = os.getcwd()
    if options.command == 'generate-config' or options.command == 'generate':
        generate_file(os.path.join(cwd, 'bwa-alignment.config'), generate_config)
    if options.command == 'generate-manifest' or options.command == 'generate':
        generate_file(os.path.join(cwd, 'bwa-alignment-manifest.csv'), generate_manifest)
    # Pipeline execution
    elif options.command == 'run':
        require(os.path.exists(options.config), '{} not found. Please run '
                                                'generate-config'.format(options.config))
        if not options.sample:
            options.sample = None
            require(os.path.exists(options.manifest), '{} not found and no sample provided. Please '
                                                       'run "generate-manifest"'.format(options.manifest))
        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(options.config).read()).iteritems()}
        inputs = argparse.Namespace(**parsed_config)
        if options.manifest:
            inputs.manifest = options.manifest

    # Launch Pipeline
    Job.Runner.startToil(Job.wrapJobFn(download_shared_files, inputs, options.sample, options.output_dir), options)


if __name__ == "__main__":
    main()
