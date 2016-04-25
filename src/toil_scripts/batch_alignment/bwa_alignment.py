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
import multiprocessing
import os

from toil.job import Job

from toil_scripts.lib.files import move_files
from toil_scripts.lib.jobs import map_job
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url_job, s3am_upload


def download_shared_files(job, inputs):
    """
    Downloads shared files that are used by all samples for alignment

    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param Namespace inputs: Input arguments (see main)
    """
    job.fileStore.logToMaster('Downloading shared files for aligment.')

    shared_files = [inputs.ref, inputs.amb, inputs.ann, inputs.bwt, inputs.pac, inputs.sa, inputs.fai]
    if inputs.alt:
        shared_files.append(inputs.alt)

    shared_ids = []
    for url in shared_files:
        shared_ids.append(job.addChildJobFn(download_url_job, url).rv())
    job.addFollowOnJobFn(parse_config, shared_ids, inputs)


def parse_config(job, shared_ids, inputs):
    """
    Stores the UUID and urls associated with the input files to be retrieved.
    Configuration file has one sample per line, with the following format:  UUID,1st_url,2nd_url

    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param dict shared_ids: stores FileStoreIDs
    :param Namespace inputs: Input arguments
    """
    samples = []

    with open(inputs.config, 'r') as f_in:
        for line in f_in:
            if line.strip():
                line = line.strip().split(',')
                assert len(line) == 3, 'Improper formatting. Expected UUID,URl1,URL2. Received: {}'.format(line)
                uuid = line[0]
                urls = line[1:]
                mock_bam = '.'.join(line[1].split('.')[:-2])[:-2] + ".bam"
                samples.append((uuid, urls, mock_bam))
    inputs.cores = multiprocessing.cpu_count()
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

    # Add read group line
    rg = "@RG\\tID:{0}\\tLB:{1}\\tPL:{2}\\tPU:{3}\\tSM:{0}".format(inputs.uuid, inputs.library,
                                                                   inputs.platform, inputs.program_unit)

    # BWA Options
    opt_args = []
    if not inputs.skip_sort:
        opt_args.append('-s')
    if inputs.trim:
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
                parameters=parameters, inputs=file_names, outputs=outputs, work_dir=work_dir, sudo=inputs.sudo)

    # BWA insists on adding an `*.aln.sam` suffix, so rename the output file
    output_file = os.path.join(work_dir, '{}.bam'.format(inputs.uuid))
    os.rename(os.path.join(work_dir, 'aligned.aln.bam'),
              output_file)

    # Either write file to local output directory or upload to S3 cloud storage
    job.fileStore.logToMaster('Aligned sample: {}'.format(inputs.uuid))
    if inputs.output_dir:
        move_files([output_file], inputs.output_dir)
    if inputs.s3_dir:
        s3am_upload(output_file, inputs.s3_dir, s3_key_path=inputs.ssec)


def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, add_help=True)
    parser.add_argument('-c', '--config', required=True, help='configuration file. One sample per line: uuid,url,url')
    parser.add_argument('--library', required=True, help='the LB (library) entry to go in the BAM header')
    parser.add_argument('--program-unit', default='12345', help='Program Unit for BAM header')
    parser.add_argument('--platform', default='ILLUMINA', help='Platform to put in the read group')
    parser.add_argument('-r', '--ref', required=True, help='Reference fasta file')
    parser.add_argument('-m', '--amb', required=True, help='Reference fasta file (amb)')
    parser.add_argument('-n', '--ann', required=True, help='Reference fasta file (ann)')
    parser.add_argument('-b', '--bwt', required=True, help='Reference fasta file (bwt)')
    parser.add_argument('-p', '--pac', required=True, help='Reference fasta file (pac)')
    parser.add_argument('-a', '--sa', required=True, help='Reference fasta file (sa)')
    parser.add_argument('-f', '--fai', required=True, help='Reference fasta file (fai)')
    parser.add_argument('-t', '--alt', required=False,
                        help='Alternate file for reference build (alt). Necessary for alt aware alignment')
    parser.add_argument('--ssec', default=None, help='Path to Key File for SSE-C Encryption')
    parser.add_argument('--output-dir', default=None, help='full path where final results will be output')
    parser.add_argument('--sudo', dest='sudo', default=False, action='store_true',
                        help='Docker usually needs sudo to execute locally, but not when running Mesos '
                             'or when a member of a Docker group.')
    parser.add_argument('--s3-dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                       'cgl-driver-projects/ckcc/rna-seq-samples/')
    parser.add_argument('--file-size', default='50G', help='Approximate input file size. Should be given as %d[TGMK], '
                                                           'e.g., for a 100 gigabyte file, use --file_size 100G')
    parser.add_argument('--trim', action='store_true', help='Trim adapters.')
    parser.add_argument('--skip-sort', action='store_true', default=False, help='Skip sorting.')
    return parser


def main():
    """
    This is a Toil pipeline used to perform alignment of fastqs.
    """
    # Define Parser object and add to Toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # Launch Pipeline
    Job.Runner.startToil(Job.wrapJobFn(download_shared_files, args), args)


if __name__ == "__main__":
    main()
