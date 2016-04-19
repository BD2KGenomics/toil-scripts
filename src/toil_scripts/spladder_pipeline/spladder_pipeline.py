#!/usr/bin/env python2.7
"""
UCSC Computational Genomics Lab Spladder Pipeline
Author: John Vivian
Affiliation: UC Santa Cruz Genomics Institute

Structure of Spladder Pipeline (per sample)

         0 ---> 4 --> 5
         |
         1
        / \
       2   3

0   -   Download Sample
1   -   STAR
2   -   QC and Variant Calling
3   -   Alternative Splice Calling
4   -   Consolidate Output
5   -   Upload results to S3
"""
import argparse
import fnmatch
import multiprocessing
import os
import shutil
import subprocess
import tarfile
from contextlib import closing
from glob import glob

from toil.job import Job

from toil_scripts.lib.files import mkdir_p
from toil_scripts.lib.files import tarball_files
from toil_scripts.lib.jobs import map_job
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.programs import which
from toil_scripts.lib.urls import s3am_upload_job
from toil_scripts.lib.urls import download_url_job
from toil_scripts.lib.urls import download_url


def parse_input_samples(job, inputs):
    """
    Parses config file to pull sample information.
    Stores samples as tuples of (uuid, URL)

    :param JobFunctionWrappingJob job: passed by Toil automatically
    :param Namespace inputs: Stores input arguments (see main)
    """
    job.fileStore.logToMaster('Parsing input samples and batching jobs')
    samples = []
    if inputs.config:
        with open(inputs.config, 'r') as f:
            for line in f.readlines():
                if not line.isspace():
                    sample = line.strip().split(',')
                    assert len(sample) == 2, 'Error: Config file is inappropriately formatted.'
                    samples.append(sample)
    job.addChildJobFn(map_job, download_sample, samples, inputs)


def download_sample(job, sample, inputs):
    """
    Download the input sample

    :param JobFunctionWrappingJob job: passed by Toil automatically
    :param tuple sample: Tuple containing (UUID,URL) of a sample
    :param Namespace inputs: Stores input arguments (see main)
    """
    uuid, url = sample
    job.fileStore.logToMaster('Downloading sample: {}'.format(uuid))
    # Download sample
    tar_id = job.addChildJobFn(download_url_job, url, s3_key_path=inputs.ssec, disk='30G').rv()
    # Create copy of inputs for each sample
    sample_inputs = argparse.Namespace(**vars(inputs))
    sample_inputs.uuid = uuid
    sample_inputs.cores = multiprocessing.cpu_count()
    # Call children and follow-on jobs
    job.addFollowOnJobFn(process_sample, sample_inputs, tar_id, cores=2, disk='60G')


def process_sample(job, inputs, tar_id):
    """
    Converts sample.tar(.gz) into two fastq files.
    Due to edge conditions... BEWARE: HERE BE DRAGONS

    :param JobFunctionWrappingJob job: passed by Toil automatically
    :param Namespace inputs: Stores input arguments (see main)
    :param str tar_id: FileStore ID of sample tar
    """
    job.fileStore.logToMaster('Processing sample into read pairs: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    # I/O
    tar_path = job.fileStore.readGlobalFile(tar_id, os.path.join(work_dir, 'sample.tar'))
    # Untar File and concat
    subprocess.check_call(['tar', '-xvf', tar_path, '-C', work_dir])
    os.remove(os.path.join(work_dir, 'sample.tar'))
    # Grab files from tarball
    fastqs = []
    for root, subdir, files in os.walk(work_dir):
        fastqs.extend([os.path.join(root, x) for x in files])
    # Check for read 1 and read 2 files
    r1 = sorted([x for x in fastqs if 'R1' in x])
    r2 = sorted([x for x in fastqs if 'R2' in x])
    if not r1 or not r2:
        # Check if using a different standard
        r1 = sorted([x for x in fastqs if '_1' in x])
        r2 = sorted([x for x in fastqs if '_2' in x])
    # Prune file name matches from each list
    if len(r1) > len(r2):
        r1 = [x for x in r1 if x not in r2]
    elif len(r2) > len(r1):
        r2 = [x for x in r2 if x not in r1]
    # Flag if data is single-ended
    assert r1 and r2, 'This pipeline does not support single-ended data. R1: {}\nR2:{}'.format(r1, r2)
    command = 'zcat' if r1[0].endswith('gz') and r2[0].endswith('gz') else 'cat'
    with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
        p1 = subprocess.Popen([command] + r1, stdout=f1)
    with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
        p2 = subprocess.Popen([command] + r2, stdout=f2)
    p1.wait()
    p2.wait()
    # Write to fileStore
    r1_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1.fastq'))
    r2_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2.fastq'))
    job.fileStore.deleteGlobalFile(tar_id)
    # Start cutadapt step
    job.addChildJobFn(cutadapt, inputs, r1_id, r2_id, disk='60G').rv()


def cutadapt(job, inputs, r1_id, r2_id):
    """
    Filters out adapters that may be left in the RNA-seq files

    :param JobFunctionWrappingJob job: passed by Toil automatically
    :param Namespace inputs: Stores input arguments (see main)
    :param str r1_id: FileStore ID of read 1 fastq
    :param str r2_id: FileStore ID of read 2 fastq
    """
    job.fileStore.logToMaster('Running CutAdapt: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    inputs.improper_pair = None
    # Retrieve files
    job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
    job.fileStore.readGlobalFile(r2_id, os.path.join(work_dir, 'R2.fastq'))
    # Cutadapt parameters
    parameters = ['-a', inputs.fwd_3pr_adapter,
                  '-m', '35',
                  '-A', inputs.rev_3pr_adapter,
                  '-o', '/data/R1_cutadapt.fastq',
                  '-p', '/data/R2_cutadapt.fastq',
                  '/data/R1.fastq', '/data/R2.fastq']
    # Call: CutAdapt
    base_docker_call = 'docker run --log-driver=none --rm -v {}:/data'.format(work_dir).split()
    if inputs.sudo:
        base_docker_call = ['sudo'] + base_docker_call
    tool = 'quay.io/ucsc_cgl/cutadapt:1.9--6bd44edd2b8f8f17e25c5a268fedaab65fa851d2'
    p = subprocess.Popen(base_docker_call + [tool] + parameters, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        if 'improperly paired' in stderr:
            inputs.improper_pair = True
            shutil.move(os.path.join(work_dir, 'R1.fastq'), os.path.join(work_dir, 'R1_cutadapt.fastq'))
            shutil.move(os.path.join(work_dir, 'R2.fastq'), os.path.join(work_dir, 'R2_cutadapt.fastq'))
    # Write to fileStore
    if inputs.improper_pair:
        r1_cutadapt = r1_id
        r2_cutadapt = r2_id
    else:
        r1_cutadapt = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1_cutadapt.fastq'))
        r2_cutadapt = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2_cutadapt.fastq'))
        job.fileStore.deleteGlobalFile(r1_id)
        job.fileStore.deleteGlobalFile(r2_id)
    # start STAR
    cores = min(inputs.cores, 16)
    job.addChildJobFn(star, inputs, r1_cutadapt, r2_cutadapt, cores=cores, disk='100G', memory='40G').rv()


def star(job, inputs, r1_cutadapt, r2_cutadapt):
    """
    Performs alignment of fastqs to BAM via STAR

    :param JobFunctionWrappingJob job: passed by Toil automatically
    :param Namespace inputs: Stores input arguments (see main)
    :param str r1_cutadapt: FileStore ID of read 1 fastq
    :param str r2_cutadapt: FileStore ID of read 2 fastq
    """
    job.fileStore.logToMaster('Aligning with STAR: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    cores = min(inputs.cores, 16)
    # Retrieve files
    job.fileStore.readGlobalFile(r1_cutadapt, os.path.join(work_dir, 'R1_cutadapt.fastq'))
    job.fileStore.readGlobalFile(r2_cutadapt, os.path.join(work_dir, 'R2_cutadapt.fastq'))
    # Get starIndex
    download_url(inputs.star_index, work_dir, 'starIndex.tar.gz')
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'starIndex.tar.gz'), '-C', work_dir])
    # Parameters
    parameters = ['--runThreadN', str(cores),
                  '--genomeDir', '/data/starIndex',
                  '--outFileNamePrefix', 'rna',
                  '--outSAMtype', 'BAM', 'SortedByCoordinate',
                  '--outSAMunmapped', 'Within',
                  '--quantMode', 'TranscriptomeSAM',
                  '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                  '--outFilterType', 'BySJout',
                  '--outFilterMultimapNmax', '20',
                  '--outFilterMismatchNmax', '999',
                  '--outFilterMismatchNoverReadLmax', '0.04',
                  '--alignIntronMin', '20',
                  '--alignIntronMax', '1000000',
                  '--alignMatesGapMax', '1000000',
                  '--alignSJoverhangMin', '8',
                  '--alignSJDBoverhangMin', '1',
                  '--sjdbScore', '1',
                  '--readFilesIn', '/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq']
    # Call: STAR Map
    docker_call(tool='quay.io/ucsc_cgl/star:2.4.2a--bcbd5122b69ff6ac4ef61958e47bde94001cfe80',
                work_dir=work_dir, parameters=parameters, sudo=inputs.sudo)
    # Call Samtools Index
    index_command = ['index', '/data/rnaAligned.sortedByCoord.out.bam']
    docker_call(work_dir=work_dir, parameters=index_command, sudo=inputs.sudo,
                tool='quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c')
    # fileStore
    bam_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.sortedByCoord.out.bam'))
    bai_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.sortedByCoord.out.bam.bai'))
    job.fileStore.deleteGlobalFile(r1_cutadapt)
    job.fileStore.deleteGlobalFile(r2_cutadapt)
    # Launch children and follow-on
    vcqc_id = job.addChildJobFn(variant_calling_and_qc, inputs, bam_id, bai_id, cores=2, disk='30G').rv()
    spladder_id = job.addChildJobFn(spladder, inputs, bam_id, bai_id, disk='30G').rv()
    job.addFollowOnJobFn(consolidate_output_tarballs, inputs, vcqc_id, spladder_id, disk='30G')


def variant_calling_and_qc(job, inputs, bam_id, bai_id):
    """
    Perform variant calling with samtools nad QC with CheckBias

    :param JobFunctionWrappingJob job: passed by Toil automatically
    :param Namespace inputs: Stores input arguments (see main)
    :param str bam_id: FileStore ID of bam
    :param str bai_id: FileStore ID of bam index file
    :return: FileStore ID of qc tarball
    :rtype: str
    """
    job.fileStore.logToMaster('Variant calling and QC: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    # Pull in alignment.bam from fileStore
    job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'alignment.bam'))
    job.fileStore.readGlobalFile(bai_id, os.path.join(work_dir, 'alignment.bam.bai'))
    # Download input files
    input_info = [(inputs.genome, 'genome.fa'), (inputs.positions, 'positions.tsv'),
                  (inputs.genome_index, 'genome.fa.fai'), (inputs.gtf, 'annotation.gtf'),
                  (inputs.gtf_m53, 'annotation.m53')]
    for url, fname in input_info:
        download_url(url, work_dir=work_dir, name=fname)

    # Part 1: Variant Calling
    variant_command = ['mpileup',
                       '-f', 'genome.fa',
                       '-l', 'positions.tsv',
                       '-v', 'alignment.bam',
                       '-t', 'DP,SP,INFO/AD,INFO/ADF,INFO/ADR,INFO/DPR,SP',
                       '-o', '/data/output.vcf.gz']
    docker_call(work_dir=work_dir, parameters=variant_command, sudo=inputs.sudo,
                tool='quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c')

    # Part 2: QC
    qc_command = ['-o', 'qc',
                  '-n', 'alignment.bam',
                  '-a', 'annotation.gtf',
                  '-m', 'annotation.m53']
    docker_call(work_dir=work_dir, parameters=qc_command, sudo=inputs.sudo,
                tool='jvivian/checkbias:612f129--b08a1fb6526a620bbb0304b08356f2ae7c3c0ec3')
    # Write output to fileStore and return ids
    output_tsv = glob(os.path.join(work_dir, '*counts.tsv*'))[0]
    output_vcf = os.path.join(work_dir, 'output.vcf.gz')
    tarball_files('vcqc.tar.gz', file_paths=[output_tsv, output_vcf], output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'vcqc.tar.gz'))


def spladder(job, inputs, bam_id, bai_id):
    """
    Run SplAdder to detect and quantify alternative splicing events

    :param JobFunctionWrappingJob job: passed by Toil automatically
    :param Namespace inputs: Stores input arguments (see main)
    :param str bam_id: FileStore ID of bam
    :param str bai_id: FileStore ID of bam index file
    :return: FileStore ID of SplAdder tarball
    :rtype: str
    """
    job.fileStore.logToMaster('SplAdder: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    # Pull in alignment.bam from fileStore
    job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'alignment.bam'))
    job.fileStore.readGlobalFile(bai_id, os.path.join(work_dir, 'alignment.bam.bai'))
    # Download input file
    download_url(url=inputs.gtf, work_dir=work_dir, name='annotation.gtf')
    download_url(url=inputs.gtf_pickle, work_dir=work_dir, name='annotation.gtf.pickle')
    # Call Spladder
    command = ['--insert_ir=y',
               '--insert_es=y',
               '--insert_ni=y',
               '--remove_se=n',
               '--validate_sg=n',
               '-b', 'alignment.bam',
               '-o ', '/data',
               '-a', 'annotation.gtf',
               '-v', 'y',
               '-c', '3',
               '-M', 'single',
               '-T', 'n',
               '-n', '50',
               '-P', 'y',
               '-p', 'n',
               '--sparse_bam', 'y']
    docker_call(work_dir=work_dir, parameters=command, sudo=inputs.sudo, tool='jvivian/spladder:1.0')
    # Write output to fileStore and return ids
    output_pickle = os.path.join(work_dir, ' ', 'spladder', 'genes_graph_conf3.alignment.pickle')
    if not os.path.exists(output_pickle):
        matches = []
        for root, dirnames, filenames in os.walk(work_dir):
            for filename in fnmatch.filter(filenames, '*genes_graph*'):
                matches.append(os.path.join(root, filename))
        if matches:
            output_pickle = matches[0]
        else:
            raise RuntimeError("Couldn't find genes file!")
    output_filt = os.path.join(work_dir, 'alignment.filt.hdf5')
    output = os.path.join(work_dir, 'alignment.hdf5')
    print os.listdir(work_dir)
    tarball_files('spladder.tar.gz', file_paths=[output_pickle, output_filt, output], output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'spladder.tar.gz'))


def consolidate_output_tarballs(job, inputs, vcqc_id, spladder_id):
    """
    Combine the contents of separate tarballs into one.

    :param JobFunctionWrappingJob job: passed by Toil automatically
    :param Namespace inputs: Stores input arguments (see main)
    :param str vcqc_id: FileStore ID of variant calling and QC tarball
    :param str spladder_id: FileStore ID of spladder tarball
    """
    job.fileStore.logToMaster('Consolidating files and uploading: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve IDs
    uuid = inputs.uuid
    # Unpack IDs
    # Retrieve output file paths to consolidate
    vcqc_tar = job.fileStore.readGlobalFile(vcqc_id, os.path.join(work_dir, 'vcqc.tar.gz'))
    spladder_tar = job.fileStore.readGlobalFile(spladder_id, os.path.join(work_dir, 'spladder.tar.gz'))
    # I/O
    fname = uuid + '.tar.gz' if not inputs.improper_pair else 'IMPROPER_PAIR' + uuid + '.tar.gz'
    out_tar = os.path.join(work_dir, fname)
    # Consolidate separate tarballs into one
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar in [vcqc_tar, spladder_tar]:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == vcqc_tar:
                            tarinfo.name = os.path.join(uuid, 'variants_and_qc', os.path.basename(tarinfo.name))
                        else:
                            tarinfo.name = os.path.join(uuid, 'spladder', os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
    # Move to output directory
    if inputs.output_dir:
        mkdir_p(inputs.output_dir)
        shutil.copy(out_tar, os.path.join(inputs.output_dir, os.path.basename(out_tar)))
    # Upload to S3
    if inputs.output_s3_dir:
        out_id = job.fileStore.writeGlobalFile(out_tar)
        job.addChildJobFn(s3am_upload_job, file_id=out_id, s3_dir=inputs.output_s3_dir,
                          num_cores=inputs.cores, file_name=fname, key_path=inputs.ssec, cores=inputs.cores)


def main():
    """
    This Toil pipeline aligns reads and performs alternative splicing analysis.

    Please read the README.md located in the same directory for run instructions.
    """
    # Define Parser object and add to toil
    url_prefix = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/'
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--config', required=True,
                        help='Path to configuration file for samples, one per line. UUID,URL_to_bamfile. '
                             'The URL may be a standard "http://", a "file://<abs_path>", or "s3://<bucket>/<key>"')
    parser.add_argument('--gtf', help='URL to annotation GTF file',
                        default=url_prefix + 'rnaseq_cgl/gencode.v23.annotation.gtf')
    parser.add_argument('--gtf-pickle', help='Pickled GTF file',
                        default=url_prefix + 'spladder/gencode.v23.annotation.gtf.pickle')
    parser.add_argument('--gtf-m53', help='M53 preprocessing annotation table',
                        default=url_prefix + 'spladder/gencode.v23.annotation.gtf.m53')
    parser.add_argument('--positions', help='URL to SNP positions over genes file (TSV)',
                        default=url_prefix + 'spladder/positions_fixed.tsv')
    parser.add_argument('--genome', help='URL to Genome fasta',
                        default=url_prefix + 'rnaseq_cgl/hg38_no_alt.fa')
    parser.add_argument('--genome-index', help='Index file (fai) of genome',
                        default=url_prefix + 'spladder/hg38_no_alt.fa.fai')
    parser.add_argument('--ssec', default=None, help='Path to master key used for downloading encrypted files.')
    parser.add_argument('--output-s3-dir', default=None, help='S3 Directory of the form: s3://bucket/directory')
    parser.add_argument('--output-dir', default=None, help='full path where final results will be output')
    parser.add_argument('--sudo', action='store_true', default=False,
                        help='Set flag if sudo is required to run Docker.')
    parser.add_argument('--star-index', help='URL to download STAR Index built from HG38/gencodev23 annotation.',
                        default=url_prefix + 'rnaseq_cgl/starIndex_hg38_no_alt.tar.gz')
    parser.add_argument('--fwd-3pr-adapter', help="Sequence for the FWD 3' Read Adapter.", default='AGATCGGAAGAG')
    parser.add_argument('--rev-3pr-adapter', help="Sequence for the REV 3' Read Adapter.", default='AGATCGGAAGAG')
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # Sanity Checks
    if args.config:
        assert os.path.isfile(args.config), 'Config not found at: {}'.format(args.config)
    if args.ssec:
        assert os.path.isfile(args.ssec), 'Encryption key not found at: {}'.format(args.config)
    if args.output_s3_dir:
        assert args.output_s3_dir.startswith('s3://'), 'Wrong format for output s3 directory'
    # Program checks
    for program in ['curl', 'docker']:
        assert which(program), 'Program "{}" must be installed on every node.'.format(program)

    Job.Runner.startToil(Job.wrapJobFn(parse_input_samples, args), args)
