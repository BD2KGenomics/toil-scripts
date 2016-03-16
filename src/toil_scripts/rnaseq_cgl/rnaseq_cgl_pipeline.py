#!/usr/bin/env python2.7
"""
UCSC Computational Genomics Lab RNA-seq Pipeline
Author: John Vivian
Affiliation: UC Santa Cruz Genomics Institute

Please see the README.md in the same directory

Structure of RNA-Seq Pipeline (per sample)

                 0
                 |
                 1
                 |
                 2 - - - - > 7 -> 8
                / \
               3   6
               |
               4
               |
               5

0 = Download sample
1 = Unpack/Merge fastqs
2 = CutAdapt
3 = STAR Alignment
4 = RSEM Quantification
5 = RSEM Post-processing
6 = Kallisto
7 = Consoliate output (into a tarball)
8 = upload results to S3 (optional)
=======================================
Dependencies
Curl:       apt-get install curl
Docker:     wget -qO- https://get.docker.com/ | sh
Toil:       pip install toil
Boto:       pip install boto
"""
import argparse
import glob
import logging
import multiprocessing
import os
import shutil
import subprocess
import tarfile
from contextlib import closing
from subprocess import PIPE
from urlparse import urlparse

from toil.job import Job

from toil_scripts.lib import copy_to_output_dir, flatten
from toil_scripts.lib.files import mkdir_p
from toil_scripts.lib.files import tarball_files
from toil_scripts.lib.jobs import map_job
from toil_scripts.lib.programs import docker_call, which
from toil_scripts.lib.urls import download_url_job, s3am_upload, s3am_upload_job, download_url

logging.basicConfig(level=logging.INFO)
_log = logging.getLogger(__name__)


def parse_input_samples(job, inputs):
    """
    Parse the various input formats then launch job batcher

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Contains argument inputs
    """
    job.fileStore.logToMaster('Parsing input Samples')
    samples = None
    if inputs.config:
        samples = parse_config(inputs.config)
    elif inputs.sample_urls:
        samples = [[os.path.basename(sample).split('.')[0], sample] for sample in inputs.sample_urls]
    elif inputs.genetorrent:
        samples = parse_config(inputs.genetorrent, genetorrent=True)
    elif inputs.dir:
        files = os.listdir(inputs.dir)
        samples = [[os.path.splitext(os.path.basename(f))[0],
                    'file://' + os.path.join(inputs.dir, os.path.basename(f))] for f in files]
    # Pass to batcher to spawn tree of jobs
    job.addChildJobFn(map_job, download_sample, samples, inputs)


def download_sample(job, sample, inputs):
    """
    Downloads sample and stores unique attributes of that sample

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param tuple sample:
    :param Namespace inputs:
    """
    sample_input = argparse.Namespace(**vars(inputs))
    uuid, url = None, None
    tar_id, r1_id, r2_id = None, None, None
    if len(sample) == 2:
        uuid, url = sample
        sample_input.tar = url
    if len(sample) == 3:
        uuid = sample[0]
        url = sample[1:]
        sample_input.R1 = url[0]
        sample_input.R2 = url[1]
    assert uuid and url, 'Issue with sample configuration retrieval: {}'.format(sample)
    job.fileStore.logToMaster('Downloading sample: {}'.format(uuid))
    # Update values unique to sample
    sample_input.uuid = uuid
    sample_input.output_dir = os.path.join(inputs.output_dir, uuid) if inputs.output_dir else None
    sample_input.cores = multiprocessing.cpu_count()
    # Download or locate local file and place in the jobStore
    if inputs.genetorrent:
        cghub_key_path = inputs.genetorrent_key
        tar_id = job.addChildJobFn(download_url_job, url, cghub_key_path=cghub_key_path, disk='20G').rv()
    elif type(url) is list and len(url) == 2:
        if urlparse(url[0]) == 'file':
            r1_id = job.fileStore.writeGlobalFile(urlparse(url[0]).path)
            r2_id = job.fileStore.writeGlobalFile(urlparse(url[1]).path)
        else:
            if inputs.ssec:
                r1_id = job.addChildJobFn(download_url_job, url[0], s3_key_path=inputs.ssec, disk='20G').rv()
                r2_id = job.addChildJobFn(download_url_job, url[1], s3_key_path=inputs.ssec, disk='20G').rv()
            else:
                r1_id = job.addChildJobFn(download_url_job, url[0], disk='20G').rv()
                r2_id = job.addChildJobFn(download_url_job, url[1], disk='20G').rv()
    elif urlparse(url).scheme == 'file':
        tar_id = job.fileStore.writeGlobalFile(urlparse(url).path)
    else:
        if inputs.ssec:
            tar_id = job.addChildJobFn(download_url_job, url, s3_key_path=inputs.ssec, disk='20G').rv()
        else:
            tar_id = job.addChildJobFn(download_url_job, url, disk='20G').rv()
    job.addFollowOnJobFn(static_dag_launchpoint, sample_input, tar_id, r1_id, r2_id)


def static_dag_launchpoint(job, inputs, tar_id, r1_id, r2_id):
    """
    Statically define rest of pipeline

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Argparse Namespace object containing argument inputs
    :param str tar_id:
    :param str r1_id:
    :param str r2_id:
    """
    job.fileStore.logToMaster('Defining static DAG')
    disk = '2G' if inputs.ci_test else '100G'
    if tar_id:
        a = job.wrapJobFn(process_sample_tar, inputs, tar_id, disk=disk).encapsulate()
    else:
        a = job.wrapJobFn(cutadapt, inputs, r1_id, r2_id, disk=disk).encapsulate()
    b = job.wrapJobFn(consolidate_output, inputs, a.rv(), disk='2G')
    # Take advantage of "encapsulate" to simplify pipeline wiring
    job.addChild(a)
    a.addChild(b)


def process_sample_tar(job, inputs, tar_id):
    """
    Converts sample.tar(.gz) into a fastq pair or single fastq if single-ended.
    Due to different naming standards and edge conditions... BEWARE: HERE BE DRAGONS

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Argparse Namespace object containing argument inputs
    :param str tar_id: fileStoreID of the tarball
    :return str: FileStoreID from Cutadapt
    """
    job.fileStore.logToMaster('Processing sample: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    r_id, r1_id, r2_id = None, None, None
    # I/O
    job.fileStore.readGlobalFile(tar_id, os.path.join(work_dir, 'sample.tar'))
    tar_path = os.path.join(work_dir, 'sample.tar')
    # Untar File and concat
    p = subprocess.Popen(['tar', '-xvf', tar_path, '-C', work_dir], stderr=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        # Handle error if tar archive is corrupt
        if 'EOF' in stderr:
            error_path = os.path.join(work_dir, inputs.uuid + '.error.txt')
            with open(error_path, 'w') as f:
                f.write(stderr)
                f.write(stdout)
            if inputs.s3_dir:
                s3am_upload(error_path, inputs.s3_dir)
        else:
            raise subprocess.CalledProcessError
    else:
        os.remove(os.path.join(work_dir, 'sample.tar'))
        # Grab files from tarball
        fastqs = []
        for root, subdir, files in os.walk(work_dir):
            fastqs.extend([os.path.join(root, x) for x in files])
        # Check for read 1 and read 2 files
        r1 = sorted([x for x in fastqs if '_1' in x])
        r2 = sorted([x for x in fastqs if '_2' in x])
        if not r1 or not r2:
            # Check if using a different standard
            r1 = sorted([x for x in fastqs if 'R1' in x])
            r2 = sorted([x for x in fastqs if 'R2' in x])
        # Prune file name matches from each list
        if len(r1) > len(r2):
            r1 = [x for x in r1 if x not in r2]
        elif len(r2) > len(r1):
            r2 = [x for x in r2 if x not in r1]
        if not r1 or not r2:
            # Sample is assumed to be single-ended
            if fastqs[0].endswith('.gz'):
                with open(os.path.join(work_dir, 'R.fastq'), 'w') as f:
                    subprocess.check_call(['zcat'] + fastqs, stdout=f)
            elif len(fastqs) > 1:
                with open(os.path.join(work_dir, 'R.fastq'), 'w') as f:
                    subprocess.check_call(['cat'] + fastqs, stdout=f)
            else:
                shutil.move(fastqs[0], os.path.join(work_dir, 'R.fastq'))
            r_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R.fastq'))
        else:
            # Sample is assumed to be paired end
            if r1[0].endswith('.gz') and r2[0].endswith('.gz'):
                with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
                    p1 = subprocess.Popen(['zcat'] + r1, stdout=f1)
                with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
                    p2 = subprocess.Popen(['zcat'] + r2, stdout=f2)
                p1.wait()
                p2.wait()
            elif len(r1) > 1 and len(r2) > 1:
                with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
                    p1 = subprocess.Popen(['cat'] + r1, stdout=f1)
                with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
                    p2 = subprocess.Popen(['cat'] + r2, stdout=f2)
                p1.wait()
                p2.wait()
            else:
                shutil.move(r1[0], os.path.join(work_dir, 'R1.fastq'))
                shutil.move(r2[0], os.path.join(work_dir, 'R2.fastq'))
            r1_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1.fastq'))
            r2_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2.fastq'))
        job.fileStore.deleteGlobalFile(tar_id)
        # Start cutadapt step
        disk = '2G' if inputs.ci_test else '125G'
        return job.addChildJobFn(cutadapt, inputs, r_id, r1_id, r2_id, disk=disk).rv()


def cutadapt(job, inputs, r_id, r1_id, r2_id):
    """
    Adapter triming for RNA-seq data

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Argparse Namespace object containing argument inputs
    :param str r_id: FileStoreID of single-end r.fastq
    :param str r1_id: FileStoreID of fastq (pair 1)
    :param str r2_id: FileStoreID of fastq (pair 2)
    :return tuple: FileStoreIDs and sample information
    """
    job.fileStore.logToMaster('CutAdapt: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    inputs.single_end, inputs.improper_pair = None, None
    # Retrieve files
    parameters = ['-a', inputs.fwd_3pr_adapter,
                  '-m', '35']
    if r_id:
        inputs.single_end = True
        job.fileStore.readGlobalFile(r_id, os.path.join(work_dir, 'R.fastq'))
        parameters.extend(['-o', '/data/R_cutadapt.fastq', '/data/R.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
        job.fileStore.readGlobalFile(r2_id, os.path.join(work_dir, 'R2.fastq'))
        parameters.extend(['-A', inputs.rev_3pr_adapter,
                           '-o', '/data/R1_cutadapt.fastq',
                           '-p', '/data/R2_cutadapt.fastq',
                           '/data/R1.fastq', '/data/R2.fastq'])
    # Call: CutAdapt
    base_docker_call = 'docker run --log-driver=none --rm -v {}:/data'.format(work_dir).split()
    if inputs.sudo:
        base_docker_call = ['sudo'] + base_docker_call
    tool = 'quay.io/ucsc_cgl/cutadapt:1.9--6bd44edd2b8f8f17e25c5a268fedaab65fa851d2'
    p = subprocess.Popen(base_docker_call + [tool] + parameters, stderr=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        if 'improperly paired' in stderr:
            inputs.improper_pair = True
            if r_id:
                shutil.move(os.path.join(work_dir, 'R.fastq'), os.path.join(work_dir, 'R_cutadapt.fastq'))
            else:
                shutil.move(os.path.join(work_dir, 'R1.fastq'), os.path.join(work_dir, 'R1_cutadapt.fastq'))
                shutil.move(os.path.join(work_dir, 'R2.fastq'), os.path.join(work_dir, 'R2_cutadapt.fastq'))
        else:
            logging.error('Stdout: {}\n\nStderr: {}'.format(stdout, stderr))
            raise subprocess.CalledProcessError(p.returncode, parameters, stderr)
    # Write to fileStore
    if r_id:
        r_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R_cutadapt.fastq'))
        r1_cut_id = None
        r2_cut_id = None
    else:
        r_cut_id = None
        r1_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1_cutadapt.fastq'))
        r2_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2_cutadapt.fastq'))
    # start STAR and Kallisto steps
    disk = '2G' if inputs.ci_test else '100G'
    memory = '6G' if inputs.ci_test else '40G'
    rsem_output = job.addChildJobFn(star, inputs, r_cut_id, r1_cut_id, r2_cut_id,
                                    cores=inputs.cores, disk=disk, memory=memory).rv()
    cores = min(inputs.cores, 16)
    kallisto_output = job.addChildJobFn(kallisto, inputs, r_cut_id, r1_cut_id, r2_cut_id, cores=cores, disk=disk).rv()
    return rsem_output, kallisto_output


def kallisto(job, inputs, r_cut_id, r1_cut_id, r2_cut_id):
    """
    RNA quantification via Kallisto

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Argparse Namespace object containing argument inputs
    :param str r_cut_id: FileStoreID of single-end r.fastq
    :param str r1_cut_id: FileStoreID of fastq (pair 1)
    :param str r2_cut_id: FileStoreID of fastq (pair 2)
    :return str: FileStoreID from Kallisto
    """
    job.fileStore.logToMaster('Kallisto: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    download_url(url=inputs.kallisto_index, name='kallisto_hg38.idx', work_dir=work_dir)
    # Retrieve files
    parameters = ['quant',
                  '-i', '/data/kallisto_hg38.idx',
                  '-t', str(inputs.cores),
                  '-o', '/data/',
                  '-b', '100']
    if r_cut_id:
        job.fileStore.readGlobalFile(r_cut_id, os.path.join(work_dir, 'R_cutadapt.fastq'))
        parameters.extend(['--single', '-l', '200', '-s', '15', '/data/R_cutadapt.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_cut_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        job.fileStore.readGlobalFile(r2_cut_id, os.path.join(work_dir, 'R2_cutadapt.fastq'))
        parameters.extend(['/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq'])
    # Call: Kallisto
    docker_call(tool='quay.io/ucsc_cgl/kallisto:0.42.4--35ac87df5b21a8e8e8d159f26864ac1e1db8cf86',
                work_dir=work_dir, parameters=parameters, sudo=inputs.sudo)
    # Tar output files together and store in fileStore
    output_files = [os.path.join(work_dir, x) for x in ['run_info.json', 'abundance.tsv', 'abundance.h5']]
    tarball_files(tar_name='kallisto.tar.gz', file_paths=output_files, output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'kallisto.tar.gz'))


def star(job, inputs, r_cut_id, r1_cut_id, r2_cut_id):
    """
    Performs alignment of fastqs to bam via STAR

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Argparse Namespace object containing argument inputs
    :param str r_cut_id: FileStoreID of single-end r.fastq
    :param str r1_cut_id: FileStoreID of fastq (pair 1)
    :param str r2_cut_id: FileStoreID of fastq (pair 2)
    :return str: FileStoreID from RSEM
    """
    job.fileStore.logToMaster('STAR: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    download_url(url=inputs.star_index, name='starIndex.tar.gz', work_dir=work_dir)
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'starIndex.tar.gz'), '-C', work_dir])
    # Parameters and input retrieval
    parameters = ['--runThreadN', str(inputs.cores),
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
                  '--sjdbScore', '1']
    if inputs.wiggle:
        parameters.extend(['--outWigType', 'bedGraph',
                           '--outWigStrand', 'Unstranded',
                           '--outWigReferencesPrefix', 'chr'])
    if r_cut_id:
        job.fileStore.readGlobalFile(r_cut_id, os.path.join(work_dir, 'R_cutadapt.fastq'))
        parameters.extend(['--readFilesIn', '/data/R_cutadapt.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_cut_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        job.fileStore.readGlobalFile(r2_cut_id, os.path.join(work_dir, 'R2_cutadapt.fastq'))
        parameters.extend(['--readFilesIn', '/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq'])
    # Call: STAR Map
    docker_call(tool='quay.io/ucsc_cgl/star:2.4.2a--bcbd5122b69ff6ac4ef61958e47bde94001cfe80',
                work_dir=work_dir, parameters=parameters, sudo=inputs.sudo)
    # Write to fileStore
    bam_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.toTranscriptome.out.bam'))
    # Save Wiggle File
    if inputs.wiggle and inputs.s3_dir:
        wiggles = [os.path.basename(x) for x in glob.glob(os.path.join(work_dir, '*.bg'))]
        # Rename extension
        for wiggle in wiggles:
            shutil.move(os.path.join(work_dir, wiggle),
                        os.path.join(work_dir, os.path.splitext(wiggle)[0] + '.bedGraph'))
        wiggles = [os.path.join(work_dir, x) for x in [os.path.splitext(x)[0] + '.bedGraph' for x in wiggles]]
        tarball_files('wiggle.tar.gz', file_paths=wiggles, output_dir=work_dir)
        wiggle_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'wiggle.tar.gz'))
        job.addChildJobFn(s3am_upload_job, file_id=wiggle_id, file_name='wiggle.tar.gz',
                          s3_dir=inputs.s3_dir, num_cores=inputs.cores)
    if inputs.save_bam and inputs.s3_dir:
        sorted_bam_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.sortedByCoord.out.bam'))
        job.addChildJobFn(s3am_upload_job, file_id=sorted_bam_id, file_name=inputs.uuid + '.sorted.bam',
                          s3_dir=inputs.s3_dir, num_cores=inputs.cores, s3_key_path=inputs.ssec)
    # RSEM doesn't tend to use more than 16 cores
    cores = min(inputs.cores, 16)
    disk = '2G' if inputs.ci_test else '50G'
    return job.addChildJobFn(rsem, inputs, bam_id, cores=cores, disk=disk).rv()


def rsem(job, inputs, bam_id):
    """
    RNA quantification with RSEM

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Argparse Namespace object containing argument inputs
    :param str bam_id: FileStoreID of transcriptome bam for quantification
    :return str: FileStoreID from RSEM postprocess
    """
    job.fileStore.logToMaster('RSEM: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    cores = 16 if inputs.cores >= 16 else inputs.cores
    download_url(url=inputs.rsem_ref, name='rsem_ref_hg38.tar.gz', work_dir=work_dir)
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'rsem_ref_hg38.tar.gz'), '-C', work_dir])
    # I/O
    job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'transcriptome.bam'))
    prefix = 'rsem'
    # Call: RSEM
    parameters = ['--quiet',
                  '--no-qualities',
                  '-p', str(cores),
                  '--forward-prob', '0.5',
                  '--seed-length', '25',
                  '--fragment-length-mean', '-1.0',
                  '--bam', '/data/transcriptome.bam',
                  '/data/rsem_ref_hg38/hg38',
                  prefix]
    if not inputs.single_end:
        parameters = ['--paired-end'] + parameters
    docker_call(tool='quay.io/ucsc_cgl/rsem:1.2.25--d4275175cc8df36967db460b06337a14f40d2f21',
                parameters=parameters, work_dir=work_dir, sudo=inputs.sudo)
    os.rename(os.path.join(work_dir, prefix + '.genes.results'), os.path.join(work_dir, 'rsem_gene.tab'))
    os.rename(os.path.join(work_dir, prefix + '.isoforms.results'), os.path.join(work_dir, 'rsem_isoform.tab'))
    # Write to FileStore
    rsem_gene_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_gene.tab'))
    rsem_isoform_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_isoform.tab'))
    job.fileStore.deleteGlobalFile(bam_id)
    # Run child jobs
    return job.addChildJobFn(rsem_postprocess, inputs, rsem_gene_id, rsem_isoform_id).rv()


def rsem_postprocess(job, inputs, rsem_gene_id, rsem_isoform_id):
    """
    Parses RSEMs output to produce the separate .tab files (TPM, FPKM, counts) for both gene and isoform.
    These are two-column files: Genes and Quantifications.
    HUGO files are also provided that have been mapped from Gencode/ENSEMBLE names.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Argparse Namespace object containing argument inputs
    :param str rsem_gene_id: FileStoreID of rsem_gene_ids
    :param str rsem_isoform_id: FileStoreID of rsem_isoform_ids
    :return str: FileStoreID from RSEM post process tarball
    """
    job.fileStore.logToMaster('RSEM Postprocess: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    # I/O
    job.fileStore.readGlobalFile(rsem_gene_id, os.path.join(work_dir, 'rsem_gene.tab'))
    job.fileStore.readGlobalFile(rsem_isoform_id, os.path.join(work_dir, 'rsem_isoform.tab'))
    # Convert RSEM files into individual .tab files.
    docker_call(tool='jvivian/rsem_postprocess', parameters=[inputs.uuid], work_dir=work_dir, sudo=inputs.sudo)
    os.rename(os.path.join(work_dir, 'rsem_gene.tab'), os.path.join(work_dir, 'rsem_genes.results'))
    os.rename(os.path.join(work_dir, 'rsem_isoform.tab'), os.path.join(work_dir, 'rsem_isoforms.results'))
    output_files = ['rsem.genes.norm_counts.tab', 'rsem.genes.raw_counts.tab', 'rsem.genes.norm_fpkm.tab',
                    'rsem.genes.norm_tpm.tab', 'rsem.isoform.norm_counts.tab', 'rsem.isoform.raw_counts.tab',
                    'rsem.isoform.norm_fpkm.tab', 'rsem.isoform.norm_tpm.tab', 'rsem_genes.results',
                    'rsem_isoforms.results']
    # Perform HUGO gene / isoform name mapping
    genes = [x for x in output_files if 'gene' in x]
    isoforms = [x for x in output_files if 'isoform' in x]
    command = ['-g'] + genes + ['-i'] + isoforms
    docker_call(tool='jvivian/gencode_hugo_mapping', parameters=command, work_dir=work_dir, sudo=inputs.sudo)
    hugo_files = [os.path.splitext(x)[0] + '.hugo' + os.path.splitext(x)[1] for x in output_files]
    # Create tarballs for outputs
    tarball_files('rsem.tar.gz', file_paths=[os.path.join(work_dir, x) for x in output_files], output_dir=work_dir)
    tarball_files('rsem_hugo.tar.gz', [os.path.join(work_dir, x) for x in hugo_files], output_dir=work_dir)
    rsem_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem.tar.gz'))
    hugo_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_hugo.tar.gz'))
    return rsem_id, hugo_id, inputs.improper_pair, inputs.single_end


def consolidate_output(job, inputs, output_ids_and_info):
    """
    Combines the contents of the outputs into one tarball and places in output directory or s3

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace inputs: Argparse Namespace object containing argument inputs
    :param list output_ids_and_info: Array of output IDs to consolidate
    """
    job.fileStore.logToMaster('Consolidating input: {}'.format(inputs.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve IDs
    rsem_id, hugo_id, improper_pair, single_end, kallisto_id = flatten(output_ids_and_info)
    # Retrieve output file paths to consolidate
    rsem_tar = job.fileStore.readGlobalFile(rsem_id, os.path.join(work_dir, 'rsem.tar.gz'))
    kallisto_tar = job.fileStore.readGlobalFile(kallisto_id, os.path.join(work_dir, 'kallisto.tar.gz'))
    hugo_tar = job.fileStore.readGlobalFile(hugo_id, os.path.join(work_dir, 'rsem_hugo.tar.gz'))
    # I/O
    if improper_pair:
        inputs.uuid = 'IMPROPERLY_PAIRED.{}'.format(inputs.uuid)
    if single_end:
        inputs.uuid = 'SINGLE-END.{}'.format(inputs.uuid)
    out_tar = os.path.join(work_dir, inputs.uuid + '.tar.gz')
    # Consolidate separate tarballs into one as streams (avoids unnecessary untaring)
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar in [rsem_tar, kallisto_tar, hugo_tar]:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == rsem_tar:
                            tarinfo.name = os.path.join(inputs.uuid, 'RSEM', os.path.basename(tarinfo.name))
                        elif tar == hugo_tar:
                            tarinfo.name = os.path.join(inputs.uuid, 'RSEM', 'Hugo', os.path.basename(tarinfo.name))
                        else:
                            tarinfo.name = os.path.join(inputs.uuid, 'Kallisto', os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
        if improper_pair:
            with open(os.path.join(work_dir, 'WARNING.txt'), 'w') as f:
                f.write('cutadapt: error: Reads are improperly paired. Uneven number of reads in fastq pair.')
            f_out.add(os.path.join(work_dir, 'WARNING.txt'))
    # Move to output directory
    if inputs.output_dir:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(inputs.uuid, inputs.output_dir))
        mkdir_p(inputs.output_dir)
        copy_to_output_dir(inputs.output_dir, fpaths=[os.path.join(work_dir, inputs.uuid + '.tar.gz')])
    # Upload to S3
    if inputs.s3_dir:
        job.fileStore.logToMaster('Uploading {} to S3: {}'.format(inputs.uuid, inputs.s3_dir))
        s3am_upload(fpath=out_tar, s3_dir=inputs.s3_dir, num_cores=inputs.cores)


# Parser, convenience functions, and custom exceptions
def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    prefix = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rnaseq_cgl/'
    group.add_argument('-c', '--config', default=None,
                       help='Path to CSV. One sample per line with the format: uuid,url (see example config). '
                            'URLs follow the format of: https://www.address.com/file.tar or file:///full/path/file.tar.'
                            ' Samples must be tarfiles that contain fastq files. \n\nAlternatively, if your data is'
                            'not tarred, you can create a CSV of one sample per line: uuid,url1,url2 .  Where url1'
                            'and url2 refer to urls to the pair of fastq files.')
    group.add_argument('-d', '--dir', default=None,
                       help='Path to directory of samples. Samples must be tarfiles that contain fastq files. '
                            'The UUID for the sample will be derived from the file basename.')
    group.add_argument('-s', '--sample-urls', nargs='+', default=None, type=str,
                       help='Sample URLs (any number). Samples must be tarfiles that contain fastq files. '
                            'URLs follow the format of: https://www.address.com/file.tar or file:///full/path/file.tar.'
                            'The UUID for the sample will be derived from the file basename.')
    group.add_argument('-g', '--genetorrent', default=None,
                       help='Path to a file with one analysis ID per line for data hosted on CGHub.')
    parser.add_argument('-k', '--genetorrent-key', default=None,
                        help='Path to a CGHub key that has access to the TCGA data being requested. An exception will'
                             'be thrown if "-g" is set but not this argument.')
    parser.add_argument('--star-index', help='URL to download STAR Index built from HG38/gencodev23 annotation.',
                        default=prefix + 'starIndex_hg38_no_alt.tar.gz')
    parser.add_argument('--kallisto-index', help='URL to download Kallisto Index built from HG38 transcriptome.',
                        default=prefix + 'kallisto_hg38.idx')
    parser.add_argument('--rsem-ref', help='URL to download RSEM reference built from HG38/gencodev23.',
                        default=prefix + 'rsem_ref_hg38_no_alt.tar.gz')
    parser.add_argument('--wiggle', default=None, action='store_true', help='Uploads a wiggle from STAR to S3.')
    parser.add_argument('--save-bam', default=None, action='store_true', help='Uploads aligned BAM to S3.')
    parser.add_argument('--fwd-3pr-adapter', help="Sequence for the FWD 3' Read Adapter.", default='AGATCGGAAGAG')
    parser.add_argument('--rev-3pr-adapter', help="Sequence for the REV 3' Read Adapter.", default='AGATCGGAAGAG')
    parser.add_argument('--ssec', default=None, help='Path to Key File for SSE-C Encryption')
    parser.add_argument('--output-dir', default=None, help='full path where final results will be output')
    parser.add_argument('--s3-dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                       's3://cgl-driver-projects/ckcc/rna-seq-samples')
    parser.add_argument('--sudo', dest='sudo', action='store_true', default=False,
                        help='Docker usually needs sudo to execute locally, but not when running Mesos or when '
                             'the user is a member of a Docker group.')
    parser.add_argument('--ci-test', action='store_true', default=False,
                        help='Pass flag indicating continuous integration testing for resource requirements')
    return parser


def parse_config(path_to_config, genetorrent=None):
    """
    Parses config file, which should follow format specified in --help.

    :param bool genetorrent: If True, parses config in genetorrent format
    :param str path_to_config: Path to configuration file
    :return list[list]: Samples containing the uuid and url of each sample
    """
    samples = []
    with open(path_to_config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                sample = [line.strip()] * 2 if genetorrent else line.strip().split(',')
                if len(sample) < 2 or len(sample) > 3:
                    raise IllegalArgumentException('Error: Bad Format! Check documentation. {}'.format(sample))
                samples.append(sample)
    return samples


class IllegalArgumentException(Exception):
    pass


def main():
    """
    This is a Toil pipeline for the UCSC CGL's RNA-seq pipeline.
    RNA-seq fastqs are combined, aligned, and quantified with 2 different methods.

    Please read the README.md located in the same directory for run instructions.
    """
    # Define Parser object and add to toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # Sanity Checks
    if args.ssec:
        assert os.path.isfile(args.ssec)
    if args.config:
        assert os.path.isfile(args.config)
    if args.dir:
        assert os.path.isdir(args.dir)
        assert os.listdir(args.dir) != []
    if args.sample_urls:
        assert args.sample_urls != []
    if args.genetorrent:
        assert os.path.isfile(args.genetorrent)
    if args.genetorrent_key:
        assert os.path.isfile(args.genetorrent_key)
    # If inputs are files, ensure path is accessible
    for url in [args.star_index, args.kallisto_index, args.rsem_ref]:
        if urlparse(url).scheme == 'file':
            assert os.path.isfile(urlparse(url).path)
    # Genetorrent key must be provided along with genetorrent option
    if args.genetorrent and not args.genetorrent_key:
        raise RuntimeError("Cannot supply -genetorrent without -genetorrent_key")
    # Program checks
    for program in ['curl', 'docker']:
        assert which(program), 'Program "{}" must be installed on every node.'.format(program)

    # Start Pipeline
    Job.Runner.startToil(Job.wrapJobFn(parse_input_samples, args), args)


if __name__ == '__main__':
    main()
