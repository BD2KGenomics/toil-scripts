#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import glob
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import tarfile
import textwrap
from ast import literal_eval
from contextlib import closing
from subprocess import PIPE

from toil.job import Job

from toil_scripts.lib import flatten, require, UserError
from toil_scripts.lib.files import mkdir_p
from toil_scripts.lib.files import tarball_files, move_files
from toil_scripts.lib.jobs import map_job
from toil_scripts.lib.programs import docker_call, which
from toil_scripts.lib.urls import download_url_job, s3am_upload, s3am_upload_job, download_url


# Start of pipeline
def download_sample(job, sample, config):
    """
    Download sample and store unique attributes

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param tuple sample:
    :param Namespace config:
    """
    # Create copy of config that is sample specific
    config = argparse.Namespace(**vars(config))
    config.file_type, config.paired, config.uuid, config.url = sample
    config.paired = True if config.paired == 'paired' else False
    config.single = True if config.paired == 'single' else False
    config.cores = multiprocessing.cpu_count()
    job.fileStore.logToMaster('Downloading sample: {}'.format(config.uuid))
    # Download or locate local file and place in the jobStore
    tar_id, r1_id, r2_id = None, None, None
    if config.file_type == 'tar':
        tar_id = job.addChildJobFn(download_url_job, config.url, cghub_key_path=config.gt_key,
                                   s3_key_path=config.ssec, disk='20G').rv()
    else:
        if config.paired:
            require(len(config.url.split(',')==2), 'Fastq pairs must have 2 URLS separated by comma')
            r1_url, r2_url = config.url.split(',')
            r1_id = job.addChildJobFn(download_url_job, r1_url, cghub_key_path=config.gt_key,
                                      s3_key_path=config.ssec, disk='20G').rv()
            r2_id = job.addChildJobFn(download_url_job, r2_url, cghub_key_path=config.gt_key,
                                      s3_key_path=config.ssec, disk='20G').rv()
        else:
            r1_id = job.addChildJobFn(download_url_job, config.url, cghub_key_path=config.gt_key,
                                      s3_key_path=config.ssec, disk='20G').rv()
    job.addFollowOnJobFn(statically_define_pipeline, config, tar_id, r1_id, r2_id)


def statically_define_pipeline(job, config, tar_id, r1_id, r2_id):
    """
    Statically define remainder of pipeline

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str tar_id: FileStoreID of sample tar
    :param str r1_id: FileStoreID of sample read 1
    :param str r2_id: FileStoreID of sample read 2
    """
    job.fileStore.logToMaster('Defining static DAG')
    disk = '2G' if config.ci_test else '100G'
    if tar_id:
        a = job.wrapJobFn(process_sample_tar, config, tar_id, disk=disk).encapsulate()
    else:
        a = job.wrapJobFn(cutadapt, config, r1_id, r2_id, disk=disk).encapsulate()
    b = job.wrapJobFn(consolidate_output, config, a.rv(), disk='2G')
    # Take advantage of "encapsulate" to simplify pipeline wiring
    job.addChild(a)
    a.addChild(b)


def process_sample_tar(job, config, tar_id):
    """
    Converts sample.tar(.gz) into a fastq pair or single fastq if single-ended.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str tar_id: fileStoreID of the tarball
    :return: FileStoreID from Cutadapt
    :rtype: str
    """
    job.fileStore.logToMaster('Processing sample: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    r1_id, r2_id = None, None
    # I/O
    job.fileStore.readGlobalFile(tar_id, os.path.join(work_dir, 'sample.tar'))
    tar_path = os.path.join(work_dir, 'sample.tar')
    # Untar File and concat
    subprocess.check_call(['tar', '-xvf', tar_path, '-C', work_dir], stderr=PIPE, stdout=PIPE)
    os.remove(os.path.join(work_dir, 'sample.tar'))
    fastqs = []
    for root, subdir, files in os.walk(work_dir):
        fastqs.extend([os.path.join(root, x) for x in files])
    if config.paired:
        r1 = sorted([x for x in fastqs if 'R1' in x])
        r2 = sorted([x for x in fastqs if 'R2' in x])
        if not r1 or not r2:
            r1 = sorted([x for x in fastqs if '_1' in x])
            r2 = sorted([x for x in fastqs if '_2' in x])
        require(len(r1) == len(r2), 'Check fastq naming, uneven number of pairs found: r1: {}, r2: {}'.format(r1, r2))
        # Concatenate fastqs
        command = 'zcat' if r1[0].endswith('.gz') and r2[0].endswith('.gz') else 'cat'
        with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
            p1 = subprocess.Popen([command] + r1, stdout=f1)
        with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
            p2 = subprocess.Popen([command] + r2, stdout=f2)
        p1.wait()
        p2.wait()
        r1_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1.fastq'))
        r2_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2.fastq'))
    elif config.single:
        command = 'zcat' if fastqs[0].endswith('.gz') else 'cat'
        with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f:
            subprocess.check_call([command] + fastqs, stdout=f)
        r1_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1.fastq'))
    job.fileStore.deleteGlobalFile(tar_id)
    # Start cutadapt step
    disk = '2G' if config.ci_test else '125G'
    return job.addChildJobFn(cutadapt, config, r1_id, r2_id, disk=disk).rv()


def cutadapt(job, config, r1_id, r2_id):
    """
    Adapter triming for RNA-seq data

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str r1_id: FileStoreID of fastq read 1
    :param str r2_id: FileStoreID of fastq read 2
    :return: FileStoreIDs and sample information
    :rtype: tuple
    """
    job.fileStore.logToMaster('CutAdapt: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    config.improper_pair = None
    # Retrieve files
    parameters = ['-a', config.fwd_3pr_adapter,
                  '-m', '35']
    if config.single:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
        parameters.extend(['-o', '/data/R1_cutadapt.fastq', '/data/R1.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
        job.fileStore.readGlobalFile(r2_id, os.path.join(work_dir, 'R2.fastq'))
        parameters.extend(['-A', config.rev_3pr_adapter,
                           '-o', '/data/R1_cutadapt.fastq',
                           '-p', '/data/R2_cutadapt.fastq',
                           '/data/R1.fastq', '/data/R2.fastq'])
    # Call: CutAdapt
    base_docker_call = 'docker run --log-driver=none --rm -v {}:/data'.format(work_dir).split()
    if config.sudo:
        base_docker_call = ['sudo'] + base_docker_call
    tool = 'quay.io/ucsc_cgl/cutadapt:1.9--6bd44edd2b8f8f17e25c5a268fedaab65fa851d2'
    p = subprocess.Popen(base_docker_call + [tool] + parameters, stderr=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        if 'improperly paired' in stderr:
            config.improper_pair = True
            if config.paired == 'single':
                shutil.move(os.path.join(work_dir, 'R1.fastq'), os.path.join(work_dir, 'R1_cutadapt.fastq'))
            else:
                shutil.move(os.path.join(work_dir, 'R1.fastq'), os.path.join(work_dir, 'R1_cutadapt.fastq'))
                shutil.move(os.path.join(work_dir, 'R2.fastq'), os.path.join(work_dir, 'R2_cutadapt.fastq'))
        else:
            logging.error('Stdout: {}\n\nStderr: {}'.format(stdout, stderr))
            raise subprocess.CalledProcessError(p.returncode, parameters, stderr)
    # Write to fileStore
    if config.single:
        r1_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1_cutadapt.fastq'))
        r2_cut_id = None
    else:
        r1_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1_cutadapt.fastq'))
        r2_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2_cutadapt.fastq'))
    # start STAR and Kallisto steps
    disk = '2G' if config.ci_test else '100G'
    memory = '6G' if config.ci_test else '40G'
    # Based on inputs, run Kallisto and/or STAR/RSEM
    rsem_output, kallisto_output = None, None
    if config.rsem_ref:
        rsem_output = job.addChildJobFn(star, config, r1_cut_id, r2_cut_id,
                                        cores=config.cores, disk=disk, memory=memory).rv()
    if config.kallisto_index:
        cores = min(config.cores, 16)
        kallisto_output = job.addChildJobFn(kallisto, config, r1_cut_id, r2_cut_id, cores=cores, disk=disk).rv()
    return rsem_output, kallisto_output, config.single, config.improper_pair


def kallisto(job, config, r1_cut_id, r2_cut_id):
    """
    RNA quantification via Kallisto

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str r1_cut_id: FileStoreID of fastq (pair 1)
    :param str r2_cut_id: FileStoreID of fastq (pair 2)
    :return: FileStoreID from Kallisto
    :rtype: str
    """
    job.fileStore.logToMaster('Kallisto: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    download_url(url=config.kallisto_index, name='kallisto_hg38.idx', work_dir=work_dir)
    # Retrieve files
    parameters = ['quant',
                  '-i', '/data/kallisto_hg38.idx',
                  '-t', str(config.cores),
                  '-o', '/data/',
                  '-b', '100']
    if config.single:
        job.fileStore.readGlobalFile(r1_cut_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        parameters.extend(['--single', '-l', '200', '-s', '15', '/data/R1_cutadapt.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_cut_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        job.fileStore.readGlobalFile(r2_cut_id, os.path.join(work_dir, 'R2_cutadapt.fastq'))
        parameters.extend(['/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq'])
    # Call: Kallisto
    docker_call(tool='quay.io/ucsc_cgl/kallisto:0.42.4--35ac87df5b21a8e8e8d159f26864ac1e1db8cf86',
                work_dir=work_dir, parameters=parameters, sudo=config.sudo)
    # Tar output files together and store in fileStore
    output_files = [os.path.join(work_dir, x) for x in ['run_info.json', 'abundance.tsv', 'abundance.h5']]
    tarball_files(tar_name='kallisto.tar.gz', file_paths=output_files, output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'kallisto.tar.gz'))


def star(job, config, r1_cut_id, r2_cut_id):
    """
    Performs alignment of fastqs to bam via STAR

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str r1_cut_id: FileStoreID of fastq (pair 1)
    :param str r2_cut_id: FileStoreID of fastq (pair 2)
    :return: FileStoreID from RSEM
    :rtype: str
    """
    job.fileStore.logToMaster('STAR: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    download_url(url=config.star_index, name='starIndex.tar.gz', work_dir=work_dir)
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'starIndex.tar.gz'), '-C', work_dir])
    os.remove(os.path.join(work_dir, 'starIndex.tar.gz'))
    # Determine tarball structure - star index contains are either in a subdir or in the tarball itself
    star_index = os.path.join('/data', os.listdir(work_dir)[0]) if len(os.listdir(work_dir)) == 1 else '/data'
    # Parameters and input retrieval
    parameters = ['--runThreadN', str(config.cores),
                  '--genomeDir', star_index,
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
    if config.wiggle:
        parameters.extend(['--outWigType', 'bedGraph',
                           '--outWigStrand', 'Unstranded',
                           '--outWigReferencesPrefix', 'chr'])
    if config.single:
        job.fileStore.readGlobalFile(r1_cut_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        parameters.extend(['--readFilesIn', '/data/R1_cutadapt.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_cut_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        job.fileStore.readGlobalFile(r2_cut_id, os.path.join(work_dir, 'R2_cutadapt.fastq'))
        parameters.extend(['--readFilesIn', '/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq'])
    # Call: STAR Map
    docker_call(tool='quay.io/ucsc_cgl/star:2.4.2a--bcbd5122b69ff6ac4ef61958e47bde94001cfe80',
                work_dir=work_dir, parameters=parameters, sudo=config.sudo)
    # Write to fileStore
    bam_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.toTranscriptome.out.bam'))
    # Save Wiggle File
    if config.wiggle and config.s3_dir:
        wiggles = [os.path.basename(x) for x in glob.glob(os.path.join(work_dir, '*.bg'))]
        # Rename extension
        for wiggle in wiggles:
            shutil.move(os.path.join(work_dir, wiggle),
                        os.path.join(work_dir, os.path.splitext(wiggle)[0] + '.bedGraph'))
        wiggles = [os.path.join(work_dir, x) for x in [os.path.splitext(x)[0] + '.bedGraph' for x in wiggles]]
        tarball_files('wiggle.tar.gz', file_paths=wiggles, output_dir=work_dir)
        wiggle_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'wiggle.tar.gz'))
        job.addChildJobFn(s3am_upload_job, file_id=wiggle_id, file_name='wiggle.tar.gz',
                          s3_dir=config.s3_dir, num_cores=config.cores)
    if config.save_bam and config.s3_dir:
        sorted_bam_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.sortedByCoord.out.bam'))
        job.addChildJobFn(s3am_upload_job, file_id=sorted_bam_id, file_name=config.uuid + '.sorted.bam',
                          s3_dir=config.s3_dir, num_cores=config.cores, s3_key_path=config.ssec)
    # RSEM doesn't tend to use more than 16 cores
    cores = min(config.cores, 16)
    disk = '2G' if config.ci_test else '50G'
    return job.addChildJobFn(rsem, config, bam_id, cores=cores, disk=disk).rv()


def rsem(job, config, bam_id):
    """
    RNA quantification with RSEM

    :param JobFunctionWrappingJob job: Passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str bam_id: FileStoreID of transcriptome bam for quantification
    :return: FileStoreID from RSEM postprocess
    :rtype: str
    """
    job.fileStore.logToMaster('RSEM: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    cores = 16 if config.cores >= 16 else config.cores
    download_url(url=config.rsem_ref, name='rsem_ref.tar.gz', work_dir=work_dir)
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'rsem_ref.tar.gz'), '-C', work_dir])
    os.remove(os.path.join(work_dir, 'rsem_ref.tar.gz'))
    # Determine tarball structure - based on it, ascertain folder name and rsem reference prefix
    rsem_files = []
    for root, directories, files in os.walk(work_dir):
        rsem_files.extend([os.path.join(root, x) for x in files])
    # "grp" is a required RSEM extension that should exist in the RSEM reference
    ref_prefix = [os.path.basename(os.path.splitext(x)[0]) for x in rsem_files if 'grp' in x][0]
    ref_folder = os.path.join('/data', os.listdir(work_dir)[0]) if len(os.listdir(work_dir)) == 1 else '/data'
    # I/O
    job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'transcriptome.bam'))
    output_prefix = 'rsem'
    # Call: RSEM
    parameters = ['--quiet',
                  '--no-qualities',
                  '-p', str(cores),
                  '--forward-prob', '0.5',
                  '--seed-length', '25',
                  '--fragment-length-mean', '-1.0',
                  '--bam', '/data/transcriptome.bam',
                  os.path.join(ref_folder, ref_prefix),
                  output_prefix]
    if not config.single:
        parameters = ['--paired-end'] + parameters
    docker_call(tool='quay.io/ucsc_cgl/rsem:1.2.25--d4275175cc8df36967db460b06337a14f40d2f21',
                parameters=parameters, work_dir=work_dir, sudo=config.sudo)
    os.rename(os.path.join(work_dir, output_prefix + '.genes.results'), os.path.join(work_dir, 'rsem_gene.tab'))
    os.rename(os.path.join(work_dir, output_prefix + '.isoforms.results'), os.path.join(work_dir, 'rsem_isoform.tab'))
    # Write to FileStore
    rsem_gene_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_gene.tab'))
    rsem_isoform_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_isoform.tab'))
    job.fileStore.deleteGlobalFile(bam_id)
    # Run child jobs
    return job.addChildJobFn(rsem_postprocess, config, rsem_gene_id, rsem_isoform_id).rv()


def rsem_postprocess(job, config, rsem_gene_id, rsem_isoform_id):
    """
    Parses RSEMs output to produce the separate .tab files (TPM, FPKM, counts) for both gene and isoform.
    These are two-column files: Genes and Quantifications.
    HUGO files are also provided that have been mapped from Gencode/ENSEMBLE names.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param str rsem_gene_id: FileStoreID of rsem_gene_ids
    :param str rsem_isoform_id: FileStoreID of rsem_isoform_ids
    :return: FileStoreID from RSEM post process tarball
    :rytpe: str
    """
    job.fileStore.logToMaster('RSEM Postprocess: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    # I/O
    job.fileStore.readGlobalFile(rsem_gene_id, os.path.join(work_dir, 'rsem_gene.tab'))
    job.fileStore.readGlobalFile(rsem_isoform_id, os.path.join(work_dir, 'rsem_isoform.tab'))
    # Convert RSEM files into individual .tab files.
    docker_call(tool='jvivian/rsem_postprocess', parameters=[config.uuid], work_dir=work_dir, sudo=config.sudo)
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
    docker_call(tool='jvivian/gencode_hugo_mapping', parameters=command, work_dir=work_dir, sudo=config.sudo)
    hugo_files = [os.path.splitext(x)[0] + '.hugo' + os.path.splitext(x)[1] for x in output_files]
    # Create tarballs for outputs
    tarball_files('rsem.tar.gz', file_paths=[os.path.join(work_dir, x) for x in output_files], output_dir=work_dir)
    tarball_files('rsem_hugo.tar.gz', [os.path.join(work_dir, x) for x in hugo_files], output_dir=work_dir)
    rsem_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem.tar.gz'))
    hugo_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_hugo.tar.gz'))
    return rsem_id, hugo_id


def consolidate_output(job, config, output_ids_and_info):
    """
    Combines the contents of the outputs into one tarball and places in output directory or s3

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param list output_ids_and_info: Array of output IDs to consolidate
    """
    job.fileStore.logToMaster('Consolidating input: {}'.format(config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve IDs
    rsem_id, hugo_id, kallisto_id, single_end, improper_pair = flatten(output_ids_and_info)
    # Retrieve output file paths to consolidate
    rsem_tar, hugo_tar, kallisto_tar = None, None, None
    if rsem_id:
        rsem_tar = job.fileStore.readGlobalFile(rsem_id, os.path.join(work_dir, 'rsem.tar.gz'))
        hugo_tar = job.fileStore.readGlobalFile(hugo_id, os.path.join(work_dir, 'rsem_hugo.tar.gz'))
    if kallisto_id:
        kallisto_tar = job.fileStore.readGlobalFile(kallisto_id, os.path.join(work_dir, 'kallisto.tar.gz'))
    # I/O
    if improper_pair:
        config.uuid = 'IMPROPERLY_PAIRED.{}'.format(config.uuid)
    if single_end:
        config.uuid = 'SINGLE-END.{}'.format(config.uuid)
    out_tar = os.path.join(work_dir, config.uuid + '.tar.gz')
    # Consolidate separate tarballs into one as streams (avoids unnecessary untaring)
    tar_list = [x for x in [rsem_tar, hugo_tar, kallisto_tar] if x is not None]
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar in tar_list:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == rsem_tar:
                            tarinfo.name = os.path.join(config.uuid, 'RSEM', os.path.basename(tarinfo.name))
                        elif tar == hugo_tar:
                            tarinfo.name = os.path.join(config.uuid, 'RSEM', 'Hugo', os.path.basename(tarinfo.name))
                        else:
                            tarinfo.name = os.path.join(config.uuid, 'Kallisto', os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
        if improper_pair:
            with open(os.path.join(work_dir, 'WARNING.txt'), 'w') as f:
                f.write('cutadapt: error: Reads are improperly paired. Uneven number of reads in fastq pair.')
            f_out.add(os.path.join(work_dir, 'WARNING.txt'))
    # Move to output directory
    if config.output_dir:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(config.uuid, config.output_dir))
        mkdir_p(config.output_dir)
        move_files(file_paths=[os.path.join(work_dir, config.uuid + '.tar.gz')], output_dir=config.output_dir)
    # Upload to S3
    if config.s3_dir:
        job.fileStore.logToMaster('Uploading {} to S3: {}'.format(config.uuid, config.s3_dir))
        s3am_upload(fpath=out_tar, s3_dir=config.s3_dir, num_cores=config.cores)


# Pipeline specific functions
def parse_samples(path_to_manifest=None, sample_urls=None):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :param list sample_urls: Sample URLs
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    samples = []
    if sample_urls:
        for url in sample_urls:
            samples.append(['tar', 'paired', os.path.basename(url.split('.')[0]), url])
    elif path_to_manifest:
        with open(path_to_manifest, 'r') as f:
            for line in f.readlines():
                if not line.isspace() and not line.startswith('#'):
                    sample = line.strip().split('\t')
                    require(len(sample) == 4, 'Bad manifest format! '
                                              'Expected 4 tab separated columns, got: {}'.format(sample))
                    file_type, paired, uuid, url = sample
                    require(file_type == 'tar' or file_type == 'fq',
                            '1st column must be "tar" or "fq": {}'.format(sample[0]))
                    require(paired == 'paired' or paired == 'single',
                            '2nd column must be "paired" or "single": {}'.format(sample[1]))
                    if file_type == 'fq' and paired == 'paired':
                        require(len(url.split(',')) == 2, 'Fastq pair requires two URLs separated'
                                                          ' by a comma: {}'.format(url))
                    samples.append(sample)
    return samples


def generate_config():
    return textwrap.dedent("""
    {
        # RNA-seq CGL Pipeline configuration file
        # This config file is formatted as a Python dictionary { key: value, key: value, ... }
        # Edit the values in this configuration file and then rerun the pipeline: "toil-rnaseq run"
        # Just Kallisto or STAR/RSEM can be run by supplying only the inputs to those tools
        # Comments (beginning with #) do not need to be removed
        ###############################################################################################
        'star-index': '',       # Required: URL (http, file, s3) to index tarball used by STAR
                                # Example: 's3://cgl-pipeline-inputs/rnaseq_cgl/starIndex_hg38_no_alt.tar.gz'
        'kallisto-index': '',   # Required: URL (http, file, s3) to kallisto index file.
                                # Example: 's3://cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx'
        'rsem-ref': '',         # Required: URL (http, file, s3) to reference tarball used by RSEM
                                # Example: 's3://cgl-pipeline-inputs/rnaseq_cgl/rsem_ref_hg38_no_alt.tar.gz',
        'output-dir': None,     # Optional: Provide a full path to where results will appear
        's3-dir': None,         # Optional: Provide an s3 path (s3://bucket/dir) where results will appear
        'ssec': None,           # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
        'gt-key': None,         # Optional: Provide a full path to a CGHub Key used to access GNOS hosted data
        'wiggle': None,         # Optional: If True, saves a "wiggle" file produced by STAR
        'save-bam': None ,      # Optional: If True, saves the aligned bam (by coordinate) produced by STAR
        'sudo': None,           # Optional: If True, will execute Docker with 'sudo' prepended.
        'ci-test': None,        # Optional: If True, uses resource requirements appropriate for continuous integration
        'fwd-3pr-adapter': 'AGATCGGAAGAG',  # Adapter sequence to trim. Defaults set for Illumina
        'rev-3pr-adapter': 'AGATCGGAAGAG'   # Adapter sequence to trim. Defaults set for Illumina
    }
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 4 tab-separated columns: filetype, paired/unpaired, UUID, URL(s) to sample
        #
        #   filetype    Filetype of the sample. Options: "tar" or "fq", for tarball/tarfile or fastq/fastq.gz
        #   paired      Indicates whether the data is paired or single-ended. Options:  "paired" or "single"
        #   UUID        This should be a unique identifier for the sample to be processed
        #   URL         A URL (http://, ftp://, file://, s3://, gnos://) pointing to the sample
        #
        #   If sample is being submitted as a fastq pair, provide two URLs separated by a comma.
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   tar paired  UUID_1  file:///path/to/sample.tar
        #   fq  paired  UUID_2  file:///path/to/R1.fq.gz,file:///path/to/R2.fq.gz
        #   tar single  UUID_3  http://sample-depot.com/single-end-sample.tar
        #   tar paired  UUID_4  s3://my-bucket-name/directory/paired-sample.tar.gz
        #   fq  single  UUID_5  s3://my-bucket-name/directory/single-end-file.fq
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
    Toil RNA-seq pipeline

    RNA-seq fastqs are combined, aligned, and quantified with 2 different methods (RSEM and Kallisto)

    General usage:
    1. Type "toil-rnaseq generate" to create an editable manifest and config in the current working directory.
    2. Parameterize the pipeline by editing the config.
    3. Fill in the manifest with information pertaining to your samples.
    4. Type "toil-rnaseq run [jobStore]" to execute the pipeline.

    Please read the README.md located in the source directory or at:
    https://github.com/BD2KGenomics/toil-scripts/tree/master/src/toil_scripts/rnaseq_cgl

    Structure of RNA-Seq Pipeline (per sample)

                  3 -- 4 -- 5
                 /          |
      0 -- 1 -- 2 ---- 6 ---7 -- 8

    0 = Download sample
    1 = Unpack/Merge fastqs
    2 = CutAdapt (adapter trimming)
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
    Boto:       pip install boto (OPTIONAL)
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the RNA-seq pipeline')
    group = parser_run.add_mutually_exclusive_group(required=True)
    parser_run.add_argument('--config', default='toil-rnaseq.config', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    group.add_argument('--manifest', default='toil-rnaseq-manifest.tsv', type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')
    group.add_argument('--samples', default=None, nargs='+', type=str,
                       help='Space delimited sample URLs (any number). Samples must be tarfiles/tarballs that contain '
                            'fastq files. URLs follow the format: http://foo.com/sample.tar, '
                            'file:///full/path/to/file.tar. The UUID for the sample will be derived from the file.'
                            'Samples passed in this way will be assumed to be paired end, if using single-end data, '
                            'please use the manifest option.')
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
        if not args.samples:
            require(os.path.exists(args.manifest), '{} not found and no samples provided. Please '
                                                   'run "toil-rnaseq generate-manifest"'.format(args.manifest))
            samples = parse_samples(path_to_manifest=args.manifest)
        else:
            samples = parse_samples(sample_urls=args.samples)
        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in literal_eval(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        # Config sanity checks
        require(config.kallisto_index or config.star_index,
                 'URLs not provided for kallisto or star, so there is nothing to do')
        if config.star_index or config.rsem_ref:
            require(config.star_index and config.rsem_ref, 'Input provided for STAR or RSEM but not both. STAR: '
                                                           '{}, RSEM: {}'.format(config.star_index, config.rsem_ref))
        # Program checks
        for program in ['curl', 'docker']:
            require(which(program), program + ' must be installed on every node.'.format(program))

        # Start the workflow by using map_job() to run the pipeline for each sample
        Job.Runner.startToil(Job.wrapJobFn(map_job, download_sample, samples, config), args)


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
