import os

import subprocess

from toil_scripts.lib.files import tarball_files
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url


def run_kallisto(job, cores, r1_id, r2_id, kallisto_index_url):
    """
    RNA quantification via Kallisto

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Number of cores to run Kallisto with
    :param str r1_id: FileStoreID of fastq (pair 1)
    :param str r2_id: FileStoreID of fastq (pair 2 if applicable, otherwise pass None for single-end)
    :param str kallisto_index_url: FileStoreID for Kallisto index file
    :return: FileStoreID from Kallisto output
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    download_url(url=kallisto_index_url, name='kallisto_hg38.idx', work_dir=work_dir)
    # Retrieve files
    parameters = ['quant',
                  '-i', '/data/kallisto_hg38.idx',
                  '-t', str(cores),
                  '-o', '/data/',
                  '-b', '100']
    if r1_id and r2_id:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        job.fileStore.readGlobalFile(r2_id, os.path.join(work_dir, 'R2_cutadapt.fastq'))
        parameters.extend(['/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        parameters.extend(['--single', '-l', '200', '-s', '15', '/data/R1_cutadapt.fastq'])

    # Call: Kallisto
    docker_call(tool='quay.io/ucsc_cgl/kallisto:0.42.4--35ac87df5b21a8e8e8d159f26864ac1e1db8cf86',
                work_dir=work_dir, parameters=parameters)
    # Tar output files together and store in fileStore
    output_files = [os.path.join(work_dir, x) for x in ['run_info.json', 'abundance.tsv', 'abundance.h5']]
    tarball_files(tar_name='kallisto.tar.gz', file_paths=output_files, output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'kallisto.tar.gz'))


def run_rsem(job, cores, bam_id, rsem_ref_url, paired=True):
    """
    RNA quantification with RSEM

    :param JobFunctionWrappingJob job: Passed automatically by Toil
    :param int cores: Number of cores to run RSEM with
    :param str bam_id: FileStoreID of transcriptome bam for quantification
    :param str rsem_ref_url: URL of RSEM reference (tarball)
    :param bool paired: If True, uses parameters for paired end data
    :return: FileStoreIDs for RSEM's gene and isoform output
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    download_url(url=rsem_ref_url, name='rsem_ref.tar.gz', work_dir=work_dir)
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
    if paired:
        parameters = ['--paired-end'] + parameters
    docker_call(tool='quay.io/ucsc_cgl/rsem:1.2.25--d4275175cc8df36967db460b06337a14f40d2f21',
                parameters=parameters, work_dir=work_dir)
    os.rename(os.path.join(work_dir, output_prefix + '.genes.results'), os.path.join(work_dir, 'rsem_gene.tab'))
    os.rename(os.path.join(work_dir, output_prefix + '.isoforms.results'), os.path.join(work_dir, 'rsem_isoform.tab'))
    # Write to FileStore
    gene_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_gene.tab'))
    isoform_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_isoform.tab'))
    return gene_id, isoform_id


def run_rsem_postprocess(job, uuid, rsem_gene_id, rsem_isoform_id):
    """
    Parses RSEMs output to produce the separate .tab files (TPM, FPKM, counts) for both gene and isoform.
    These are two-column files: Genes and Quantifications.
    HUGO files are also provided that have been mapped from Gencode/ENSEMBLE names.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: UUID to mark the samples with
    :param str rsem_gene_id: FileStoreID of rsem_gene_ids
    :param str rsem_isoform_id: FileStoreID of rsem_isoform_ids
    :return: FileStoreID from RSEM post process tarball
    :rytpe: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    # I/O
    job.fileStore.readGlobalFile(rsem_gene_id, os.path.join(work_dir, 'rsem_gene.tab'), mutable=True)
    job.fileStore.readGlobalFile(rsem_isoform_id, os.path.join(work_dir, 'rsem_isoform.tab'), mutable=True)
    # Convert RSEM files into individual .tab files.
    docker_call(tool='jvivian/rsem_postprocess', parameters=[uuid], work_dir=work_dir)
    os.rename(os.path.join(work_dir, 'rsem_gene.tab'), os.path.join(work_dir, 'rsem_genes.results'))
    os.rename(os.path.join(work_dir, 'rsem_isoform.tab'), os.path.join(work_dir, 'rsem_isoforms.results'))
    output_files = ['rsem.genes.norm_counts.tab', 'rsem.genes.raw_counts.tab', 'rsem.isoform.norm_counts.tab',
                    'rsem.isoform.raw_counts.tab', 'rsem_genes.results', 'rsem_isoforms.results']
    # Perform HUGO gene / isoform name mapping
    genes = [x for x in output_files if 'rsem.genes' in x]
    isoforms = [x for x in output_files if 'rsem.isoform' in x]
    command = ['-g'] + genes + ['-i'] + isoforms
    docker_call(tool='jvivian/gencode_hugo_mapping', parameters=command, work_dir=work_dir)
    hugo_files = [os.path.splitext(x)[0] + '.hugo' + os.path.splitext(x)[1] for x in genes + isoforms]
    # Create tarballs for outputs
    tarball_files('rsem.tar.gz', file_paths=[os.path.join(work_dir, x) for x in output_files], output_dir=work_dir)
    tarball_files('rsem_hugo.tar.gz', [os.path.join(work_dir, x) for x in hugo_files], output_dir=work_dir)
    rsem_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem.tar.gz'))
    hugo_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_hugo.tar.gz'))
    return rsem_id, hugo_id
