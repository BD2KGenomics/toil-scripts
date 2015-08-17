#!/usr/bin/env python2.7
# John Vivian
# 8/15
"""
RNA-Seq Pipeline

Dependencies:
Curl:       apt-get install curl
Docker:     apt-get install docker.io # linux, o.w. just docker
Samtools:   apt-get install samtools
Unzip:      apt-get install unzip
Toil:       pip install git+https://github.com/BD2KGenomics/toil.git
"""

import argparse
import base64
from collections import OrderedDict
import glob
import hashlib
import os
import subprocess
import errno
import multiprocessing
import shutil
from toil.job import Job


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', required=True, help='Path to config')
    parser.add_argument('-u', '--unc', required=True, help='URL to unc_hg19.bed')
    parser.add_argument('-f', '--fasta', required=True, help='URL to hg19_M_rCRS_ref.transcripts.fa')
    parser.add_argument('-x', '--composite_exons', required=True, help='URL to composite_exons.bed')
    parser.add_argument('-n', '--normalize', required=True, help='URL to normalizeBedToolsExonQuant.pl')
    parser.add_argument('-r', '--rsem_ref', required=True, help='RSEM_REF URL')
    parser.add_argument('-chr', '--chromosomes', required=True, help='Chromosomes Directory')
    parser.add_argument('-e', '--ebwt', required=True, help='EBWT Directory')
    parser.add_argument('-s', '--ssec', help='Path to Key File for SSE-C Encryption')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to output final results')
    return parser


# Convenience functions used in the pipeline
def mkdir_p(path):
    """
    It is Easier to Ask for Forgiveness than Permission
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def download_encrypted_file(job, input_args, ids, name):
    """
    Downloads encrypted files
    """
    work_dir = job.fileStore.getLocalTempDir()
    key_path = input_args['ssec']
    file_path = os.path.join(work_dir, name)
    url = input_args[name]

    with open(key_path, 'r') as f:
        key = f.read()
    if len(key) != 32:
        raise RuntimeError('Invalid Key! Must be 32 bytes: {}'.format(key))
    encoded_key = base64.b64encode(key)
    encoded_key_md5 = base64.b64encode( hashlib.md5(key).digest() )
    h1 = 'x-amz-server-side-encryption-customer-algorithm:AES256'
    h2 = 'x-amz-server-side-encryption-customer-key:{}'.format(encoded_key)
    h3 = 'x-amz-server-side-encryption-customer-key-md5:{}'.format(encoded_key_md5)
    try:
        subprocess.check_call(['curl', '-fs', '--retry', '5', '-H', h1, '-H', h2, '-H', h3, url, '-o', file_path])
    except OSError:
        raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    job.fileStore.updateGlobalFile(ids[name], file_path)


def download_from_url(job, input_args, ids, name):
    """
    Downloads a URL that was supplied as an argument to running this script in LocalTempDir.
    After downloading the file, it is stored in the FileStore.
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, name)
    url = input_args[name]
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    job.fileStore.updateGlobalFile(ids[name], file_path)


def return_input_paths(job, work_dir, ids, *args):
    """
    Returns the paths of files from the FileStore if they are not present.
    """
    paths = OrderedDict()
    for name in args:
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = job.fileStore.readGlobalFile(ids[name], os.path.join(work_dir, name))
        else:
            file_path = os.path.join(work_dir, name)
        paths[name] = file_path
        if len(args) == 1:
            return file_path

    return paths.values()


def docker_path(filepath):
    """
    Given a path, returns that files path inside the docker mount directory (/data).
    """
    return os.path.join('/data', os.path.basename(filepath))


def docker_call(tool, tool_parameters, work_dir):
    """
    Makes subprocess call of a command to a docker container.
    work_dir MUST BE AN ABSOLUTE PATH or the call will fail.
    """
    # base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    base_docker_call = 'docker run -v {}:/data'.format(work_dir)
    call = base_docker_call.split() + [tool] + tool_parameters
    try:
        subprocess.check_call(call)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status for cmd {}'.format(call))
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


def move_to_output_dir(work_dir, output_dir, uuid=None, files=list()):
    """
    A list of files to move from work_dir to output_dir.
    """
    for fname in files:
        if uuid is None:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


# Start of job Functions
def batch_start(job, input_args):
    """
    Creates FileStoreIDs for files that are shared as input by every sample.
    Downloads and stores these files in the FileStore.
    """
    shared_files = ['unc.bed', 'hg19.transcripts.fa','composite_exons.bed', 'normalize.pl', 'rsem_ref.zip',
                    'ebwt.zip', 'chromosomes.zip']
    shared_ids = {x: job.fileStore.getEmptyFileStoreID() for x in shared_files}
    # Download shared files used by all samples in the pipeline
    for file in shared_files:
        job.addChildJobFn(download_from_url, input_args, shared_ids, file)
    job.addFollowOnJobFn(spawn_batch_jobs, shared_ids, input_args)


def spawn_batch_jobs(job, shared_ids, input_args):
    """
    Launches pipeline for each sample.
    """
    samples = []
    config = input_args['config']
    with open(config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                samples.append(line.strip().split(','))
    for sample in samples:
        job.addChildJobFn(start_node, shared_ids, input_args, sample)


def start_node(job, shared_ids, input_args, sample):
    """
    Defines variables unique to a sample that are used in the rest of the pipeline.
    """
    uuid, sample_url = sample
    # Symbolic names for sample specific files
    symbolic_names = ['sample.zip', 'R1.fastq.gz', 'R2.fastq.gz', 'R1.fastq', 'R2.fastq', 'alignments.bam',
                      'rg_alignments.bam', 'sorted.bam', 'sorted.bam.bai', 'sort_by_ref.bam', 'transcriptome.bam',
                      'filtered.bam', 'rsem_gene.tab', 'rsem_isoform.tab', 'rsem.genes.norm_counts.tab', 'rsem.genes.raw_counts.tab',
                      'rsem.genes.norm_fpkm.tab', 'rsem.genes.norm_tpm.tab', 'rsem.isoform.raw_counts.tab',
                      'rsem.isoform.norm_counts.tab', 'rsem.isoform.norm_fpkm.tab', 'rsem.isoform.norm_tpm.tab',
                      'gene_raw_count.tab', 'gene_norm_count.tab', 'gene_fpkm.tab', 'gene_tpm.tab',
                      'isoform_raw_count.tab', 'isoform_norm_count.tab', 'isoform_fpkm.tab', 'isoform_tpm.tab',
                      'stats.txt', 'stats2.txt', 'stats_all.txt', 'mapping.tab', 'exon_quant', 'exon_quant.bed']
    ids = shared_ids.copy()
    ids.update( {x: job.fileStore.getEmptyFileStoreID() for x in symbolic_names} )
    sample_input = dict(input_args)
    # Update input
    sample_input['uuid'] = uuid
    sample_input['sample.zip'] = sample_url
    sample_input['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    job_vars = (sample_input, ids)

    # Download sample fastqs and launch pipeline
    job.addChildJobFn(download_encrypted_file, sample_input, ids, 'sample.zip', disk='12 G')
    job.addFollowOnJobFn(unzip, job_vars, disk='50 G')


def unzip(job, job_vars):
    """
    Unzips input sample and concats the Read1 and Read2 groups together.
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    sample = return_input_paths(job, work_dir, ids, 'sample.zip')
    # Unzip File
    subprocess.check_call(['unzip', sample, '-d', work_dir])
    # Remove large files before creating concat versions.
    os.remove(os.path.join(work_dir, 'sample.zip'))
    # Zcat files in parallel
    # TODO: FIX THIS! Files should follow R1 and R2 format. NOT left/right -- change interleaved script.
    # R1_files = glob.glob(os.path.join(work_dir, '*R1*'))
    # R2_files = glob.glob(os.path.join(work_dir, '*R2*'))
    R1_files = glob.glob(os.path.join(work_dir, '*left*'))
    R2_files = glob.glob(os.path.join(work_dir, '*right*'))
    with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
        p1 = subprocess.Popen(['zcat'] + R1_files, stdout=f1)
    with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
        p2 = subprocess.Popen(['zcat'] + R2_files, stdout=f2)
    p1.wait()
    p2.wait()
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['R1.fastq'], os.path.join(work_dir, 'R1.fastq'))
    job.fileStore.updateGlobalFile(ids['R2.fastq'], os.path.join(work_dir, 'R2.fastq'))
    job.fileStore.deleteGlobalFile(ids['sample.zip'])
    # Run children and follow-on
    job.addFollowOnJobFn(mapsplice, job_vars, cpu=36, memory='30 G', disk='100 G')


def mapsplice(job, job_vars):
    """
    Maps RNA-Seq reads to a reference genome.
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cpus = input_args['cpu_count']

    # I/O
    r1, r2, ebwt, chromosomes = return_input_paths(job, work_dir, ids, 'R1.fastq', 'R2.fastq', 'ebwt.zip', 'chromosomes.zip')
    for fname in ['chromosomes.zip', 'ebwt.zip']:
        subprocess.check_call(['unzip', '-o', os.path.join(work_dir, fname), '-d', work_dir])
    # Command and call
    parameters = ['-p', str(cpus),
                  '-s', '25',
                  '--bam',
                  '--min-map-len', '50',
                  '-x', docker_path('ebwt'),
                  '-c', docker_path('chromosomes'),
                  '-1', docker_path(r1),
                  '-2', docker_path(r2),
                  '-o', '/data']
    docker_call(tool='computationalgenomicslab/mapsplice', tool_parameters=parameters, work_dir=work_dir)
    # Update FileStore
    for fname in ['alignments.bam', 'stats.txt']:
        job.fileStore.updateGlobalFile(ids[fname], os.path.join(work_dir, fname))
    for fname in ['R1.fastq', 'R2.fastq']:
        job.fileStore.deleteGlobalFile(ids[fname])
    # Run child job
    job.addChildJobFn(add_read_groups, job_vars, disk='30 G')
    job.addChildJobFn(mapping_stats, job_vars)


def add_read_groups(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    alignments = return_input_paths(job, work_dir, ids, 'alignments.bam')
    output = os.path.join(work_dir, 'rg_alignments.bam')
    # Command and call
    parameter = ['AddOrReplaceReadGroups',
                 'INPUT={}'.format(docker_path(alignments)),
                 'OUTPUT={}'.format(docker_path(output)),
                 'RGSM={}'.format(input_args['uuid']),
                 'RGID={}'.format(input_args['uuid']),
                 'RGLB=TruSeq',
                 'RGPL=illumina',
                 'RGPU=barcode',
                 'VALIDATION_STRINGENCY=SILENT']
    docker_call(tool='computationalgenomicslab/picardtools', tool_parameters=parameter, work_dir=work_dir)
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['rg_alignments.bam'], output)
    # Run child job
    job.addChildJobFn(bamsort_and_index, job_vars, disk='30 G')


def bamsort_and_index(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    rg_alignments = return_input_paths(job, work_dir, ids, 'rg_alignments.bam')
    output = os.path.join(work_dir, 'sorted.bam')
    # Command -- second argument is "Output Prefix"
    cmd1 = ['sort', docker_path(rg_alignments), docker_path('sorted')]
    cmd2 = ['index', docker_path(output)]
    docker_call(tool='computationalgenomicslab/samtools', tool_parameters=cmd1, work_dir=work_dir)
    docker_call(tool='computationalgenomicslab/samtools', tool_parameters=cmd2, work_dir=work_dir)
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['sorted.bam'], output)
    job.fileStore.updateGlobalFile(ids['sorted.bam.bai'], os.path.join(work_dir, 'sorted.bam.bai'))
    # Run child job
    job.addChildJobFn(sort_bam_by_reference, job_vars, disk='50 G')
    job.addChildJobFn(rseq_qc, job_vars, disk='20 G')


def rseq_qc(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = os.path.join(input_args['output_dir'], 'rseq_qc')
    mkdir_p(output_dir)
    uuid = input_args['uuid']

    #I/O
    sorted_bam, sorted_bai = return_input_paths(job, work_dir, ids, 'sorted.bam', 'sorted.bam.bai')
    # Command
    docker_call(tool='jvivian/qc', tool_parameters=['/opt/cgl-docker-lib/RseqQC_v2.sh', docker_path(sorted_bam), uuid], work_dir=work_dir)
    # Update FileStore
    output_files = [f for f in glob.glob(os.path.join(work_dir, '*')) if 'sorted.bam' not in f]
    move_to_output_dir(work_dir, output_dir, uuid=None, files=[os.path.basename(x) for x in output_files])


def sort_bam_by_reference(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    sorted_bam, sorted_bai = return_input_paths(job, work_dir, ids, 'sorted.bam', 'sorted.bam.bai')
    output = os.path.join(work_dir, 'sort_by_ref.bam')
    # Command
    ref_seqs = []
    handle = subprocess.Popen(["samtools", "view", "-H", sorted_bam], stdout=subprocess.PIPE).stdout
    for line in handle:
        if line.startswith("@SQ"):
            tmp = line.split("\t")
            chrom = tmp[1].split(":")[1]
            ref_seqs.append(chrom)
    handle.close()
    # Iterate through chromosomes to create mini-bams
    for chrom in ref_seqs:
        # job.addChildJobFn(sbbr_child, chrom, os.path.join(work_dir, chrom), sorted_bam)
        cmd_view = ["samtools", "view", "-b", sorted_bam, chrom]
        cmd_sort = ["samtools", "sort", "-m", "3000000000", "-n", "-", os.path.join(work_dir, chrom)]
        p1 = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
        subprocess.check_call(cmd_sort, stdin=p1.stdout)
    sorted_files = [os.path.join(work_dir, chrom) + '.bam' for chrom in ref_seqs]
    cmd = ["samtools", "cat", "-o", output] + sorted_files
    subprocess.check_call(cmd)
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['sort_by_ref.bam'], output)
    job.addChildJobFn(transcriptome, job_vars, disk='30 G')
    job.addChildJobFn(exon_count, job_vars, disk='30 G')


def exon_count(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)

    # I/O
    sort_by_ref, normalize_pl, composite_bed = return_input_paths(job, work_dir, ids, 'sort_by_ref.bam',
                                                                  'normalize.pl', 'composite_exons.bed')
    # Command
    tool = 'jvivian/bedtools'
    cmd_1 = ['coverage',
               '-split',
              '-abam', docker_path(sort_by_ref),
              '-b', docker_path(composite_bed)]
    cmd_2 = ['perl',
              os.path.join(work_dir, 'normalize.pl'),
              sort_by_ref,
              composite_bed]

    p = subprocess.Popen(['sudo', 'docker', 'run', '-v', '{}:/data'.format(work_dir), tool] + cmd_1, stdout=subprocess.PIPE)
    with open(os.path.join(work_dir, 'exon_quant'), 'w') as f:
        subprocess.check_call(cmd_2, stdin=p.stdout, stdout=f)

    p1 = subprocess.Popen(['cat', os.path.join(work_dir, 'exon_quant')], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['tr', '":"', '"\t"'], stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(['tr', '"-"', '"\t"'], stdin=p2.stdout, stdout=subprocess.PIPE)
    with open(os.path.join(work_dir, 'exon_quant.bed'), 'w') as f:
        subprocess.check_call(['cut', '-f1-4'], stdin=p3.stdout, stdout=f)

    job.fileStore.updateGlobalFile(ids['exon_quant'], os.path.join(work_dir, 'exon_quant'))
    job.fileStore.updateGlobalFile(ids['exon_quant.bed'], os.path.join(work_dir, 'exon_quant.bed'))

    move_to_output_dir(work_dir, output_dir, uuid=input_args['uuid'],
                       files=['exon_quant.bed', 'exon_quant'])


def transcriptome(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    sort_by_ref, bed, hg19_fa = return_input_paths(job, work_dir, ids, 'sort_by_ref.bam',
                                                   'unc.bed', 'hg19.transcripts.fa')
    output = os.path.join(work_dir, 'transcriptome.bam')
    # Command
    parameters = ['sam-xlate',
                  '--bed', docker_path(bed),
                  '--in', docker_path(sort_by_ref),
                  '--order', docker_path(hg19_fa),
                  '--out', docker_path(output),
                  '--xgtag',
                  '--reverse']
    docker_call(tool='computationalgenomicslab/ubu', tool_parameters=parameters, work_dir=os.path.dirname(output))
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['transcriptome.bam'], output)
    # Run child job
    job.addChildJobFn(filter_bam, job_vars, disk='30 G')


def filter_bam(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    transcriptome_bam = return_input_paths(job, work_dir, ids, 'transcriptome.bam')
    output = os.path.join(work_dir, 'filtered.bam')
    # Command
    parameters = ['sam-filter',
                  '--strip-indels',
                  '--max-insert', '1000',
                  '--mapq', '1',
                  '--in', docker_path(transcriptome_bam),
                  '--out', docker_path(output)]
    docker_call(tool='computationalgenomicslab/ubu', tool_parameters=parameters, work_dir=os.path.dirname(output))
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['filtered.bam'], output)
    # Run child job
    job.addChildJobFn(rsem, job_vars, cpu=32, disk='30 G')


def rsem(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cpus = input_args['cpu_count']

    # I/O
    filtered_bam, rsem_ref = return_input_paths(job, work_dir, ids, 'filtered.bam', 'rsem_ref.zip')
    subprocess.check_call(['unzip', '-o', os.path.join(work_dir, 'rsem_ref.zip'), '-d', work_dir])
    output_prefix = 'rsem'
    # Make tool call to Docker
    parameters = ['--quiet',
                  '--no-qualities',
                  '--paired-end',
                  '-p', str(cpus),
                  '--forward-prob', '0.5',
                  '--seed-length', '25',
                  '--fragment-length-mean', '-1.0',
                  '--bam', docker_path(filtered_bam),
                  '/data/rsem_ref/hg19_M_rCRS_ref',
                  output_prefix]

    # TODO: CHANGE TOOL ONCE RSEM IS PUSHED TO CGL
    docker_call(tool='jvivian/rsem', tool_parameters=parameters, work_dir=work_dir)
    os.rename(os.path.join(work_dir, output_prefix+'.genes.results'), os.path.join(work_dir, 'rsem_gene.tab'))
    os.rename(os.path.join(work_dir, output_prefix+'.isoforms.results'), os.path.join(work_dir, 'rsem_isoform.tab'))
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['rsem_gene.tab'], os.path.join(work_dir, 'rsem_gene.tab'))
    job.fileStore.updateGlobalFile(ids['rsem_isoform.tab'], os.path.join(work_dir, 'rsem_isoform.tab'))
    # Run child jobs
    job.addChildJobFn(rsem_postprocess, job_vars)


def rsem_postprocess(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    return_input_paths(job, work_dir, ids, 'rsem_gene.tab', 'rsem_isoform.tab')
    # Command
    sample = input_args['uuid']
    docker_call(tool='jvivian/rsem_postprocess', tool_parameters=[sample], work_dir=work_dir)
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['rsem.genes.norm_counts.tab'], os.path.join(work_dir, 'rsem.genes.norm_counts.tab'))
    job.fileStore.updateGlobalFile(ids['rsem.genes.raw_counts.tab'], os.path.join(work_dir, 'rsem.genes.raw_counts.tab'))
    job.fileStore.updateGlobalFile(ids['rsem.genes.norm_fpkm.tab'], os.path.join(work_dir, 'rsem.genes.norm_fpkm.tab'))
    job.fileStore.updateGlobalFile(ids['rsem.genes.norm_tpm.tab'], os.path.join(work_dir, 'rsem.genes.norm_tpm.tab'))
    job.fileStore.updateGlobalFile(ids['rsem.isoform.norm_counts.tab'], os.path.join(work_dir, 'rsem.isoform.norm_counts.tab'))
    job.fileStore.updateGlobalFile(ids['rsem.isoform.raw_counts.tab'], os.path.join(work_dir, 'rsem.isoform.raw_counts.tab'))
    job.fileStore.updateGlobalFile(ids['rsem.isoform.norm_fpkm.tab'], os.path.join(work_dir, 'rsem.isoform.norm_fpkm.tab'))
    job.fileStore.updateGlobalFile(ids['rsem.isoform.norm_tpm.tab'], os.path.join(work_dir, 'rsem.isoform.norm_tpm.tab'))
    # Run child job
    job.addChildJobFn(rsem_reduce, job_vars)


def rsem_reduce(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    uuid = input_args['uuid']
    output_dir = input_args['output_dir']
    return_input_paths(job, work_dir, ids, 'rsem.genes.norm_counts.tab', 'rsem.genes.raw_counts.tab',
                       'rsem.genes.norm_fpkm.tab', 'rsem.genes.norm_tpm.tab', 'rsem.isoform.norm_counts.tab',
                       'rsem.isoform.raw_counts.tab', 'rsem.isoform.norm_fpkm.tab', 'rsem.isoform.norm_tpm.tab')

    # Command
    docker_call(tool='jvivian/reduce', tool_parameters=[uuid], work_dir=work_dir)
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['gene_raw_count.tab'], os.path.join(work_dir, '{}.gene_raw_count.tab'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['gene_norm_count.tab'], os.path.join(work_dir, '{}.gene_norm_count.tab'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['gene_fpkm.tab'], os.path.join(work_dir, '{}.gene_fpkm.tab'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['gene_tpm.tab'], os.path.join(work_dir, '{}.gene_tpm.tab'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['isoform_raw_count.tab'], os.path.join(work_dir, '{}.isoform_raw_count.tab'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['isoform_norm_count.tab'], os.path.join(work_dir, '{}.isoform_norm_count.tab'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['isoform_fpkm.tab'], os.path.join(work_dir, '{}.isoform_fpkm.tab'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['isoform_tpm.tab'], os.path.join(work_dir, '{}.isoform_tpm.tab'.format(uuid)))
    # Move files to output
    move_to_output_dir(work_dir, output_dir, uuid=None, files=['{}.gene_raw_count.tab'.format(uuid),
                                                               '{}.gene_norm_count.tab'.format(uuid),
                                                               '{}.gene_fpkm.tab'.format(uuid),
                                                               '{}.gene_tpm.tab'.format(uuid),
                                                               '{}.isoform_raw_count.tab'.format(uuid),
                                                               '{}.isoform_norm_count.tab'.format(uuid),
                                                               '{}.isoform_fpkm.tab'.format(uuid),
                                                               '{}.isoform_tpm.tab'.format(uuid)])


def mapping_stats(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    return_input_paths(job, work_dir, ids, 'stats.txt')
    uuid = input_args['uuid']
    output_dir = input_args['output_dir']
    # Command
    docker_call(tool='jvivian/mapping_stats', tool_parameters=[uuid], work_dir=work_dir)
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['stats2.txt'], os.path.join(work_dir, '{}_stats2.txt'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['stats_all.txt'], os.path.join(work_dir, '{}_stats_all.txt'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['mapping.tab'], os.path.join(work_dir, '{}_mapping.tab'.format(uuid)))
    # Move files to output_dir
    move_to_output_dir(work_dir, output_dir, uuid=None, files=['{}_stats2.txt'.format(uuid),
                                                               '{}_stats_all.txt'.format(uuid),
                                                               '{}_mapping.tab'.format(uuid)])


if __name__ == "__main__":
    # Define Parser object and add to toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    # Store input_URLs for downloading
    inputs = {'config': args.config,
              'unc.bed': args.unc,
              'hg19.transcripts.fa': args.fasta,
              'composite_exons.bed': args.composite_exons,
              'normalize.pl': args.normalize,
              'output_dir': args.output_dir,
              'rsem_ref.zip': args.rsem_ref,
              'chromosomes.zip': args.chromosomes,
              'ebwt.zip': args.ebwt,
              'uuid': None,
              'samples.zip': None,
              'ssec': args.ssec,
              'cpu_count': multiprocessing.cpu_count()}

    # Launch jobs
    Job.Runner.startToil(Job.wrapJobFn(batch_start, inputs), args)
    Job.Runner.cleanup(args)
