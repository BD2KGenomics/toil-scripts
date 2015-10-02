#!/usr/bin/env python2.7
# John Vivian
# 8/15
"""
Tree Structure of RNA-Seq Pipeline (per sample)

    0---> 2
    |     |
    1     3
         / \
       *4   5
            |
            6
           / \
         *7   8
             / \
           *9   10
                |
                11
                |
                12
                |
               *13

0 = Start Node
1 = Download Sample
2 = Unzip
3 = Mapsplice
4 = Mapping Stats
5 = Add Read Groups
6 = Bamsort and Index
7 = Rseq-QC
8 = Sort Bam by Reference
9 = Exon Quantification
10 = Transcriptome
11 = Filter
12 = RSEM
13 = RSEM Post-Process

4,7,9,13 contribute to producing the final output

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
from contextlib import closing
import glob
import hashlib
import os
import subprocess
import errno
import multiprocessing
import shutil
import tarfile
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
    parser.add_argument('-o', '--output_dir', default=None, help='full path where final results will be output')
    parser.add_argument('-3', '--s3_dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                             'cgl-driver-projects/ckcc/rna-seq-samples/')
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


def generate_unique_key(master_key_path, url):
    """
    Input1: Path to the BD2K Master Key (for S3 Encryption)
    Input2: S3 URL (e.g. https://s3-us-west-2.amazonaws.com/cgl-driver-projects-encrypted/wcdt/exome_bams/DTB-111-N.bam)

    Returns: 32-byte unique key generated for that URL
    """
    with open(master_key_path, 'r') as f:
        master_key = f.read()
    assert len(master_key) == 32, 'Invalid Key! Must be 32 characters. ' \
                                  'Key: {}, Length: {}'.format(master_key, len(master_key))
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is invalid and is not 32 characters: {}'.format(new_key)
    return new_key


def download_encrypted_file(job, input_args, ids, name):
    """
    Downloads encrypted files from S3 via header injection

    Input1: input arguments defined in main()
    Input2: dictionary of jobStore IDs
    Input3: symbolic name associated with file
    """
    work_dir = job.fileStore.getLocalTempDir()
    key_path = input_args['ssec']
    file_path = os.path.join(work_dir, name)
    url = input_args[name]

    with open(key_path, 'r') as f:
        key = f.read()
    if len(key) != 32:
        raise RuntimeError('Invalid Key! Must be 32 bytes: {}'.format(key))

    key = generate_unique_key(key_path, url)

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

    Input1: input arguments defined in main()
    Input2: dictionary of jobStore IDs
    Input3: symbolic name associated with file
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
    This function should be unpacked for every item being returned, unless none
    of the paths are needed in which case it should be unassigned.

    Input1: working directory files should be placed in
    Input2: dictionary of jobStore IDs
    Input3: symoblic filenames that are being retrieved
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

    Input1: Docker tool to be pulled and run (repo/tool_name)
    Input2: parameters to the Docker tool being called
    Input3: working directory where input files are located
    """
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
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


def tarball_files(work_dir, zip_name, uuid=None, files=None):
    with tarfile.open(os.path.join(work_dir, zip_name), 'w:gz') as f_out:
        for fname in files:
            if uuid:
                f_out.add(os.path.join(work_dir, fname), arcname=uuid + '.' + fname)
            else:
                f_out.add(os.path.join(work_dir, fname), arcname=fname)


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
                      'stats.txt', 'stats2.txt', 'stats_all.txt', 'mapping.tab', 'exon_quant', 'exon_quant.bed',
                      'map.tar.gz', 'exon.tar.gz', 'qc.tar.gz', 'rsem.tar.gz', 'uuid.tar.gz']
    ids = shared_ids.copy()
    ids.update( {x: job.fileStore.getEmptyFileStoreID() for x in symbolic_names} )
    sample_input = dict(input_args)
    # Update input
    sample_input['uuid'] = uuid
    sample_input['sample.zip'] = sample_url
    sample_input['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    sample_input['cpu_count'] = multiprocessing.cpu_count()
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
    cores = input_args['cpu_count']

    # I/O
    sample = return_input_paths(job, work_dir, ids, 'sample.zip')
    # Unzip File
    subprocess.check_call(['unzip', sample, '-d', work_dir])
    # Remove large files before creating concat versions.
    os.remove(os.path.join(work_dir, 'sample.zip'))
    # Zcat files in parallel
    R1_files = sorted(glob.glob(os.path.join(work_dir, '*R1*')))
    R2_files = sorted(glob.glob(os.path.join(work_dir, '*R2*')))
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
    job.addChildJobFn(mapsplice, job_vars, cores=cores, memory='30 G', disk='100 G')


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


def mapping_stats(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # I/O
    return_input_paths(job, work_dir, ids, 'stats.txt')
    uuid = input_args['uuid']
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)
    # Command
    docker_call(tool='jvivian/mapping_stats', tool_parameters=[uuid], work_dir=work_dir)
    # Update FileStore
    job.fileStore.updateGlobalFile(ids['stats2.txt'], os.path.join(work_dir, '{}_stats2.txt'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['stats_all.txt'], os.path.join(work_dir, '{}_stats_all.txt'.format(uuid)))
    job.fileStore.updateGlobalFile(ids['mapping.tab'], os.path.join(work_dir, '{}_mapping.tab'.format(uuid)))
    # Zip output files and store
    output_files = ['{}_stats2.txt'.format(uuid), '{}_stats_all.txt'.format(uuid), '{}_mapping.tab'.format(uuid)]
    tarball_files(work_dir, zip_name='map.tar.gz', files=output_files)
    job.fileStore.updateGlobalFile(ids['map.tar.gz'], os.path.join(work_dir, 'map.tar.gz'))


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
    tarball_files(work_dir, zip_name='qc.tar.gz', uuid=None, files=output_files)
    job.fileStore.updateGlobalFile(ids['qc.tar.gz'], os.path.join(work_dir, 'qc.tar.gz'))


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
    job.addChildJobFn(transcriptome, job_vars, disk='30 G', memory='30 G')
    job.addChildJobFn(exon_count, job_vars, disk='30 G')


def exon_count(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    uuid = input_args['uuid']
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
    # Create zip, upload to fileStore, and move to output_dir as a backup
    output_files = ['exon_quant.bed', 'exon_quant']
    tarball_files(work_dir, zip_name='exon.tar.gz', uuid=uuid, files=output_files)
    job.fileStore.updateGlobalFile(ids['exon.tar.gz'], os.path.join(work_dir, 'exon.tar.gz'))


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
    cores = input_args['cpu_count']

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
    job.addChildJobFn(rsem, job_vars, cores=cores, disk='30 G')


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
    uuid = input_args['uuid']
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)

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

    output_files = ['rsem.genes.norm_counts.tab', 'rsem.genes.raw_counts.tab', 'rsem.genes.norm_fpkm.tab',
                    'rsem.genes.norm_tpm.tab', 'rsem.isoform.norm_counts.tab', 'rsem.isoform.raw_counts.tab',
                    'rsem.isoform.norm_fpkm.tab', 'rsem.isoform.norm_tpm.tab']
    tarball_files(work_dir, zip_name='rsem.tar.gz', uuid=uuid, files=output_files)
    job.fileStore.updateGlobalFile(ids['rsem.tar.gz'], os.path.join(work_dir, 'rsem.tar.gz'))
    # Run children
    job.addChildJobFn(consolidate_output, job_vars)



def consolidate_output(job, job_vars):
    """
    Combine the contents of separate zipped outputs into one via streaming
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)

    # I/O
    rsem_tar, map_tar, exon_tar, qc_tar =  return_input_paths(job, work_dir, ids, 'rsem.tar.gz', 'map.tar.gz',
                                                              'exon.tar.gz', 'qc.tar.gz')
    out_tar = os.path.join(work_dir, uuid + '.tar.gz')
    # Consolidate separate tarballs
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar in [rsem_tar, map_tar, exon_tar, qc_tar]:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == qc_tar:
                            tarinfo.name = os.path.join(uuid, 'rseq_qc', os.path.basename(tarinfo.name))
                        else:
                            tarinfo.name = os.path.join(uuid, os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
    # Update file store
    job.fileStore.updateGlobalFile(ids['uuid.tar.gz'], out_tar)
    if input_args['output_dir']:
        move_to_output_dir(work_dir, output_dir, uuid=None, files=[uuid + '.tar.gz'])
    # If S3 bucket argument specified, upload to S3
    if input_args['s3_dir']:
        job.addChildJobFn(upload_to_s3, job_vars)


def upload_to_s3(job, job_vars):
    """
    If s3_dir is specified in arguments, file will be uploaded to S3 using boto.
    WARNING: ~/.boto credentials are necessary for this to succeed!
    """
    import boto
    from boto.s3.key import Key

    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    # Parse s3_dir
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.split('/')[0]
    bucket_dir = '/'.join(s3_dir.split('/')[1:])
    # I/O
    uuid_tar = return_input_paths(job, work_dir, ids, 'uuid.tar.gz')
    # Upload to S3 via boto
    conn = boto.connect_s3()
    bucket = conn.get_bucket(bucket_name)
    k = Key(bucket)
    k.key = os.path.join(bucket_dir, uuid + '.tar.gz')
    k.set_contents_from_filename(uuid_tar)


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
              'ssec': args.ssec,
              's3_dir': args.s3_dir,
              'uuid': None,
              'samples.zip': None,
              'cpu_count': None}

    # Launch jobs
    Job.Runner.startToil(Job.wrapJobFn(batch_start, inputs), args)