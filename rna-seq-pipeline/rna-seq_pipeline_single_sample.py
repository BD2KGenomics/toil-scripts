#!/usr/bin/env python2.7
# John Vivian
# 8/15
"""
RNA-Seq Pipeline

Dependencies:
Boto:       pip install boto
Curl:       apt-get install curl
Docker:     apt-get install docker.io # linux, o.w. just docker
Samtools:   apt-get install samtools
"""

import argparse
from collections import OrderedDict
import os
import subprocess
import errno
import multiprocessing
from toil.job import Job

import boto
from boto.exception import S3ResponseError


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r1', '--sample_1', required=True, help='URL to Sample: Read 1')
    parser.add_argument('-r2', '--sample_2', required=True, help='URL to Sample: Read 2')
    parser.add_argument('-u', '--unc', required=True, help='URL to unc_hg19.bed')
    parser.add_argument('-f', '--fasta', required=True, help='URL to hg19_M_rCRS_ref.transcripts.fa')
    parser.add_argument('-c', '--composite_exons', required=True, help='URL to composite_exons.bed')
    parser.add_argument('-n', '--normalize', required=True, help='URL to normalizeBedToolsExonQuant.pl')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to output final results')
    parser.add_argument('-i', '--uuid', required=True, help='UUID for these samples')
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


def download_from_URL(job, input_args, ids, name):
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, name)
    url = input_args[name]
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
        except subprocess.CalledProcessError:
            raise RuntimeError('\nNecessary file could not be acquired: {}. Check input URL'.format(url))
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    job.fileStore.updateGlobalFile(ids[name], file_path)


def return_input_paths(job, work_dir, ids, *args):
    """
    Given one or more strings representing file_names, return the paths to those files.
    """
    # for every file in *args, place in work_dir via the FileStore and then return the mounted docker path.
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
    return os.path.join('/data', os.path.basename(filepath))


def docker_call(tool, tool_parameters, work_dir):
    """
    Makes subprocess call of a command to a docker container.
    work_dir MUST BE AN ABSOLUTE PATH or the call will fail.
    """
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    call = base_docker_call.split() + [tool] + tool_parameters
    try:
        subprocess.check_call(call)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status for cmd {}'.format(call))
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


def download_s3_directory(bucket_name, dir_name):
    """
    This requires a "~/.boto" configuration file
    that has access to the UCSC AWS Account
    """
    conn = boto.connect_s3()
    try:
        bucket = conn.get_bucket(bucket_name)
    except S3ResponseError:
        raise RuntimeError('Bucket does not exist')
    for key in bucket.list(os.path.basename(dir_name)):
        try:
            key.get_contents_to_filename(os.path.join(os.path.dirname(dir_name), key.name))
        except:
            raise RuntimeError('Necessary input file could not be acquired: {}. to Dir: {}'.format(key, dir_name))


def move_to_output_dir(work_dir, output_dir, uuid=None, files=list()):
    """
    A list of files to move from work_dir to output_dir. Must contain a UUID
    """
    for fname in files:
        if uuid is None:
            os.rename(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            os.rename(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


# Start of job Functions
def start_node(job, input_args):
    """
    Launchpoint for the script. Defines variables used in the rest of the pipeline.
    """
    symbolic_names = input_args.keys() + ['R1.fastq', 'R2.fastq', 'alignments.bam', 'stats.txt', 'rg_alignments.bam',
                                          'sorted.bam', 'sorted.bam.bai', 'sort_by_ref.bam', 'transcriptome.bam',
                                          'filtered.bam', 'rsem_gene.tab', 'rsem_isoform.tab',
                                          'rsem.genes.norm_counts.tab', 'rsem.genes.raw_counts.tab',
                                          'rsem.genes.norm_fpkm.tab', 'rsem.genes.norm_tpm.tab',
                                          'rsem.isoform.raw_counts.tab', 'rsem.isoform.norm_counts.tab',
                                          'rsem.isoform.norm_fpkm.tab', 'rsem.isoform.norm_tpm.tab',
                                          'gene_raw_count.tab', 'gene_norm_count.tab', 'gene_fpkm.tab', 'gene_tpm.tab',
                                          'isoform_raw_count.tab', 'isoform_norm_count.tab', 'isoform_fpkm.tab',
                                          'isoform_tpm.tab','stats.txt', 'stats2.txt', 'stats_all.txt', 'mapping.tab',
                                          'exon_quant', 'exon_quant.bed']
    ids = {x: job.fileStore.getEmptyFileStoreID() for x in symbolic_names}
    job_vars = (input_args, ids)

    # Download initial sample fastqs and launch pipeline
    job.addChildJobFn(download_from_URL, input_args, ids, 'sample_1.fastq.gz')
    job.addChildJobFn(download_from_URL, input_args, ids, 'sample_2.fastq.gz')
    job.addFollowOnJobFn(unzip, job_vars)


def unzip(job, job_vars):
    """
    Unzips input sample files
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()

    # Run children and follow-on
    job.addChildJobFn(unzip_worker_fn, job_vars, 'sample_1.fastq.gz', 'R1.fastq')
    job.addChildJobFn(unzip_worker_fn, job_vars, 'sample_2.fastq.gz', 'R2.fastq')
    job.addFollowOnJobFn(mapsplice, job_vars)


def unzip_worker_fn(job, job_vars, symbolic_name, output_name):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # I/O
    sample = return_input_paths(job, work_dir, ids, symbolic_name)
    # Command
    command = ['zcat', sample]
    with open(os.path.join(work_dir, output_name), 'w') as f:
        subprocess.check_call(command, stdout=f)
    # Update FileStore
    job.fileStore.updateGlobalFile(ids[output_name], os.path.join(work_dir, output_name))
    job.fileStore.deleteGlobalFile(ids[os.path.basename(sample)])


def mapsplice(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    ebwt_dir = os.path.join(work_dir, 'ebwt')
    chromsome_dir = os.path.join(work_dir, 'chromosomes')
    mkdir_p(ebwt_dir)
    mkdir_p(chromsome_dir)
    cpus = input_args['cpu_count']

    # I/O
    r1, r2 = return_input_paths(job, work_dir, ids, 'R1.fastq', 'R2.fastq')
    # Get Bowtie Index Dir and Chromosome Dir
    download_s3_directory(bucket_name='rna-seq-pipeline', dir_name=ebwt_dir)
    download_s3_directory(bucket_name='rna-seq-pipeline', dir_name=chromsome_dir)
    # Command and call
    parameters = ['-p', str(cpus),
                  '-s', '25',
                  '--bam',
                  '--min-map-len', '50',
                  '-x', docker_path(ebwt_dir),
                  '-c', docker_path(chromsome_dir),
                  '-1', docker_path(r1),
                  '-2', docker_path(r2),
                  '-o', '/data']
    docker_call(tool='computationalgenomicslab/mapsplice', tool_parameters=parameters, work_dir=work_dir)
    # TODO: TEMPORARY REPLACEMENT TO TEST OTHER FUNCTIONS IN PIPELINE
    # align_url = "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/alignments.bam"
    # stats_url = "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/stats.txt"
    # subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', align_url, '-o', os.path.join(work_dir, 'alignments.bam')])
    # subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', stats_url, '-o', os.path.join(work_dir, 'stats.txt')])

    # Update FileStore
    job.fileStore.updateGlobalFile(ids['alignments.bam'], os.path.join(work_dir, 'alignments.bam'))
    job.fileStore.updateGlobalFile(ids['stats.txt'], os.path.join(work_dir, 'stats.txt'))
    # Run child job
    job.addChildJobFn(add_read_groups, job_vars)
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
    job.addChildJobFn(bamsort_and_index, job_vars)


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
    job.addChildJobFn(sort_bam_by_reference, job_vars)


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
    job.addChildJobFn(download_transcriptome_files, job_vars)
    job.addChildJobFn(download_exon_count_files, job_vars)


def download_exon_count_files(job, job_vars):
    input_args, ids = job_vars
    job.addChildJobFn(download_from_URL, input_args, ids, 'composite_exons.bed')
    job.addChildJobFn(download_from_URL, input_args, ids,  'normalize.pl')
    job.addFollowOnJobFn(exon_count, job_vars)


def exon_count(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']

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


def download_transcriptome_files(job, job_vars):
    input_args, ids = job_vars
    job.addChildJobFn(download_from_URL, input_args, ids, 'unc.bed')
    job.addChildJobFn(download_from_URL, input_args, ids, 'hg19.transcripts.fa')
    job.addFollowOnJobFn(transcriptome, job_vars)


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
    job.addChildJobFn(filter_bam, job_vars)


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
    job.addChildJobFn(rsem, job_vars)


def rsem(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cpus = input_args['cpu_count']

    # I/O
    filtered_bam = return_input_paths(job, work_dir, ids, 'filtered.bam')
    rsem_ref_dir = os.path.join(work_dir, 'rsem_ref')
    mkdir_p(rsem_ref_dir)
    download_s3_directory(bucket_name='rna-seq-pipeline', dir_name=rsem_ref_dir)
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
    gene, isoform = return_input_paths(job, work_dir, ids, 'rsem_gene.tab', 'rsem_isoform.tab')
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
                       'rsem.isoform.raw_counts.tab', 'rsem.isoform.norm_fpkm.tab','rsem.isoform.norm_tpm.tab')

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
    inputs = {'sample_1.fastq.gz': args.sample_1,
              'sample_2.fastq.gz': args.sample_2,
              'unc.bed': args.unc,
              'hg19.transcripts.fa': args.fasta,
              'composite_exons.bed': args.composite_exons,
              'normalize.pl': args.normalize,
              'output_dir': args.output_dir,
              'uuid': args.uuid,
              'cpu_count': multiprocessing.cpu_count()}

    # Launch jobs
    Job.Runner.startToil(Job.wrapJobFn(start_node, inputs), args)
    Job.Runner.cleanup(args)
