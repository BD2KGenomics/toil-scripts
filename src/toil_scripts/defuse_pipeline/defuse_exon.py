"""
Module for Arjun's exon filter
"""

from __future__ import print_function

import os
import csv
import time
import multiprocessing
import defuse_lib as lib

from operator import itemgetter
from collections import defaultdict
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.files import copy_files_to

def exon_filter_pipeline(job, fastq1, fastq2, defuse_job, sample_options, tool_options):
    star = job.wrapJobFn(run_star, fastq1, fastq2, sample_options, tool_options['star'])
    bedtools = job.wrapJobFn(run_bedtools_coverage, star.rv('rnaAligned.sortedByCoord.out.bam'),
                             tool_options['gencode'],
                             sample_options, tool_options['bedtools'])
    rsem = job.wrapJobFn(lib.run_rsem, star.rv('rnaAligned.toTranscriptome.out.bam'),
                         sample_options, tool_options['rsem'])

    filter = job.wrapJobFn(run_exon_filter, bedtools.rv(), rsem.rv(), defuse_job.rv(),
                           sample_options, tool_options)

    job.addChild(defuse_job)
    defuse_job.addChild(star)
    star.addFollowOn(rsem)
    rsem.addChild(bedtools)
    bedtools.addFollowOn(filter)



def prepare_gtf(job, tool_options):
    if tool_options['filter'] == 'exon':
        try:
            gencode_gtf = tool_options['gencode']['gencode_gtf']
            tool_options['exon_gtf'] = job.addChildJobFn(get_exon_gtf, gencode_gtf).rv()
        except KeyError:
            raise ValueError('Requires gencode gtf for exon filtering')
    return tool_options


def get_exon_gtf(job, gencode_id):
    work_dir = job.fileStore.getLocalTempDir()

    input_files = {'gencode.gtf': gencode_id}
    lib.get_files_from_filestore(job, input_files, work_dir)

    gtf_path = os.path.join(work_dir, 'gencode.gtf')
    output_path = os.path.join(work_dir, 'gencode.exons.gtf')

    reader = csv.reader(open(gtf_path, 'rb'), delimiter='\t')
    line = reader.next()
    headers = []
    while line[0].startswith('##'):
        headers.append(line)
        line = reader.next()

    outfile = open(output_path, 'w')
    for header in headers:
        outfile.write('\t'.join(header) + '\n')

    if line[2] == 'exon':
        outfile.write('\t'.join(line) + '\n')

    for line in reader:
        if line[2] == 'exon':
            outfile.write('\t'.join(line) + '\n')
    outfile.close()
    exon_gtf = job.fileStore.writeGlobalFile(output_path)
    return job.addChildJobFn(sort_gtf, exon_gtf).rv()

def sort_gtf(job, gtf_id):
    job.fileStore.logToMaster('Sorting GTF')
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {'gencode.gtf': gtf_id}
    lib.get_files_from_filestore(job, input_files, work_dir)

    gtf_path = os.path.join(work_dir, 'gencode.gtf')
    output_path = os.path.join(work_dir, 'gencode.sorted.gtf')

    reader = csv.reader(open(gtf_path, 'rb'), delimiter='\t')
    line = reader.next()
    headers = []
    while line[0].startswith('##'):
        headers.append(line)
        line = reader.next()

    outfile = open(output_path, 'w')
    for header in headers:
        outfile.write('\t'.join(header) + '\n')

    data = defaultdict(list)

    chromosome_format = 'chr'
    chromosomes = ['{}{}'.format(chromosome_format, num) for num in range(1, 23) + ['X', 'Y', 'M']]

    chromosome = 0
    start_pos = 3

    while True:
        if len(line) == 9:
            data[line[chromosome]].append((line, int(line[start_pos])))
        else:
            raise ValueError("GTF file must have 9 features")
        try:
            line = reader.next()
        except StopIteration:
            break

    for chrom in chromosomes:
        gtfs = sorted(data[chrom], key=itemgetter(1))
        lines = ['\t'.join(gtf) for gtf, _ in gtfs]
        for line in lines:
            outfile.write(line + '\n')
    return job.fileStore.writeGlobalFile(output_path)

def run_star(job, fastq1, fastq2, sample_options, star_options):
    """
    This module uses STAR to align the RNA fastqs to the reference

    ARGUMENTS
    1. fastqs: REFER RETURN VALUE of run_cutadapt()
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
              +- 'dockerhub': <dockerhub to use>
    3. star_options: Dict of parameters specific to STAR/home/ja/home/jacob/code/protectcob/code/protect
         star_options
             |- 'index_tar': <JSid for the STAR index tarball>
             +- 'n': <number of threads to allocate>
    RETURN VALUES
    1. output_files: Dict of aligned bams
         output_files
             |- 'rnaAligned.toTranscriptome.out.bam': <JSid>
             +- 'rnaAligned.sortedByCoord.out.bam': Dict of genome bam + bai
                                |- 'rna_fix_pg_sorted.bam': <JSid>
                                +- 'rna_fix_pg_sorted.bam.bai': <JSid>

    This module corresponds to node 9 on the tree
    """
    job.fileStore.logToMaster('Running STAR on %s' % sample_options['patient_id'])

    assert star_options['type'] in ('star', 'starlong')
    work_dir = job.fileStore.getLocalTempDir()

    input_files = {
        'rna_cutadapt_1.fastq': fastq1,
        'rna_cutadapt_2.fastq': fastq2,
        'star_index': star_options['index']}

    input_files = lib.get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    parameters = ['--runThreadN', str(star_options['ncores']),
                  '--genomeDir', input_files['star_index'],
                  '--outFileNamePrefix', 'rna',
                  '--readFilesIn',
                  input_files['rna_cutadapt_1.fastq'],
                  input_files['rna_cutadapt_2.fastq'],
                  '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                  '--outSAMtype', 'BAM', 'SortedByCoordinate',
                  '--quantMode', 'TranscriptomeSAM',
                  '--outSAMunmapped', 'Within']

    result_files = ['rnaAligned.toTranscriptome.out.bam',
                    'rnaAligned.sortedByCoord.out.bam',
                    'rnaAligned.sortedByCoord.out.bam.bai']

    result_urls = {key: star_options[key] for key in result_files}

    if star_options['type'] == 'star':
        docker_call(tool='star', parameters=parameters, work_dir=work_dir, mock=True, outputs=result_urls)
    else:
        docker_call(tool='starlong', parameters=parameters, work_dir=work_dir, mock=True, outputs=result_urls)
    output_files = defaultdict()
    for bam_file in ['rnaAligned.toTranscriptome.out.bam',
                     'rnaAligned.sortedByCoord.out.bam']:
        output_files[bam_file] = job.fileStore.writeGlobalFile('/'.join([
            work_dir, bam_file]))
    job.fileStore.deleteGlobalFile(fastq1)
    job.fileStore.deleteGlobalFile(fastq2)
    return output_files

    # index_star = job.wrapJobFn(index_bamfile,
    #                            output_files['rnaAligned.sortedByCoord.out.bam'],
    #                            'rna', sample_options, disk='120G')
    # job.addChild(index_star)
    # output_files['rnaAligned.sortedByCoord.out.bam'] = index_star.rv()
    # return output_files


def run_bedtools_coverage(job, bam, gencode_options, sample_options, bedtools_options):
    job.fileStore.logToMaster('Running deFuse on %s' % sample_options['patient_id'])

    work_dir = job.fileStore.getLocalTempDir()

    input_files = {'bedtools.bam': bam,
                   'bedtools.gtf': gencode_options['gencode_gtf']}

    input_files = lib.get_files_from_filestore(job, input_files, work_dir, docker=True)

    cores = multiprocessing.cpu_count()
    parameters = ['coverage',
                  '-a', '/data/bedtools.gtf',
                  '-b', '/data/bedtools.bam',
                  '--1fastq', '/data/rna_1.fastq',
                  '--2fastq', '/data/rna_2.fastq',
                  '--name', sample_options['patient_id'],
                  '-split', '-s', '-sorted']

    results_url = bedtools_options['output']

    docker_call('aarjunrao/bedtools2', parameters=parameters, work_dir=work_dir,
                mock=True, outputs={'bedtools.cov': results_url})

    results_path = os.path.join(work_dir, 'bedtools.cov')
    return job.fileStore.writeGlobalFile(results_path)


def run_exon_filter(job, exon_cov, rsem_out, defuse_out, sample_options, tool_options):
    job.fileStore.logToMaster('Running exon filter on %s' % sample_options['patient_id'])

    work_dir = job.fileStore.getLocalTempDir()

    input_files = {'results.tsv': defuse_out,
                   'exon_coverage': exon_cov,
                   'rsem.isoforms': rsem_out}

    input_files = lib.get_files_from_filestore(job, input_files, work_dir)

    fusion_file = open(input_files['results.tsv'], 'r')

    output_path = os.path.join(work_dir, 'exon_debugger')

    output_file = open(output_path, 'w')

    fusion_reader = csv.reader(fusion_file, delimiter='\t')

    for line in fusion_reader:
        job.fileStore.logToMaster(repr(line))





