#! /usr/bin/env python2.7
# Jacob Pfeil
# Spring 2016
"""
Toil workflow for running deFuse

Tree structure

            0
            |
            1
            |
            2
          /   \
         3    4

0 = Start node
1 = Download references
2 = Run deFuse
3 = Awk filter
4 = Exon filter
"""


from __future__ import print_function
import sys
import os
import re
import csv
import argparse
import shutil
import yaml
import multiprocessing
import defuse_lib as lib

from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url_job
from toil.job import Job


def parse_config(config_path):
    writeToDebug("Parse_config")
    if not os.path.exists(config_path):
        raise ValueError("Could not locate config file")
    if os.stat(config_path).st_size == 0:
        raise ValueError('Config file is empty')
    with open(config_path, 'r') as stream:
        return yaml.load(stream)


def get_shared_files(job, config, univ_options):
    shared_files = set(['defuse_index'])
    work_dir = job.fileStore.getLocalTempDir()

    if len(shared_files.intersection(set(config.keys()))) == 0:
        raise ValueError('Need defuse_index in config file')
    for key, value in config.iteritems():
        univ_options[key] = job.addChildJobFn(download_url_job, value, name=key, s3_key_path=univ_options['ssec']).rv()
    return univ_options


def start_pipeline(job, manifest, univ_options):
    with open(manifest, 'r') as m:
        for line in m:
            try:
                sample_name, rna_url = line.strip().split()
            except ValueError:
                raise ValueError('Manifest file is not properly formatted')
            job.addChildJobFn(defuse_pipeline, sample_name, rna_url, univ_options.copy())


def defuse_pipeline(job, sample_name, url, univ_options):
    """
    Runs defuse pipeline

    :param job:
    :param str sample_name: name for sample
    :param str url:
    :param univ_options:
    :return:
    """
    univ_options['patient_id'] = sample_name
    if url.endswith('.gz'):
        univ_options['gzipped'] = True

    get_sample = job.wrapJobFn(download_sample, url)
    defuse = job.wrapJobFn(run_defuse, get_sample.rv(0), get_sample.rv(1), univ_options)
    defuse_out = job.wrapJobFn(move_file, defuse.rv(), '{}.results.tsv'.format(sample_name), univ_options)
    awk = job.wrapJobFn(run_simple_filter, defuse.rv())
    awk_out = job.wrapJobFn(move_file, awk.rv(), '{}.results.awk.tsv'.format(sample_name), univ_options)


    job.addFollowOn(get_sample)
    get_sample.addFollowOn(defuse)
    if univ_options['filter']:
        defuse.addFollowOn(awk)
        awk.addFollowOn(awk_out)
    else:
        defuse.addFollowOn(defuse_out)


def move_file(job, fileStoreID, filename, univ_options):
    work_dir = job.fileStore.getLocalTempDir()
    filepath = os.path.join(work_dir, filename)
    job.fileStore.readGlobalFile(fileStoreID, filepath)
    shutil.copy(filepath, univ_options['output_dir'])


def download_sample(job, url):

    fq, gz = '', ''

    if url.endswith('.gz'):
        basname_fq, gz = os.path.splitext(url)
        basename, fq = os.path.splitext(basname_fq)
    else:
        basename, fq = os.path.splitext(url)

    if basename.endswith('1'):
        basename2 = re.sub('1$', '2', basename)
        url2 = ''.join([basename2, fq, gz])

    else:
        raise ValueError('Could not find paired sample. Does RNA filename end with 1?')

    file_id1 = job.addChildJobFn(download_url_job, url).rv()
    file_id2 = job.addChildJobFn(download_url_job, url2).rv()
    return file_id1, file_id2


def run_defuse(job, fastq1, fastq2, univ_options):
    cores = multiprocessing.cpu_count()

    work_dir = job.fileStore.getLocalTempDir()

    fq_extn = '.gz' if univ_options['gzipped'] else ''

    input_files = {'rna_1.fastq' + fq_extn: fastq1,
                   'rna_2.fastq' + fq_extn: fastq2,
                   'defuse_index': univ_options['defuse_index']}

    input_files = lib.get_files_from_filestore(job, input_files, work_dir, docker=True)

    parameters = ['--config', '/data/defuse_index',
                  '--output', '/data',
                  '--1fastq', '/data/rna_1.fastq',
                  '--2fastq', '/data/rna_2.fastq',
                  '--name', univ_options['patient_id'],
                  '--local', '/data/jobDir',
                  '--parallel', str(cores)]

    results_url = 'File:///home/jacob/munge/defuse/defuse_test/results.tsv'

    docker_call('jpfeil/defuse:0.6.2', parameters=parameters, work_dir=work_dir,
                mock=True, outputs={'results.tsv': results_url})

    results_path = os.path.join(work_dir, 'results.tsv')
    if not os.path.exists(results_path):
        raise RuntimeError('Defuse did not create a results.tsv file!')
    return job.fileStore.writeGlobalFile(results_path)


def run_simple_filter(job, tsvID):

    work_dir = job.fileStore.getLocalTempDir()

    input_files = {'results.tsv': tsvID}

    lib.get_files_from_filestore(job, input_files, work_dir)

    tsvFile = os.path.join(work_dir, 'results.tsv')

    outpath = os.path.join(work_dir, 'results.awk.tsv')

    features = {'splitr_count': lambda x: x > 1,
                'span_count': lambda x: x > 10,
                'orf': lambda x: x == 'Y',
                'adjacent': lambda x: x == 'N',
                'altsplice': lambda x: x == 'N',
                'min_map_count': lambda x: x >= 1,
                'gene_chromosome1': lambda x: x != 'MT',
                'gene_chromosome2': lambda x: x != 'MT'}

    featureIndex = {}
    with open(tsvFile, 'rb') as f, open(outpath, 'wb') as g:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(g, delimiter='\t')
        header = reader.next()
        writer.writerow(header)
        for feature in features:
            try:
                featureIndex[feature] = header.index(feature)
            except ValueError:
                raise ValueError('deFuse results.tsv file is in a different format. Use version 0.6.2')

        for line in reader:
            save = True
            for feature, filter in features.iteritems():
                index = featureIndex[feature]
                if not filter(line[index]):
                    save = False
                    break
            if save:
                writer.writerow(line)
    # The header is size 990
    if os.stat(outpath).st_size <= 990:
        raise ValueError('No gene fusions passed the simple filter')
    return job.fileStore.writeGlobalFile(outpath)


def writeToDebug(msgs, outfile='/tmp/toil_debug'):
    if isinstance(msgs, str):
        msgs = [msgs]
    with open(outfile, 'a') as f:
        for msg in msgs:
            f.write(msg + '\n')

def main():
    """
    This is a Toil pipeline used to perform variant analysis (usually on exomes) from Tumor/Normal BAMs.
    All samples are co-cleaned (GATK Indel Realignment (IR) and Base Quality Score Recalibration (BQSR))
    before variant analysis is performed by MuTect.  The final output of this pipeline is a tarball
    containing the output of MuTect (.vcf, .cov, .out).

    Please see the associated README.md for an overview and quickstart walkthrough.
    """
    # Define Parser object and add to jobTree
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--config', required=True, help='Config for deFuse Toil workflow')
    parser.add_argument('-m', '--manifest', required=True, help='Manifest of sample URLs')
    parser.add_argument('-f', '--filter', action='store_true', help='Filter deFuse output')
    parser.add_argument('-o', '--output_dir', default=None, help='Full path to final output dir')
    parser.add_argument('-s', '--ssec', help='A key that can be used to fetch encrypted data')
    parser.add_argument('-3', '--s3_dir', default=None, help='S3 Directory, starting with bucket name.')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    writeToDebug('Debug log')

    Job.Runner.addToilOptions(parser)

    params = parser.parse_args()

    global_options = {'ssec': params.ssec,
                      'output_dir': params.output_dir,
                      's3_dir': params.s3_dir,
                      'filter': params.filter}

    root = Job.wrapFn(parse_config, os.path.abspath(params.config))
    download_refs = Job.wrapJobFn(get_shared_files, root.rv(), global_options)
    start = Job.wrapJobFn(start_pipeline, os.path.abspath(params.manifest), download_refs.rv())

    root.addFollowOn(download_refs)
    download_refs.addFollowOn(start)


    Job.Runner.startToil(root, params)

if __name__ == '__main__':
    main()
