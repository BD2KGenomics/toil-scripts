#! /usr/bin/env python2.7
# Jacob Pfeil
# Spring 2016
"""
Toil workflow for running deFuse

Tree structure

                  6 -- 7 -- 8
                /
0 -- 1 -- 2 -- 3 -- 4 -- 5

0 = Start node
1 = Download references
2 = Download Sample
3 = Run cutadapt
4 = Run deFuse
5 = Awk filter
6 = Run STAR
7 = Run RSEM
8 = Exon filter
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
import defuse_exon as exon_lib


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


def get_shared_files(job, config, tool_options):
    important_parameters = {'defuse': {'index_url': ('required', 'download', 'index'),
                                       'output': ('optional', None, 'output')},
                            'gencode': {'url': ('optional', 'download', 'gencode_gtf')},
                            'bedtools': {'output': ('optional', None, 'output')},
                            'star': {'type': ('optional', None, 'type'),
                                     'ncores': ('optional', None, 'ncores'),
                                     'index_url': ('optional', 'download', 'index'),
                                     'rnaAligned.toTranscriptome.out.bam': ('optional', None,
                                                                            'rnaAligned.toTranscriptome.out.bam'),
                                     'rnaAligned.sortedByCoord.out.bam': ('optional', None,
                                                                          'rnaAligned.sortedByCoord.out.bam'),
                                     'rnaAligned.sortedByCoord.out.bam.bai': ('optional', None,
                                                                              'rnaAligned.sortedByCoord.out.bam.bai')},
                            'rsem': {'index_url': ('optional', 'download', 'index'),
                                     'ncores': ('optional', None, 'ncores'),
                                     'output': ('optional', None, 'output')}
                            }

    for tool, parameters in important_parameters.iteritems():
        for parameter, (status, action, name) in parameters.iteritems():
            if tool not in tool_options:
                tool_options[tool] = {}
            try:
                value = config[tool][parameter]

                if action == 'download':
                    tool_options[tool][name] = job.addChildJobFn(download_url_job, value, name=name,
                                                                 s3_key_path=tool_options['ssec']).rv()
                else:
                    tool_options[tool][name] = value

            except KeyError:
                if status == 'required':
                    raise KeyError('Tool {} requires parameter {}'.format(tool, parameter))
                elif status == 'optional':
                    msg = 'Depending on configuration, tool {} may require {}'.format(tool, parameter)
                    job.fileStore.logToMaster(msg)
    return tool_options


def start_pipeline(job, manifest, univ_options):
    important_paramters = ['url', 'adapter3', 'adapter5']
    with open(manifest, 'r') as m:
        for sample_name, sample_options in yaml.load(m).iteritems():
            if all(param in sample_options for param in important_paramters):
                sample_options['patient_id'] = sample_name
                job.addChildJobFn(defuse_pipeline, sample_options, univ_options.copy())
            else:
                msg = "Sample {} does not have all required parameters: {}".format(sample_name,
                                                                                   str(important_paramters))
                raise ValueError(msg)


def defuse_pipeline(job, sample_options, tool_options):
    """
    Runs defuse pipeline

    :param job:
    :param dict sample_options: sample features
    :param dict tool_options: paramters for bioinformatics tools
    :return:
    """
    if sample_options['url'].endswith('.gz'):
        sample_options['gzipped'] = True
    else:
        sample_options['gzipped'] = False

    sample = job.wrapJobFn(download_sample, sample_options)
    cutadapt = job.wrapJobFn(lib.run_cutadapt, sample.rv(0), sample.rv(1), sample_options)
    defuse = job.wrapJobFn(run_defuse, cutadapt.rv(0), cutadapt.rv(1), sample_options, tool_options['defuse'])
    defuse_out = job.wrapJobFn(move_file, defuse.rv(), '{}.results.tsv'.format(sample_options['patient_id']),
                               tool_options)
    exon_filter = job.wrapJobFn(exon_lib.exon_filter_pipeline, cutadapt.rv(0), cutadapt.rv(1),
                                defuse, sample_options, tool_options)

    awk = job.wrapJobFn(run_simple_filter, defuse.rv())
    awk_out = job.wrapJobFn(move_file, awk.rv(), '{}.results.awk.tsv'.format(sample_options['patient_id']),
                            tool_options)


    job.addFollowOn(sample)
    sample.addFollowOn(cutadapt)

    if tool_options['filter'] == 'awk':
        cutadapt.addFollowOn(defuse)
        defuse.addFollowOn(awk)
        awk.addFollowOn(awk_out)
    elif tool_options['filter'] == 'exon':
        cutadapt.addFollowOn(exon_filter)
    else:
        cutadapt.addFollowOn(defuse)
        defuse.addFollowOn(defuse_out)


def move_file(job, fileStoreID, filename, univ_options):
    work_dir = job.fileStore.getLocalTempDir()

    src = os.path.join(work_dir, filename)
    job.fileStore.readGlobalFile(fileStoreID, src)

    dest = os.path.join(univ_options['output_dir'], filename)

    if os.path.exists(dest):
        job.fileStore.logToMaster('File {} already exists'.format(dest))
    else:
        shutil.copy(src, dest)


def download_sample(job, sample_options):

    fq, gz = '', ''

    url = sample_options['url']

    if url.endswith('.gz'):
        basname_fq, gz = os.path.splitext(url)
        basename, fq = os.path.splitext(basname_fq)
    else:
        basename, fq = os.path.splitext(url)

    if basename.endswith('1'):
        basename2 = re.sub('1$', '2', basename)
        url2 = ''.join([basename2, fq, gz])

    elif 'url2' in sample_options.keys():
        url2 = sample_options['url2']

    else:
        raise ValueError('Could not find paired sample. Use url2 parameter to specify paired reads')

    file_id1 = job.addChildJobFn(download_url_job, url).rv()
    file_id2 = job.addChildJobFn(download_url_job, url2).rv()
    return file_id1, file_id2


def run_defuse(job, fastq1, fastq2, sample_options, defuse_options):
    job.fileStore.logToMaster('Running deFuse on %s' % sample_options['patient_id'])

    work_dir = job.fileStore.getLocalTempDir()

    fq_extn = '.gz' if sample_options['gzipped'] else ''

    input_files = {'rna_1.fastq' + fq_extn: fastq1,
                   'rna_2.fastq' + fq_extn: fastq2,
                   'defuse_index': defuse_options['index']}

    input_files = lib.get_files_from_filestore(job, input_files, work_dir, docker=True)

    cores = multiprocessing.cpu_count()
    parameters = ['--config', '/data/defuse_index',
                  '--output', '/data',
                  '--1fastq', '/data/rna_1.fastq',
                  '--2fastq', '/data/rna_2.fastq',
                  '--name', sample_options['patient_id'],
                  '--local', '/data/jobDir',
                  '--parallel', str(cores)]

    results_url = defuse_options['output']

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
        print(header)
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
    parser.add_argument('-o', '--output_dir', required=True, help='Full path to final output dir')
    parser.add_argument('-f', '--filter', default=None, help='Filter deFuse output')
    parser.add_argument('-s', '--ssec', default=None, help='A key that can be used to fetch encrypted data')
    parser.add_argument('-3', '--s3_dir', default=None, help='S3 Directory, starting with bucket name.')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    writeToDebug('Debug log')

    Job.Runner.addToilOptions(parser)

    params = parser.parse_args()

    if params.filter not in {None, 'awk', 'exon'}:
        raise ValueError('Filter must be either None, awk, or exon')

    global_options = {'ssec': params.ssec,
                      'output_dir': os.path.abspath(params.output_dir),
                      's3_dir': params.s3_dir,
                      'filter': params.filter}

    root = Job.wrapFn(parse_config, os.path.abspath(params.config))
    download_refs = Job.wrapJobFn(get_shared_files, root.rv(), global_options)
    update_parameters = Job.wrapJobFn(exon_lib.prepare_gtf, download_refs.rv())
    start = Job.wrapJobFn(start_pipeline, os.path.abspath(params.manifest), update_parameters.rv())

    root.addFollowOn(download_refs)
    download_refs.addFollowOn(update_parameters)
    update_parameters.addFollowOn(start)


    Job.Runner.startToil(root, params)

if __name__ == '__main__':
    main()
