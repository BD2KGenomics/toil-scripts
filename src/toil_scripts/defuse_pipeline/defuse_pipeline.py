#! /usr/bin/env python2.7
# Jacob Pfeil
# Spring 2016
"""
Toil workflow for running deFuse

0 = Start node
1 = Download references
2 = Download Sample
3 = Run cutadapt
4 = Run deFuse
5 = Optional ilter
"""

from __future__ import print_function
import sys
import os
import re
import csv
import argparse
import shutil
import yaml
import textwrap
import multiprocessing
import defuse_lib as lib

from toil.job import Job
from toil_scripts.lib import require
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url_job
from bd2k.util.processes import which
from collections import defaultdict


def download_shared_files(job, config):
    required_files = {'defuse': [('index_url', 'index')]}
    for tool, url_keys in required_files.iteritems():
        for (key, name) in url_keys:
            if name is None:
                name = key
            config[tool][name] = job.addChildJobFn(download_url_job, config[tool][key], s3_key_path=config['ssec']).rv()
    return config


def defuse_pipeline(job, uuid, sample_options, config):
    """
    Runs defuse pipeline

    :param job:
    :param dict sample_options: sample features
    :param dict tool_options: paramters for bioinformatics tools
    :return:
    """
    job.fileStore.logToMaster('Started deFuse pipeline for sample: {}'.format(uuid))

    # Unpack sample parameters
    url = sample_options['url']
    url2 = sample_options.get('url2', None)
    adapter3 = sample_options['adapter3']
    adapter5 = sample_options['adapter5']
    is_gzipped = False
    if url.endswith('.gz'):
        is_gzipped = True

    sample = job.wrapJobFn(download_sample, url, url2, is_gzipped)
    cutadapt = job.wrapJobFn(lib.run_cutadapt, uuid, sample.rv(0), sample.rv(1), adapter3, adapter5,
                             is_gzipped=is_gzipped, mock=config['mock'])
    defuse = job.wrapJobFn(run_defuse, uuid, cutadapt.rv(0), cutadapt.rv(1), config['defuse'],
                           is_gzipped=is_gzipped, mock=config['mock'])
    defuse_out = job.wrapJobFn(move_file, defuse.rv(), '{}.results.tsv'.format(uuid),
                               config['output_dir'])
    filter = job.wrapJobFn(run_simple_filter, uuid, defuse.rv(), mock=config['mock'])
    filter_out = job.wrapJobFn(move_file, filter.rv(), '{}.results.filter.tsv'.format(uuid), config['output_dir'])

    job.addFollowOn(sample)
    sample.addFollowOn(cutadapt)
    cutadapt.addFollowOn(defuse)
    defuse.addFollowOn(defuse_out)

    if config['filter']:
        defuse.addFollowOn(filter)
        filter.addFollowOn(filter_out)


def move_file(job, fileStoreID, filename, output_dir):
    work_dir = job.fileStore.getLocalTempDir()
    src = os.path.join(work_dir, filename)
    job.fileStore.readGlobalFile(fileStoreID, src)

    if output_dir is None:
        output_dir = os.getcwd()
    dest = os.path.join(output_dir, filename)

    if os.path.exists(dest):
        job.fileStore.logToMaster('File {} already exists'.format(dest))
    else:
        shutil.copy(src, dest)


def download_sample(job, url, url2=None, is_gzipped=False):
    job.fileStore.logToMaster("Downloading url: {}".format(url))
    fq, gz = '', ''

    # If url2 is not provided, try to find it using the provided URL
    if url2 is None:
        if is_gzipped:
            basname_fq, gz = os.path.splitext(url)
            basename, fq = os.path.splitext(basname_fq)
        else:
            basename, fq = os.path.splitext(url)

        if basename.endswith('1'):
            basename2 = re.sub('1$', '2', basename)
            url2 = ''.join([basename2, fq, gz])
        else:
            raise ValueError('Could not find paired sample. Use url2 parameter to specify paired reads')

    file_id1 = job.addChildJobFn(download_url_job, url).rv()
    file_id2 = job.addChildJobFn(download_url_job, url2).rv()
    return file_id1, file_id2


def run_defuse(job, uuid, fastq1, fastq2, tool_options, mock=False, is_gzipped=False):
    job.fileStore.logToMaster('Running deFuse on %s' % uuid)
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()

    fq_extn = '.gz' if is_gzipped else ''
    input = {'rna_1.fastq' + fq_extn: fastq1,
             'rna_2.fastq' + fq_extn: fastq2,
             'defuse_index': tool_options['index']}
    input = lib.get_files_from_filestore(job, input, work_dir, docker=True)
    parameters = ['--config', '/data/defuse_index',
                  '--output', '/data',
                  '--1fastq', '/data/rna_1.fastq',
                  '--2fastq', '/data/rna_2.fastq',
                  '--name', uuid,
                  '--local', '/data/jobDir',
                  '--parallel', str(cores)]

    results_url = tool_options.get('mock_output', None)
    output = {'results.tsv': results_url}
    docker_call('jpfeil/defuse:0.6.2', parameters=parameters, work_dir=work_dir,
                mock=mock, inputs=input.keys(), outputs=output)
    results_path = os.path.join(work_dir, 'results.tsv')
    return job.fileStore.writeGlobalFile(results_path)


def run_simple_filter(job, uuid, tsvID, mock=False):
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {'results.tsv': tsvID}
    lib.get_files_from_filestore(job, input_files, work_dir)
    tsvFile = os.path.join(work_dir, 'results.tsv')
    outpath = os.path.join(work_dir, 'results.filter.tsv')

    if mock:
        with open(outpath, 'w') as f:
            f.write('content')
        return job.fileStore.writeGlobalFile(outpath)

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
    # The header is around 990 bytes
    if os.stat(outpath).st_size <= 990:
        job.fileStore.logToMaster('No gene fusions passed the filter for sample: {}'.format(uuid))
    return job.fileStore.writeGlobalFile(outpath)

def generate_config():
    return textwrap.dedent("""
    # CGL Germline Variant Pipeline configuration file
    # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
    # Edit the values in this configuration file and then rerun the pipeline: "toil-variant run"
    # URLs can take the form: http://, file://, s3://, gnos://.
    # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank
    ####################################################################################################################
    # Required: URL to defuse index tar file
    defuse:
     index_url:
     mock_output:

    # Optional: If True, filter defuse output
    filter:

    # Optional: Provide a full path to where results will appear
    output-dir:

    # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
    ssec:

    # Optional: If true, run in mock mode
    mock:
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample pair to be run.
        #   Manifest is written in yaml format
        #
        #   UUID            This should be a unique identifier for the sample to be processed
        #   url:      A URL (http://, file://, s3://, gnos://) pointing to the RNA-Seq R1 fastq
        #   url2:      A URL (http://, file://, s3://, gnos://) pointing to the RNA-Seq R2 fastq
        #   adapter1:  3' adapter sequence
        #   adapter2:  5' adapter sequence
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1:
        #    url:  file:///path/to/sample_R1.fq
        #    url2: file:///path/to/sample_R2.fq
        #    adapter3: AGATCGGAAGAG
        #    adapter5: AGATCGGAAGAG
        #   UUID_2:
        #    url:  http://sample-depot.com/sample_R1.fq
        #    url2: http://sample-depot.com/sample_R2.fq
        #    adapter3: AGATCGGAAGAG
        #    adapter5: AGATCGGAAGAG
        #   UUID_3:
        #    url:  s3://my-bucket-name/directory/sample_R1.fq
        #    url2: s3://my-bucket-name/directory/sample_R2.fq
        #    adapter3: AGATCGGAAGAG
        #    adapter5: AGATCGGAAGAG
        #
        #   Place your samples below, one per line.
        """[1:])

def writeToDebug(msgs, outfile='/tmp/toil_debug'):
    if isinstance(msgs, str):
        msgs = [msgs]
    with open(outfile, 'a') as f:
        for msg in msgs:
            f.write(msg + '\n')


def generate_file(file_path, generate_func):
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))

def check_for_required_parameters(config):
    """
    Parses config dictionary and checks for missing parameters
    :param config:
    :return:
    """
    missing_params = []
    required_tool_params = {'defuse': ['index_url']}
    for tool, params in required_tool_params.iteritems():
        for param in params:
            try:
                config[tool][param]
            except KeyError:
                missing_params.append((tool, param))
    if missing_params:
        raise ValueError("Missing following parameters in config file:\n{}".format('\n'.join(missing_params)))


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
    subparsers = parser.add_subparsers(dest='command')
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the CGL defuse pipeline')
    parser_run.add_argument('--config', default='config-toil-defuse.yaml', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--manifest', default='manifest-toil-defuse.tsv', type=str,
                            help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--fq', default=None, type=str,
                            help='URL for the sample BAM. URLs can take the form: http://, file://, s3://, '
                                 'and gnos://. The UUID for the sample must be given with the "--uuid" flag.')
    parser_run.add_argument('--uuid', default=None, type=str, help='Provide the UUID of a sample when using the'
                                                                   '"--bam" option')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    writeToDebug('Debug log')

    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()

    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'config-toil-defuse.yaml'), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-toil-defuse.tsv'), generate_manifest)
    if 'generate' in args.command:
        sys.exit()
    if args.command == 'run':
        # Read in config yaml file and set the default value to None
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        check_for_required_parameters(config)

        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program)), program + ' must be installed on every node.'.format(program))

        if args.fq or args.uuid:
            require(args.fq and args.uuid, '"--fq" and "--uuid" must supplied')
            samples = [[args.uuid, args.fq]]
        else:
            samples = yaml.load(open(args.manifest))

        download_config_files = Job.wrapJobFn(download_shared_files, config)

        required_params = {'url', 'adapter3', 'adapter5'}
        for uuid, parameters in samples.iteritems():
            given_params = set(parameters)
            parameter_diff = required_params - given_params
            require(given_params >= required_params,
                    'Sample {} does not have required parameters\n{}'.format(uuid, parameter_diff))
            download_config_files.addFollowOnJobFn(defuse_pipeline, uuid, parameters, download_config_files.rv())

    Job.Runner.startToil(download_config_files, args)

if __name__ == '__main__':
    main()
