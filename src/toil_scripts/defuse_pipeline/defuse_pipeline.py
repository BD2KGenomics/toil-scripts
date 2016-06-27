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
import textwrap
import multiprocessing
import defuse_lib as lib
import defuse_exon as exon_lib

from toil.job import Job
from toil_scripts.lib import require
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url_job
from bd2k.util.processes import which
from urlparse import urlparse


def download_shared_files(job, config):
    print(config, file=sys.stderr)



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


def defuse_pipeline(job, uuid, sample_parameters, tool_parameters):
    """
    Runs defuse pipeline

    :param job:
    :param dict sample_options: sample features
    :param dict tool_options: paramters for bioinformatics tools
    :return:
    """
    job.fileStore.logToMaster('Started deFuse pipeline for {}'.format(uuid))

    if sample_parameters['url'].endswith('.gz'):
        sample_parameters['is_gzipped'] = True
    else:
        sample_parameters['is_gzipped'] = False

    sample = job.wrapJobFn(download_sample, sample_parameters['url'], url2=sample_parameters['url2'],
                           is_gzipped=sample_parameters['is_gzipped'])
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

def move_file(job, fileStoreID, filename, univ_options):
    work_dir = job.fileStore.getLocalTempDir()

    src = os.path.join(work_dir, filename)
    job.fileStore.readGlobalFile(fileStoreID, src)

    dest = os.path.join(univ_options['output_dir'], filename)

    if os.path.exists(dest):
        job.fileStore.logToMaster('File {} already exists'.format(dest))
    else:
        shutil.copy(src, dest)


def download_sample(job, url, url2=None, is_gzipped=False):

    fq, gz = '', ''

    if is_gzipped:
        basname_fq, gz = os.path.splitext(url)
        basename, fq = os.path.splitext(basname_fq)
    else:
        basename, fq = os.path.splitext(url)

    if basename.endswith('1'):
        basename2 = re.sub('1$', '2', basename)
        url2 = ''.join([basename2, fq, gz])

    else:
        if url2 is None:
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
                mock=defuse_options['mock'], outputs={'results.tsv': results_url})

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

def generate_config():
    return textwrap.dedent("""
    # CGL Germline Variant Pipeline configuration file
    # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
    # Edit the values in this configuration file and then rerun the pipeline: "toil-variant run"
    # URLs can take the form: http://, file://, s3://, gnos://.
    # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank
    ####################################################################################################################
    # Required: URL to reference genome
    index: https://s3-us-west-2.amazonaws.com/pimmuno-test-data/CI_test_input/references/defuse_input_gencode74_with_config.tgz

    # Optional: Provide a full path to where results will appear
    output-dir:

    # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
    ssec:

    # Optional: Add True to run in mock mode
    mock:
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample pair to be run.
        #   There are 2tab-separated columns: UUID, Sample BAM URL
        #
        #   UUID            This should be a unique identifier for the sample to be processed
        #   Sample URL      A URL (http://, file://, s3://, gnos://) pointing to the normal bam
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1:
        #    url1: file:///path/to/sample_R1.fq
        #    url2: file:///path/to/sample_R2.fq
        #    adapter1: AGATCGGAAGAG
        #    adapter2: AGATCGGAAGAG
        #   UUID_2:
        #    url1: http://sample-depot.com/sample_R1.fq
        #    url2: http://sample-depot.com/sample_R2.fq
        #    adapter1: AGATCGGAAGAG
        #    adapter2: AGATCGGAAGAG
        #   UUID_3:
        #    url1: s3://my-bucket-name/directory/sample_R1.fq
        #    url2: s3://my-bucket-name/directory/sample_R2.fq
        #    adapter1: AGATCGGAAGAG
        #    adapter2: AGATCGGAAGAG
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
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}

        require(set(config) > {'index'},
                'deFuse pipeline is missing reference data. Check config file')

        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program)), program + ' must be installed on every node.'.format(program))

        if args.fq or args.uuid:
            require(args.fq and args.uuid, '"--fq" and "--uuid" must supplied')
            samples = [[args.uuid, args.fq,]]
        else:
            samples = yaml.load(open(args.manifest))

        shared_files = Job.wrapJobFn(download_shared_files, config)

        required_params = {'url1', 'adapter1', 'adapter2'}
        for uuid, parameters in samples.iteritems():
            given_params = set(parameters)
            parameter_diff = required_params - given_params
            require(given_params >= required_params,
                    'Sample {} does not have all required parameters\n{}'.format(uuid, parameter_diff))
            shared_files.addFollowOnJobFn(defuse_pipeline, uuid, parameters, shared_files.rv())

    Job.Runner.startToil(shared_files, args)

if __name__ == '__main__':
    main()
