#!/usr/bin/env python2.7

""" 
GATK HaplotypeCaller Genotyping Mode

            Tree Structure of GATK Pipeline
            0 --> 1 --> 2 --> 3 --> 4
                                    |
                                    5
                                   / \
                                  6   7
                                 /     \
                                8       9
0 = Start Node
1 = Download Reference
2 = Index Reference
3 = Reference Dictionary
4 = Index Samples
5 = HaplotypeCaller SNP & Indel
6 = VariantRecalibrator SNPs
7 = VariantRecalibrator Indels
8 = ApplyRecalibration SNPs
9 = ApplyRecalibration Indels

===================================================================
:Dependencies:
curl            - apt-get install curl
docker          - apt-get install docker (or 'docker.io' for linux)
toil            - pip install --pre toil
"""
from __future__ import print_function
import argparse
import errno
import shutil
import os
import subprocess
import multiprocessing
from collections import OrderedDict
from toil.job import Job


def build_parser():
    """
    Create parser object containing necessary input files
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', required=True, help="Reference Genome URL")
    parser.add_argument('-f', '--config', required=True, help="Each line contains (CSV): UUID,Normal_URL,Tumor_URL")
    parser.add_argument('-p', '--phase', required=True, help='1000G_phase1.indels.b37.vcf URL')
    parser.add_argument('-m', '--mills', required=True, help='Mills_and_1000G_gold_standard.indels.b37.vcf URL')
    parser.add_argument('-d', '--dbsnp', required=True, help='dbsnp_137.b37.vcf URL')
    parser.add_argument('-n', '--omni', required=True, help='1000G_omni.5.b37.vcf URL')
    parser.add_argument('-t', '--hapmap', required=True, help='hapmap_3.3.b37.vcf URL')
    parser.add_argument('-o', '--output_dir', default="./data", help='Full path to final output dir')
    return parser


# Convenience functions used in the pipeline
def mkdir_p(path):
    """
    Makes an output directory if it does not already exist.

    :param path: path to output directory
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def docker_call(work_dir, tool_parameters, tool, input_files=None, output_files=None):
    """
    Runs docker call and checks that input and output files exist.

    :param work_dir: working directory
    :param tool_parameters: parameters for docker container tool, list
    :param tool: docker container tool, str
    :param input_files: list of input files
    :param output_files: list of output files

    :returns: None
    """
    for input_file in input_files:
        try:
            path = os.path.join(work_dir, input_file)
            assert os.path.exists(path)
        except AssertionError:
            assert os.path.exists(input_file)
    base_docker_call = 'docker run -v {work_dir}:/data'.format(work_dir=work_dir)
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_parameters)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')
    for output_file in output_files:
        try:
            path = os.path.join(work_dir, output_file)
            assert os.path.exists(path)
        except AssertionError:
            assert os.path.exists(output_file)


def download_from_url(job, url, filename):
    """
    Downloads a file from a URL and places it in the jobStore

    :param job: Job instance
    :param url: data url, str
    :param filename: name given to downloaded data, str

    :return: fileStore promise
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, filename)
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
        except subprocess.CalledProcessError:
            raise RuntimeError(
                '\nNecessary file could not be acquired: {}. Check input URL'.format(url))
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    return job.fileStore.writeGlobalFile(file_path)


def return_input_paths(job, work_dir, ids, *filenames):
    """
    Given one or more filenames, reads the file from the fileStore and
    copies it (or makes a hard link) into the working directory.

    :param job: Job instance
    :param work_dir: working directory, str
    :param ids: dictionary of shared file ids
    :param filenames: remaining arguments are filenames

    :returns: list of paths to files
    """
    paths = OrderedDict()
    for filename in filenames:
        if not os.path.exists(os.path.join(work_dir, filename)):
            file_path = job.fileStore.readGlobalFile(ids[filename], os.path.join(work_dir, filename))
        else:
            file_path = os.path.join(work_dir, filename)
        paths[filename] = file_path
        if len(args) == 1:
            return file_path
    return paths.values()


def write_to_filestore(job, work_dir, ids, *filenames):
    """
    Given one or more file names in working directory, write
    files to filestore and stores the filestore promise in a dictionary

    :param job: Job instance
    :param work_dir: working directory
    :param ids: shared files promises, dict
    :param filenames: remaining arguments are keys for ids

    :returns: shared ids dictionary
    """
    for filename in filenames:
        ids[filename] = job.fileStore.writeGlobalFile(os.path.join(work_dir, filename))
    return ids


def read_from_filestore(job, work_dir, ids, *filenames):
    """
    Given one or more file_names, move file in working directory.

    :param job: Job instance
    :param work_dir: working directory
    :param ids: shared file promises, dict
    :param filenames: remaining arguments are filenames
    """
    for filename in filenames:
        if not os.path.exists(os.path.join(work_dir, filename)):
            job.fileStore.readGlobalFile(ids[filename], os.path.join(work_dir, filename))


def move_to_output_dir(work_dir, output_dir, *filenames):
    """`
    A list of files to move from the working directory to the output directory.

    :param work_dir: the working directory
    :param output_dir: the output directory
    :param filenames: remaining arguments are filenames
    """
    for filename in filenames:
        origin = os.path.join(work_dir, filename)
        dest = os.path.join(output_dir, filename)
        shutil.move(origin, dest)


def batch_start(job, input_args):
    """ 
    Downloads reference data and stores fileStore promises,
    spawns next step in pipeline - reference indexing

    :param job: Job instance
    :param input_args: command line arguments, dict

    :returns: None
    """
    shared_files = ['ref.fa', 'phase.vcf', 'omni.vcf', 'dbsnp.vcf', 'hapmap.vcf', 'mills.vcf']
    shared_ids = {}
    for file_name in shared_files:
        url = input_args[file_name]
        shared_ids[file_name] = job.addChildJobFn(download_from_url, url, file_name).rv()
    job.addFollowOnJobFn(create_reference_index, shared_ids, input_args)


def create_reference_index(job, shared_ids, input_args):
    """
    Uses samtools to create reference index file in working directory,
    spawns next job in pipeline - create reference dictionary

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    # Unpack convenience variables for job
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = return_input_paths(job, work_dir, shared_ids, 'ref.fa')
    faidx_output = os.path.join(work_dir, 'ref.fa.fai')
    # Call: Samtools
    faidx_command = ['faidx', 'ref.fa']
    docker_call(work_dir, faidx_command, 'quay.io/ucsc_cgl/samtools', [ref_path], [faidx_output])
    # Update fileStore for output
    shared_ids['ref.fa.fai'] = job.fileStore.writeGlobalFile(faidx_output)
    job.addChildJobFn(create_reference_dict, shared_ids, input_args)


def create_reference_dict(job, shared_ids, input_args):
    """
    Uses Picardtools to create sequence dictionary for reference genome

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    # Unpack convenience variables for job
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    ref_path = return_input_paths(job, work_dir, shared_ids, 'ref.fa')
    # Call: picardtools
    picard_output = os.path.join(work_dir, 'ref.dict')
    command = ['CreateSequenceDictionary', 'R=ref.fa', 'O=ref.dict']
    docker_call(work_dir, command, 'quay.io/ucsc_cgl/picardtools', [ref_path], [picard_output])
    # Update fileStore for output
    shared_ids['ref.dict'] = job.fileStore.writeGlobalFile(picard_output)
    job.addChildJobFn(spawn_batch_jobs, shared_ids, input_args)


def spawn_batch_jobs(job, shared_ids, input_args):
    """
    Reads config file and dynamically starts a job for each sample.

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    # Names for every input file used in the pipeline by each sample
    samples = []
    config = input_args['config']
    with open(config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                samples.append(line.strip().split(','))
    for sample in samples:
        job.addChildJobFn(start, shared_ids, input_args, sample)


def start(job, shared_ids, input_args, sample):
    """
    Configures parameters for sample and starts GATK pipeline.

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    :param sample: tuple with uuid and file url
    """
    uuid, url = sample
    ids = shared_ids.copy()
    # Update input
    input_args['uuid'] = uuid
    # Sample bam file holds a url?
    input_args['bam_url'] = url
    input_args['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    if input_args['ssec'] is None:
        ids['toil.bam'] = job.addChildJobFn(download_from_url, url, 'toil.bam').rv()
    else:
        pass
    job.addFollowOnJobFn(index, ids, input_args)


def index(job, shared_ids, input_args):
    """
    Index sample bam using samtools, calls haplotypeCaller.

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    bam_path = return_input_paths(job, work_dir, shared_ids, 'toil.bam')
    output_path = os.path.join(work_dir, 'toil.bam.bai')
    # Call: index the normal.bam
    parameters = ['index', 'toil.bam']
    docker_call(work_dir, parameters, 'quay.io/ucsc_cgl/samtools', [bam_path], [output_path])
    # Update FileStore and call child
    shared_ids['toil.bam.bai'] = job.fileStore.writeGlobalFile(output_path)
    job.addChildJobFn(haplotype_caller, shared_ids, input_args)


def haplotype_caller(job, shared_ids, input_args):
    """
    Uses GATK HaplotypeCaller to identify SNPs and Indels

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    input_files = ['ref.fa', 'ref.fa.fai', 'ref.dict', 'toil.bam', 'toil.bam.bai']
    read_from_filestore(job, work_dir, shared_ids, *input_files)
    output = 'unified.raw.BOTH.gatk.vcf'
    # Call GATK -- HaplotypeCaller
    command = ['-nct', input_args['cpu_count'],
               '-R', 'ref.fa',
               '-T', 'HaplotypeCaller',
               '--genotyping_mode', 'Discovery',
               '--output_mode', 'EMIT_VARIANTS_ONLY',
               '-I', 'toil.bam',
               '-o', 'unified.raw.BOTH.gatk.vcf',
               '-stand_emit_conf', '10.0',
               '-stand_call_conf', '30.0']
    docker_call(work_dir, command, 'quay.io/ucsc_cgl/gatk', inputs, [output])
    # Update fileStore and spawn child job
    shared_ids['unified.raw.BOTH.gatk.vcf'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, output))
    job.addChildJobFn(vqsr_snp, shared_ids, input_args)
    job.addChildJobFn(vqsr_indel, shared_ids, input_args)


def vqsr_snp(job, shared_ids, input_args):
    """
    Variant quality score recalibration for SNP variants, calls SNP recalibration

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    input_files = ['ref.fa', 'ref.fa.fai', 'ref.dict', 'unified.raw.BOTH.gatk.vcf',
                   'hapmap.vcf', 'omni.vcf', 'dbsnp.vcf', 'phase.vcf']
    read_from_filestore(job, work_dir, shared_ids, *input_files)
    outputs = ['HAPSNP.recal', 'HAPSNP.tranches', 'HAPSNP.plots']
    command = ['-T', 'VariantRecalibrator',
               '-R', 'ref.fa',
               '-input', 'unified.raw.BOTH.gatk.vcf',
               '-nt', input_args['cpu_count'],
               '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', 'hapmap.vcf',
               '-resource:omni,known=false,training=true,truth=false,prior=12.0', 'omni.vcf',
               '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0', 'dbsnp.vcf',
               '-resource:1000G,known=false,training=true,truth=false,prior=10.0', 'phase.vcf',
               '-an', 'QD', '-an', 'DP', '-an', 'FS', '-an', 'ReadPosRankSum',
               '-mode', 'SNP', '-minNumBad', '1000',
               '-recalFile', 'HAPSNP.recal',
               '-tranchesFile', 'HAPSNP.tranches',
               '-rscriptFile', 'HAPSNP.plots']
    docker_call(work_dir, command, 'quay.io/ucsc_cgl/gatk', inputs, outputs)
    shared_ids = write_to_filestore(job, work_dir, shared_ids, *outputs)
    job.addChildJobFn(apply_vqsr_snp, shared_ids, input_args)


def apply_vqsr_snp(job, shared_ids, input_args):
    """
    Apply variant quality score recalibration for SNP variants.
    Writes vcf file to output directory

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)
    uuid = input_args['uuid']
    inputs = ['ref.fa', 'ref.fa.fai', 'ref.dict', 'unified.raw.BOTH.gatk.vcf',
              'HAPSNP.tranches', 'HAPSNP.recal']
    read_from_filestore(job, work_dir, shared_ids, *inputs)
    output = '{}.HAPSNP.vqsr.SNP.vcf'.format(uuid)
    command = ['-T', 'ApplyRecalibration',
               '-input', 'unified.raw.BOTH.gatk.vcf',
               '-o', output,
               '-R', 'ref.fa',
               '-nt', '1',
               '-ts_filter_level', '99.0',
               '-tranchesFile', 'HAPSNP.tranches',
               '-recalFile', 'HAPSNP.recal',
               '-mode', 'SNP']
    docker_call(work_dir, command, 'quay.io/ucsc_cgl/gatk', inputs, [output])
    move_to_output_dir(work_dir, output_dir, output)


# TODO Figure out if it is okay to use maxGaussian flag
#Indel Recalibration
def vqsr_indel(job, shared_ids, input_args):
    """
    Variant quality score recalibration for Indel variants, calls Indel recalibration

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    input_files = ['ref.fa', 'ref.fa.fai', 'ref.dict', 'unified.raw.BOTH.gatk.vcf', 'mills.vcf']
    read_from_filestore(job, work_dir, shared_ids, *input_files)
    outputs = ['HAPINDEL.recal', 'HAPINDEL.tranches', 'HAPINDEL.plots']
    command = ['-T', 'VariantRecalibrator',
               '-R', 'ref.fa',
               '-input', 'unified.raw.BOTH.gatk.vcf',
               '-nt', input_args['cpu_count'],
               '-resource:mills,known=true,training=true,truth=true,prior=12.0', 'mills.vcf',
               '-an', 'DP', '-an', 'FS', '-an', 'ReadPosRankSum',
               '-mode', 'INDEL',
               '-minNumBad', '1000',
               '-recalFile', 'HAPINDEL.recal',
               '-tranchesFile', 'HAPINDEL.tranches',
               '-rscriptFile', 'HAPINDEL.plots',
               '--maxGaussians', '4']
    docker_call(work_dir, command, 'quay.io/ucsc_cgl/gatk', inputs, outputs)
    shared_ids = write_to_filestore(job, work_dir, shared_ids, *outputs)
    job.addChildJobFn(apply_vqsr_indel, shared_ids, input_args)


def apply_vqsr_indel(job, shared_ids, input_args):
    """
    Apply variant quality score recalibration for Indel variants.
    Writes vcf file to output directory

    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)
    uuid = input_args['uuid']
    inputs = ['ref.fa', 'ref.fa.fai', 'ref.dict', 'unified.raw.BOTH.gatk.vcf',
              'HAPINDEL.recal', 'HAPINDEL.tranches', 'HAPINDEL.plots']
    read_from_filestore(job, work_dir, shared_ids, *inputs)
    output = '{}.HAPSNP.vqsr.INDEL.vcf'.format(uuid)
    command = ['-T', 'ApplyRecalibration',
               '-input', 'unified.raw.BOTH.gatk.vcf',
               '-o', output,
               '-R', 'ref.fa',
               '-nt', '1',
               '-ts_filter_level', '99.0',
               '-tranchesFile', 'HAPINDEL.tranches',
               '-recalFile', 'HAPINDEL.recal',
               '-mode', 'INDEL']
    docker_call(work_dir, command, 'quay.io/ucsc_cgl/gatk', inputs, [output])
    move_to_output_dir(work_dir, output_dir, output)

if __name__ == '__main__':
    args_parser = build_parser()
    Job.Runner.addToilOptions(args_parser)
    args = args_parser.parse_args()

    inputs = {'ref.fa': args.reference,
              'config': args.config,
              'phase.vcf': args.phase,
              'mills.vcf': args.mills,
              'dbsnp.vcf': args.dbsnp,
              'hapmap.vcf': args.hapmap,
              'omni.vcf': args.omni,
              'output_dir': args.output_dir,
              'uuid': None,
              'cpu_count': str(multiprocessing.cpu_count()),
              'ssec':None}
    
    Job.Runner.startToil(Job.wrapJobFn(batch_start, inputs), args)
