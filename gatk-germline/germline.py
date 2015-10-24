#!/usr/bin/env python2.7

""" 
GATK Haplotype VARIANT Callling

    Tree Structure of GATK Pipeline (per sample)
                         0 --> 1 --> 2 --> 3
                                          / \  
                                         4   5
                                        /     \
                                       6       7
0 = Start node
1 = Download references
2 = Index Samples
3 = HaplotypeCaller SNP & Indel
4 = VariantRecalibrator SNPs
5 = VariantRecalibrator Indels
6 = ApplyRecalibration SNPs
7 = ApplyRecalibration Indels

===================================================================
:Dependencies:
curl            - apt-get install curl
docker          - apt-get install docker (or 'docker.io' for linux)
jobTree         - https://github.com/benedictpaten/jobTree
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

debug = False

def build_parser():
    """
    Contains arguments for the all of necessary input files
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
    It is Easier to Ask for Forgiveness than Permission
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def docker_path(file_path):
    """
    Returns file names to the perceived path inside of the tools' Docker containers
    """
    return os.path.join('/data', os.path.basename(file_path))


# TODO Assert files exist
def docker_call(work_dir, tool_parameters, tool, input_files=None, output_files=None):
    """
    Takes parameters to run a docker container and checks input and output files
    """
    for input_file in input_files:
        assert os.path.exists(input_file)
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    if debug:
        base_docker_call = 'echo {}'.format(base_docker_call)
        for output_file in output_files:
            f = open(output_file, 'w')
            f.write('debug')
            f.close()
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_parameters)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')
    for output_file in output_files:
        assert os.path.exists(output_file)


def download_from_url(job, url, fname):
    """
    Downloads a file from a URL and places it in the jobStore
    Input1: Toil job instance
    Input2: Input arguments
    Input3: jobstore id dictionary
    Input4: Name of key used to access url in input_args
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, fname)
    if not os.path.exists(file_path):
        try:
            if debug:
                f = open(file_path, 'w')
                f.write('debug')
                f.close()
            else:
                subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
        except subprocess.CalledProcessError:
            raise RuntimeError(
                '\nNecessary file could not be acquired: {}. Check input URL'.format(url))
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    return job.fileStore.writeGlobalFile(file_path)


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


def move_to_output_dir(work_dir, output_dir, uuid=None, filenames=None):
    """`
    A list of files to move from work_dir to output_dir.
    """
    for fname in filenames:
        if uuid is None:
            origin = os.path.join(work_dir, fname)
            dest = os.path.join(output_dir, fname)
            debug_log.write(origin)
            debug_log.write(dest)
            shutil.move(origin, dest)
        else:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


# TODO I can download the sequences but it crashes after that...
# Added childjob instead of follow on job
def batch_start(job, input_args):
    """ 
    Downloads shared files that are used by all samples for alignment and calibration and places 
    them in job store
    """
    shared_files = ['ref.fa', 'phase.vcf', 'omni.vcf', 'dbsnp.vcf', 'hapmap.vcf', 'mills.vcf']
    shared_ids = {}
    for file_name in shared_files:
        url = input_args[file_name]
        shared_ids[file_name] = job.addChildJobFn(download_from_url, url, file_name).rv()
    job.addFollowOnJobFn(spawn_batch_jobs, shared_ids, input_args)


def spawn_batch_jobs(job, shared_ids, input_args):
    """
    Reads in the input files with uuid and url information and starts a job
    for each sample
    """
    #Names for every input file used in the pipeline by each sample
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
    Takes a sample and creates a copy of shared data, then
    starts the index function
    """
    uuid, url = sample
    ids = shared_ids.copy()
    #Update input
    input_args['uuid'] = uuid
    #Sample bam file holds a url?
    input_args['bam_url'] = url
    input_args['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    if input_args['ssec'] is None:
        ids['bam'] = job.addChildJobFn(download_from_url, url, 'bam').rv()
    else:
        pass
    job.addChildJobFn(index, ids, input_args)


# TODO remove debug path
def index(job, ids, input_args):
    work_dir = job.fileStore.getLocalTempDir()
    #Retrieve file path
    bam_path = return_input_paths(job, work_dir, ids, 'bam')
    #Call: index the normal.bam
    parameters = ['index', bam_path]
    output_path = '{}.bai'.format(bam_path)
    docker_call(work_dir, parameters, 'computationalgenomicslab/samtools', [bam_path], [output_path])
    #Update FileStore and call child
    ids['bai'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, output_path))
    job.addChildJobFn(haplotype_caller, ids, input_args)


##TODO There is something wrong with this function ... it crashes every time
def haplotype_caller(job, ids, input_args):
    """ snps & indels together """
    work_dir = job.fileStore.getLocalTempDir()
    ref_fasta, bam, bai = return_input_paths(job, work_dir, ids, 'ref.fa', 'bam', 'bai')
    output = os.path.join(work_dir, 'unified.raw.BOTH.gatk.vcf')
    #Call GATK -- HaplotypeCaller
    command = ['-nct', input_args['cpu_count'],
               '-R', ref_fasta,
               '-T', 'HaplotypeCaller',
               '--genotyping_mode', 'Discovery',
               '--output_mode', 'EMIT_VARIANTS_ONLY',
               '-I', bam,
               '-o', output,
               '-stand_emit_conf', '10.0',
               '-stand_call_conf', '30.0']
    docker_call(work_dir, command, 'computationalgenomicslab/gatk', [ref_fasta, bam, bai], [output])
    #Update fileStore and spawn child job
    ids['unified.raw.BOTH.gatk.vcf'] = job.fileStore.writeGlobalFile(output)
    job.addChildJobFn(vqsr_snp, ids, input_args)
    job.addChildJobFn(vqsr_indel, ids, input_args)


def vqsr_snp(job, ids, input_args):
    work_dir = job.fileStore.getLocalTempDir()
    ref_fasta, raw_vcf, hapmap, omni, \
    dbsnp, phase = return_input_paths(job, work_dir, ids, 'ref.fa', 'unified.raw.BOTH.gatk.vcf',
                                      'hapmap.vcf', 'omni.vcf', 'dbsnp.vcf', 'phase.vcf')
    inputs = [ref_fasta, raw_vcf, hapmap, omni, dbsnp, phase]
    recalFile = os.path.join(work_dir, 'HAPSNP.recalFile')
    tranches = os.path.join(work_dir, 'HAPSNP.tranches')
    rscriptFile = os.path.join(work_dir, 'HAPSNP.plots')
    outputs = [recalFile, tranches, rscriptFile]
    command = ['-T', 'VariantRecalibrator',
               '-R', ref_fasta,
               '-input', raw_vcf,
               '-nt', input_args['cpu_count'],
               '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', hapmap,
               '-resource:omni,known=false,training=true,truth=false,prior=12.0', omni,
               '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0', dbsnp,
               '-resource:1000G,known=false,training=true,truth=false,prior=10.0', phase,
               '-an QD', '-an DP', '-an FS', '-an ReadPosRankSum', '-mode SNP', '-minNumBad 1000',
               '-recalFile', recalFile,
               '-tranchesFile', tranches,
               '-rscriptFile', rscriptFile]
    docker_call(work_dir, command, 'computationalgenomicslab/gatk', inputs, outputs)
    ids['HAPSNP.recal'] = job.fileStore.writeGlobalFile(recalFile)
    ids['HAPSNP.tranches'] = job.fileStore.writeGlobalFile(tranches)
    ids['HAPSNP.plots'] = job.fileStore.writeGlobalFile(rscriptFile)
    job.addChildJobFn(apply_vqsr_snp, ids, input_args)


#TODO Figure out out to move to output_dir
#Apply Snp Recalibration
def apply_vqsr_snp(job, ids, input_args):
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)
    uuid = input_args['uuid']
    output_name = '{}.HAPSNP.vqsr.SNP.vcf'.format(uuid)
    output_path = os.path.join(work_dir, output_name)
    ref_fasta, raw_vcf, tranches, recal = return_input_paths(job, work_dir, ids, 'ref.fa', 'unified.raw.BOTH.gatk.vcf',
                                                             'HAPSNP.tranches', 'HAPSNP.recal')
    inputs = [ref_fasta, raw_vcf, tranches, recal]
    command = ['-T', 'ApplyRecalibration',
               '-input', raw_vcf,
               '-o', output_path, '-R', ref_fasta,
               '-nt', input_args['cpu_count'],
               '-ts_filter_level 99.0',
               '-tranchesFile', tranches,
               '-recalFile', recal,
               '-mode', 'SNP']
    docker_call(work_dir, command, 'computationalgenomicslab/gatk', inputs, [output_path])
    move_to_output_dir(work_dir, output_dir, uuid=None, filenames=[output_name])


#Indel Recalibration
def vqsr_indel(job, ids, input_args):
    work_dir = job.fileStore.getLocalTempDir()
    ref_fasta, raw_vcf, mills = return_input_paths(job, work_dir, ids, 'ref.fa', 'unified.raw.BOTH.gatk.vcf',
                                                   'mills.vcf')
    inputs = [ref_fasta, raw_vcf, mills]
    recalFile = os.path.join(work_dir, 'HAPINDEL.recalFile')
    tranches = os.path.join(work_dir, 'HAPINDEL.tranches')
    rscriptFile = os.path.join(work_dir, 'HAPINDEL.plots')
    outputs = [recalFile, tranches, rscriptFile]
    command = ['-T', 'VariantRecalibrator',
               '-R', ref_fasta,
               '-input', raw_vcf,
               '-nt', input_args['cpu_count'],
               '-resource:mills,known=true,training=true,truth=true,prior=12.0', mills,
               '-an DP', '-an FS', '-an ReadPosRankSum', '-mode INDEL',
               '-minNumBad 1000',
               '-recalFile', recalFile,
               '-tranchesFile', tranches,
               '-rscriptFile', rscriptFile]
    docker_call(work_dir, command, 'computationalgenomicslab/gatk', inputs, outputs)
    ids['HAPINDEL.recal'] = job.fileStore.writeGlobalFile(recalFile)
    ids['HAPINDEL.tranches'] = job.fileStore.writeGlobalFile(tranches)
    ids['HAPINDEL.plots'] = job.fileStore.writeGlobalFile(rscriptFile)
    job.addChildJobFn(apply_vqsr_indel, ids, input_args)

#Apply Indel Recalibration
def apply_vqsr_indel(job, ids, input_args):
    work_dir = job.fileStore.getLocalTempDir()
    output_dir = input_args['output_dir']
    mkdir_p(output_dir)
    uuid = input_args['uuid']
    output_name = '{}.HAPSNP.vqsr.INDEL.vcf'.format(uuid)
    output_path = os.path.join(work_dir, output_name)
    ref_fasta, raw_vcf, tranches, recal = return_input_paths(job, work_dir, ids,
                                                             'ref.fa',
                                                             'unified.raw.BOTH.gatk.vcf',
                                                             'HAPINDEL.tranches',
                                                             'HAPINDEL.recal')
    inputs = [ref_fasta, raw_vcf, tranches, recal]
    command = ['-T', 'ApplyRecalibration',
               '-input', raw_vcf,
               '-o', output_path,
               '-R', ref_fasta,
               '-nt', input_args['cpu_count'],
               '-ts_filter_level', '99.0',
               '-tranchesFile', tranches,
               '-recalFile', recal,
               '-mode', 'INDEL']
    docker_call(work_dir, command, 'computationalgenmicslab/gatk', inputs, [output_path])
    move_to_output_dir(work_dir, output_dir, uuid=None, filenames=[output_name])


if __name__ == '__main__':
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()


    input_args = {'ref.fa': args.reference,
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
    
    Job.Runner.startToil(Job.wrapJobFn(batch_start, input_args), args)
