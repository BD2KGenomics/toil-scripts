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


def docker_call(work_dir, tool_parameters, tool):
    """
    Makes subprocess call of a command to a docker container.
    :parameter tool: name of the Docker image to be used (e.g. computationalgenomicslab/samtools)
    """
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_parameters)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')

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
            file_path = docker_path(job.fileStore.readGlobalFile(ids[name], os.path.join(work_dir, name)))
        else:
            file_path = docker_path(name)
        paths[name] = file_path
        if len(args) == 1:
            paths = file_path
    return paths


def move_to_output_dir(work_dir, output_dir, uuid=None, files=None):
    """
    A list of files to move from work_dir to output_dir.
    """
    for fname in files:
        if uuid is None:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.move(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


# TODO I can download the sequences but it crashes after that...
def batch_start(job, input_args, debug=False):
    """ 
    Downloads shared files that are used by all samples for alignment and calibration and places 
    them in job store
    """
    shared_files = ['ref.fa', 'phase.vcf', 'omni.vcf', 'dbsnp.vcf', 'hapmap.vcf', 'mills.vcf']
    shared_ids = {}
    if debug:
        shared_files = ['debug']
    for file_name in shared_files:
        url = input_args[file_name]
        shared_ids[file_name] = job.addChildJobFn(download_from_url, url, file_name).rv()
    job.addFollowOnJobFn(spawn_batch_jobs, shared_ids, input_args, debug=True)

def spawn_batch_jobs(job, shared_ids, input_args, debug=False):
    # Names for every input file used in the pipeline by each sample
    if not debug:
        samples = []
        config = input_args['config']
        with open(config, 'r') as f:
            for line in f.readlines():
                if not line.isspace():
                    samples.append(line.strip().split(','))
        for sample in samples:
            job.addChildJobFn(start, shared_ids, input_args, sample)
    else:
        job.addChildJobFn(start, shared_ids, input_args, 'sample', debug=True)

def start(job, shared_ids, input_args, sample, debug=False):
    if not debug:
        uuid, url = sample
        ids = shared_ids.copy()
        # Update input
        input_args['uuid'] = uuid
        # Sample bam file holds a url?
        input_args['bam'] = url
        input_args['output_dir'] = os.path.join(input_args['output_dir'], uuid)

        if input_args['ssec'] is None:
            ids['bam'] = job.addChildJobFn(download_from_url, input_args, ids, sample)
        else:
            pass
        job.addChildJobFn(index, ids, input_args, 'sample', cores=1, memory='1 G', disk='8 G')
    else:
        ids = shared_ids.copy()
        ids['bam'] = job.addChildJobFn(index, ids, input_args, sample, debug=True)

def index(job, ids, input_args, sample, debug=False, **kwargs):
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve file path
    bam = '{}.bam'.format(sample)
    if not debug:
        path = return_input_paths(job, work_dir, ids, bam)
        # Call: index the normal.bam
        parameters = ['index', '{}'.format(docker_path(path))]
        docker_call(work_dir, parameters, 'computationalgenomicslab/samtools')
        # Update FileStore and call child
        ids['bai'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, bam) + '.bai')
        job.addChildJobFn(haplotype_caller, ids, input_args, sample, cores=int(input_args['cpu_count']), memory='10 G', disk='15 G')
    else:
        fileid = job.fileStore.readGlobalFile(ids['debug'])
        ids['debug' + '_index'] = job.fileStore.writeGlobalFile(fileid)
        # job.addChildJobFn(haplotype_caller, ids, input_args, sample, cores=None, debug=True)

# TODO There is something wrong with this function ... it crashes every time
def haplotype_caller(job, ids, input_args, sample, debug=False, **kwargs):
   """ snps & indels together """
   work_dir = job.fileStore.getLocalTempDir()
   if not debug:
       pass
       ref_fasta, bam = return_input_paths(job, work_dir, ids, 'ref.fa', 'bam')
       output = os.path.join(work_dir, 'unified.raw.BOTH.gatk.vcf')
       # Call GATK -- HaplotypeCaller
       command = ['-nct', multiprocessing.cpu_count,
                  '-R', docker_path(ref_fasta),
                  '-T', 'HaplotypeCaller',
                  '--genotyping_mode', 'Discovery',
                  '--output_mode', 'EMIT_VARIANTS_ONLY',
                  '-I', docker_path(bam),
                  '-o', docker_path(output),
                  '-stand_emit_conf', '10.0',
                  '-stand_call_conf', '30.0']
       docker_call(work_dir, command, 'computationalgenomicslab/gatk')
      #Update fileStore and spawn child job
      # ids['unified.raw.BOTH.gatk.vcf'] = job.fileStore.writeGlobalFile(output)
      # job.addChildJobFn(vqsr_snp, ids, input_args, sample)
      # job.addChildJobFn(vqsr_indel, ids, input_args, sample)
   else:
       fileid = job.fileStore.readGlobalFile(ids['debug'])
       ids['debug_haplotype_caller'] = job.fileStore.writeGlobalFile(fileid)
       job.addChildJobFn(vqsr_snp, ids, input_args, sample, debug=True)

def vqsr_snp(job, ids, input_args, sample, debug=False):
    if not debug:
        work_dir = job.fileStore.getLocalTempDir()
        ref_fasta, raw_vcf, hapmap, omni, \
        dbsnp, phase = return_input_paths(job, work_dir, ids, 'ref.fa', 'unified.raw.BOTH.gatk.vcf',
                                        'hapmap.vcf', 'omni.vcf', 'dbsnp.vcf', 'phase.vcf')

        output_recalFile = os.path.join(work_dir, 'HAPSNP.recalFile')
        output_tranches = os.path.join(work_dir, 'HAPSNP.tranches')
        output_rscriptFile = os.path.join(work_dir, 'HAPSNP.plots')
        command = ['-T', 'VariantRecalibrator',
                   '-R', ref_fasta,
                   '-input', raw_vcf,
                   '-nt', input_args['cpu_count'],
                   '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', hapmap,
                   '-resource:omni,known=false,training=true,truth=false,prior=12.0', omni,
                   '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0', dbsnp,
                   '-resource:1000G,known=false,training=true,truth=false,prior=10.0', phase,
                   '-an QD', '-an DP', '-an FS', '-an ReadPosRankSum', '-mode SNP', '-minNumBad 1000',
                   '-recalFile', docker_path(output_recalFile),
                   '-tranchesFile', docker_path(output_tranches),
                   '-rscriptFile', docker_path(output_rscriptFile)]
        docker_call(work_dir, command, tool='computationalgenomicslab/gatk')
        ids['HAPSNP.recal'] = job.fileStore.writeGlobalFile(output_recalFile)
        ids['HAPSNP.tranches'] = job.fileStore.writeGlobalFile(output_tranches)
        ids['HAPSNP.plots'] = job.fileStore.writeGlobalFile(output_rscriptFile)
        job.addChildJobFn(apply_vqsr_snp, ids, input_args, sample, cores=int(input_args['cpu_count']), memory='10 G', disk='15 G')
    else:
        fileid = job.fileStore.readGlobalFile(ids['debug_haplotype_caller'])
        ids['debug_vqsr_snp'] = job.fileStore.writeGlobalFile(fileid)
        job.addChildJobFn(apply_vqsr_snp, ids, input_args, debug=True)

#TODO Figure out out to move to output_dir
#Apply Snp Recalibration
def apply_vqsr_snp(job, ids, input_args, debug=False):
    if not debug:
        work_dir = job.fileStore.getLocalTempDir()
        output_dir = input_args['output_dir']
        mkdir_p(output_dir)
        uuid = input_args['uuid']
        output_name = '{}.HAPSNP.vqsr.SNP.vcf'.format(uuid)
        output_path = docker_path(os.path.join(output_dir, output_name))
        ref_fasta, raw_vcf, tranches, recal = return_input_paths(job, work_dir, ids, 'ref.fa', 'unified.raw.BOTH.gatk.vcf',
                                                                  'HAPSNP.tranches', 'HAPSNP.recal')
        command = ['-T', 'ApplyRecalibration',
                   '-input', raw_vcf,
                   '-o', output_path, '-R', ref_fasta,
                   '-nt', multiprocessing.cpu_count(),
                   '-ts_filter_level 99.0',
                   '-tranchesFile', tranches,
                   '-recalFile', recal,
                   '-mode', 'SNP']
        docker_call(work_dir, command, tool='computationalgenomicslab/gatk')
        move_to_output_dir(work_dir, output_dir, uuid=None, files=[os.path.join(work_dir, output_name)])
    else:
        fileid = job.fileStore.readGlobalFile(ids['debug_vqsr_snp'])
        ids['debug_apply_vqsr_snp'] = job.fileStore.writeGlobalFile(fileid)



#Indel Recalibration
def vqsr_indel(job, ids, input_args, sample, debug=False):
    if not debug:
        work_dir = job.fileStore.getLocalTempDir()
        ref_fasta, raw_vcf, mills = return_input_paths(job, work_dir, ids, 'ref.fa', 'unified.raw.BOTH.gatk.vcf',
                                                       'mills.vcf')
        output_recalFile = os.path.join(work_dir, 'HAPINDEL.recalFile')
        output_tranches = os.path.join(work_dir, 'HAPINDEL.tranches')
        output_rscriptFile = os.path.join(work_dir, 'HAPINDEL.plots')
        command = ['-T', 'VariantRecalibrator',
                   '-R', ref_fasta,
                   '-input', raw_vcf,
                   '-nt', multiprocessing.cpu_count(),
                   '-resource:mills,known=true,training=true,truth=true,prior=12.0', mills,
                   '-an DP', '-an FS', '-an ReadPosRankSum', '-mode INDEL',
                   '-minNumBad 1000',
                   '-recalFile', docker_path(output_recalFile),
                   '-tranchesFile', docker_path(output_tranches),
                   '-rscriptFile', docker_path(output_rscriptFile)]
        docker_call(work_dir, command, tool='computationalgenomicslab/gatk')
        ids['HAPINDEL.recal'] = job.fileStore.writeGlobalFile(output_recalFile)
        ids['HAPINDEL.tranches'] = job.fileStore.writeGlobalFile(output_tranches)
        ids['HAPINDEL.plots'] = job.fileStore.writeGlobalFile(output_rscriptFile)
        # job.addChildJobFn(apply_vqsr_indel, ids, input_args, sample, cores=int(input_args['cpu_count']), memory='10 G', disk='15 G')
    else:
        fileid = job.fileStore.readGlobalFile(ids['debug_haplotype_caller'])
        ids['debug_apply_vqsr_indel'] = job.fileStore.writeGlobalFile(fileid)
        job.addChildJobFn(apply_vqsr_indel, ids, input_args, debug=True)

# Apply Indel Recalibration
def apply_vqsr_indel(job, ids, input_args, sample, debug=False, **kwargs):
    if not debug:
        work_dir = job.fileStore.getLocalTempDir()
        output_dir = input_args['output_dir']
        mkdir_p(output_dir)
        uuid = input_args['uuid']
        output_name = '{}.HAPSNP.vqsr.SNP.vcf'.format(uuid)
        output_path = docker_path(os.path.join(work_dir, output_name))
        ref_fasta, raw_vcf, tranches, recal = return_input_paths(job, work_dir, ids,
                                                                 'ref.fa',
                                                                 'unified.raw.BOTH.gatk.vcf',
                                                                 'HAPINDEL.tranches',
                                                                 'HAPINDEL.recal')
        command = ['-T', 'ApplyRecalibration',
                   '-input', ids['{}.unified.raw.BOTH.gatk.vcf'.format(sample)],
                   '-o', output_path,
                   '-R', ref_fasta,
                   '-nt', multiprocessing.cpu_count(),
                   '-ts_filter_level', 99.0,
                   '-tranchesFile', tranches,
                   '-recalFile', recal,
                   '-mode', 'INDEL']
        docker_call(work_dir, command, tool='computationalgenmicslab/gatk')
        move_to_output_dir(work_dir, output_dir, uuid=None, files=[os.path.join(work_dir, output_name)])
    else:
        fileid = job.fileStore.readGlobalFile(ids['debug_vqsr_indel'])
        ids['debug_apply_vqsr_indel'] = job.fileStore.writeGlobalFile(fileid)



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
                  'cpu_count': None}

    debug_args = {'debug': args.hapmap}
    
    Job.Runner.startToil(Job.wrapJobFn(batch_start, debug_args, debug=True), args)
