#!/usr/bin/env python2.7
"""
@author Frank Austin Nothaft fnothaft@berkeley.edu
@date 12/30/2015

Pipeline to go from FASTQ to VCF using both the ADAM+HaplotypeCaller pipeline
as well as the GATK best practices pipeline.

This doesn't contain the GATK best practices pipeline yet though. More work TBD.

        0 --> 1 --> 2 --> 3 --> 4 --> 5
                                      |++(6)
                                      7 --> 9 --> 10 --> 11 --> 12 <snip1>
                                       ++(8)

        <snip1 /> 12 --> 13 --> 14 --> 15 --> 16 --> 17
                                                     |
                                                     18
                                                    /  \
                                                  19    20
                                                 /        \
                                               21          22

0   bwa alignment to a reference
1   samtools sam to bam conversion (no sort)
2   Fix header
3   Add read groups
4   Upload to S3
5   Start master
6   Master Service
7   Start Workers
8   Worker Service
9   Download Data
10  ADAM Convert
11  ADAM Transform
12  Upload Data
13  Start GATK box
14  Download reference
15  Index reference
16  Build reference dictionary
17  Index samples
18  Run HaplotypeCaller
19  Run VQSR on SNPs
20  Run VQSR on INDELs
21  Apply VQSR model to SNPs
22  Apply VQSR model to INDELs

However, the pipeline in this file is actually just three encapsulated jobs:

        A --> B --> C

A  Run BWA (jobs 0-4)
B  Run ADAM (jobs 5-12)
C  Run GATK (jobs 13-22)

Those hypens should be en dashes. Alas, TIL (2/1/16) that Python does not like
non-ASCII text in comments.

===================================================================
:Dependencies:
curl            - apt-get install curl
Toil            - pip install --pre toil
Docker          - http://docs.docker.com/engine/installation/

Optional:
S3AM            - pip install --s3am (requires ~/.boto config file)
"""

# import from python system libraries
import argparse
import multiprocessing
import os

# import toil features
from toil.job import Job

# import job steps from other toil pipelines
from toil_scripts.adam_pipeline.spark_toil_script import *
from toil_scripts.batch_alignment.bwa_alignment import *
from toil_scripts.gatk_germline.germline import *

def build_parser():

    parser = argparse.ArgumentParser()

    # add sample uuid
    parser.add_argument('-U', '--uuid', required = True,
                        help = 'Sample UUID.')

    # add bucket args
    parser.add_argument('-3', '--s3_bucket', required = True,
                        help = 'S3 Bucket URI')
    parser.add_argument('-3r', '--bucket_region', default = "us-east-1",
                        help = 'Region of the S3 bucket. Defaults to us-east-1.')
    parser.add_argument('-y', '--aws_access_key', required = True,
                        help = 'Amazon web services access key')
    parser.add_argument('-S', '--aws_secret_key', required = True,
                        help = 'Amazon web services secret key')

    # add file size argument
    parser.add_argument('-se', '--file_size', default = '100G',
                        help = 'Approximate input file size. Should be given as %d[TGMK], e.g., for a 100 gigabyte file, use --file_size 100G')

    # add bwa args
    parser.add_argument('-r', '--ref', required = True,
                        help = 'Reference fasta file')
    parser.add_argument('-m', '--amb', required = True,
                        help = 'Reference fasta file (amb)')
    parser.add_argument('-n', '--ann', required = True,
                        help = 'Reference fasta file (ann)')
    parser.add_argument('-b', '--bwt', required = True,
                        help = 'Reference fasta file (bwt)')
    parser.add_argument('-p', '--pac', required = True,
                        help = 'Reference fasta file (pac)')
    parser.add_argument('-a', '--sa', required = True,
                        help = 'Reference fasta file (sa)')
    parser.add_argument('-f', '--fai', required = True,
                        help = 'Reference fasta file (fai)')
    parser.add_argument('-u', '--sudo', dest = 'sudo', action = 'store_true',
                        help = 'Docker usually needs sudo to execute '
                        'locally, but not''when running Mesos '
                        'or when a member of a Docker group.')
    parser.add_argument('-k', '--use_bwakit', action='store_true', help='Use bwakit instead of the binary build of bwa')
    parser.add_argument('-t', '--alt', required=False, help='Alternate file for reference build (alt). Necessary for alt aware alignment.')
    parser.set_defaults(sudo = False)

    # add ADAM args
    parser.add_argument('-N', '--num_nodes', type = int, required = True,
                        help = 'Number of nodes to use')
    parser.add_argument('-d', '--driver_memory', required = True,
                        help = 'Amount of memory to allocate for Spark Driver.')
    parser.add_argument('-q', '--executor_memory', required = True,
                        help = 'Amount of memory to allocate per Spark Executor.')

    # add GATK args
    parser.add_argument('-P', '--phase', required = True,
                        help = '1000G_phase1.indels.b37.vcf URL')
    parser.add_argument('-M', '--mills', required = True,
                        help = 'Mills_and_1000G_gold_standard.indels.b37.vcf URL')
    parser.add_argument('-s', '--dbsnp', required = True,
                        help = 'dbsnp_137.b37.vcf URL')
    parser.add_argument('-O', '--omni', required = True,
                        help = '1000G_omni.5.b37.vcf URL')
    parser.add_argument('-H', '--hapmap', required = True,
                        help = 'hapmap_3.3.b37.vcf URL')
    
    # return built parser
    return parser

def static_dag(job, s3_bucket, bucket_region, uuid, bwa_inputs, adam_inputs, gatk_inputs):
    """
    Prefer this here as it allows us to pull the job functions from other jobs
    without rewrapping the job functions back together.

    bwa_inputs: Input arguments to be passed to BWA.
    adam_inputs: Input arguments to be passed to ADAM.
    gatk_inputs: Input arguments to be passed to the GATK.
    """

    # get work directory
    work_dir = job.fileStore.getLocalTempDir()

    # what region is our bucket in?
    if bucket_region == "us-east-1":
        bucket_region = ""
    else:
        bucket_region = "-%s" % bucket_region

    # does the work directory exist?
    if not os.path.exists(work_dir):
        os.mkdirs(work_dir)

    # write config for bwa
    bwa_config_path = os.path.join(work_dir, "%s_bwa_config.csv" % uuid)
    bwafp = open(bwa_config_path, "w")
    print >> bwafp, "%s,https://s3%s.amazonaws.com/%s/sequence/%s_1.fastq.gz,https://s3%s.amazonaws.com/%s/sequence/%s_2.fastq.gz" % (uuid, bucket_region, s3_bucket, uuid, bucket_region, s3_bucket, uuid)
    bwafp.flush()
    bwafp.close()
    bwa_inputs['config'] = job.fileStore.writeGlobalFile(bwa_config_path)

    # write config for gatk
    gatk_config_path = os.path.join(work_dir, "%s_gatk_config.csv" % uuid)
    gatkfp = open(gatk_config_path, "w")
    print >> gatkfp, "%s,https://s3%s.amazonaws.com/%s/analysis/%s.bam" % (uuid, bucket_region, s3_bucket, uuid)
    gatkfp.flush()
    gatkfp.close()
    gatk_inputs['config'] = job.fileStore.writeGlobalFile(gatk_config_path)
    
    # get head bwa job function and encapsulate it
    bwa = job.wrapJobFn(download_shared_files,
                        bwa_inputs).encapsulate()

    # get head ADAM job function and encapsulate it
    adam = job.wrapJobFn(start_master,
                         adam_inputs).encapsulate()
    
    # get head GATK job function and encapsulate it
    gatk = job.wrapJobFn(batch_start,
                         gatk_inputs).encapsulate()

    # wire up dag
    job.addChild(bwa)
    bwa.addChild(adam)
    adam.addChild(gatk)

if __name__ == '__main__':
    
    args_parser = build_parser()
    Job.Runner.addToilOptions(args_parser)
    args = args_parser.parse_args()

    bwa_inputs = {'ref.fa': args.ref,
                  'ref.fa.amb': args.amb,
                  'ref.fa.ann': args.ann,
                  'ref.fa.bwt': args.bwt,
                  'ref.fa.pac': args.pac,
                  'ref.fa.sa': args.sa,
                  'ref.fa.fai': args.fai,
                  'ref.fa.alt': args.alt,
                  'ssec': None,
                  'output_dir': None,
                  'sudo': args.sudo,
                  's3_dir': "%s/alignment" % args.s3_bucket,
                  'lb': args.uuid,
                  'uuid': args.uuid,
                  'cpu_count': None,
                  'file_size': args.file_size,
                  'use_bwakit': args.use_bwakit,
                  'aws_access_key':  args.aws_access_key,
                  'aws_secret_key':  args.aws_secret_key}
    
    adam_inputs = {'numWorkers': args.num_nodes - 1,
                   'outDir':     's3://%s/analysis/%s.bam' % (args.s3_bucket, args.uuid),
                   'knownSNPs':  args.dbsnp.replace("https://s3-us-west-2.amazonaws.com/", "s3://"),
                   'accessKey':  args.aws_access_key,
                   'secretKey':  args.aws_secret_key,
                   'driverMemory': args.driver_memory,
                   'executorMemory': args.executor_memory,
                   'sudo': args.sudo,
                   'bamName': 's3://%s/alignment/%s.bam' % (args.s3_bucket, args.uuid)}

    gatk_inputs = {'ref.fa': args.ref,
                   'phase.vcf': args.phase,
                   'mills.vcf': args.mills,
                   'dbsnp.vcf': args.dbsnp,
                   'hapmap.vcf': args.hapmap,
                   'omni.vcf': args.omni,
                   'output_dir': None,
                   'uuid': None,
                   'cpu_count': str(multiprocessing.cpu_count()),
                   'ssec': None,
                   's3_dir': "%s/%s/analysis" % (args.s3_bucket, args.uuid),
                   'file_size': args.file_size,
                   'aws_access_key':  args.aws_access_key,
                   'aws_secret_key':  args.aws_secret_key}

    Job.Runner.startToil(Job.wrapJobFn(static_dag,
                                       args.s3_bucket,
                                       args.bucket_region,
                                       args.uuid,
                                       bwa_inputs,
                                       adam_inputs,
                                       gatk_inputs), args)
