#!/usr/bin/env python2.7
"""
@author Frank Austin Nothaft fnothaft@berkeley.edu
@date 12/30/2015

Pipeline to go from FASTQ to VCF using both the ADAM+HaplotypeCaller pipeline
as well as the GATK best practices pipeline.

  0 --> ... --> 4 --> 5
                |     |++(6)
                |     7 --> 9 --> ... --> 12 --> 13 --> ... --> 17
                |     ++(8)                                     |
                |                                               18
                |                                              /  \
                |                                            19    20
                |                                           /        \
                |                                         21          22
                |
                |
                + --> 23 --> ... --> 34 --> 35 --> ... --> 39
                                                           |
                                                           40
                                                          /  \
                                                        41    42
                                                       /        \
                                                      43         44


BWA alignment

0   bwa alignment to a reference
1   samtools sam to bam conversion (no sort)
2   Fix header
3   Add read groups
4   Upload to S3

ADAM preprocessing

5   Start master
6   Master Service
7   Start Workers
8   Worker Service
9   Download Data
10  ADAM Convert
11  ADAM Transform
12  Upload Data

GATK haplotype caller

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

GATK preprocessing

23  Download shared data
24  Reference preprocessing
25  Download sample 
26  Index
27  Sort
28  Mark duplicates
29  Index
30  Realigner target 
31  Indel realignment
32  Index
33  Base recalibration
34  Output BQSR file

GATK haplotype caller

35  Start GATK box
36  Download reference
37  Index reference
38  Build reference dictionary
39  Index samples
40  Run HaplotypeCaller
41  Run VQSR on SNPs
42  Run VQSR on INDELs
43  Apply VQSR model to SNPs
44  Apply VQSR model to INDELs


However, the pipeline in this file is actually just five encapsulated jobs:

        A
       / \
      B   D
      |   |
      C   E

A  Run BWA alignment (jobs 0-4)
B  Run ADAM preprocessing (jobs 5-12)
C  Run GATK haplotype caller (jobs 13-22)
D  Run GATK preprocessing (jobs 23-34)
E  Run GATK haplotype caller (jobs 35-44)

===================================================================
:Dependencies:
curl            - apt-get install curl
Toil            - pip install --pre toil
Docker          - http://docs.docker.com/engine/installation/

Optional:
S3AM            - pip install --s3am (requires ~/.boto config file)
"""

# import from python system libraries
import copy

# import boto sdb connection for autoscaling
import boto.sdb

# import toil features
from toil.job import Job

# import job steps from other toil pipelines
from toil_scripts.adam_pipeline.adam_preprocessing import *
from toil_scripts.batch_alignment.bwa_alignment import *
from toil_scripts.gatk_germline.germline import *
from toil_scripts.gatk_processing.gatk_preprocessing import *

# these don't seem necessary! but, must be imported here due to a serialization issue
from toil_scripts.spark_utils.spawn_cluster import *

# import autoscaling tools
from toil_scripts.adam_uberscript.automated_scaling import Samples

def build_parser():

    parser = argparse.ArgumentParser()

    # add sample uuid
    parser.add_argument('-U', '--uuid_manifest', required = False,
                        help = 'Sample UUID.')

    # optionally add suffix
    parser.add_argument('--dir_suffix',
                        help = 'Optional suffix to add to output directory names.',
                        default = '')

    # what pipeline are we running
    parser.add_argument('-PR', '--pipeline_to_run',
                        help = "Whether we should run 'adam', 'gatk', or 'both'. Default is 'both'.",
                        default = 'both')
    parser.add_argument('-SA', '--skip_alignment',
                        help = "Skip alignment and start running from preprocessing.",
                        action = 'store_true')
    parser.add_argument('-SP', '--skip_preprocessing',
                        help = "Skip preprocessing and start running from variant calling. Implies --skip_alignment.",
                        action = 'store_true')

    # what directory are our sequence files in?
    parser.add_argument('-SD', '--sequence_dir',
                        help = 'Directory where raw sequences are.',
                        default = 'sequence')

    # are we automatically scaling the cluster?
    parser.add_argument('--autoscale_cluster', action='store_true', help = "Scales cluster during pipeline run")

    # add bucket args
    parser.add_argument('-3', '--s3_bucket', required = not mock_mode(),
                        help = 'S3 Bucket URI')
    parser.add_argument('-3r', '--bucket_region', default = "us-west-2",
                        help = 'Region of the S3 bucket. Defaults to us-east-1.')

    # add file size argument
    parser.add_argument('-se', '--file_size', default = '100G',
                        help = 'Approximate input file size. Should be given as %d[TGMK], e.g., '
                               'for a 100 gigabyte file, use --file_size 100G')

    # add bwa args
    parser.add_argument('-r', '--ref', required = not mock_mode(),
                        help = 'Reference fasta file')
    parser.add_argument('-m', '--amb', required = not mock_mode(),
                        help = 'Reference fasta file (amb)')
    parser.add_argument('-n', '--ann', required = not mock_mode(),
                        help = 'Reference fasta file (ann)')
    parser.add_argument('-b', '--bwt', required = not mock_mode(),
                        help = 'Reference fasta file (bwt)')
    parser.add_argument('-p', '--pac', required = not mock_mode(),
                        help = 'Reference fasta file (pac)')
    parser.add_argument('-a', '--sa', required = not mock_mode(),
                        help = 'Reference fasta file (sa)')
    parser.add_argument('-f', '--fai', required = not mock_mode(),
                        help = 'Reference fasta file (fai)')
    parser.add_argument('-u', '--sudo', dest = 'sudo', action = 'store_true',
                        help = 'Docker usually needs sudo to execute '
                        'locally, but not''when running Mesos '
                        'or when a member of a Docker group.')
    parser.add_argument('-k', '--use_bwakit', action='store_true', help='Use bwakit instead of the binary build of bwa')
    parser.add_argument('-t', '--alt', required=False, help='Alternate file for reference build (alt). Necessary for alt aware alignment.')
    parser.set_defaults(sudo = False) 
    parser.add_argument('--trim', action='store_true', help='Trim adapters during alignment.')

    # add ADAM args
    parser.add_argument('-N', '--num_nodes', type = int, required = False, default = None,
                        help = 'Number of Toil nodes to allocate per Spark subcluster. Exclusive of --master_ip.')
    parser.add_argument('-MI', '--master_ip', required = False, default = None,
                        help ="IP or hostname of host running for Spark master and HDFS namenode. Should be provided "
                              "if pointing at a static (external or standalone) Spark cluster. The special value "
                              "'auto' indicates the master of standalone cluster, i.e. one that is managed by the "
                              "uberscript.")
    parser.add_argument('-d', '--driver_memory', required = not mock_mode(),
                        help = 'Amount of memory to allocate for Spark Driver.')
    parser.add_argument('-q', '--executor_memory', required = not mock_mode(),
                        help = 'Amount of memory to allocate per Spark Executor.')

    # add GATK args
    parser.add_argument('-P', '--phase', required = not mock_mode(),
                        help = '1000G_phase1.indels.b37.vcf URL')
    parser.add_argument('-M', '--mills', required = not mock_mode(),
                        help = 'Mills_and_1000G_gold_standard.indels.b37.vcf URL')
    parser.add_argument('-s', '--dbsnp', required = not mock_mode(),
                        help = 'dbsnp_137.b37.vcf URL')
    parser.add_argument('-O', '--omni', required = not mock_mode(),
                        help = '1000G_omni.5.b37.vcf URL')
    parser.add_argument('-H', '--hapmap', required = not mock_mode(),
                        help = 'hapmap_3.3.b37.vcf URL')

    # return built parser
    return parser


def sample_loop(job,
                bucket_region,
                s3_bucket,
                uuid_list,
                bwa_inputs,
                adam_inputs,
                gatk_preprocess_inputs,
                gatk_adam_call_inputs,
                gatk_gatk_call_inputs,
                pipeline_to_run,
                skip_alignment,
                skip_preprocessing,
                autoscale_cluster,
                sequence_dir,
                dir_suffix):
  """
  Loops over the sample_ids (uuids) in the manifest, creating child jobs to process each
  """

  for uuid_rg in uuid_list:

    uuid_items = uuid_rg.split(',')
    uuid = uuid_items[0]
    rg_line = None
    if len(uuid_items) > 1:
        rg_line = uuid_items[1]

    uuid_bwa_inputs = copy.deepcopy(bwa_inputs)
    uuid_adam_inputs = copy.deepcopy(adam_inputs)
    uuid_gatk_preprocess_inputs = copy.deepcopy(gatk_preprocess_inputs)
    uuid_gatk_adam_call_inputs = copy.deepcopy(gatk_adam_call_inputs)
    uuid_gatk_gatk_call_inputs = copy.deepcopy(gatk_gatk_call_inputs)

    ## set uuid inputs
    uuid_bwa_inputs['lb'] = uuid
    uuid_bwa_inputs['uuid'] = uuid
    uuid_bwa_inputs['rg_line'] = rg_line
    uuid_adam_inputs['outDir'] = 's3://{s3_bucket}/analysis{dir_suffix}/{uuid}'.format(**locals())
    uuid_adam_inputs['bamName'] = 's3://{s3_bucket}/alignment{dir_suffix}/{uuid}.bam'.format(**locals())
    uuid_gatk_preprocess_inputs['s3_dir'] = '{s3_bucket}/analysis{dir_suffix}/{uuid}'.format(**locals())
    uuid_gatk_adam_call_inputs['s3_dir'] = '{s3_bucket}/analysis{dir_suffix}/{uuid}'.format(**locals())
    uuid_gatk_gatk_call_inputs['s3_dir'] = '{s3_bucket}/analysis{dir_suffix}/{uuid}'.format(**locals())

    # are we autoscaling? if so, bump the cluster size by one now
    if autoscale_cluster:
        Samples.increase_nodes(uuid, 1)

    job.addChildJobFn(static_dag,
                      bucket_region,
                      s3_bucket,
                      uuid,
                      uuid_bwa_inputs,
                      uuid_adam_inputs,
                      uuid_gatk_preprocess_inputs,
                      uuid_gatk_adam_call_inputs,
                      uuid_gatk_gatk_call_inputs,
                      pipeline_to_run,
                      skip_alignment,
                      skip_preprocessing,
                      autoscale_cluster,
                      sequence_dir,
                      dir_suffix)


def static_dag(job,
               bucket_region,
               s3_bucket,
               uuid,
               bwa_inputs,
               adam_inputs,
               gatk_preprocess_inputs,
               gatk_adam_call_inputs,
               gatk_gatk_call_inputs,
               pipeline_to_run,
               skip_alignment,
               skip_preprocessing,
               autoscale_cluster,
               sequence_dir,
               dir_suffix):
    """
    Prefer this here as it allows us to pull the job functions from other jobs
    without rewrapping the job functions back together.

    bwa_inputs: Input arguments to be passed to BWA.
    adam_inputs: Input arguments to be passed to ADAM.
    gatk_preprocess_inputs: Input arguments to be passed to GATK preprocessing.
    gatk_adam_call_inputs: Input arguments to be passed to GATK haplotype caller for the result of ADAM preprocessing.
    gatk_gatk_call_inputs: Input arguments to be passed to GATK haplotype caller for the result of GATK preprocessing.
    """

    # get work directory
    work_dir = job.fileStore.getLocalTempDir()

    # what region is our bucket in?
    if bucket_region == "us-east-1":
        bucket_region = ""
    else:
        bucket_region = "-%s" % bucket_region

    # write config for bwa
    bwa_config_path = os.path.join(work_dir, '{uuid}_bwa_config.csv'.format(**locals()))
    bwafp = open(bwa_config_path, "w")
    print >> bwafp, ('{uuid}'
                     ',s3://{s3_bucket}/{sequence_dir}/{uuid}_1.fastq.gz'
                     ',s3://{s3_bucket}/{sequence_dir}/{uuid}_2.fastq.gz'.format(**locals()))
    bwafp.flush()
    bwafp.close()
    bwa_inputs['config'] = job.fileStore.writeGlobalFile(bwa_config_path)

    # write config for GATK preprocessing
    gatk_preprocess_config_path = os.path.join(work_dir, '{uuid}_gatk_preprocess_config.csv'.format(**locals()))
    gatk_preprocess_fp = open(gatk_preprocess_config_path, "w")
    print >> gatk_preprocess_fp, '{uuid},s3://{s3_bucket}/alignment{dir_suffix}/{uuid}.bam'.format(**locals())
    gatk_preprocess_fp.flush()
    gatk_preprocess_fp.close()
    gatk_preprocess_inputs['cpu_count'] = multiprocessing.cpu_count()
    gatk_preprocess_inputs['config'] = job.fileStore.writeGlobalFile(gatk_preprocess_config_path)

    # write config for GATK haplotype caller for the result of ADAM preprocessing
    gatk_adam_call_config_path = os.path.join(work_dir, '{uuid}_gatk_adam_call_config.csv'.format(**locals()))
    gatk_adam_call_fp = open(gatk_adam_call_config_path, "w")
    print >> gatk_adam_call_fp, '{uuid},s3://{s3_bucket}/analysis{dir_suffix}/{uuid}/{uuid}.adam.bam'.format(**locals())
    gatk_adam_call_fp.flush()
    gatk_adam_call_fp.close()
    gatk_adam_call_inputs['cpu_count'] = multiprocessing.cpu_count()
    gatk_adam_call_inputs['config'] = job.fileStore.writeGlobalFile(gatk_adam_call_config_path)

    # write config for GATK haplotype caller for the result of GATK preprocessing
    gatk_gatk_call_config_path = os.path.join(work_dir, '{uuid}_gatk_gatk_call_config.csv'.format(**locals()))
    gatk_gatk_call_fp = open(gatk_gatk_call_config_path, "w")
    print >> gatk_gatk_call_fp, '{uuid},s3://{s3_bucket}/analysis{dir_suffix}/{uuid}/{uuid}.gatk.bam'.format(**locals())
    gatk_gatk_call_fp.flush()
    gatk_gatk_call_fp.close()
    gatk_gatk_call_inputs['cpu_count'] = multiprocessing.cpu_count()
    gatk_gatk_call_inputs['config'] = job.fileStore.writeGlobalFile(gatk_gatk_call_config_path)

    # get head BWA alignment job function and encapsulate it
    bwa = job.wrapJobFn(download_shared_files,
                        bwa_inputs).encapsulate()

    # get head ADAM preprocessing job function and encapsulate it
    adam_preprocess = job.wrapJobFn(static_adam_preprocessing_dag,
                                    adam_inputs).encapsulate()

    # get head GATK preprocessing job function and encapsulate it
    gatk_preprocess = job.wrapJobFn(download_gatk_files,
                                    gatk_preprocess_inputs).encapsulate()

    # get head GATK haplotype caller job function for the result of ADAM preprocessing and encapsulate it
    gatk_adam_call = job.wrapJobFn(batch_start,
                                   gatk_adam_call_inputs).encapsulate()

    # get head GATK haplotype caller job function for the result of GATK preprocessing and encapsulate it
    gatk_gatk_call = job.wrapJobFn(batch_start,
                                   gatk_gatk_call_inputs).encapsulate()

    # add code to bump the number of jobs after alignment
    # start with -1 since we already have a single node for our sample
    nodes_needed_after_alignment = -1

    # adam needs:
    # - a spark driver
    # - a spark master/hdfs namenode
    # - _n_ spark workers/hdfs datanodes
    if (pipeline_to_run == "adam" or
        pipeline_to_run == "both"):
        if adam_inputs['masterIP']:
            nodes_needed_after_alignment += 1
        else:
            nodes_needed_after_alignment += (adam_inputs['numWorkers'] + 2)

    # gatk needs one node
    if (pipeline_to_run == "gatk" or
        pipeline_to_run == "both"):
        nodes_needed_after_alignment += 1

    # create a job that runs after alignment and increases the number of nodes in the system
    increase_nodes_after_alignment = job.wrapJobFn(increase_node_count,
                                                   nodes_needed_after_alignment - 1,
                                                   uuid)

    # since we evaluate this conditional repeatedly, just calculate it once
    # we should only schedule the job that increases the node count if we
    # are using autoscaling _and_ we are increasing the node count
    autoscale_after_alignment = (nodes_needed_after_alignment > 1) and autoscale_cluster

    # create a job that runs after ADAM's preprocessing to decrease the number of nodes
    decrease_nodes_after_adam_preprocess = job.wrapJobFn(decrease_node_count,
                                                         adam_inputs['numWorkers'] + 1,
                                                         uuid)

    # wire up dag
    if not skip_alignment:
        job.addChild(bwa)

        if autoscale_after_alignment:
            bwa.addChild(increase_nodes_after_alignment)
    elif autoscale_after_alignment:
        job.addChild(increase_nodes_after_alignment)

    if (pipeline_to_run == "adam" or
        pipeline_to_run == "both"):

        if skip_preprocessing:

            # TODO: currently, we don't support autoscaling if you're just
            # running variant calling
            job.addChild(gatk_adam_call)
        else:
            if autoscale_after_alignment:
                increase_nodes_after_alignment.addChild(adam_preprocess)
            elif skip_alignment:
                job.addChild(adam_preprocess)
            else:
                bwa.addChild(adam_preprocess)

            # if we are running autoscaling, then we should decrease the cluster
            # size after adam completes
            if autoscale_cluster and not adam_inputs['masterIP']:
                adam_preprocess.addChild(decrease_nodes_after_adam_preprocess)
                decrease_nodes_after_adam_preprocess.addChild(gatk_adam_call)
            else:
                adam_preprocess.addChild(gatk_adam_call)

    if (pipeline_to_run == "gatk" or
        pipeline_to_run == "both"):

        if skip_preprocessing:

            # TODO: currently, we don't support autoscaling if you're just
            # running variant calling
            job.addChild(gatk_gatk_call)
        else:
            if autoscale_after_alignment:
                increase_nodes_after_alignment.addChild(gatk_preprocess)
            elif skip_alignment:
                job.addChild(gatk_preprocess)
            else:
                bwa.addChild(gatk_preprocess)

            gatk_preprocess.addChild(gatk_gatk_call)


def increase_node_count(job, nodes_to_add, uuid):

    Samples.increase_nodes(uuid, nodes_to_add)


def decrease_node_count(job, nodes_to_add, uuid):

    Samples.decrease_nodes(uuid, nodes_to_add)


if __name__ == '__main__':

    args_parser = build_parser()
    Job.Runner.addToilOptions(args_parser)
    args = args_parser.parse_args()

    if mock_mode():
        from mock_inputs import mock_inputs
        for k,v in mock_inputs.iteritems():
            setattr(args, k, v)

    ## Parse manifest file
    uuid_list = []
    with open(args.uuid_manifest) as f_manifest:
      for uuid in f_manifest:
        uuid_list.append(uuid.strip())

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
                  's3_dir': "%s/alignment%s" % (args.s3_bucket, args.dir_suffix),
                  'cpu_count': None,
                  'file_size': args.file_size,
                  'use_bwakit': args.use_bwakit,
                  'sort': False,
                  'trim': args.trim}

    if bool(args.master_ip) == bool(args.num_nodes):
        raise ValueError("Exactly one of --master_ip (%s) and --num_nodes (%d) must be provided." %
                         (args.master_ip, args.num_nodes))

    if args.num_nodes <= 1 and not args.master_ip:
        raise ValueError("--num_nodes allocates one Spark/HDFS master and n-1 workers, and thus must be greater than 1. %d was passed." %
                         args.num_nodes)

    adam_inputs = {'numWorkers': args.num_nodes - 1,
                   'knownSNPs':  args.dbsnp.replace("https://s3-us-west-2.amazonaws.com/", "s3://"),
                   'driverMemory': args.driver_memory,
                   'executorMemory': args.executor_memory,
                   'sudo': args.sudo,
                   'suffix': '.adam',
                   'masterIP': args.master_ip}

    gatk_preprocess_inputs = {'ref.fa': args.ref,
                              'phase.vcf': args.phase,
                              'mills.vcf': args.mills,
                              'dbsnp.vcf': args.dbsnp,
                              'output_dir': None,
                              'sudo': args.sudo,
                              'ssec': None,
                              'cpu_count': None,
                              'suffix': '.gatk',
                              'memory': args.executor_memory}

    gatk_adam_call_inputs = {'ref.fa': args.ref,
                             'phase.vcf': args.phase,
                             'mills.vcf': args.mills,
                             'dbsnp.vcf': args.dbsnp,
                             'hapmap.vcf': args.hapmap,
                             'omni.vcf': args.omni,
                             'output_dir': None,
                             'uuid': None,
                             'cpu_count': None,
                             'ssec': None,
                             'file_size': args.file_size,
                             'suffix': '.adam',
                             'indexed': False,
                             'sudo': args.sudo,
                             'memory': args.executor_memory}

    gatk_gatk_call_inputs = {'ref.fa': args.ref,
                             'phase.vcf': args.phase,
                             'mills.vcf': args.mills,
                             'dbsnp.vcf': args.dbsnp,
                             'hapmap.vcf': args.hapmap,
                             'omni.vcf': args.omni,
                             'output_dir': None,
                             'uuid': None,
                             'cpu_count': None,
                             'ssec': None,
                             'file_size': args.file_size,
                             'suffix': '.gatk',
                             'indexed': True,
                             'sudo': args.sudo,
                             'memory': args.executor_memory,
                             'sudo': args.sudo}

    if (args.pipeline_to_run != "adam" and
        args.pipeline_to_run != "gatk" and
        args.pipeline_to_run != "both"):
        raise ValueError("--pipeline_to_run must be either 'adam', 'gatk', or 'both'. %s was passed." % args.pipeline_to_run)

    Job.Runner.startToil(Job.wrapJobFn(sample_loop,
                                       args.bucket_region,
                                       args.s3_bucket,
                                       uuid_list,
                                       bwa_inputs,
                                       adam_inputs,
                                       gatk_preprocess_inputs,
                                       gatk_adam_call_inputs,
                                       gatk_gatk_call_inputs,
                                       args.pipeline_to_run,
                                       args.skip_alignment,
                                       args.skip_preprocessing,
                                       args.autoscale_cluster,
                                       args.sequence_dir,
                                       args.dir_suffix), args)
