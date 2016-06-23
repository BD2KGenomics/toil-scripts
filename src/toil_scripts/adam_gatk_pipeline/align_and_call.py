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
import argparse
import copy
import textwrap
from multiprocessing import cpu_count

import yaml
# import toil features
from toil.job import Job
# these don't seem necessary! but, must be imported here due to a serialization issue
from toil.lib.spark import spawn_spark_cluster

# import job steps from other toil pipelines
from toil_scripts.adam_pipeline.adam_preprocessing import * #static_adam_preprocessing_dag
from toil_scripts.bwa_alignment.bwa_alignment import * #download_shared_files
from toil_scripts.gatk_germline.germline import * #batch_start
from toil_scripts.gatk_processing.gatk_preprocessing import * #download_gatk_files
from toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline import generate_file

from toil_scripts.lib.programs import mock_mode

# import autoscaling tools
from toil_scripts.adam_uberscript.automated_scaling import Samples

def sample_loop(job, uuid_list, inputs):
  """
  Loops over the sample_ids (uuids) in the manifest, creating child jobs to process each
  """

  for uuid_rg in uuid_list:

    uuid_items = uuid_rg.split(',')
    uuid = uuid_items[0]
    rg_line = None
    if len(uuid_items) > 1:
        rg_line = uuid_items[1]

    # are we autoscaling? if so, bump the cluster size by one now
    if inputs.autoscale_cluster:
        Samples.increase_nodes(uuid, 1)

    job.addChildJobFn(static_dag, uuid, rg_line, inputs)


def static_dag(job, uuid, rg_line, inputs):
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

    inputs.cpu_count = cpu_count()
    inputs.maxCores = sys.maxint
    args = {'uuid': uuid,
            's3_bucket': inputs.s3_bucket,
            'sequence_dir': inputs.sequence_dir,
            'dir_suffix': inputs.dir_suffix}

    # get head BWA alignment job function and encapsulate it
    inputs.rg_line = rg_line
    inputs.output_dir = 's3://{s3_bucket}/alignment{dir_suffix}'.format(**args)
    bwa = job.wrapJobFn(download_reference_files,
                        inputs,
                        [[uuid,
                         ['s3://{s3_bucket}/{sequence_dir}/{uuid}_1.fastq.gz'.format(**args),
                          's3://{s3_bucket}/{sequence_dir}/{uuid}_2.fastq.gz'.format(**args)]]]).encapsulate()

    # get head ADAM preprocessing job function and encapsulate it
    adam_preprocess = job.wrapJobFn(static_adam_preprocessing_dag,
                                    inputs,
                                    's3://{s3_bucket}/alignment{dir_suffix}/{uuid}.bam'.format(**args),
                                    's3://{s3_bucket}/analysis{dir_suffix}/{uuid}'.format(**args),
                                    suffix='.adam').encapsulate()

    # get head GATK preprocessing job function and encapsulate it
    gatk_preprocess = job.wrapJobFn(download_gatk_files,
                                    inputs,
                                    [uuid,'s3://{s3_bucket}/alignment{dir_suffix}/{uuid}.bam'.format(**args)],
                                    's3://{s3_bucket}/analysis{dir_suffix}/{uuid}'.format(**args),
                                    suffix='.gatk').encapsulate()

    adam_call_inputs = inputs
    gatk_call_inputs = copy.deepcopy(inputs)
    adam_call_inputs.indexed = False
    gatk_call_inputs.indexed = True

    # get head GATK haplotype caller job function for the result of ADAM preprocessing and encapsulate it
    gatk_adam_call = job.wrapJobFn(batch_start,
                                   adam_call_inputs,
                                   [uuid,'s3://{s3_bucket}/analysis{dir_suffix}/{uuid}/{uuid}.adam.bam'.format(**args)],
                                   's3://{s3_bucket}/analysis{dir_suffix}/{uuid}'.format(**args),
                                   suffix='.adam').encapsulate()

    # get head GATK haplotype caller job function for the result of GATK preprocessing and encapsulate it
    gatk_gatk_call = job.wrapJobFn(batch_start,
                                   gatk_call_inputs,
                                   [uuid,'s3://{s3_bucket}/analysis{dir_suffix}/{uuid}/{uuid}.gatk.bam'.format(**args)],
                                   's3://{s3_bucket}/analysis{dir_suffix}/{uuid}'.format(**args),
                                   suffix='.gatk').encapsulate()

    # add code to bump the number of jobs after alignment
    # start with -1 since we already have a single node for our sample
    nodes_needed_after_alignment = -1

    # adam needs:
    # - a spark driver
    # - a spark master/hdfs namenode
    # - _n_ spark workers/hdfs datanodes
    if (inputs.pipeline_to_run == "adam" or
        inputs.pipeline_to_run == "both"):
        if inputs.master_ip:
            nodes_needed_after_alignment += 1
        else:
            nodes_needed_after_alignment += (inputs.num_nodes + 1)

    # gatk needs one node
    if (inputs.pipeline_to_run == "gatk" or
        inputs.pipeline_to_run == "both"):
        nodes_needed_after_alignment += 1

    # since we evaluate this conditional repeatedly, just calculate it once
    # we should only schedule the job that increases the node count if we
    # are using autoscaling _and_ we are increasing the node count
    autoscale_after_alignment = (nodes_needed_after_alignment > 1) and inputs.autoscale_cluster

    if autoscale_after_alignment:
        # create a job that runs after alignment and increases the number of nodes in the system
        increase_nodes_after_alignment = job.wrapJobFn(increase_node_count,
                                                       nodes_needed_after_alignment - 1,
                                                       uuid)

        # create a job that runs after ADAM's preprocessing to decrease the number of nodes
        decrease_nodes_after_adam_preprocess = job.wrapJobFn(decrease_node_count,
                                                             inputs.num_nodes,
                                                             uuid)

    # wire up dag
    if not inputs.skip_alignment:
        job.addChild(bwa)

        if autoscale_after_alignment:
            bwa.addChild(increase_nodes_after_alignment)
    elif autoscale_after_alignment:
        job.addChild(increase_nodes_after_alignment)

    if (inputs.pipeline_to_run == "adam" or
        inputs.pipeline_to_run == "both"):

        if inputs.skip_preprocessing:

            # TODO: currently, we don't support autoscaling if you're just
            # running variant calling
            job.addChild(gatk_adam_call)
        else:
            if autoscale_after_alignment:
                increase_nodes_after_alignment.addChild(adam_preprocess)
            elif inputs.skip_alignment:
                job.addChild(adam_preprocess)
            else:
                bwa.addChild(adam_preprocess)

            # if we are running autoscaling, then we should decrease the cluster
            # size after adam completes
            if inputs.autoscale_cluster and not inputs.master_ip:
                adam_preprocess.addChild(decrease_nodes_after_adam_preprocess)
                decrease_nodes_after_adam_preprocess.addChild(gatk_adam_call)
            else:
                adam_preprocess.addChild(gatk_adam_call)

    if (inputs.pipeline_to_run == "gatk" or
        inputs.pipeline_to_run == "both"):

        if inputs.skip_preprocessing:

            # TODO: currently, we don't support autoscaling if you're just
            # running variant calling
            job.addChild(gatk_gatk_call)
        else:
            if autoscale_after_alignment:
                increase_nodes_after_alignment.addChild(gatk_preprocess)
            elif inputs.skip_alignment:
                job.addChild(gatk_preprocess)
            else:
                bwa.addChild(gatk_preprocess)

            gatk_preprocess.addChild(gatk_gatk_call)


def increase_node_count(job, nodes_to_add, uuid):

    Samples.increase_nodes(uuid, nodes_to_add)


def decrease_node_count(job, nodes_to_add, uuid):

    Samples.decrease_nodes(uuid, nodes_to_add)

def generate_mock_config():

    return textwrap.dedent(""" 
        # ADAM/GATK Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline.
        # This configuration file is pre-filled for use in MOCK MODE
        ##############################################################################################################  
        # MOCK INPUTS
        pipeline-to-run: both
        skip-alignment: False
        skip-preprocessing: False
        sequence-dir: sequence
        autoscale-cluster: False
        s3-bucket: adam-gatk-pipeline-mock-files
        cpu-count:
        program-unit: 12345
        platform: ILLUMINA
        ref: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_ref.fa
        amb: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_ref.fa.amb
        ann: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_ref.fa.ann
        bwt: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_ref.fa.bwt
        pac: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_ref.fa.pac
        sa: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_ref.fa.sa
        fai: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_ref.fa.fai
        alt: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_ref.fa.alt
        phase: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_phase.vcf
        mills: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_mills.vcf
        dbsnp: s3://adam-gatk-pipeline-mock-files/mock-pipeline-inputs/bqsr1.vcf
        omni: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_omni.vcf
        hapmap: https://s3-us-west-2.amazonaws.com/adam-gatk-pipeline-mock-files/mock-pipeline-inputs/mock_hapmap.vcf
        trim-adapters: False
        file-size: 10M
        s3-bucket: adam-gatk-pipeline-mock-files
        memory: 2
        dir-suffix: /mock
        num-nodes: #3
        master-ip: spark-master
        ssec:
    """[1:])


def generate_config():
    if mock_mode():
        return generate_mock_config()

    return textwrap.dedent("""
        # ADAM/GATK Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        pipeline-to-run: both     #
        skip-alignment: False     #
        skip-preprocessing: False #
        sequence-dir: sequence    #
        autoscale-cluster: False  #
        s3-bucket:                # S3 Bucket URI
        cpu-count:                # Optional:
        program-unit: 12345       #
        platform: ILLUMINA        #
        ref:                      # Required: Reference fasta file
        amb:                      # Required: Reference fasta file (amb)
        ann:                      # Required: Reference fasta file (ann)
        bwt:                      # Required: Reference fasta file (bwt)
        pac:                      # Required: Reference fasta file (pac)
        sa:                       # Required: Reference fasta file (sa)
        fai:                      # Required: Reference fasta file (fai)
        alt:                      # Optional: Alternate file for reference build (alt). Necessary for alt aware alignment.
        phase:                    # Required: URL (1000G_phase1.indels.hg19.sites.fixed.vcf)
        mills:                    # Required: URL (Mills_and_1000G_gold_standard.indels.hg19.sites.vcf)
        dbsnp:                    # Required: URL (dbsnp_132_b37.leftAligned.vcf)
        hapmap:                   # Required: URL (hapmap_3.3.b37.vcf)
        omni:                     # Required: URL (1000G_omni.5.b37.vcf)
        trim-adapters: False      # Trim adapters.
        num-nodes: 9              # Number of nodes to use. Do not set if providing master_ip.
        master-ip:                # Optional: IP or hostname of host running for Spark master and HDFS namenode.
                                  # Should be provided instead of num-nodes if pointing at a static (external or
                                  # standalone) Spark cluster. The special value 'auto' indicates the master of
                                  # an externally autoscaled cgcloud spark cluster, i.e. one that is managed by
                                  # the uberscript.
        file-size: 100G           # Approximate input file size. Should be given as %d[TGMK], e.g.,
                                  # for a 100 gigabyte file, use file_size: '100G'
        ssec:                     # Optional: (string) Path to Key File for SSE-C Encryption
        dir-suffix:               # Optional: suffix to add to output directory names.
        memory:                   # Required: Amount of available memory on each worker node.                                   
    """[1:])


def generate_mock_manifest():
    return textwrap.dedent("""
        # This manifest was generated for use in MOCK MODE
        mouse_chrM_a
        mouse_chrM_b
        """[1:])

def generate_manifest():
    if mock_mode():
        return generate_mock_manifest()
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There is a single column: UUID
        #
        #   UUID        This should be a unique identifier for the sample to be processed that corresponds to 
        #               the prefix of the filenames of the input fastq files.
        #   
        #   Example:
        #   If your input fastq file pairs were input_file_name_1.illumina_1.fastq.gz, input_file_name_1.illumina_2.fastq.gz and 
        #   input_file_name_2.illumina_1.fastq.gz, input_file_name_2.illumina_2.fastq.gz, the manifest would be:
        #
        #   input_file_name_1.illumina
        #   input_file_name_2.illumina   
        #
        #   Input fastq files MUST be named according to the filename_1.fastq.gz, filename_2.fastq.gz convention
        #
        #   Place your samples below, one per line.
        """[1:])


def main():
    """
    This is a Toil pipeline used to perform alignment of fastqs.
    """
    # Define Parser object and add to Toil
    if mock_mode():
        usage_msg = 'You have the TOIL_SCRIPTS_MOCK_MODE environment variable set, so this pipeline ' \
                    'will run in mock mode. To disable mock mode, set TOIL_SCRIPTS_MOCK_MODE=0'
    else:
        usage_msg = None

    parser = argparse.ArgumentParser(usage=usage_msg)
    subparsers = parser.add_subparsers(dest='command')
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')
    # Run subparser                                                                                                              
    parser_run = subparsers.add_parser('run', help='Runs the ADAM/GATK pipeline')
    default_config = 'adam-gatk-mock.config' if mock_mode() else 'adam-gatk.config'
    default_manifest = 'adam-gatk-mock-manifest.csv' if mock_mode() else 'adam-gatk-manifest.csv'
    parser_run.add_argument('--config', default=default_config, type=str,
                            help='Path to the (filled in) config file, generated with "generate-config".')
    parser_run.add_argument('--manifest', default=default_manifest,
                            type=str, help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                                           '\nDefault value: "%(default)s".')
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()

    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, default_config), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, default_manifest), generate_manifest)
    # Pipeline execution
    elif args.command == 'run':
        require(os.path.exists(args.config), '{} not found. Please run '
                                             'generate-config'.format(args.config))
        if not hasattr(args, 'sample'):
            require(os.path.exists(args.manifest), '{} not found and no samples provided. Please '
                                                   'run "generate-manifest"'.format(args.manifest))
        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        inputs = argparse.Namespace(**parsed_config)

        # Parse manifest file
        uuid_list = []
        with open(args.manifest) as f_manifest:
            for line in f_manifest:
                if not line.isspace() and not line.startswith('#'):
                    uuid_list.append(line.strip())

        inputs.sort = False
        if not inputs.dir_suffix:
            inputs.dir_suffix = ''
        if not inputs.s3_bucket:
            inputs.s3_bucket = ''

        if inputs.master_ip and inputs.num_nodes:
            raise ValueError("Exactly one of master_ip (%s) and num_nodes (%d) must be provided." %
                             (inputs.master_ip, inputs.num_nodes))

        if not hasattr(inputs, 'master_ip') and inputs.num_nodes <= 1:
            raise ValueError('num_nodes allocates one Spark/HDFS master and n-1 workers, and thus must be greater '
                             'than 1. %d was passed.' % inputs.num_nodes)

        if (inputs.pipeline_to_run != "adam" and
            inputs.pipeline_to_run != "gatk" and
            inputs.pipeline_to_run != "both"):
            raise ValueError("pipeline_to_run must be either 'adam', 'gatk', or 'both'. %s was passed." % inputs.pipeline_to_run)

        Job.Runner.startToil(Job.wrapJobFn(sample_loop, uuid_list, inputs), args)

if __name__=="__main__":
    main()
