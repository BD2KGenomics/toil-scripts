#!/usr/bin/env python2.7
"""
UCSC Computational Genomics Lab ADAM/Spark pipeline

@author Audrey Musselman-Brown, almussel@ucsc.edu
@author Frank Austin Nothaft, fnothaft@berkeley.edu

Please see the README.md in the same directory

Toil pipeline for ADAM preprocessing

                        Tree structure of ADAM pipeline
                                       0
                                       |++(1)
                                       2 --> 4 --> 5 --> 6 -> 7
                                        ++(3)

0 = Start Master
1 = Master Service
2 = Start Workers
3 = Worker Service
4 = Do All The Things (Download from s3, convert to ADAM, preprocess, upload to s3)

================================================================================
:Dependencies
docker          - apt-get install docker (or 'docker.io' for linux)
toil            - pip install toil
"""

import argparse
import logging
import multiprocessing
import os
import sys
import textwrap
from subprocess import check_call, check_output

import yaml
from toil.job import Job
from toil.lib.spark import spawn_spark_cluster

from toil_lib import require
from toil_lib.files import copy_files, move_files
from toil_lib.programs import docker_call, mock_mode
from toil_lib.tools.spark_tools import call_adam, call_conductor, MasterAddress, HDFS_MASTER_PORT, SPARK_MASTER_PORT

from toil_lib.files import generate_file

log = logging.getLogger(__name__)


def remove_file(master_ip, filename, spark_on_toil):
    """
    Remove the given file from hdfs with master at the given IP address

    :type masterIP: MasterAddress
    """
    master_ip = master_ip.actual

    ssh_call = ['ssh', '-o', 'StrictHostKeyChecking=no', master_ip]

    if spark_on_toil:
        output = check_output(ssh_call + ['docker', 'ps'])
        container_id = next(line.split()[0] for line in output.splitlines() if 'apache-hadoop-master' in line)
        ssh_call += ['docker', 'exec', container_id]

    try:
        check_call(ssh_call + ['hdfs', 'dfs', '-rm', '-r', '/' + filename])
    except:
        pass


def truncate_file(master_ip, filename, spark_on_toil):
    """
    Truncate the given hdfs file to 10 bytes with master at the given IP address

    :type masterIP: MasterAddress
    """
    master_ip = master_ip.actual

    ssh_call = ['ssh', '-o', 'StrictHostKeyChecking=no', master_ip]
    hdfs = ['hdfs']

    if spark_on_toil:
        output = check_output(ssh_call + ['docker', 'ps'])
        container_id = next(line.split()[0] for line in output.splitlines() if 'apache-hadoop-master' in line)
        ssh_call += ['docker', 'exec', container_id]
        hdfs = ['/opt/apache-hadoop/bin/hdfs']

    try:
        check_call(ssh_call + hdfs + ['dfs', '-truncate', '-w', '10', '/' + filename])
    except:
        pass


def download_data(master_ip, inputs, known_snps, bam, hdfs_snps, hdfs_bam):
    """
    Downloads input data files from S3.

    :type masterIP: MasterAddress
    """

    log.info("Downloading known sites file %s to %s.", known_snps, hdfs_snps)
    call_conductor(master_ip, known_snps, hdfs_snps, memory=inputs.memory)

    log.info("Downloading input BAM %s to %s.", bam, hdfs_bam)
    call_conductor(master_ip, bam, hdfs_bam, memory=inputs.memory)


def adam_convert(master_ip, inputs, in_file, in_snps, adam_file, adam_snps, spark_on_toil):
    """
    Convert input sam/bam file and known SNPs file into ADAM format
    """

    log.info("Converting input BAM to ADAM.")
    call_adam(master_ip,
              ["transform", in_file, adam_file],
              memory=inputs.memory,
              run_local=inputs.run_local,
              native_adam_path=inputs.native_adam_path)

    in_file_name = in_file.split("/")[-1]
    remove_file(master_ip, in_file_name, spark_on_toil)

    log.info("Converting known sites VCF to ADAM.")

    call_adam(master_ip,
              ["vcf2adam", "-only_variants", in_snps, adam_snps],
              memory=inputs.memory,
              run_local=inputs.run_local,
              native_adam_path=inputs.native_adam_path)

    in_snps_name = in_snps.split("/")[-1]
    remove_file(master_ip, in_snps_name, spark_on_toil)


def adam_transform(master_ip, inputs, in_file, snp_file, hdfs_dir, out_file, spark_on_toil):
    """
    Preprocess in_file with known SNPs snp_file:
        - mark duplicates
        - realign indels
        - recalibrate base quality scores
    """

    log.info("Marking duplicate reads.")
    call_adam(master_ip,
              ["transform",
               in_file,  hdfs_dir + "/mkdups.adam",
               "-aligned_read_predicate",
               "-limit_projection",
               "-mark_duplicate_reads"],
              memory=inputs.memory,
              run_local=inputs.run_local,
              native_adam_path=inputs.native_adam_path)

    #FIXME
    in_file_name = in_file.split("/")[-1]
    remove_file(master_ip, in_file_name + "*", spark_on_toil)

    log.info("Realigning INDELs.")
    call_adam(master_ip,
              ["transform",
               hdfs_dir + "/mkdups.adam",
               hdfs_dir + "/ri.adam",
               "-realign_indels"],
              memory=inputs.memory,
              run_local=inputs.run_local,
              native_adam_path=inputs.native_adam_path)

    remove_file(master_ip, hdfs_dir + "/mkdups.adam*", spark_on_toil)

    log.info("Recalibrating base quality scores.")
    call_adam(master_ip,
              ["transform",
               hdfs_dir + "/ri.adam",
               hdfs_dir + "/bqsr.adam",
               "-recalibrate_base_qualities",
               "-known_snps", snp_file],
              memory=inputs.memory,
              run_local=inputs.run_local,
              native_adam_path=inputs.native_adam_path)

    remove_file(master_ip, "ri.adam*", spark_on_toil)

    log.info("Sorting reads and saving a single BAM file.")
    call_adam(master_ip,
              ["transform",
               hdfs_dir + "/bqsr.adam",
               out_file,
               "-sort_reads", "-single"],
              memory=inputs.memory,
              run_local=inputs.run_local,
              native_adam_path=inputs.native_adam_path)

    remove_file(master_ip, "bqsr.adam*", spark_on_toil)

    return out_file


def upload_data(master_ip, inputs, hdfs_name, upload_name, spark_on_toil):
    """
    Upload file hdfsName from hdfs to s3
    """

    if mock_mode():
        truncate_file(master_ip, hdfs_name, spark_on_toil)

    log.info("Uploading output BAM %s to %s.", hdfs_name, upload_name)
    call_conductor(master_ip, hdfs_name, upload_name, memory=inputs.memory)
    remove_file(master_ip, hdfs_name, spark_on_toil)


def download_run_and_upload(job, master_ip, inputs, spark_on_toil):
    """
    Monolithic job that calls data download, conversion, transform, upload.
    Previously, this was not monolithic; change came in due to #126/#134.
    """
    master_ip = MasterAddress(master_ip)

    bam_name = inputs.sample.split('://')[-1].split('/')[-1]
    sample_name = ".".join(os.path.splitext(bam_name)[:-1])

    hdfs_subdir = sample_name + "-dir"

    if inputs.run_local:
        inputs.local_dir = job.fileStore.getLocalTempDir()
        if inputs.native_adam_path is None:
            hdfs_dir = "/data/"
        else:
            hdfs_dir = inputs.local_dir
    else:
        inputs.local_dir = None
        hdfs_dir = "hdfs://{0}:{1}/{2}".format(master_ip, HDFS_MASTER_PORT, hdfs_subdir)

    try:
        hdfs_prefix = hdfs_dir + "/" + sample_name
        hdfs_bam = hdfs_dir + "/" + bam_name

        hdfs_snps = hdfs_dir + "/" + inputs.dbsnp.split('://')[-1].split('/')[-1]

        if not inputs.run_local:
            download_data(master_ip, inputs, inputs.dbsnp, inputs.sample, hdfs_snps, hdfs_bam)
        else:
            copy_files([inputs.sample, inputs.dbsnp], inputs.local_dir)

        adam_input = hdfs_prefix + ".adam"
        adam_snps = hdfs_dir + "/snps.var.adam"
        adam_convert(master_ip, inputs, hdfs_bam, hdfs_snps, adam_input, adam_snps, spark_on_toil)

        adam_output = hdfs_prefix + ".processed.bam"
        adam_transform(master_ip, inputs, adam_input, adam_snps, hdfs_dir, adam_output, spark_on_toil)

        out_file = inputs.output_dir + "/" + sample_name + inputs.suffix + ".bam"

        if not inputs.run_local:
            upload_data(master_ip, inputs, adam_output, out_file, spark_on_toil)
        else:
            local_adam_output = "%s/%s.processed.bam" % (inputs.local_dir, sample_name)
            move_files([local_adam_output], inputs.output_dir)

        remove_file(master_ip, hdfs_subdir, spark_on_toil)
    except:
        remove_file(master_ip, hdfs_subdir, spark_on_toil)
        raise


def static_adam_preprocessing_dag(job, inputs, sample, output_dir, suffix=''):
    """
    A Toil job function performing ADAM preprocessing on a single sample
    """
    inputs.sample = sample
    inputs.output_dir = output_dir
    inputs.suffix = suffix

    if inputs.master_ip is not None or inputs.run_local:
        if not inputs.run_local and inputs.master_ip == 'auto':
            # Static, standalone Spark cluster managed by uberscript
            spark_on_toil = False
            scale_up = job.wrapJobFn(scale_external_spark_cluster, 1)
            job.addChild(scale_up)
            spark_work = job.wrapJobFn(download_run_and_upload,
                                       inputs.master_ip, inputs, spark_on_toil)
            scale_up.addChild(spark_work)
            scale_down = job.wrapJobFn(scale_external_spark_cluster, -1)
            spark_work.addChild(scale_down)
        else:
            # Static, external Spark cluster
            spark_on_toil = False
            spark_work = job.wrapJobFn(download_run_and_upload,
                                       inputs.master_ip, inputs, spark_on_toil)
            job.addChild(spark_work)
    else:
        # Dynamic subclusters, i.e. Spark-on-Toil
        spark_on_toil = True
        cores = multiprocessing.cpu_count()
        master_ip = spawn_spark_cluster(job,
                                        False, # Sudo
                                        inputs.num_nodes-1,
                                        cores=cores,
                                        memory=inputs.memory)
        spark_work = job.wrapJobFn(download_run_and_upload,
                                   master_ip, inputs, spark_on_toil)
        job.addChild(spark_work)


def scale_external_spark_cluster(num_samples=1):
    from toil_scripts.adam_uberscript.adam_uberscript import standalone_spark_semaphore_name
    from toil_scripts.adam_uberscript.automated_scaling import Semaphore
    sem = Semaphore.load_in_toil(standalone_spark_semaphore_name)
    if num_samples > 0:
        sem.acquire(delta=num_samples)
    else:
        sem.release(delta=-num_samples)


def generate_config():
    return textwrap.dedent("""
        # ADAM Preprocessing Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        num-nodes: 9              # Optional: Number of nodes to use. Do not set if providing master_ip.
        master-ip:                # Optional: IP or hostname of host running for Spark master and HDFS namenode.
                                  # Should be provided instead of num-nodes if pointing at a static (external or
                                  # standalone) Spark cluster. The special value 'auto' indicates the master of
                                  # an externally autoscaled cgcloud spark cluster, i.e. one that is managed by
                                  # the uberscript.
        dbsnp:                    # Required: The full s3 url of a VCF file of known snps
        memory:                   # Required: Amount of memory to allocate for Spark Driver and executor.
                                  # This should be equal to the available memory on each worker node.
        run-local:                # Optional: If true, runs ADAM locally and doesn't connect to a cluster.
        local-dir:                # Required if run-local is true. Sets the local directory to use for input.
        native-adam-path:         # Optional: If set, runs ADAM using the local build of ADAM at this path.
    """[1:])


def main():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the ADAM preprocessing pipeline')
    parser_run.add_argument('--config', default='adam_preprocessing.config', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--sample', help='The S3 URL or local path to the input SAM or BAM file.'
                            'NOTE: unlike other pipelines, we do not support ftp://, gnos://, etc. schemes.')
    parser_run.add_argument('--output-dir', required=True, default=None,
                            help='full path where final results will be output')
    parser_run.add_argument('-s', '--suffix', default='',
                            help='Additional suffix to add to the names of the output files')

    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()
    cwd = os.getcwd()
    if args.command == 'generate-config':
        generate_file(os.path.join(cwd, 'adam-preprocessing.config'), generate_config)
    # Pipeline execution
    elif args.command == 'run':
        require(os.path.exists(args.config), '{} not found. Please run '
                                             'generate-config'.format(args.config))
        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        inputs = argparse.Namespace(**parsed_config)

        require(not (inputs.master_ip and inputs.num_nodes),
            'Only one of master_ip and num_nodes can be provided.')

        if not hasattr(inputs, 'master_ip'):
            require(inputs.num_nodes > 1,
                'num_nodes allocates one Spark/HDFS master and n-1 workers, and '
                'thus must be greater than 1. %d was passed.' % inputs.num_nodes)

        for arg in [inputs.dbsnp, inputs.memory]:
            require(arg, 'Required argument {} missing from config'.format(arg))

        Job.Runner.startToil(Job.wrapJobFn(static_adam_preprocessing_dag, inputs,
                                           args.sample, args.output_dir), args)

if __name__ == "__main__":
    main()
