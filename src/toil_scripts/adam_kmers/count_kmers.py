"""
Spark demo application.
"""

# imports from python core
import argparse
import logging
import shlex

# imports from outside of toil_scripts and python core
from toil.lib.spark import spawn_spark_cluster

# imports from toil_scripts
from toil_lib import require
from toil_lib.tools.spark_tools import call_adam, call_conductor, \
    MasterAddress, HDFS_MASTER_PORT, SPARK_MASTER_PORT


_log = logging.getLogger(__name__)


def kmer_dag(job,
             input_file,
             output_path,
             kmer_length,
             spark_conf,
             workers,
             cores,
             memory,
             sudo):
    '''
    Optionally launches a Spark cluster and then runs ADAM to count k-mers on an
    input file.

    :param job: Toil job
    :param input_file: URL/path to input file to count k-mers on
    :param output_path: URL/path to save k-mer counts at
    :param kmer_length: The length of k-mer substrings to count.
    :param spark_conf: Optional Spark configuration. If set, workers should \
    not be set.
    :param workers: Optional number of Spark workers to launch. If set, \
    spark_conf should not be set, and cores and memory should be set.
    :param cores: Number of cores per Spark worker. Must be set if workers is \
    set.
    :param memory: Amount of memory to provided to Spark workers. Must be set \
    if workers is set.
    :param sudo: Whether or not to run Spark containers with sudo.

    :type job: toil.Job
    :type input_file: string
    :type output_path: string
    :type kmer_length: int or string
    :type spark_conf: string or None
    :type workers: int or None
    :type cores: int or None
    :type memory: int or None
    :type sudo: boolean
    '''

    require((spark_conf is not None and workers is None) or
            (workers is not None and cores is not None and memory is not None and spark_conf is not None),
            "Either worker count (--workers) must be defined or user must pass in Spark configuration (--spark-conf).")

    # if we do not have a spark configuration, then we must spawn a cluster
    if spark_conf is None:
        master_hostname = spawn_spark_cluster(job,
                                              sudo,
                                              workers,
                                              cores)
    else:
        spark_conf = shlex.split(spark_conf)

    job.addChildJobFn(download_count_upload,
                      masterHostname,
                      input_file, output_file, kmer_length,
                      spark_conf, memory, sudo)

def download_count_upload(job,
                          master_ip,
                          input_file,
                          output_file,
                          kmer_length,
                          spark_conf,
                          memory,
                          sudo):
    '''
    Runs k-mer counting.

    1. If the input file is located in S3, the file is copied into HDFS.
    2. If the input file is not in Parquet format, the file is converted into Parquet.
    3. The k-mers are counted and saved as text.
    4. If the output path is an S3 URL, the file is copied back to S3.

    :param job: Toil job
    :param input_file: URL/path to input file to count k-mers on
    :param output_file: URL/path to save k-mer counts at
    :param kmer_length: The length of k-mer substrings to count.
    :param spark_conf: Optional Spark configuration. If set, memory should \
    not be set.
    :param memory: Amount of memory to provided to Spark workers. Must be set \
    if spark_conf is not set.
    :param sudo: Whether or not to run Spark containers with sudo.

    :type job: toil.Job
    :type input_file: string
    :type output_file: string
    :type kmer_length: int or string
    :type spark_conf: list of string or None
    :type memory: int or None
    :type sudo: boolean
    '''

    if master_ip is not None:
        hdfs_dir = "hdfs://{0}:{1}/".format(master_ip, HDFS_MASTER_PORT)
    else:
        _log.warn('Master IP is not set. If default filesystem is not set, jobs may fail.')
        hdfs_dir = ""

    # if the file isn't already in hdfs, copy it in
    hdfs_input_file = hdfs_dir
    if input_file.startswith("s3://"):

        # append the s3 file name to our hdfs path
        hdfs_input_file += input_file.split("/")[-1]

        # run the download
        _log.info("Downloading input file %s to %s.", input_file, hdfs_input_file)
        call_conductor(master_ip, input_file, hdfs_input_file,
                       memory=memory, override_parameters=spark_conf)

    else:
        if not input_file.startswith("hdfs://"):
            _log.warn("If not in S3, input file (%s) expected to be in HDFS (%s).",
                      input_file, hdfs_dir)

    # where are we writing the output to? is it going to a location in hdfs or not?
    run_upload = True
    hdfs_output_file = hdfs_dir + "kmer_output.txt"
    if output_file.startswith(hdfs_dir):
        run_upload = False
        hdfs_output_file = output_file
    
    # do we need to convert to adam?
    if (hdfs_input_file.endswith('.bam') or
        hdfs_input_file.endswith('.sam') or
        hdfs_input_file.endswith('.fq') or
        hdfs_input_file.endswith('.fastq')):
        
        hdfs_tmp_file = hdfs_input_file

        # change the file extension to adam
        hdfs_input_file = '.'.join(hdfs_input_file.split('.')[:-1].append('adam'))

        # convert the file
        _log.info('Converting %s into ADAM format at %s.', hdfs_tmp_file, hdfs_input_file)
        call_adam(master_ip,
                  ['transform',
                   hdfs_tmp_file, hdfs_input_file],
                  memory=memory, override_parameters=spark_conf)
        
    # run k-mer counting
    _log.info('Counting %d-mers in %s, and saving to %s.',
              kmer_length, hdfs_input_file, hdfs_output_file)
    call_adam(master_ip,
              ['count_kmers',
               hdfs_input_file, hdfs_output_file,
               str(kmer_length)],
              memory=memory, override_parameters=spark_conf)

    # do we need to upload the file back? if so, run upload
    if run_upload:
        _log.info("Uploading output file %s to %s.", hdfs_output_file, output_file)
        call_conductor(master_ip, hdfs_output_file, output_file,
                       memory=memory, override_parameters=spark_conf)
        

def main():
    '''
    Sets up command line parser for Toil/ADAM based k-mer counter, and launches
    k-mer counter with optional Spark cluster.
    '''

    parser = argparse.ArgumentParser()

    # add parser arguments
    parser.add_argument('--input_path',
                        help='The full path to the input SAM/BAM/ADAM/FASTQ file.')
    parser.add_argument('--output-path',
                        help='full path where final results will be output.')
    parser.add_argument('--kmer-length',
                        help='Length to use for k-mer counting. Defaults to 20.',
                        default=20,
                        type=int)
    parser.add_argument('--spark-conf',
                        help='Optional configuration to pass to Spark commands. Either this or --workers must be specified.',
                        default=None)
    parser.add_argument('--memory',
                        help='Optional memory configuration for Spark workers/driver. This must be specified if --workers is specified.',
                        default=None,
                        type=int)
    parser.add_argument('--cores',
                        help='Optional core configuration for Spark workers/driver. This must be specified if --workers is specified.',
                        default=None,
                        type=int)
    parser.add_argument('--workers',
                        help='Number of workers to spin up in Toil. Either this or --spark-conf must be specified. If this is specified, --memory and --cores must be specified.',
                        default=None,
                        type=int)
    parser.add_argument('--sudo',
                        help='Run docker containers with sudo. Defaults to False.',
                        default=False,
                        action='store_true')

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    Job.Runner.startToil(Job.wrapJobFn(kmer_dag,
                                       args.kmer_length,
                                       args.input_path,
                                       args.output_path,
                                       args.spark_conf,
                                       args.workers,
                                       args.cores,
                                       args.memory,
                                       args.sudo,
                                       checkpoint=True), args)
    
if __name__ == "__main__":
    main()
