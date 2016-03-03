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
4 = Do All The Things (Download from S3, convert to ADAM, preprocess, upload to S3)

================================================================================
:Dependencies
docker          - apt-get install docker (or 'docker.io' for linux)
toil            - pip install --pre toil
"""

# python core libraries
import argparse
import logging
import multiprocessing
import os
from subprocess import call, check_call, check_output
import sys
import time

# toil imports
from toil.job import Job

# toil scripts imports
from toil_scripts.batch_alignment.bwa_alignment import docker_call
from toil_scripts.spark_utils.spawn_cluster import *

SPARK_MASTER_PORT = "7077"
HDFS_MASTER_PORT = "8020"
log = logging.getLogger(__name__)

def call_conductor(masterIP, inputs, src, dst):
    """
    Invokes the conductor container.
    """
    docker_call(no_rm = True,
                work_dir = os.getcwd(),
                tool = "quay.io/ucsc_cgl/conductor",
                docker_parameters = ["--net=host"],
                tool_parameters = ["--master", "spark://"+masterIP+":"+SPARK_MASTER_PORT,
                 "--conf", "spark.driver.memory=%sg" % inputs["driverMemory"],
                 "--conf", "spark.executor.memory=%sg" % inputs["executorMemory"],
                 "--", "-C", src, dst],
                sudo = inputs['sudo'])


def call_adam(masterIP, inputs, arguments):

    default_params = ["--master", ("spark://%s:%s" % (masterIP, SPARK_MASTER_PORT)), 
                      "--conf", ("spark.driver.memory=%sg" % inputs["driverMemory"]),
                      "--conf", ("spark.executor.memory=%sg" % inputs["executorMemory"]),
                      "--conf", ("spark.hadoop.fs.default.name=hdfs://%s:%s" % (masterIP, HDFS_MASTER_PORT)),
                      "--"]

    docker_call(no_rm = True,
                work_dir = os.getcwd(),
                tool = "quay.io/ucsc_cgl/adam:962-ehf--6e7085f8cac4b9a927dc9fb06b48007957256b80",
                docker_parameters = ["--net=host"],
                tool_parameters = default_params + arguments,
                sudo = inputs['sudo'])


def remove_file(masterIP, filename, sparkOnToil):
    """
    Remove the given file from hdfs with master at the given IP address
    """
    if sparkOnToil:
        containerID = check_output(["ssh", "-o", "StrictHostKeyChecking=no", masterIP, "docker", "ps", \
                                "|", "grep", "apache-hadoop-master", "|", "awk", "'{print $1}'"])[:-1]
        check_call(["ssh", "-o", "StrictHostKeyChecking=no", masterIP, "docker", "exec", containerID, \
                "/opt/apache-hadoop/bin/hdfs", "dfs", "-rm", "-r", "/"+filename])
    else:
        log.warning("Cannot remove file %s. Can only remove files when running Spark-on-Toil", filename)


def download_data(masterIP, inputs, sparkOnToil):
    """
    Downloads input data files from s3.
    """

    
    snpFileSystem, snpPath = inputs['knownSNPs'].split('://')
    snpName = snpPath.split('/')[-1]
    hdfsSNPs = "hdfs://"+masterIP+":"+HDFS_MASTER_PORT+"/"+snpName
    
    log.info("Downloading known sites file %s to %s.", inputs['knownSNPs'], hdfsSNPs)
    call_conductor(masterIP, inputs, inputs['knownSNPs'], hdfsSNPs)
        
    bamFileSystem, bamPath = inputs['bamName'].split('://')
    bamName = bamPath.split('/')[-1]
    hdfsBAM = "hdfs://"+masterIP+":"+HDFS_MASTER_PORT+"/"+bamName

    log.info("Downloading input BAM %s to %s.", inputs['bamName'], hdfsBAM)
    call_conductor(masterIP, inputs, inputs['bamName'], hdfsBAM)
     
    return (hdfsBAM, hdfsSNPs)


def adam_convert(masterIP, inFile, snpFile, inputs, sparkOnToil):
    """
    Convert input sam/bam file and known SNPs file into ADAM format
    """

    log.info("Converting input BAM to ADAM.")
    adamFile = ".".join(os.path.splitext(inFile)[:-1])+".adam"
    
    call_adam(masterIP,
              inputs,
              ["transform", 
               inFile, adamFile])
              
    inFileName = inFile.split("/")[-1]
    remove_file(masterIP, inFileName, sparkOnToil)

    log.info("Converting known sites VCF to ADAM.")
    adamSnpFile = ".".join(os.path.splitext(snpFile)[:-1])+".var.adam"

    call_adam(masterIP,
              inputs,
              ["vcf2adam", 
               "-only_variants", 
               snpFile, adamSnpFile])

    snpFileName = snpFile.split("/")[-1]
    remove_file(masterIP, snpFileName, sparkOnToil)
 
    return (adamFile, adamSnpFile)

def adam_transform(masterIP, inFile, snpFile, inputs, sparkOnToil):
    """
    Preprocess inFile with known SNPs snpFile:
        - mark duplicates
        - realign indels
        - recalibrate base quality scores
    """

    outFile = ".".join(os.path.splitext(inFile)[:-1])+".processed.bam"

    log.info("Marking duplicate reads.")
    call_adam(masterIP,
              inputs,
              ["transform", 
               inFile,  "hdfs://%s:%s/mkdups.adam" % (masterIP, HDFS_MASTER_PORT),
               "-aligned_read_predicate",
               "-limit_projection",
               "-mark_duplicate_reads"])

    inFileName = inFile.split("/")[-1]
    remove_file(masterIP, inFileName+"*", sparkOnToil)

    log.info("Realigning INDELs.")
    call_adam(masterIP,
              inputs,
              ["transform", 
               "hdfs://%s:%s/mkdups.adam" % (masterIP, HDFS_MASTER_PORT),
               "hdfs://%s:%s/ri.adam" % (masterIP, HDFS_MASTER_PORT),
               "-realign_indels"])

    remove_file(masterIP, "mkdups.adam*", sparkOnToil)

    log.info("Recalibrating base quality scores.")
    call_adam(masterIP,
              inputs,
              ["transform", 
               "hdfs://%s:%s/ri.adam" % (masterIP, HDFS_MASTER_PORT),
               "hdfs://%s:%s/bqsr.adam" % (masterIP, HDFS_MASTER_PORT),
               "-recalibrate_base_qualities", 
               "-known_snps", snpFile])
              
    remove_file(masterIP, "ri.adam*", sparkOnToil)

    log.info("Sorting reads and saving a single BAM file.")
    call_adam(masterIP,
              inputs,
              ["transform", 
               "hdfs://%s:%s/bqsr.adam" % (masterIP, HDFS_MASTER_PORT), 
               outFile,
               "-sort_reads", "-single"])

    remove_file(masterIP, "bqsr.adam*", sparkOnToil)

    return outFile


def upload_data(masterIP, hdfsName, inputs, sparkOnToil):
    """
    Upload file hdfsName from hdfs to s3
    """

    fileSystem, path = hdfsName.split('://')
    nameOnly = path.split('/')[-1]
    
    uploadName = "%s/%s" % (inputs['outDir'], nameOnly.replace('.processed', ''))
    if inputs['suffix']:
        uploadName = uploadName.replace('.bam', '%s.bam' % inputs['suffix'])

    log.info("Uploading output BAM %s to %s.", hdfsName, uploadName)
    call_conductor(masterIP, inputs, hdfsName, uploadName)


# toil jobs

def download_run_and_upload(job, masterIP, inputs, sparkOnToil):
    """
    Monolithic job that calls data download, conversion, transform, upload.
    Previously, this was not monolithic; change came in due to #126/#134.
    """
    try:
        bam, snps = download_data(masterIP, inputs, sparkOnToil)
        adamInput, adamSnps = adam_convert(masterIP, bam, snps, inputs, sparkOnToil)
        adamFile = adam_transform(masterIP, adamInput, adamSnps, inputs, sparkOnToil)
        upload_data(masterIP, adamFile, inputs, sparkOnToil)
    except:
        if sparkOnToil:
            # if a stage failed, we must clean up HDFS for the pipeline to succeed on retry
            remove_file(masterIP, '*', sparkOnToil)
        else:
            log.warning("Jobs failed, but cannot clean up files on a non-toil cluster.")

        raise

def static_adam_preprocessing_dag(job, inputs):

    masterIP = inputs['masterIP']
    followJob = job
    sparkOnToil = not masterIP

    # if the master IP string was not passed in, then we need to start a spark cluster
    if sparkOnToil:
        startMaster = job.wrapJobFn(start_spark_hdfs_master,
                                    inputs['numWorkers'],
                                    inputs['executorMemory'],
                                    inputs['sudo'])
        job.addChild(startMaster)
        masterIP = startMaster.rv
        
        startWorkers = job.wrapJobFn(start_spark_hdfs_workers,
                                     masterIP,
                                     inputs['numWorkers'],
                                     inputs['executorMemory'],
                                     inputs['sudo'])
        startMaster.addChild(startWorkers)

        followJob = startWorkers

    # transformations should follow spark cluster creation if spark cluster created,
    # otherwise, should run first
    sparkWork = job.wrapJobFn(download_run_and_upload, masterIP, inputs, sparkOnToil)
    followJob.addChild(sparkWork)

def build_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file_name', required = True,
                        help = "The full s3 url of the input SAM or BAM file")
    parser.add_argument('-n', '--num_nodes', type = int, required = False,
                        default = None,
                        help = 'Number of nodes to use. Do not set if providing --master_ip.')
    parser.add_argument('-o', '--output_directory', required = True,
                        help = 's3 directory url in which to place output files')
    parser.add_argument('-k', '--known_SNPs', required = True,
                        help = 'The full s3 url of a VCF file of known snps')
    parser.add_argument('-d', '--driver_memory', required = True,
                        help = 'Amount of memory to allocate for Spark Driver.')
    parser.add_argument('-q', '--executor_memory', required = True,
                        help = 'Amount of memory to allocate per Spark Executor.')
    parser.add_argument('-j', '--jobstore', required = True,
                        help = 'Name of the jobstore')
    parser.add_argument('-m', '--master_ip', required = False, default = None,
                        help = 'IP for the Spark master/HDFS Namenode, if not spawning a cluster.')
    parser.add_argument('-u', '--sudo',
                        dest='sudo', action='store_true',
                        help='Docker usually needs sudo to execute '
                        'locally, but not''when running Mesos '
                        'or when a member of a Docker group.')

    return parser


def main(args):
    
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    options = parser.parse_args()

    if not ((options.master_ip and not options.num_nodes) or
            (not options.master_ip and options.num_nodes)):
        raise ValueError("Only one of --master_ip (%s) and --num_nodes (%d) can be provided." % 
                         (options.master_ip, options.num_nodes))

    if options.num_nodes <= 1:
        raise ValueError("--num_nodes allocates one Spark/HDFS master and n-1 workers, and thus must be greater than 1. %d was passed." %
                         options.num_nodes)

    inputs = {'numWorkers': options.num_nodes - 1,
              'outDir':     options.output_directory,
              'bamName':    options.input_file_name,
              'knownSNPs':  options.known_SNPs,
              'driverMemory': options.driver_memory,
              'executorMemory': options.executor_memory,
              'sudo': options.sudo,
              'suffix': None,
              'masterIP': options.master_ip}

    Job.Runner.startToil(Job.wrapJobFn(start_master, inputs), options)

if __name__=="__main__":
    sys.exit(main(sys.argv[1:]))
