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
toil            - pip install toil
"""

import argparse
import logging
import multiprocessing
import os
from subprocess import check_call, check_output
import sys

from toil.job import Job

from toil_scripts.adam_uberscript.automated_scaling import SparkMasterAddress
from toil_scripts.lib.programs import docker_call, mock_mode
from toil_scripts.spark_utils.spawn_cluster import start_spark_hdfs_cluster

SPARK_MASTER_PORT = "7077"
HDFS_MASTER_PORT = "8020"
log = logging.getLogger(__name__)


class MasterAddress(str):
    """
    A string containing the hostname or IP of the Spark/HDFS master. The Spark master expects its own address to
    match what the client uses to connect to it. For example, if the master is configured with a host name,
    the driver can't use an IP address to connect to it, and vice versa. This class works around by distinguishing
    between the notional master address (self) and the actual one (self.actual) and adds support for the special
    master address "auto" in order to implement auto-discovery of the master of a standalone.

    >>> foo = MasterAddress('foo')
    >>> foo == 'foo'
    True
    >>> foo.actual == 'foo'
    True
    >>> foo.actual == foo
    True
    """
    def __init__(self, master_ip):
        # TBD: this could do more tricks like always mapping an IP address to 'spark-master' etc.
        if master_ip == 'auto':
            super(MasterAddress, self).__init__('spark-master')
            self.actual = SparkMasterAddress.load_in_toil().value
        else:
            super(MasterAddress, self).__init__(master_ip)
            self.actual = self

    def docker_parameters(self, docker_parameters=None):
        """
        Augment a list of "docker run" arguments with those needed to map the notional Spark master address to the
        real one, if they are different.
        """
        if self != self.actual:
            add_host_option = '--add-host=spark-master:' + self.actual
            if docker_parameters is None:
                docker_parameters = [add_host_option]
            else:
                docker_parameters.append(add_host_option)
        return docker_parameters


def call_conductor(masterIP, inputs, src, dst):
    """
    Invokes the Conductor container to copy files between S3 and HDFS

    :type masterIP: MasterAddress
    """
    docker_call(rm = False,
                tool = "quay.io/ucsc_cgl/conductor",
                docker_parameters = masterIP.docker_parameters(["--net=host"]),
                parameters = ["--master", "spark://"+masterIP+":"+SPARK_MASTER_PORT,
                 "--conf", "spark.driver.memory=%sg" % inputs["driverMemory"],
                 "--conf", "spark.executor.memory=%sg" % inputs["executorMemory"],
                 "--", "-C", src, dst],
                sudo = inputs['sudo'],
                mock=False)


def call_adam(masterIP, inputs, arguments):
    """
    Invokes the ADAM container

    :type masterIP: MasterAddress
    """
    default_params = ["--master", ("spark://%s:%s" % (masterIP, SPARK_MASTER_PORT)),
                      "--conf", ("spark.driver.memory=%sg" % inputs["driverMemory"]),
                      "--conf", ("spark.executor.memory=%sg" % inputs["executorMemory"]),
                      "--conf", ("spark.hadoop.fs.default.name=hdfs://%s:%s" % (masterIP, HDFS_MASTER_PORT)),
                      "--conf", "spark.driver.maxResultSize=0", # set max result size to unlimited, see #177
                      "--"]
    docker_call(rm = False,
                tool = "quay.io/ucsc_cgl/adam:962-ehf--6e7085f8cac4b9a927dc9fb06b48007957256b80",
                docker_parameters = masterIP.docker_parameters(["--net=host"]),
                parameters = default_params + arguments,
                sudo = inputs['sudo'],
                mock=False)


def remove_file(masterIP, filename, sparkOnToil):
    """
    Remove the given file from hdfs with master at the given IP address

    :type masterIP: MasterAddress
    """
    masterIP = masterIP.actual
    
    ssh_call = ['ssh', '-o', 'StrictHostKeyChecking=no', masterIP]

    if sparkOnToil:
        output = check_output(ssh_call +  ['docker', 'ps'])
        containerID = next(line.split()[0] for line in output.splitlines() if 'apache-hadoop-master' in line)
        ssh_call += ['docker', 'exec', containerID]

    try:
        check_call(ssh_call + ['hdfs', 'dfs', '-rm', '-r', '/' + filename])
    except:
        pass


def truncate_file(masterIP, filename, sparkOnToil):
    """
    Truncate the given hdfs file to 10 bytes with master at the given IP address

    :type masterIP: MasterAddress
    """
    masterIP = masterIP.actual

    ssh_call = ['ssh', '-o', 'StrictHostKeyChecking=no', masterIP]

    if sparkOnToil:
        output = check_output(ssh_call +  ['docker', 'ps'])
        containerID = next(line.split()[0] for line in output.splitlines() if 'apache-hadoop-master' in line)
        ssh_call += ['docker', 'exec', containerID]

    try:
        check_call(ssh_call + ['hdfs', 'dfs', '-truncate', '-w', '10', '/' + filename])
    except:
        pass


def download_data(masterIP, inputs, knownSNPs, bam, hdfsSNPs, hdfsBAM):
    """
    Downloads input data files from S3.

    :type masterIP: MasterAddress
    """

    log.info("Downloading known sites file %s to %s.", knownSNPs, hdfsSNPs)
    call_conductor(masterIP, inputs, knownSNPs, hdfsSNPs)

    log.info("Downloading input BAM %s to %s.", bam, hdfsBAM)
    call_conductor(masterIP, inputs, bam, hdfsBAM)


def adam_convert(masterIP, inputs, inFile, inSnps, adamFile, adamSnps, sparkOnToil):
    """
    Convert input sam/bam file and known SNPs file into ADAM format
    """

    log.info("Converting input BAM to ADAM.")
    call_adam(masterIP,
              inputs,
              ["transform",
               inFile, adamFile])

    inFileName = inFile.split("/")[-1]
    remove_file(masterIP, inFileName, sparkOnToil)

    log.info("Converting known sites VCF to ADAM.")

    call_adam(masterIP,
              inputs,
              ["vcf2adam",
               "-only_variants",
               inSnps, adamSnps])

    inSnpsName = inSnps.split("/")[-1]
    remove_file(masterIP, inSnpsName, sparkOnToil)


def adam_transform(masterIP, inputs, inFile, snpFile, hdfsDir, outFile, sparkOnToil):
    """
    Preprocess inFile with known SNPs snpFile:
        - mark duplicates
        - realign indels
        - recalibrate base quality scores
    """

    log.info("Marking duplicate reads.")
    call_adam(masterIP,
              inputs,
              ["transform",
               inFile,  hdfsDir + "/mkdups.adam",
               "-aligned_read_predicate",
               "-limit_projection",
               "-mark_duplicate_reads"])

    #FIXME
    inFileName = inFile.split("/")[-1]
    remove_file(masterIP, inFileName+"*", sparkOnToil)

    log.info("Realigning INDELs.")
    call_adam(masterIP,
              inputs,
              ["transform",
               hdfsDir + "/mkdups.adam",
               hdfsDir + "/ri.adam",
               "-realign_indels"])

    remove_file(masterIP, hdfsDir + "/mkdups.adam*", sparkOnToil)

    log.info("Recalibrating base quality scores.")
    call_adam(masterIP,
              inputs,
              ["transform",
               hdfsDir + "/ri.adam",
               hdfsDir + "/bqsr.adam",
               "-recalibrate_base_qualities",
               "-known_snps", snpFile])

    remove_file(masterIP, "ri.adam*", sparkOnToil)

    log.info("Sorting reads and saving a single BAM file.")
    call_adam(masterIP,
              inputs,
              ["transform",
               hdfsDir + "/bqsr.adam",
               outFile,
               "-sort_reads", "-single"])

    remove_file(masterIP, "bqsr.adam*", sparkOnToil)

    return outFile


def upload_data(masterIP, inputs, hdfsName, uploadName, sparkOnToil):
    """
    Upload file hdfsName from hdfs to s3
    """

    if mock_mode():
        truncate_file(masterIP, hdfsName, sparkOnToil)

    log.info("Uploading output BAM %s to %s.", hdfsName, uploadName)
    call_conductor(masterIP, inputs, hdfsName, uploadName)


def download_run_and_upload(job, masterIP, inputs, sparkOnToil):
    """
    Monolithic job that calls data download, conversion, transform, upload.
    Previously, this was not monolithic; change came in due to #126/#134.
    """
    masterIP = MasterAddress(masterIP)

    bamName = inputs['bamName'].split('://')[-1].split('/')[-1]
    sampleName = ".".join(os.path.splitext(bamName)[:-1])
    hdfsSubdir = sampleName + "-dir"
    hdfsDir = "hdfs://{0}:{1}/{2}".format(masterIP, HDFS_MASTER_PORT, hdfsSubdir)

    try:
        hdfsPrefix = hdfsDir + "/" + sampleName
        hdfsBAM = hdfsDir + "/" + bamName
        
        hdfsSNPs = hdfsDir + "/" + inputs['knownSNPs'].split('://')[-1].split('/')[-1]

        download_data(masterIP, inputs, inputs['knownSNPs'], inputs['bamName'], hdfsSNPs, hdfsBAM)

        adamInput = hdfsPrefix + ".adam"
        adamSNPs = hdfsDir + "/snps.var.adam"
        adam_convert(masterIP, inputs, hdfsBAM, hdfsSNPs, adamInput, adamSNPs, sparkOnToil)

        adamOutput = hdfsPrefix + ".processed.adam" 
        adam_transform(masterIP, inputs, adamInput, adamSNPs, hdfsDir, adamOutput, sparkOnToil)

        suffix = inputs['suffix'] if inputs['suffix'] else ''
        outFile = inputs['outDir'] + "/" + sampleName + suffix + ".bam"

        upload_data(masterIP, inputs, adamOutput, outFile, sparkOnToil)

    except:
        if sparkOnToil:
            # if a stage failed, we must clean up HDFS for the pipeline to succeed on retry
            remove_file(masterIP, hdfsSubdir, sparkOnToil)
        else:
            log.warning("Jobs failed, but cannot clean up files on a non-toil cluster.")

        raise


def static_adam_preprocessing_dag(job, inputs):
    """
    A Toil job function performing ADAM preprocessing on a single sample
    """
    masterIP = inputs['masterIP']
    if not masterIP:
        # Dynamic subclusters, i.e. Spark-on-Toil
        sparkOnToil = True
        cores = multiprocessing.cpu_count()
        startCluster = job.wrapJobFn(start_spark_hdfs_cluster,
                                     inputs['numWorkers'],
                                     inputs['executorMemory'],
                                     inputs['sudo'],
                                     download_run_and_upload,
                                     jArgs = (inputs, sparkOnToil),
                                     jCores = cores,
                                     jMemory = "%s G" % inputs['driverMemory']).encapsulate()
        job.addChild(startCluster)
    elif masterIP == 'auto':
        # Static, standalone Spark cluster managed by uberscript
        sparkOnToil = False
        scaleUp = job.wrapJobFn(scale_external_spark_cluster, 1)
        job.addChild(scaleUp)
        sparkWork = job.wrapJobFn(download_run_and_upload, masterIP, inputs, sparkOnToil)
        scaleUp.addChild(sparkWork)
        scaleDown = job.wrapJobFn(scale_external_spark_cluster, -1)
        sparkWork.addChild(scaleDown)
    else:
        # Static, external Spark cluster
        sparkOnToil = False
        sparkWork = job.wrapJobFn(download_run_and_upload, masterIP, inputs, sparkOnToil)
        job.addChild(sparkWork)


def scale_external_spark_cluster(num_samples=1):
    from toil_scripts.adam_uberscript.adam_uberscript import standalone_spark_semaphore_name
    from toil_scripts.adam_uberscript.automated_scaling import Semaphore
    sem = Semaphore.load_in_toil(standalone_spark_semaphore_name)
    if num_samples > 0:
        sem.acquire(delta=num_samples)
    else:
        sem.release(delta=-num_samples)


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
    parser.add_argument('-m', '--master_ip', required = False,
                        help="IP or hostname of host running for Spark master and HDFS namenode. Should be provided "
                             "if pointing at a static (external or standalone) Spark cluster. The special value "
                             "'auto' indicates the master of standalone cluster, i.e. one that is managed by the "
                             "uberscript.")
    parser.add_argument('-u', '--sudo',
                        dest='sudo', action='store_true',
                        help='Docker usually needs sudo to execute '
                        'locally, but not''when running Mesos '
                        'or when a member of a Docker group.')

    return parser

# FIXME: unused parameter args

def main(args):

    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    options = parser.parse_args()

    if bool(options.master_ip) == bool(options.num_nodes):
        raise ValueError("Only one of --master_ip (%s) and --num_nodes (%d) can be provided." %
                         (options.master_ip, options.num_nodes))

    if options.num_nodes <= 1:
        raise ValueError("--num_nodes allocates one Spark/HDFS master and n-1 workers, and thus must be greater "
                         "than 1. %d was passed." % options.num_nodes)

    inputs = {'numWorkers': options.num_nodes - 1,
              'outDir':     options.output_directory,
              'bamName':    options.input_file_name,
              'knownSNPs':  options.known_SNPs,
              'driverMemory': options.driver_memory,
              'executorMemory': options.executor_memory,
              'sudo': options.sudo,
              'suffix': None,
              'masterIP': options.master_ip}

    Job.Runner.startToil(Job.wrapJobFn(static_adam_preprocessing_dag, inputs), options)

if __name__=="__main__":
    # FIXME: main() doesn't return anything
    sys.exit(main(sys.argv[1:]))
