#!/usr/bin/env python2.7
"""
UCSC Computational Genomics Lab RNA-seq Pipeline
Author: Audrey Musselman-Brown
Affiliation: UC Santa Cruz Genomics Institute

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

import argparse
import os
import multiprocessing
from subprocess import call, check_call, check_output
import sys
import time
from toil.job import Job
from toil_scripts.batch_alignment.bwa_alignment import docker_call

SPARK_MASTER_PORT = "7077"
HDFS_MASTER_PORT = "8020"
log = open("python.log", 'a')

# JOB FUNCTIONS

def start_master(job, inputs):
    """
    Starts the master service.
    """

    log.write("master job\n")
    log.flush()
    masterIP = job.addService(MasterService(inputs['sudo'], "%s G" % inputs['executorMemory']))
    job.addChildJobFn(start_workers, masterIP, inputs)


def start_workers(job, masterIP, inputs):
    """
    Starts the worker services.
    """
    log.write("workers job\n")
    log.flush()
    for i in range(inputs['numWorkers']):
        job.addService(WorkerService(masterIP, inputs['sudo'], "%s G" % inputs['executorMemory']))
    job.addFollowOnJobFn(download_run_and_upload, masterIP, inputs, memory = "%s G" % inputs['driverMemory'])


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


def remove_file(masterIP, filename):
    """
    Remove the given file from hdfs with master at the given IP address
    """
    containerID = check_output(["ssh", "-o", "StrictHostKeyChecking=no", masterIP, "docker", "ps", \
                                "|", "grep", "apache-hadoop-master", "|", "awk", "'{print $1}'"])[:-1]
    check_call(["ssh", "-o", "StrictHostKeyChecking=no", masterIP, "docker", "exec", containerID, \
                "/opt/apache-hadoop/bin/hdfs", "dfs", "-rm", "-r", "/"+filename])


def download_run_and_upload(job, masterIP, inputs):
    """
    Monolithic job that calls data download, conversion, transform, upload.
    Previously, this was not monolithic; change came in due to #126/#134.
    """
    try:
        bam, snps = download_data(masterIP, inputs)
        adam_convert(masterIP, bam, snps, inputs)
        adamFile = adam_transform(masterIP, bam, snps, inputs)
        upload_data(masterIP, adamFile, inputs)
    except:
        # if a stage failed, we must clean up HDFS for the pipeline to succeed on retry
        remove_file(masterIP, '*')


def download_data(masterIP, inputs):
    """
    Downloads input data files from s3.
    """
    log.write("download data\n")
    log.flush()
    
    snpFileSystem, snpPath = inputs['knownSNPs'].split('://')
    snpName = snpPath.split('/')[-1]
    hdfsSNPs = "hdfs://"+masterIP+":"+HDFS_MASTER_PORT+"/"+snpName
    
    call_conductor(masterIP, inputs, inputs['knownSNPs'], hdfsSNPs)
        
    bamFileSystem, bamPath = inputs['bamName'].split('://')
    bamName = bamPath.split('/')[-1]
    hdfsBAM = "hdfs://"+masterIP+":"+HDFS_MASTER_PORT+"/"+bamName

    call_conductor(masterIP, inputs, inputs['bamName'], hdfsBAM)
     
    return (hdfsBAM, hdfsSNPs)


def adam_convert(masterIP, inFile, snpFile, inputs):
    """
    Convert input sam/bam file and known SNPs file into ADAM format
    """
    log.write("adam convert\n")
    log.flush()

    adamFile = ".".join(os.path.splitext(inFile)[:-1])+".adam"
    
    call_adam(masterIP,
              inputs,
              ["transform", 
               inFile, adamFile])
              
    inFileName = inFile.split("/")[-1]
    remove_file(masterIP, inFileName)

    adamSnpFile = ".".join(os.path.splitext(snpFile)[:-1])+".var.adam"

    call_adam(masterIP,
              inputs,
              ["vcf2adam", 
               "-only_variants", 
               snpFile, adamSnpFile])

    snpFileName = snpFile.split("/")[-1]
    remove_file(masterIP, snpFileName)
 

def adam_transform(masterIP, inFile, snpFile, inputs):
    """
    Preprocess inFile with known SNPs snpFile:
        - mark duplicates
        - realign indels
        - recalibrate base quality scores
    """
    log.write("adam transform\n")
    log.flush()

    outFile = ".".join(os.path.splitext(inFile)[:-1])+".processed.bam"

    call_adam(masterIP,
              inputs,
              ["transform", 
               inFile,  "hdfs://%s:%s/mkdups.adam" % (masterIP, HDFS_MASTER_PORT),
               "-aligned_read_predicate",
               "-limit_projection",
               "-mark_duplicate_reads"])

    inFileName = inFile.split("/")[-1]
    remove_file(masterIP, inFileName+"*")

    call_adam(masterIP,
              inputs,
              ["transform", 
               "hdfs://%s:%s/mkdups.adam" % (masterIP, HDFS_MASTER_PORT),
               "hdfs://%s:%s/ri.adam" % (masterIP, HDFS_MASTER_PORT),
               "-realign_indels"])

    remove_file(masterIP, "mkdups.adam*")

    call_adam(masterIP,
              inputs,
              ["transform", 
               "hdfs://%s:%s/ri.adam" % (masterIP, HDFS_MASTER_PORT),
               "hdfs://%s:%s/bqsr.adam" % (masterIP, HDFS_MASTER_PORT),
               "-recalibrate_base_qualities", 
               "-known_snps", snpFile])
              
    remove_file(masterIP, "ri.adam*")

    call_adam(masterIP,
              inputs,
              ["transform", 
               "hdfs://%s:%s/bqsr.adam" % (masterIP, HDFS_MASTER_PORT), 
               outFile,
               "-sort_reads", "-single"])

    remove_file(masterIP, "bqsr.adam*")

    return outfile


def upload_data(masterIP, hdfsName, inputs):
    """
    Upload file hdfsName from hdfs to s3
    """
    log.write("write data\n")
    log.flush()

    fileSystem, path = hdfsName.split('://')
    nameOnly = path.split('/')[-1]
    
    uploadName = "%s/%s" % (inputs['outDir'], nameOnly.replace('.processed', ''))
    if inputs['suffix']:
        uploadName = uploadName.replace('.bam', '%s.bam' % inputs['suffix'])

    call_conductor(masterIP, inputs, hdfsName, uploadName)
    
    
# SERVICE CLASSES

class MasterService(Job.Service):

    def __init__(self, sudo, memory):

        self.sudo = sudo
        self.memory = memory
        self.cores = multiprocessing.cpu_count()
        Job.Service.__init__(self, memory = self.memory, cores = self.cores)

    def start(self):
        """
        Start spark and hdfs master containers
        """
        log.write("start masters\n")
        log.flush()
        
        if (os.uname()[0] == "Darwin"):
            machine = check_output(["docker-machine", "ls"]).split("\n")[1].split()[0]
            self.IP = check_output(["docker-machine", "ip", machine]).strip().rstrip()
        else:
            self.IP = check_output(["hostname", "-f",])[:-1]

        self.sparkContainerID = docker_call(no_rm = True,
                                            work_dir = os.getcwd(),
                                            tool = "quay.io/ucsc_cgl/apache-spark-master:1.5.2",
                                            docker_parameters = ["--net=host",
                                                                 "-d",
                                                                 "-v", "/mnt/ephemeral/:/ephemeral/:rw",
                                                                 "-e", "SPARK_MASTER_IP="+self.IP,
                                                                 "-e", "SPARK_LOCAL_DIRS=/ephemeral/spark/local",
                                                                 "-e", "SPARK_WORKER_DIR=/ephemeral/spark/work"],
                                            tool_parameters = [],
                                            sudo = self.sudo,
                                            check_output = True)[:-1]
        self.hdfsContainerID = docker_call(no_rm = True,
                                           work_dir = os.getcwd(),
                                           tool = "quay.io/ucsc_cgl/apache-hadoop-master:2.6.2",
                                           docker_parameters = ["--net=host",
                                                                "-d"],
                                           tool_parameters = [self.IP],
                                           sudo = self.sudo,
                                           check_output = True)[:-1]
        return self.IP

    def stop(self):
        """
        Stop and remove spark and hdfs master containers
        """
        log.write("stop masters\n")
        log.flush()
        
        sudo = []
        if self.sudo:
            sudo = ["sudo"]

        call(sudo + ["docker", "exec", self.sparkContainerID, "rm", "-r", "/ephemeral/spark"])
        call(sudo + ["docker", "stop", self.sparkContainerID])
        call(sudo + ["docker", "rm", self.sparkContainerID])
        call(sudo + ["docker", "stop", self.hdfsContainerID])
        call(sudo + ["docker", "rm", self.hdfsContainerID])

        return

                
class WorkerService(Job.Service):
    
    def __init__(self, masterIP, sudo, memory):
        self.masterIP = masterIP
        self.sudo = sudo
        self.memory = memory
        self.cores = multiprocessing.cpu_count()
        Job.Service.__init__(self, memory = self.memory, cores = self.cores)

    def start(self):
        """
        Start spark and hdfs worker containers
        """
        log.write("start workers\n")
        log.flush()

        self.sparkContainerID = docker_call(no_rm = True,
                                            work_dir = os.getcwd(),
                                            tool = "quay.io/ucsc_cgl/apache-spark-worker:1.5.2",
                                            docker_parameters = ["--net=host", 
                                                                 "-d",
                                                                 "-v", "/mnt/ephemeral/:/ephemeral/:rw",
                                                                 "-e", "\"SPARK_MASTER_IP="+self.masterIP+":"+SPARK_MASTER_PORT+"\"",
                                                                 "-e", "SPARK_LOCAL_DIRS=/ephemeral/spark/local",
                                                                 "-e", "SPARK_WORKER_DIR=/ephemeral/spark/work"],
                                            tool_parameters = [self.masterIP+":"+SPARK_MASTER_PORT],
                                            sudo = self.sudo,
                                            check_output = True)[:-1]
        
        self.hdfsContainerID = docker_call(no_rm = True,
                                           work_dir = os.getcwd(),
                                           tool = "quay.io/ucsc_cgl/apache-hadoop-worker:2.6.2",
                                           docker_parameters = ["--net=host",
                                                                "-d",
                                                                "-v", "/mnt/ephemeral/:/ephemeral/:rw"],
                                           tool_parameters = [self.masterIP],
                                           sudo = self.sudo,
                                           check_output = True)[:-1]
        
        # fake do/while to check if HDFS is up
        hdfs_down = True
        retries = 0
        while hdfs_down and (retries < 5):

            sys.stderr.write("Sleeping 30 seconds before checking HDFS startup.")
            time.sleep(30)
            clusterID = ""
            try:
                clusterID = check_output(["docker",
                                          "exec",
                                          self.hdfsContainerID,
                                          "grep",
                                          "clusterID",
                                          "-R",
                                          "/opt/apache-hadoop/logs"])
            except:
                # grep returns a non-zero exit code if the pattern is not found
                # we expect to not find the pattern, so a non-zero code is OK
                pass

            if "Incompatible" in clusterID:
                sys.stderr.write("Hadoop Datanode failed to start with: %s" % clusterID)
                sys.stderr.write("Retrying container startup, retry #%d." % retries)
                retries += 1

                sys.stderr.write("Removing ephemeral hdfs directory.")
                check_call(["docker",
                            "exec",
                            self.hdfsContainerID,
                            "rm",
                            "-rf",
                            "/ephemeral/hdfs"])

                sys.stderr.write("Killing container %s." % self.hdfsContainerID)
                check_call(["docker",
                            "kill",
                            self.hdfsContainerID])

                # todo: this is copied code. clean up!
                sys.stderr.write("Restarting datanode.")
                self.hdfsContainerID = docker_call(no_rm = True,
                                                   work_dir = os.getcwd(),
                                                   tool = "quay.io/ucsc_cgl/apache-hadoop-worker:2.6.2",
                                                   docker_parameters = ["--net=host",
                                                                        "-d",
                                                                        "-v", "/mnt/ephemeral/:/ephemeral/:rw"],
                                                   tool_parameters = [self.masterIP],
                                                   sudo = self.sudo,
                                                   check_output = True)[:-1]

            else:
                sys.stderr.write("HDFS datanode started up OK!")
                hdfs_down = False

        if retries >= 5:
            raise RuntimeError("Failed %d times trying to start HDFS datanode." % retries)
                                   
        return

    def stop(self):
        """
        Stop spark and hdfs worker containers
        """
        log.write("stop workers\n")
        log.flush()

        sudo = []
        if self.sudo:
            sudo = ['sudo']

        call(sudo + ["docker", "exec", self.sparkContainerID, "rm", "-r", "/ephemeral/spark"])
        call(sudo + ["docker", "stop", self.sparkContainerID])
        call(sudo + ["docker", "rm", self.sparkContainerID])
        call(sudo + ["docker", "exec", self.hdfsContainerID, "rm", "-r", "/ephemeral/hdfs"])
        call(sudo + ["docker", "stop", self.hdfsContainerID])
        call(sudo + ["docker", "rm", self.hdfsContainerID])

        return


def build_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file_name', required = True,
                        help = "The full s3 url of the input SAM or BAM file")
    parser.add_argument('-n', '--num_nodes', type = int, required = True,
                        help = 'Number of nodes to use')
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

    inputs = {'numWorkers': options.num_nodes - 1,
              'outDir':     options.output_directory,
              'bamName':    options.input_file_name,
              'knownSNPs':  options.known_SNPs,
              'driverMemory': options.driver_memory,
              'executorMemory': options.executor_memory,
              'sudo': options.sudo,
              'suffix': None}

    Job.Runner.startToil(Job.wrapJobFn(start_master, inputs), options)

if __name__=="__main__":
    sys.exit(main(sys.argv[1:]))
