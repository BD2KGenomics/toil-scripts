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
4 = Download Data
5 = ADAM Convert
6 = ADAM Transform
7 = Upload Data

================================================================================
:Dependencies
docker          - apt-get install docker (or 'docker.io' for linux)
toil            - pip install --pre toil
"""

import argparse
import os
from subprocess import call, check_call, check_output
import sys
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
    masterIP = job.addService(MasterService(inputs['sudo']))
    job.addChildJobFn(start_workers, masterIP, inputs)


def start_workers(job, masterIP, inputs):
    """
    Starts the worker services.
    """
    log.write("workers job\n")
    log.flush()
    for i in range(inputs['numWorkers']):
        job.addService(WorkerService(masterIP, inputs['sudo']))
    job.addFollowOnJobFn(download_data, masterIP, inputs)


def call_conductor(masterIP, inputs, src, dst):
    """
    Invokes the conductor container.
    """
    docker_call(no_rm = True,work_dir = os.getcwd(), tool = "quay.io/ucsc_cgl/conductor",
                docker_parameters = ["--net=host",
                                     "-e", "AWS_ACCESS_KEY="+inputs['accessKey'],
                                     "-e", "AWS_SECRET_KEY="+inputs['secretKey']],
                tool_parameters = ["--master", "spark://"+masterIP+":"+SPARK_MASTER_PORT,
                 "--conf", "spark.driver.memory=%s" % inputs["driverMemory"],
                 "--conf", "spark.executor.memory=%s" % inputs["executorMemory"],
                 "--", "-C", src, dst],
                sudo = inputs['sudo'])


def call_adam(inputs, masterIP, arguments):

    docker_call(no_rm = True,work_dir = os.getcwd(), tool = "quay.io/ucsc_cgl/adam:cd6ef41", 
                docker_parameter = ["--net=host"],
                tool_parameters = ["--master", "spark://"+masterIP+":"+SPARK_MASTER_PORT, 
                 "--conf", "spark.driver.memory=%s" % inputs["driverMemory"],
                 "--conf", "spark.executor.memory=%s" % inputs["executorMemory"],
                 "--conf", "spark.hadoop.fs.default.name=hdfs://%s:%s" % (masterIP, HDFS_MASTER_PORT),
                 "--"] + arguments,
                sudo = inputs['sudo'])


def remove_file(masterIP, filename):
    """
    Remove the given file from hdfs with master at the given IP address
    """
    containerID = check_output(["ssh", "-o", "StrictHostKeyChecking=no", masterIP, "docker", "ps", \
                                "|", "grep", "apache-hadoop-master", "|", "awk", "'{print $1}'"])[:-1]
    check_call(["ssh", "-o", "StrictHostKeyChecking=no", masterIP, "docker", "exec", containerID, \
                "/opt/apache-hadoop/bin/hdfs", "dfs", "-rm", "-r", "/"+filename])


def download_data(job, masterIP, inputs):
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
     
    job.addFollowOnJobFn(adam_convert, masterIP, hdfsBAM, hdfsSNPs, inputs)


def adam_convert(job, masterIP, inFile, snpFile, inputs):
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
 
    job.addFollowOnJobFn(adam_transform, masterIP, adamFile, adamSnpFile, inputs)


def adam_transform(job, masterIP, inFile, snpFile, inputs):
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

    job.addFollowOnJobFn(upload_data, masterIP, outFile, inputs)


def upload_data(job, masterIP, hdfsName, inputs):
    """
    Upload file hdfsName from hdfs to s3
    """
    log.write("write data\n")
    log.flush()

    fileSystem, path = hdfsName.split('://')
    nameOnly = path.split('/')[-1]

    call_conductor(masterIP, inputs, hdfsName, inputs['outDir']+"/"+nameOnly)
    
# SERVICE CLASSES

class MasterService(Job.Service):

    def __init__(self, sudo):

        Job.Service.__init__(self)
        self.sudo = sudo


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

        self.sparkContainerID = docker_call(no_rm = True,work_dir = os.getcwd(), tool = "quay.io/ucsc_cgl/apache-spark-master:1.5.2",
                                            docker_parameters = ["--net=host",
                                                                 "-d",
                                                                 "-v", "/mnt/ephemeral/:/ephemeral/:rw",
                                                                 "-e", "SPARK_MASTER_IP="+self.IP,
                                                                 "-e", "SPARK_LOCAL_DIRS=/ephemeral/spark/local",
                                                                 "-e", "SPARK_WORKER_DIR=/ephemeral/spark/work"],
                                            tool_parameters = [],
                                            sudo = self.sudo,
                                            check_output = True)[:-1]
        self.hdfsContainerID = docker_call(no_rm = True,work_dir = os.getcwd(), tool = "quay.io/ucsc_cgl/apache-hadoop-master:2.6.2",
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
    
    def __init__(self, masterIP, sudo):
        Job.Service.__init__(self)
        self.masterIP = masterIP
        self.sudo = sudo

    def start(self):
        """
        Start spark and hdfs worker containers
        """
        log.write("start workers\n")
        log.flush()

        self.sparkContainerID = docker_call(no_rm = True,work_dir = os.getcwd(), tool = "quay.io/ucsc_cgl/apache-spark-worker:1.5.2",
                                            docker_parameters = ["--net=host", 
                                                                 "-d",
                                                                 "-v", "/mnt/ephemeral/:/ephemeral/:rw",
                                                                 "-e", "\"SPARK_MASTER_IP="+self.masterIP+":"+SPARK_MASTER_PORT+"\"",
                                                                 "-e", "SPARK_LOCAL_DIRS=/ephemeral/spark/local",
                                                                 "-e", "SPARK_WORKER_DIR=/ephemeral/spark/work"],
                                            tool_parameters = [self.masterIP+":"+SPARK_MASTER_PORT],
                                            sudo = inputs['sudo'],
                                            check_output = True)[:-1]
        self.hdfsContainerID = docker_call(no_rm = True,work_dir = os.getcwd(), tool = "quay.io/ucsc_cgl/apache-hadoop-worker:2.6.2",
                                           docker_parameters = ["--net=host",
                                                                "-d",
                                                                "-v", "/mnt/ephemeral/:/ephemeral/:rw"],
                                           tool_parameters = [self.masterIP],
                                           sudo = inputs['sudo'],
                                           check_output = True)[:-1]
                                           
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
    parser.add_argument('-a', '--aws_access_key', required = True,
                        help = 'Amazon web services access key')
    parser.add_argument('-s', '--aws_secret_key', required = True,
                        help = 'Amazon web services secret key')
    parser.add_argument('-d', '--driver_memory', required = True,
                        help = 'Amount of memory to allocate for Spark Driver.')
    parser.add_argument('-q', '--executor_memory', required = True,
                        help = 'Amount of memory to allocate per Spark Executor.')
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
              'accessKey':  options.aws_access_key,
              'secretKey':  options.aws_secret_key,
              'driverMemory': options.driver_memory,
              'executorMemory': options.executor_memory,
              'sudo': options.sudo}

    Job.Runner.startToil(Job.wrapJobFn(start_master, inputs), options)

if __name__=="__main__":
    sys.exit(main(sys.argv[1:]))
