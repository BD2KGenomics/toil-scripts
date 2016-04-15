"""
Spawns a Spark cluster backed by HDFS.

@author Audrey Musselman-Brown, almussel@ucsc.edu
@author Frank Austin Nothaft, fnothaft@berkeley.edu
"""

import logging
import multiprocessing
import os
from subprocess import call, check_call, check_output
import time

from toil.job import Job
from toil_scripts.lib.programs import docker_call

_log = logging.getLogger(__name__)

def start_spark_hdfs_cluster(job,
                             numWorkers,
                             executorMemory,
                             sudo,
                             jFn,
                             jArgs = [],
                             jCores = None,
                             jMemory = None,
                             jDisk = None):

    # build job requirement dictionary
    jReqs = {}
    if jCores:
        jReqs['cores'] = jCores
    if jMemory:
        jReqs['memory'] = jMemory
    if jDisk:
        jReqs['disk'] = jDisk

    masterIP = start_spark_hdfs_master(job, executorMemory, sudo)
    job.addChildJobFn(start_spark_hdfs_workers,
                      masterIP,
                      numWorkers,
                      executorMemory,
                      sudo,
                      jFn,
                      jArgs,
                      jReqs)


def start_spark_hdfs_master(job, executorMemory, sudo):
    """
    Starts the master service.
    """

    _log.info("Starting Spark master and HDFS namenode.")

    masterIP = job.addService(MasterService(sudo, "%s G" % executorMemory))

    _log.info("Spark Master and HDFS Namenode started.")

    return masterIP

def start_spark_hdfs_workers(job, masterIP, numWorkers, executorMemory, sudo, jFn, jArgs, jReqs):
    """
    Starts the worker services.
    """
    _log.info("Starting %d Spark workers and HDFS Datanodes.", numWorkers)

    for i in range(numWorkers):
        job.addService(WorkerService(masterIP, sudo, "%s G" % executorMemory))

    job.addChildJobFn(jFn, masterIP, *jArgs, **jReqs)

class MasterService(Job.Service):

    def __init__(self, sudo, memory):

        self.sudo = sudo
        self.memory = memory
        self.cores = multiprocessing.cpu_count()
        Job.Service.__init__(self, memory = self.memory, cores = self.cores)

    def start(self, fileStore):
        """
        Start spark and hdfs master containers

        fileStore: Unused
        """
        
        self.IP = check_output(["hostname", "-f",])[:-1]

        _log.info("Started Spark master container.")
        self.sparkContainerID = docker_call(tool = "quay.io/ucsc_cgl/apache-spark-master:1.5.2",
                                            docker_parameters = ["--net=host",
                                                                 "-d",
                                                                 "-v", "/mnt/ephemeral/:/ephemeral/:rw",
                                                                 "-e", "SPARK_MASTER_IP="+self.IP,
                                                                 "-e", "SPARK_LOCAL_DIRS=/ephemeral/spark/local",
                                                                 "-e", "SPARK_WORKER_DIR=/ephemeral/spark/work"],
                                            rm=False,
                                            sudo = self.sudo,
                                            check_output = True,
                                            mock = False)[:-1]
        _log.info("Started HDFS Datanode.")
        self.hdfsContainerID = docker_call(tool = "quay.io/ucsc_cgl/apache-hadoop-master:2.6.2",
                                           docker_parameters = ["--net=host",
                                                                "-d"],
                                           parameters = [self.IP],
                                           rm=False,
                                           sudo = self.sudo,
                                           check_output = True,
                                           mock = False)[:-1]
        return self.IP


    def stop(self, fileStore):
        """
        Stop and remove spark and hdfs master containers

        fileStore: Unused
        """
        
        sudo = []
        if self.sudo:
            sudo = ["sudo"]

        call(sudo + ["docker", "exec", self.sparkContainerID, "rm", "-r", "/ephemeral/spark"])
        call(sudo + ["docker", "stop", self.sparkContainerID])
        call(sudo + ["docker", "rm", self.sparkContainerID])
        _log.info("Stopped Spark master.")

        call(sudo + ["docker", "stop", self.hdfsContainerID])
        call(sudo + ["docker", "rm", self.hdfsContainerID])
        _log.info("Stopped HDFS namenode.")

        return


    def check(self):
        """
        Checks to see if Spark master and HDFS namenode are still running.
        """

        containers = check_output(["docker", "ps", "-q"])

        return ((self.sparkContainerID in containers) and
                (self.hdfsContainerID in containers))


SPARK_MASTER_PORT = "7077"

class WorkerService(Job.Service):
    
    def __init__(self, masterIP, sudo, memory):
        self.masterIP = masterIP
        self.sudo = sudo
        self.memory = memory
        self.cores = multiprocessing.cpu_count()
        Job.Service.__init__(self, memory = self.memory, cores = self.cores)

    def start(self, fileStore):
        """
        Start spark and hdfs worker containers

        fileStore: Unused
        """
        # start spark and our datanode
        self.sparkContainerID = docker_call(tool = "quay.io/ucsc_cgl/apache-spark-worker:1.5.2",
                                            docker_parameters = ["--net=host", 
                                                                 "-d",
                                                                 "-v", "/mnt/ephemeral/:/ephemeral/:rw",
                                                                 "-e", "\"SPARK_MASTER_IP="+self.masterIP+":"+SPARK_MASTER_PORT+"\"",
                                                                 "-e", "SPARK_LOCAL_DIRS=/ephemeral/spark/local",
                                                                 "-e", "SPARK_WORKER_DIR=/ephemeral/spark/work"],
                                            parameters = [self.masterIP+":"+SPARK_MASTER_PORT],
                                            rm=False,
                                            sudo = self.sudo,
                                            check_output = True,
                                            mock = False)[:-1]
        self.__start_datanode()
        
        # fake do/while to check if HDFS is up
        hdfs_down = True
        retries = 0
        while hdfs_down and (retries < 5):

            _log.info("Sleeping 30 seconds before checking HDFS startup.")
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
                _log.warning("Hadoop Datanode failed to start with: %s", clusterID)
                _log.warning("Retrying container startup, retry #%d.", retries)
                retries += 1

                _log.warning("Removing ephemeral hdfs directory.")
                check_call(["docker",
                            "exec",
                            self.hdfsContainerID,
                            "rm",
                            "-rf",
                            "/ephemeral/hdfs"])

                _log.warning("Killing container %s.", self.hdfsContainerID)
                check_call(["docker",
                            "kill",
                            self.hdfsContainerID])

                # todo: this is copied code. clean up!
                _log.info("Restarting datanode.")
                self.__start_datanode()

            else:
                _log.info("HDFS datanode started up OK!")
                hdfs_down = False

        if retries >= 5:
            raise RuntimeError("Failed %d times trying to start HDFS datanode." % retries)
                                   
        return

    def __start_datanode(self):
        """
        Launches the Hadoop datanode.
        """
        self.hdfsContainerID = docker_call(tool = "quay.io/ucsc_cgl/apache-hadoop-worker:2.6.2",
                                           docker_parameters = ["--net=host",
                                                                "-d",
                                                                "-v", "/mnt/ephemeral/:/ephemeral/:rw"],
                                           parameters = [self.masterIP],
                                           rm=False,
                                           sudo = self.sudo,
                                           check_output = True,
                                           mock = False)[:-1]


    def stop(self, fileStore):
        """
        Stop spark and hdfs worker containers

        fileStore: Unused
        """

        sudo = []
        if self.sudo:
            sudo = ['sudo']

        call(sudo + ["docker", "exec", self.sparkContainerID, "rm", "-r", "/ephemeral/spark"])
        call(sudo + ["docker", "stop", self.sparkContainerID])
        call(sudo + ["docker", "rm", self.sparkContainerID])
        _log.info("Stopped Spark worker.")

        call(sudo + ["docker", "exec", self.hdfsContainerID, "rm", "-r", "/ephemeral/hdfs"])
        call(sudo + ["docker", "stop", self.hdfsContainerID])
        call(sudo + ["docker", "rm", self.hdfsContainerID])
        _log.info("Stopped HDFS datanode.")

        return


    def check(self):
        """
        Checks to see if Spark master and HDFS namenode are still running.
        """

        containers = check_output(["docker", "ps", "-q"])

        return ((self.sparkContainerID in containers) and
                (self.hdfsContainerID in containers))
