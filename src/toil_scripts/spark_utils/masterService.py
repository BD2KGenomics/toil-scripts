"""
Master service job for a Spark/HDFS cluster.

Starts up a Spark master and HDFS namenode. Return value of the service
job is the instance IP.

@author Audrey Musselman-Brown, almussel@ucsc.edu
@author Frank Austin Nothaft, fnothaft@berkeley.edu
"""

import logging
import multiprocessing
import os
import sys
import time

from toil.job import Job
from toil_scripts.batch_alignment.bwa_alignment import docker_call

class MasterService(Job.Service):

    def __init__(self, sudo, memory):

        self.sudo = sudo
        self.memory = memory
        self.cores = multiprocessing.cpu_count()
        self._log = logging.getLogger(__name__)
        Job.Service.__init__(self, memory = self.memory, cores = self.cores)

    def start(self):
        """
        Start spark and hdfs master containers
        """
        
        self.IP = check_output(["hostname", "-f",])[:-1]

        self._log.info("Started Spark master container.")
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
        self._log.info("Started HDFS Datanode.")
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
        
        sudo = []
        if self.sudo:
            sudo = ["sudo"]

        call(sudo + ["docker", "exec", self.sparkContainerID, "rm", "-r", "/ephemeral/spark"])
        call(sudo + ["docker", "stop", self.sparkContainerID])
        call(sudo + ["docker", "rm", self.sparkContainerID])
        self._log.info("Stopped Spark master.")

        call(sudo + ["docker", "stop", self.hdfsContainerID])
        call(sudo + ["docker", "rm", self.hdfsContainerID])
        self._log.info("Stopped HDFS datanode.")

        return
