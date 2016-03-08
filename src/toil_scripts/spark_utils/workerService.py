"""
Worker service job for a Spark/HDFS cluster.

Starts up a Spark worker and HDFS datanode. Associates with a provided Spark
master and HDFS namenode, which are presumed to be running at the same master
IP, and on the default Spark (7077) and HDFS (8020) master ports.

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

SPARK_MASTER_PORT = "7077"

class WorkerService(Job.Service):
    
    def __init__(self, masterIP, sudo, memory):
        self.masterIP = masterIP
        self.sudo = sudo
        self.memory = memory
        self.cores = multiprocessing.cpu_count()
        self._log = logging.getLogger(__name__)
        Job.Service.__init__(self, memory = self.memory, cores = self.cores)

    def start(self):
        """
        Start spark and hdfs worker containers
        """

        # start spark and our datanode
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
        __start_datanode()
        
        # fake do/while to check if HDFS is up
        hdfs_down = True
        retries = 0
        while hdfs_down and (retries < 5):

            self._log.info("Sleeping 30 seconds before checking HDFS startup.")
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
                self._log.warning("Hadoop Datanode failed to start with: %s", clusterID)
                self._log.warning("Retrying container startup, retry #%d.", retries)
                retries += 1

                self._log.warning("Removing ephemeral hdfs directory.")
                check_call(["docker",
                            "exec",
                            self.hdfsContainerID,
                            "rm",
                            "-rf",
                            "/ephemeral/hdfs"])

                self._log.warning("Killing container %s.", self.hdfsContainerID)
                check_call(["docker",
                            "kill",
                            self.hdfsContainerID])

                # todo: this is copied code. clean up!
                self._log.info("Restarting datanode.")
                __start_datanode()

            else:
                self._log.info("HDFS datanode started up OK!")
                hdfs_down = False

        if retries >= 5:
            raise RuntimeError("Failed %d times trying to start HDFS datanode." % retries)
                                   
        return

    def __start_datanode(self):
        """
        Launches the Hadoop datanode.
        """
        self.hdfsContainerID = docker_call(no_rm = True,
                                           work_dir = os.getcwd(),
                                           tool = "quay.io/ucsc_cgl/apache-hadoop-worker:2.6.2",
                                           docker_parameters = ["--net=host",
                                                                "-d",
                                                                "-v", "/mnt/ephemeral/:/ephemeral/:rw"],
                                           tool_parameters = [self.masterIP],
                                           sudo = self.sudo,
                                           check_output = True)[:-1]


    def stop(self):
        """
        Stop spark and hdfs worker containers
        """

        sudo = []
        if self.sudo:
            sudo = ['sudo']

        call(sudo + ["docker", "exec", self.sparkContainerID, "rm", "-r", "/ephemeral/spark"])
        call(sudo + ["docker", "stop", self.sparkContainerID])
        call(sudo + ["docker", "rm", self.sparkContainerID])
        self._log.info("Stopped Spark worker.")

        call(sudo + ["docker", "exec", self.hdfsContainerID, "rm", "-r", "/ephemeral/hdfs"])
        call(sudo + ["docker", "stop", self.hdfsContainerID])
        call(sudo + ["docker", "rm", self.hdfsContainerID])
        self._log.info("Stopped HDFS namenode.")

        return
