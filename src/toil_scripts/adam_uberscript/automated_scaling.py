# Audrey Musselman-Brown
#

import time
import boto.sdb
import boto.exception
import threading
from collections import OrderedDict
import sys

class Samples(object):

    def __init__(self, conn, dom, version, samples):
        self.conn = conn
        self.dom = dom
        self.version = version
        self.samples = samples

    @classmethod
    def load(cls, conn, dom):
        nodes_per_sample = dom.get_attributes("nodes_per_sample", consistent_read=True)
        if 'version' in nodes_per_sample:
            version = int(nodes_per_sample.pop('version'))
        else:
            version = 0
        samples = OrderedDict(map(lambda x:[x[0], map(int, x[1].split(','))],
                                  sorted(nodes_per_sample.items(),
                                         key=(lambda x: x[1])))) 
        return cls(conn, dom, version, samples)
    
    def save(self):
        attributes = dict({key: "{0},{1}".format(value[0], value[1])
                           for key, value in self.samples.items()},
                           version=str(self.version + 1))
        try:
            self.dom.put_attributes("nodes_per_sample", attributes,
                                    expected_value=('version', str(self.version) if self.version else False))

            self.samples = OrderedDict(sorted(self.samples.items(),
                                              key=(lambda x: x[1])))
            self.version += 1
            return True
        except boto.exception.SDBResponseError as e:
            if e.error_code == 'ConditionalCheckFailed':
                return False
            else:
                raise

    @classmethod
    def increase_nodes(cls, sampleID, n):
        """
        Increases the node count for sampleID by n and returns once the nodes
        become available
        This method can ONLY be called from within a toil script
        """
        
        aws, region, domain = sys.argv[1].split(':')
        conn = boto.sdb.connect_to_region(region)
        dom = conn.get_domain("{0}--files".format(domain))

        while True:
            # load nodes per sample from sdb
            nodes_per_sample = cls.load(conn, dom)

            # update or insert new nodes
            if sampleID in nodes_per_sample.samples:
                nodes_per_sample.samples[sampleID][0] = nodes_per_sample.version + 1
                nodes_per_sample.samples[sampleID][1] += n
            else:
                nodes_per_sample.samples[sampleID] = [nodes_per_sample.version + 1, n]

            # attempt to write back to sdb; retry if failed
            if nodes_per_sample.save():
                break

        while True:
            # get the current cluster size
            cluster_size = ClusterSize.load(conn, dom)

            sum = 0
            for sample, [sample_version, nodes] in nodes_per_sample.samples.items():
                sum += nodes
                if sample == sampleID and sum <= cluster_size.size:    
                    return
                if sum > cluster_size.size:
                    break
            nodes_per_sample = cls.load(conn,dom)
            time.sleep(2)

    @classmethod
    def decrease_nodes(cls, sampleID, n):
        """
        Decreases the node count for sampleID by n
        """

        aws, region, domain = sys.argv[1].split(':')
        conn = boto.sdb.connect_to_region(region)
        dom = conn.get_domain("{0}--files".format(domain))

        while True:
            nodes_per_sample = cls.load(conn, dom)
            nodes_per_sample.samples[sampleID][1] -= n
            print "decrease", nodes_per_sample.samples
            nodes_per_sample.samples[sampleID][1] = min(nodes_per_sample.samples[sampleID], 0)
            if nodes_per_sample.save():
                break


class ClusterSize(object):

    def __init__(self, conn, dom, version, size):
        self.conn = conn
        self.dom = dom
        self.version = version
        self.size = size

    @classmethod
    def load(cls, conn, dom):
        cluster_size = dom.get_attributes("cluster_size", consistent_read=True)
        if 'version' in cluster_size:
            version = int(cluster_size['version'])
        else:
            version = 0
        if 'size' in cluster_size:
            size = int(cluster_size['size'])
        else:
            size = 0
        return cls(conn, dom, version, size)


    def save(self):
        attributes = dict(size=str(self.size), version=str(self.version+1))
        try:
            self.dom.put_attributes("cluster_size", attributes,
                                    expected_value=('version', str(self.version) if self.version else False))
            return True
        except boto.exception.SDBResponseError as e:
            if e.error_code == 'ConditionalCheckFailed':
                return False
            else:
                raise

    @classmethod
    def change_size(cls, conn, dom, n):
        while True:
            cluster_size = cls.load(conn, dom)
            cluster_size.size = n
            if cluster_size.save():
                break


def test_cluster_size(conn, dom):
    
    i=0
    while i < 7:
        cluster_size = ClusterSize.load(conn, dom)
        nodes_per_sample = Samples.load(conn, dom)        
        node_requests = sum(map(lambda (v,n): n, nodes_per_sample.samples.values()))
        if cluster_size.size != node_requests:
            print "changing cluster size from", cluster_size.size, "to", node_requests
            i += 1
            cluster_size.change_size(conn, dom, node_requests)


def test_samples(conn, dom):
        
    samples_and_nodes = Samples.load(conn, dom)
    cluster_size = ClusterSize.load(conn, dom)
    print "samples and cluster size loaded"
    print "initial nodes", samples_and_nodes.samples, samples_and_nodes.version
    print "initial cluster size", cluster_size.size, cluster_size.version

    print "~~~~~~~~~~~~~~increase  sample 2 nodes by 1~~~~~~~~~~~~~~~~~~~~"
    samples_and_nodes.increase_nodes("2", 1)
    samples_and_nodes = Samples.load(conn, dom)
    print "updated:", samples_and_nodes.samples, samples_and_nodes.version

    print "~~~~~~~~~~~~~~increase  sample 1 nodes by 1~~~~~~~~~~~~~~~~~~~~"
    samples_and_nodes.increase_nodes("1", 1)
    samples_and_nodes = Samples.load(conn, dom)    
    print "updated:", samples_and_nodes.samples, samples_and_nodes.version

    print "~~~~~~~~~~~~~~increase  sample 3 nodes by 1~~~~~~~~~~~~~~~~~~~~"
    samples_and_nodes.increase_nodes("3", 1)
    samples_and_nodes = Samples.load(conn, dom)
    print "updated:", samples_and_nodes.samples, samples_and_nodes.version

    print "~~~~~~~~~~~~~~increase  sample 1 nodes by 2~~~~~~~~~~~~~~~~~~~~"
    samples_and_nodes.increase_nodes("1", 2)
    samples_and_nodes = Samples.load(conn, dom)
    print "updated:", samples_and_nodes.samples, samples_and_nodes.version

    print "~~~~~~~~~~~~~~increase  sample 2 nodes by 3~~~~~~~~~~~~~~~~~~~~"
    samples_and_nodes.increase_nodes("2", 3)
    samples_and_nodes = Samples.load(conn, dom)
    print samples_and_nodes.samples, samples_and_nodes.version

    print "~~~~~~~~~~~~~~decrease  sample 1 nodes by 3~~~~~~~~~~~~~~~~~~~~"
    samples_and_nodes.decrease_nodes("1", 3)
    samples_and_nodes = Samples.load(conn, dom)
    print samples_and_nodes.samples, samples_and_nodes.version

    print "~~~~~~~~~~~~~~decrease  sample 2 nodes by 4~~~~~~~~~~~~~~~~~~~~"
    samples_and_nodes.decrease_nodes("2", 4)
    samples_and_nodes = Samples.load(conn, dom)
    print samples_and_nodes.samples, samples_and_nodes.version


if __name__=="__main__":

    aws, region, domain = sys.argv[1].split(':')
    conn = boto.sdb.connect_to_region(region)
    dom = conn.create_domain("{0}--files".format(domain))

    dom.delete_attributes("cluster_size")
    dom.delete_attributes("nodes_per_sample")
    print "attributes deleted"

    #test_cluster_size(conn,dom)
    cluster_size_thread = threading.Thread(target=test_cluster_size, args=(conn, dom))
    samples_thread = threading.Thread(target=test_samples, args=(conn, dom))
    cluster_size_thread.start()
    samples_thread.start()
    cluster_size_thread.join()
    samples_thread.join()

    conn.delete_domain("{0}--files".format(domain))
    sys.exit()
