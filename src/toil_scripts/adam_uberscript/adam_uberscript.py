#!/usr/bin/env python2.7

import logging

# Initialize logging before the remaining imports to prevent those imported modules from snatching that one shot at 
# basicConfig.

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)-15s:%(levelname)s:%(name)s:%(message)s',
                    datefmt='%m-%d %H:%M:%S')

from StringIO import StringIO
import argparse
from collections import namedtuple
import csv
import os
import subprocess
import boto
from boto.exception import BotoServerError, EC2ResponseError
import boto.ec2.cloudwatch
import time
from uuid import uuid4
from tqdm import tqdm
import errno
from boto_lib import get_instance_ids
from datetime import datetime, timedelta
import threading
import boto.sdb
from boto.ec2 import connect_to_region

from automated_scaling import ClusterSize, Samples

metric_endtime_margin = timedelta(hours=1)
metric_initial_wait_period_in_seconds = 0
metric_collection_interval_in_seconds = 3600
metric_start_time_margin = 1800

scaling_initial_wait_period_in_seconds = 300
cluster_scaling_interval_in_seconds = 300

cluster_size_lock = threading.Lock()

aws_region = 'us-west-2'

def launch_cluster(params):
    """
    Launches a toil cluster with 1 worker, with shared dir S, of instance type I

    params: argparse.Namespace      Input arguments
    """
    log.info('Launching cluster of size: {} and type: {}'.format(1, params.instance_type))

    # if user provides a string to add to /etc/hosts, let's pass that through
    etc = []
    if params.add_to_etc_hosts:
        etc = ['-O', 'etc_hosts_entries=%s' % params.add_to_etc_hosts]

    subprocess.check_call(['cgcloud',
                           'create-cluster',
                           '--zone', '{0}a'.format(aws_region),
                           '--leader-instance-type', params.leader_type,
                           '--instance-type', params.instance_type,
                           '--num-workers', '1',
                           '--cluster-name', params.cluster_name,
                           '--leader-on-demand',
                           '--ssh-opts',
                           '"StrictHostKeyChecking=no"'] +
                          etc +
                          ['toil'])
    subprocess.check_call(['cgcloud',
                           'rsync',
                           '--zone', '{0}a'.format(aws_region),
                           '--cluster-name', params.cluster_name,
                           '--ssh-opts="-o StrictHostKeyChecking=no"',
                           'toil-leader',
                           '-a',
                           params.manifest_path, ':~/manifest'])
    subprocess.check_call(['cgcloud',
                           'rsync',
                           '--zone', '{0}a'.format(aws_region),
                           '--cluster-name', params.cluster_name,
                           '--ssh-opts="-o StrictHostKeyChecking=no"',
                           'toil-leader',
                           '-a',
                           params.share.rstrip('/'), ':'])

def place_boto_on_leader(params):
    log.info('Adding a .boto to leader to avoid credential timeouts.')
    subprocess.check_call(['cgcloud',
                           'rsync',
                           '--zone', '{0}a'.format(aws_region),
                           '--cluster-name', params.cluster_name,
                           '--ssh-opts="-o StrictHostKeyChecking=no"',
                           'toil-leader',
                           params.boto_path, ':~/.boto'])


def launch_pipeline(params):
    """
    Launches pipeline on toil-leader in a screen named the cluster run name

    params: argparse.Namespace      Input arguments
    """
    if not params.jobstore:
        jobstore = '{}-{}'.format(uuid4(), str(datetime.utcnow().date()))
    else:
        jobstore = params.jobstore
    restart = '--restart' if params.restart else ''
    log.info('Launching Pipeline and blocking. Check log.txt on leader for stderr and stdout')
    try:
        # Create screen session
        subprocess.check_call(['cgcloud',
                               'ssh',
                               '--zone', '{0}a'.format(aws_region),
                               '--cluster-name', params.cluster_name,
                               'toil-leader',
                               '-o', 'StrictHostKeyChecking=no',
                               'screen', '-dmS', params.cluster_name])

        # do we have a defined master ip?
        masterIP_arg = ''
        if params.master_ip:
            masterIP_arg = '--master_ip %s' % params.master_ip

        # Run command on screen session        
        if params.reference_genome == 'GRCh38':
            from toil_scripts.adam_uberscript.input_files import GRCh38_inputs as inputs
        elif params.reference_genome == 'hg19':
            from toil_scripts.adam_uberscript.input_files import hg19_inputs as inputs
        else:
            assert False, 'Invalid ref genome %s' % params.reference_genome

        pipeline_command = ('PYTHONPATH=$PYTHONPATH:~/toil-scripts/src python -m toil_scripts.adam_gatk_pipeline.align_and_call ' +
                            'aws:{region}:{j} ' +
                            '--autoscale_cluster ' +
                            '--sequence_dir {sequence_dir} ' +
                            '--retryCount 1 ' +
                            '--s3_bucket {b} ' +
                            '--bucket_region {region} ' +
                            '--uuid_manifest ~/manifest ' +
                            '--ref {ref} ' +
                            '--amb {amb} ' +
                            '--ann {ann} ' +
                            '--bwt {bwt} ' +
                            '--pac {pac} ' +
                            '--sa {sa} ' +
                            '--fai {fai} ')
        
        if 'alt' in inputs:
            pipeline_command += '--alt {alt} '
        
        pipeline_command += ('--use_bwakit ' +
                            '--num_nodes {s} ' +
                            '--driver_memory {m} ' +
                            '--executor_memory {m} ' +
                            '--phase {phase} ' +
                            '--mills {mills} ' +
                            '--dbsnp {dbsnp} ' +
                            '--omni {omni} ' +
                            '--hapmap {hapmap} ' +
                            '--batchSystem mesos ' +
                            '--mesosMaster $(hostname -i):5050 ' +
                            '--workDir /var/lib/toil ' +
                            '--file_size {fs} ' +
                            '--logInfo ' +
                            masterIP_arg +
                            '{r} 2>&1 | tee toil_output\n')

        pipeline_command = pipeline_command.format(j=jobstore,
                                                   b=params.bucket,
                                                   region=aws_region,
                                                   s=params.spark_nodes,
                                                   m=params.memory,
                                                   fs=params.file_size,
                                                   r=restart,
                                                   sequence_dir=params.sequence_dir,
                                                   **inputs)

        for chunk in [pipeline_command[i:i+500] for i in range(0, len(pipeline_command),500)]:
            subprocess.check_call(['cgcloud',
                                   'ssh',
                                   '--zone', '{0}a'.format(aws_region),
                                   '--cluster-name', params.cluster_name,
                                   'toil-leader',
                                   '-o', 'StrictHostKeyChecking=no',
                                   'screen', '-S', params.cluster_name,
                                   '-X', 'stuff', '"{0}"'.format(chunk)])
    
    except subprocess.CalledProcessError as e:
        log.info('Pipeline exited with non-zero status code: {}'.format(e))


def get_metric(cw, metric, instance_id, start, stop):
    """
    returns metric object associated with a paricular instance ID

    metric_name: str            Name of Metric to be Collected
    instance_id: str            Instance ID
    start: float                ISO format of UTC time start point
    stop: float                 ISO format of UTC time stop point
    :return: metric object
    """
    namespace, metric_name = metric.rsplit('/', 1)
    metric_object = cw.get_metric_statistics(namespace=namespace,
                                             metric_name=metric_name,
                                             dimensions={'InstanceId': instance_id},
                                             start_time=start,
                                             end_time=stop,
                                             period=300,
                                             statistics=['Average'])
    return metric_object


def get_cluster_size(cluster_name):
    """
    Returns the number of running toil-worker nodes
    """
    return len(list_workers(cluster_name))


def list_workers(cluster_name):
    """
    Returns list of dictionaries, each dictionary representing a worker node. Each dictinoary has the following keys:
    cluster_name, role_name, ordinal, cluster_ordinal, private_ip_address, ip_address, instance_id, instance_type,
    launch_time, state and zone
    """
    return parse_cgcloud_list_output(subprocess.check_output(['cgcloud', 'list', '-c', cluster_name, 'toil-worker']))


def parse_cgcloud_list_output(output):
    return [row for row in csv.DictReader(StringIO(output), delimiter='\t')]


def get_desired_cluster_size(conn, dom):
    
    nodes_per_sample = Samples.load(conn, dom)
    return sum(map(lambda x: x[1], nodes_per_sample.samples.values()))


def update_cluster_size(conn, dom, n):
    ClusterSize.change_size(conn, dom, n)


def grow_cluster(nodes, instance_type, cluster_name, etc):
    """
    Grow the cluster by n nodes
    """
    nodes_left = nodes
    while nodes_left > 0:
        log.info('Attempting to grow cluster by %i node(s) of type: %s', nodes_left, instance_type)

        output = ''
        cmd = (['cgcloud',
                'grow-cluster',
                '--list',
                '--instance-type', instance_type,
                '--num-workers', str(nodes_left),
                '--cluster-name', cluster_name] +
               etc +
               ['toil'])
        try:
            output = subprocess.check_output(cmd)
        except subprocess.CalledProcessError as cpe:
            log.warn('Running command %s returned with error code %d.',
                     ' '.join(cmd),
                     cpe.returncode)
            output = cpe.output

        added_nodes = len(parse_cgcloud_list_output(output))
        if added_nodes == 0:
            log.warn("Wasn't able to add any nodes. Wating 5 min.")
            time.sleep(5 * 60)
        assert added_nodes <= nodes_left
        nodes_left -= added_nodes
        log.info('Added %d node(s), %d node(s) left.', added_nodes, nodes_left)
    log.info('Successfully grew cluster by %i node(s) of type %s.', nodes, instance_type)


def manage_metrics_and_cluster_scaling(params):
    conn = boto.sdb.connect_to_region(aws_region)
    dom = conn.get_domain('{0}--files'.format(params.jobstore))
    grow_cluster_thread = threading.Thread(target=monitor_cluster_size, args=(params, conn, dom))
    metric_collection_thread = threading.Thread(target=collect_realtime_metrics, args=(params, conn, dom))
    grow_cluster_thread.start()
    metric_collection_thread.start()
    grow_cluster_thread.join()
    metric_collection_thread.join()

def monitor_cluster_size(params, conn, dom):
    """
    Monitors cluster size and grows it if the desired size is larger than the current size
    """
    
    # if user provides a string to add to /etc/hosts, let's pass that through
    etc = []
    if params.add_to_etc_hosts:
        etc = ['-O', 'etc_hosts_entries=%s' % params.add_to_etc_hosts]

    log.info('Cluster size monitor has started.')
    time.sleep(scaling_initial_wait_period_in_seconds)
    while True:
        size_check_time = time.time()
        # If cluster is too small, grow it
        cluster_size = get_cluster_size(params.cluster_name)
        desired_cluster_size = get_desired_cluster_size(conn, dom) + 1
        if cluster_size < desired_cluster_size:
            with cluster_size_lock:
                grow_cluster(desired_cluster_size - cluster_size, params.instance_type, params.cluster_name, etc)
                update_cluster_size(conn, dom, desired_cluster_size)

        # Sleep
        resize_time = time.time() - size_check_time
        log.info('Cluster is {} nodes as of {}'.format(get_cluster_size(params.cluster_name), resize_time))
        wait_time = cluster_scaling_interval_in_seconds - resize_time
        if wait_time > 0:
            time.sleep(wait_time)

# FIXME: unused parameters conn and dom

def collect_realtime_metrics(params, conn, dom, threshold=0.5, region='us-west-2'):
    """
    Collect metrics from AWS instances in 1 hour intervals.
    Instances that have gone idle (below threshold CPU value) are terminated.

    params: argparse.Namespace      Input arguments
    region: str                     AWS region metrics are being collected from
    uuid: str                       UUID of metric collection
    """
    list_of_metrics = ['AWS/EC2/CPUUtilization',
                       'CGCloud/MemUsage',
                       'CGCloud/DiskUsage_mnt_ephemeral',
                       'CGCloud/DiskUsage_root',
                       'AWS/EC2/NetworkIn',
                       'AWS/EC2/NetworkOut',
                       'AWS/EC2/DiskWriteOps',
                       'AWS/EC2/DiskReadOps']

    # Create output directory
    uuid = str(uuid4())
    date = str(datetime.utcnow().date())
    dir_path = '{}_{}_{}'.format(params.cluster_name, uuid, date)
    mkdir_p(dir_path)

    start = time.time() - metric_start_time_margin

    # Create connections to ec2 and cloudwatch
    conn = boto.ec2.connect_to_region(region)
    cw = boto.ec2.cloudwatch.connect_to_region(region)
    sdbconn = boto.sdb.connect_to_region(region)
    dom = sdbconn.get_domain('{0}--files'.format(params.jobstore))

    # Create initial variables
    start = datetime.utcfromtimestamp(start)
    DataPoint = namedtuple('datapoint', ['instance_id', 'value', 'timestamp'])
    timestamps = {}
    # Begin loop
    log.info('Metric collection has started. '
             'Waiting {} seconds before initial collection.'.format(metric_initial_wait_period_in_seconds))
    time.sleep(metric_initial_wait_period_in_seconds)
    
    while True:
        # FIXME: why doesn't filter_cluster=params.cluster_name work?
        ids = get_instance_ids(filter_name=params.namespace.strip('/').rstrip('/') + '_toil-worker')
        if not ids:
            break
        metric_collection_time = time.time()
        try:
            for instance_id in tqdm(ids):
                idle = False
                for metric in list_of_metrics:
                    datapoints = []
                    aws_start = timestamps.get(instance_id, start)
                    aws_stop = datetime.utcnow() + metric_endtime_margin
                    metric_object = get_metric(cw, metric, instance_id, aws_start, aws_stop)
                    for datum in metric_object:
                        d = DataPoint(instance_id=instance_id, value=datum['Average'], timestamp=datum['Timestamp'])
                        datapoints.append(d)
                    # Save data in local directory
                    if datapoints:
                        datapoints = sorted(datapoints, key=lambda x: x.timestamp)
                        with open(os.path.join(dir_path, '{}.csv'.format(os.path.basename(metric))), 'a') as f:
                            writer = csv.writer(f, delimiter='\t')
                            writer.writerows(datapoints)
                    # Check if instance's CPU has been idle the last 30 minutes.
                    if metric == 'AWS/EC2/CPUUtilization':
                        averages = [x.value for x in sorted(datapoints, key=lambda x: x.timestamp)][-6:]
                        # If there is at least 30 minutes of data points and max is below threshold, flag to be killed.
                        if len(averages) == 6:
                            if max(averages) < threshold:
                                idle = True
                                log.info('Flagging {} to be killed. '
                                         'Max CPU {} for last 30 minutes.'.format(instance_id, max(averages)))
                # Kill instance if idle and cluster is too large
                if idle:
                    try:
                        with cluster_size_lock:
                            cluster_size = get_cluster_size(params.cluster_name)
                            if cluster_size > get_desired_cluster_size(sdbconn, dom):
                                log.info('Terminating Instance: {}'.format(instance_id))
                                log.info('Killing instance {0}\n'.format(instance_id))
                                conn.terminate_instances(instance_ids=[instance_id])
                                update_cluster_size(sdbconn, dom, cluster_size - 1)
                    except (EC2ResponseError, BotoServerError) as e:
                        log.info('Error terminating instance: {}\n{}'.format(instance_id, e))
                # Set start point to be last collected timestamp
                timestamps[instance_id] = max(x.timestamp for x in datapoints) if datapoints else start
        except BotoServerError:
            log.error('Giving up trying to fetch metric for this interval')

        # Sleep
        collection_time = time.time() - metric_collection_time
        log.info('Metric collection took: {} seconds. Waiting one hour.'.format(collection_time))
        wait_time = metric_collection_interval_in_seconds - collection_time
        if wait_time < 0:
            log.warning('Collection time exceeded metric collection interval by: %i', -wait_time)
        else:
            time.sleep(wait_time)

    log.info('Metric collection has finished.')


def mkdir_p(path):
    """
    It is Easier to Ask for Forgiveness than Permission
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def main():
    """
    Modular script for running toil pipelines
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # Launch Cluster
    parser_cluster = subparsers.add_parser('launch-cluster', help='Launches AWS cluster via CGCloud')
    parser_cluster.add_argument('-c', '--cluster-name', required=True, help='Name of cluster.')
    parser_cluster.add_argument('-S', '--share', required=True,
                                help='Full path to directory: pipeline script, launch script, and master key.')
    parser_cluster.add_argument('-t', '--instance-type', default='r3.8xlarge',
                                help='Worker instance type. e.g.  m4.large or c3.8xlarge.')
    parser_cluster.add_argument('-T', '--leader-type', default='m3.medium', help='Sets leader instance type.')
    parser_cluster.add_argument('-b', '--boto-path', default='/home/mesosbox/.boto', type=str,
                                help='Path to local .boto file to be placed on leader.')
    parser_cluster.add_argument('-M', '--manifest-path', required=True, help='Path to manifest file.')
    parser_cluster.add_argument('-etc', '--add-to-etc-hosts', default=None,
                                required=False, help='Optional entry to add to /etc/hosts')

    # Launch Pipeline
    parser_pipeline = subparsers.add_parser('launch-pipeline', help='Launches pipeline')
    parser_pipeline.add_argument('-c', '--cluster-name', required=True, help='Name of cluster.')
    parser_pipeline.add_argument('-j', '--jobstore', default=None,
                                 help='Name of jobstore. Defaults to UUID-Date if not set')
    parser_pipeline.add_argument('--restart', default=None, action='store_true',
                                 help='Attempts to restart pipeline, requires existing jobstore.')
    parser_pipeline.add_argument('--master_ip', default=None, help = 'Spark Master IP.')
    parser_pipeline.add_argument('-B', '--bucket', help='Set destination bucket.')
    parser_pipeline.add_argument('-m', '--memory', default='200g', help='The memory per worker node in GB') 
    parser_pipeline.add_argument('-f', '--file_size', default='100G', help='Approximate size of the BAM files')
    parser_pipeline.add_argument('-s', '--spark_nodes', default='9', help='The number of nodes needed for the spark cluster')
    parser_pipeline.add_argument('-SD', '--sequence_dir',
                                 help = 'Directory where raw sequences are.',
                                 default = 'sequence')

    parser_pipeline.add_argument('-R', '--reference_genome', required=True, choices=['GRCh38','hg19'], help='Reference Genome to align and call against.  Choose between GRCh38 and hg19')

    # Launch Metric Collection
    parser_metric = subparsers.add_parser('launch-metrics', help='Launches metric collection thread')
    parser_metric.add_argument('-c', '--cluster-name', required=True, help='Name of cluster')
    parser_metric.add_argument('-j', '--jobstore', required=True, help='Name of jobstore')
    parser_metric.add_argument('--namespace', default='jtvivian', help='CGCloud NameSpace')
    parser_metric.add_argument('-t', '--instance-type', default='r3.8xlarge',
                               help='Worker instance type. e.g.  m4.large or c3.8xlarge.')
    parser_metric.add_argument('-etc', '--add-to-etc-hosts', default=None,
                               required=False, help='Optional entry to add to /etc/hosts')
    
    # Parse args
    params = parser.parse_args()

    # Modular Run Sequence
    if params.command == 'launch-cluster':
        launch_cluster(params)
        place_boto_on_leader(params)
    elif params.command == 'launch-pipeline': 
        launch_pipeline(params)
    elif params.command == 'launch-metrics':
        manage_metrics_and_cluster_scaling(params)


if __name__ == '__main__':
    main()
