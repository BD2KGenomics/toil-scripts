#!/usr/bin/env python2.7

# Initialize logging before the remaining imports to prevent those imported modules from snatching that one shot at 
# basicConfig.

import logging
from contextlib import contextmanager
from subprocess import check_call, CalledProcessError, check_output

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)-15s:%(levelname)s:%(name)s:%(message)s',
                    datefmt='%m-%d %H:%M:%S')

import argparse
import csv
import errno
import os
import threading
import time
from StringIO import StringIO
from collections import namedtuple
from datetime import datetime, timedelta
from pipes import quote
from uuid import uuid4

import boto
import boto.ec2.cloudwatch
import boto.sdb
from automated_scaling import ClusterSize, Samples, Semaphore, SparkMasterAddress
from toil_scripts.lib.programs import mock_mode
from boto.ec2 import connect_to_region
from boto.exception import BotoServerError, EC2ResponseError
from boto_lib import get_instance_ids
from tqdm import tqdm

metric_endtime_margin = timedelta(hours=1)
metric_initial_wait_period_in_seconds = 0
metric_collection_interval_in_seconds = 3600
metric_start_time_margin = 1800

scaling_initial_wait_period_in_seconds = 60
cluster_scaling_interval_in_seconds = 150

# Protects against concurrent changes to the size of the Toil cluster, by both the metric thread terminating idle
# nodes and the cluster-growing thread adding nodes.
cluster_size_lock = threading.RLock()

def launch_cluster(params):
    """
    Launches a toil cluster with 1 worker, with shared dir S, of instance type I

    :param argparse.Namespace params: parsed command line arguments and options 
    """
    log.info('Launching cluster of size: {} and type: {}'.format(1, params.instance_type))

    check_call(['cgcloud',
                'create-cluster',
                '--zone', params.zone,
                '--cluster-name', params.cluster_name,
                '--leader-instance-type', params.leader_type,
                '--instance-type', params.instance_type,
                '--num-workers', '1',
                '--ssh-opts',
                '"StrictHostKeyChecking=no"'] +
               role_options(params) +
               spot_options(params) +
               ['toil'])
    check_call(['cgcloud',
                'rsync',
                '--zone', params.zone,
                '--cluster-name', params.cluster_name,
                '--ssh-opts="-o StrictHostKeyChecking=no"',
                'toil-leader',
                '-a',
                params.manifest_path, ':~/manifest'])
    check_call(['cgcloud',
                'rsync',
                '--zone', params.zone,
                '--cluster-name', params.cluster_name,
                '--ssh-opts="-o StrictHostKeyChecking=no"',
                'toil-leader',
                '-a',
                params.share.rstrip('/'), ':'])


def place_boto_on_leader(params):
    log.info('Adding a .boto to leader to avoid credential timeouts.')
    check_call(['cgcloud',
                'rsync',
                '--zone', params.zone,
                '--cluster-name', params.cluster_name,
                '--ssh-opts="-o StrictHostKeyChecking=no"',
                'toil-leader',
                params.boto_path, ':~/.boto'])


def launch_pipeline(params):
    """
    Launches pipeline in a screen session on toil-leader. 

    :param argparse.Namespace params: parsed command line arguments and options 
    """
    if not params.jobstore:
        jobstore = '{}-{}'.format(uuid4(), str(datetime.utcnow().date()))
    else:
        jobstore = params.jobstore
    restart = '--restart' if params.restart else ''
    log.info('Launching Pipeline and blocking. Check log.txt on leader for stderr and stdout')
    try:
        # Create screen session
        check_call(['cgcloud',
                    'ssh',
                    '--zone', params.zone,
                    '--cluster-name', params.cluster_name,
                    'toil-leader',
                    '-o', 'StrictHostKeyChecking=no',
                    'screen', '-dmS', params.cluster_name])

        if params.reference_genome == 'GRCh38':
            from toil_scripts.adam_uberscript.input_files import GRCh38_inputs as inputs
        elif params.reference_genome == 'hg19':
            from toil_scripts.adam_uberscript.input_files import hg19_inputs as inputs
        else:
            assert False, 'Invalid ref genome %s' % params.reference_genome

        # Assemble pipeline command to be stuffed into a screen session
        
        pipeline_command = ['PYTHONPATH=$PYTHONPATH:~/toil-scripts/src',
                            'python -m toil_scripts.adam_gatk_pipeline.align_and_call',
                            'aws:{region}:{j}',
                            '--autoscale_cluster',
                            '--retryCount 1',
                            '--use_bwakit',
                            '--driver_memory {m}',
                            '--executor_memory {m}',
                            '--batchSystem mesos',
                            '--mesosMaster $(hostname -i):5050',
                            '--workDir /var/lib/toil',
                            '--logInfo']
        if mock_mode():
            pipeline_command = ['ADAM_GATK_MOCK_MODE=1'] + \
                               pipeline_command + \
                               ['--dir_suffix /mock']
        else:    
            pipeline_command += ['--s3_bucket {b}',
                                 '--bucket_region {region}',
                                 '--sequence_dir {sequence_dir}',
                                 '--dir_suffix /{genome}',
                                 '--uuid_manifest ~/manifest',
                                 '--ref {ref}',
                                 '--amb {amb}',
                                 '--ann {ann}',
                                 '--bwt {bwt}',
                                 '--pac {pac}',
                                 '--sa {sa}',
                                 '--fai {fai}',
                                 '--phase {phase}',
                                 '--mills {mills}',
                                 '--dbsnp {dbsnp}',
                                 '--omni {omni}',
                                 '--hapmap {hapmap}',
                                 '--file_size {fs}']

            if 'alt' in inputs:
                pipeline_command.append('--alt {alt}')

        # Do we have a defined master IP?
        if params.master_ip:
            pipeline_command.append('--master_ip %s' % params.master_ip)
        elif params.spark_nodes:
            pipeline_command.append('--num_nodes %s' % params.spark_nodes)

        pipeline_command.append('{r} 2>&1 | tee toil_output\n')

        pipeline_command = ' '.join(pipeline_command)
        pipeline_command = pipeline_command.format(j=jobstore,
                                                   b=params.bucket,
                                                   region=region_of_zone(params.zone),
                                                   m=params.memory,
                                                   fs=params.file_size,
                                                   r=restart,
                                                   sequence_dir=params.sequence_dir,
                                                   genome=params.reference_genome,
                                                   **inputs)

        chunk_size = 500
        for chunk in [pipeline_command[i:i + chunk_size] for i in range(0, len(pipeline_command), chunk_size)]:
            check_call(['cgcloud',
                        'ssh',
                        '--zone', params.zone,
                        '--cluster-name', params.cluster_name,
                        'toil-leader',
                        '-o', 'StrictHostKeyChecking=no',
                        'screen', '-S', params.cluster_name,
                        '-X', 'stuff', quote(chunk)])

    except CalledProcessError as e:
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
    return len(list_nodes(cluster_name))


def list_nodes(cluster_name, role='toil-worker'):
    """
    Returns list of dictionaries, each dictionary representing a worker node. Each dictinoary has the following keys:
    cluster_name, role_name, ordinal, cluster_ordinal, private_ip_address, ip_address, instance_id, instance_type,
    launch_time, state and zone
    """
    return parse_cgcloud_list_output(check_output(['cgcloud', 'list', '-c', cluster_name, role]))


def parse_cgcloud_list_output(output):
    return [row for row in csv.DictReader(StringIO(output), delimiter='\t')]


def get_desired_cluster_size(domain):
    nodes_per_sample = Samples.load(domain)

    sample_request = sum(map(lambda x: x[1].nodes, nodes_per_sample.value.iteritems()))
    log.info("Samples are currently requesting %d nodes: %s",
             sample_request,
             nodes_per_sample.value)

    # add one because there needs to be a single communal worker, else deadlock is upon us
    return sample_request + 1


def update_cluster_size(domain, n):
    """
    Only to be invoked by the uberscript
    """
    with cluster_size_lock:
        cluster_size = ClusterSize.load(domain)
        cluster_size.value = n
        # Holding the lock should prevent concurrent modifications
        cluster_size.save()


def grow_cluster(num_nodes, instance_type, cluster_name, cluster_type='toil', options=None):
    """
    Grow a cluster by a given number of nodes
    """
    if options is None:
        options = []

    nodes = []
    num_nodes_left = num_nodes
    while num_nodes_left > 0:
        log.info('Attempting to grow cluster by %i node(s) of type: %s', num_nodes_left, instance_type)
        cmd = (['cgcloud',
                'grow-cluster',
                '--list',
                '--instance-type', instance_type,
                '--num-workers', str(num_nodes_left),
                '--cluster-name', cluster_name] +
               options +
               [cluster_type])
        try:
            output = check_output(cmd)
        except CalledProcessError as e:
            log.warn('Command %r failed with status code %d.', cmd, e.returncode)
            output = e.output
        nodes_added = parse_cgcloud_list_output(output)
        nodes.extend(nodes_added)
        num_nodes_added = len(nodes_added)
        if num_nodes_added == 0:
            log.warn("Wasn't able to add any nodes. Wating 5 min.")
            time.sleep(5 * 60)
        assert num_nodes_added <= num_nodes_left
        num_nodes_left -= num_nodes_added
        log.info('Added %d node(s), %d node(s) left.', num_nodes_added, num_nodes_left)
    assert len(nodes) == num_nodes
    log.info('Successfully grew cluster by %i node(s) of type %s.', num_nodes, instance_type)
    return nodes


def manage_metrics_and_cluster_scaling(params):
    conn = boto.sdb.connect_to_region(region_of_zone(params.zone))
    dom = conn.get_domain('{0}--files'.format(params.jobstore))
    grow_cluster_thread = threading.Thread(target=manage_toil_cluster, args=(params, dom))
    metric_collection_thread = threading.Thread(target=collect_realtime_metrics, args=(params,))
    grow_cluster_thread.start()
    metric_collection_thread.start()
    grow_cluster_thread.join()
    metric_collection_thread.join()


def role_options(params):
    """
    :type params: argparse.Namespace
    """
    options = []
    if params.add_to_etc_hosts:
        options.extend(['-O', 'etc_hosts_entries=%s' % params.add_to_etc_hosts])
    return options


def spot_options(params):
    """
    :type params: argparse.Namespace
    """
    options = []
    if params.spot_price:
        options = ['--spot-bid', str(params.spot_price)]
    return options

@contextmanager
def throttle(interval):
    """
    A context manager for ensuring that the execution of its body takes at least a given amount of time, sleeping if
    necessary.
    """
    start = time.time()
    yield
    duration = time.time() - start
    remainder = interval - duration
    if remainder > 0:
        time.sleep(remainder)


def manage_toil_cluster(params, domain):
    log.info('Auto-scaling of Toil cluster has started.')
    time.sleep(scaling_initial_wait_period_in_seconds)
    while True:
        with throttle(cluster_scaling_interval_in_seconds):
            with cluster_size_lock:
                cluster_size = len(list_nodes(params.cluster_name))
                desired_cluster_size = get_desired_cluster_size(domain)
                if cluster_size < desired_cluster_size:
                    log.info("Cluster size (%d) is smaller than requested (%d).",
                             cluster_size,
                             desired_cluster_size)
                    num_nodes = desired_cluster_size - cluster_size
                    grow_cluster(num_nodes, params.instance_type, params.cluster_name, options=role_options(params) + spot_options(params))
                    update_cluster_size(domain, desired_cluster_size)
            cluster_size = len(list_nodes(params.cluster_name))
            log.info('Toil cluster is now at %i node(s).', cluster_size)


standalone_spark_semaphore_name = 'standalone_spark_sample_slots'


def manage_standalone_spark_cluster(params, domain):
    semaphore = Semaphore.create(domain=domain,
                                 name=standalone_spark_semaphore_name,
                                 value=params.spark_sample_slots)
    while True:
        with throttle(cluster_scaling_interval_in_seconds):
            # How many samples are being computed (or waiting to be computed) on the Spark cluster? ...
            num_unused_sample_slots = semaphore.value
            num_samples = params.spark_sample_slots - num_unused_sample_slots
            assert num_samples >= 0
            if num_samples:
                # ... At least one sample is. How many Spark workers are needed for all of them.
                spark_workers_per_sample = params.spark_nodes - 1
                num_desired_workers = num_samples * spark_workers_per_sample
                num_workers = len(list_nodes(params.cluster_name, role='spark-slave'))
                num_masters = len(list_nodes(params.cluster_name, role='spark-master'))
                if num_masters == 0:
                    # Absence of master indicates absence of cluster, so create it now.
                    assert num_workers == 0, 'Orphaned Spark workers detected'
                    output = check_output(['cgcloud',
                                           'create-cluster',
                                           '--list',
                                           '--zone', params.zone,
                                           '--cluster-name', params.cluster_name,
                                           '--leader-instance-type', params.spark_master_type or params.leader_type,
                                           '--instance-type', params.spark_instance_type or params.instance_type,
                                           '--num-workers', str(num_desired_workers),  # TODO: spot options
                                           '--ssh-opts', '"StrictHostKeyChecking=no"'] +
                                          role_options(params) +
                                          ['spark'])
                    output = parse_cgcloud_list_output(output)
                    assert len(output) == 1
                    master = output[0]
                    assert master['role_name'] == 'spark-master'
                    SparkMasterAddress.create(domain, master['private_ip_address']).save(force=True)
                elif num_masters == 1:
                    # There already is a Spark cluster, so scale it if necessary. 
                    if num_desired_workers > num_workers:
                        grow_cluster(num_nodes=num_desired_workers - num_workers,
                                     instance_type=params.instance_type,
                                     cluster_name=params.cluster_name,
                                     cluster_type='spark')
                    elif num_desired_workers < num_workers:
                        pass  # TODO: think about scaling down
                    else:
                        pass  # Cluster already has desired size
                else:
                    assert False, 'Unexpected number of master nodes: %i' % num_masters
        semaphore = Semaphore.load(domain, name=standalone_spark_semaphore_name)


def collect_realtime_metrics(params, threshold=0.5):
    """
    Collect metrics from AWS instances in 1 hour intervals. Instances that have gone idle (below threshold CPU value)
    are terminated.

    :type params: argparse.Namespace
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
    region = region_of_zone(params.zone)
    conn = boto.ec2.connect_to_region(region)
    cw = boto.ec2.cloudwatch.connect_to_region(region)
    sdbconn = boto.sdb.connect_to_region(region)
    domain = sdbconn.get_domain('{0}--files'.format(params.jobstore))

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
                            else:
                                log.info('Max CPU for {} was {} for last 30 minutes.'.format(instance_id, max(averages)))

                # Kill instance if idle and cluster is too large
                if idle:
                    try:
                        with cluster_size_lock:
                            cluster_size = get_cluster_size(params.cluster_name)
                            desired_cluster_size = get_desired_cluster_size(domain)
                            if cluster_size > desired_cluster_size:
                                log.info('Cluster size (%d) is larger than requested (%d).'
                                         'Terminating idle instance %s.',
                                         cluster_size,
                                         desired_cluster_size,
                                         instance_id)
                                
                                cmd = ['cgcloud',
                                       'terminate',
                                       '--instance-id', instance_id,
                                       '--cluster-name', params.cluster_name,
                                       'toil']

                                try:
                                    check_call(cmd)
                                    log.info("Successfully terminated instance via %s.",
                                             " ".join(cmd))
                                except:
                                    log.error("Terminating instance with %s failed.",
                                              " ".join(cmd))
                                    raise

                                update_cluster_size(domain, cluster_size - 1)
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
    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    
    # Launch cluster
    cluster_sp = subparsers.add_parser('launch-cluster',
                                       help='Launches EC2 cluster via CGCloud')
    cluster_sp.add_argument('-S', '--share', required=True,
                            help='Full path to directory containing pipeline script, launch script, and master key.')
    cluster_sp.add_argument('-T', '--leader-type', default='m3.medium',
                            help='Sets leader instance type.')
    cluster_sp.add_argument('-b', '--boto-path', default='/home/mesosbox/.boto', type=str,
                            help='Path to local .boto file to be placed on leader.')
    cluster_sp.add_argument('-M', '--manifest-path', default = None if not mock_mode() else "/home/ubuntu/toil-scripts/src/toil_scripts/adam_gatk_pipeline/mock_manifest",
                            required=not mock_mode(), help='Path to manifest file.')

    # Launch pipeline
    pipeline_sp = subparsers.add_parser('launch-pipeline',
                                        help='Launches pipeline')
    pipeline_sp.add_argument('-j', '--jobstore', default=None,
                             help='Name of jobstore. Defaults to UUID-Date if not set')
    pipeline_sp.add_argument('--restart', default=None, action='store_true',
                             help='Attempts to restart pipeline, requires existing job store.')
    pipeline_sp.add_argument('--master_ip', default=None,
                             help="The address of an external Spark master or 'auto' when using a standalone Spark "
                                  "cluster managed by this script. In that latter case you must pass the "
                                  "--max-samples-on-spark option to the launch-metric command.")
    pipeline_sp.add_argument('-B', '--bucket',
                             help='The name of the destination bucket.')
    pipeline_sp.add_argument('-m', '--memory', default='200g' if not mock_mode() else '3g',
                             help='The amount of memory per worker node in GiB. Must match what EC2 provides on '
                                  'the specified worker instance type. Defaults to 3 in mock mode')
    pipeline_sp.add_argument('-f', '--file_size', default='100G' if not mock_mode() else '10M',
                             help='Approximate size of the BAM files. Defaults to 10M in mock mode')
    pipeline_sp.add_argument('-s', '--spark_nodes', type=int, default=(8 if not mock_mode() else 2) + 1,
                             help="The number of Spark nodes, including the master, to allocate per sample. Relevant "
                                  "with separate running against a standSpark cluster managed, the master will be "
                                  "shared by all samples and the actual number of workers allocated per sample will be "
                                  "one less than specified here. Otherwise, each sample's subcluster will get its own "
                                  "master node. Default 9 in production mode, 3 in mock mode.")
    pipeline_sp.add_argument('-SD', '--sequence_dir', default='sequence',
                             help='Directory where raw sequences are.')
    pipeline_sp.add_argument('-R', '--reference_genome', default='GRCh38',
                             choices=['GRCh38', 'hg19'],
                             help='Reference genome to align and call against. Choose between GRCh38 and hg19.')

    # Launch metric collection
    metric_sp = subparsers.add_parser('launch-metrics',
                                      help='Launches metric collection thread')
    metric_sp.add_argument('-j', '--jobstore', required=True,  # differs subtly from launch-pipeline's --jobstore
                           help='Name of jobstore')
    metric_sp.add_argument('--namespace', default=os.environ.get('CGCLOUD_NAMESPACE', '/'),
                           help='CGCloud NameSpace')
    metric_sp.add_argument('--spark-sample-slots', required=False, default=0, type=int,
                           help='The maximum number of samples to be computed concurrently on a standalone Spark '
                                'cluster managed by this script. To be used in conjunction with the --master_ip=auto '
                                'option of the launch-cluster command. The default of 0 disables the standalone Spark '
                                'cluster.')

    # Common options
    cgcloud_zone = os.environ.get('CGCLOUD_ZONE')
    for sp in cluster_sp, pipeline_sp, metric_sp:
        sp.add_argument('-c', '--cluster-name', required=True,
                        help='The CGCloud cluster name for Toil leader and workers.')
        sp.add_argument('-z', '--zone', required=cgcloud_zone is None, default=cgcloud_zone,
                        help="The EC2 availability zone in which to place on-demand instances like the leaders of the "
                             "Toil and standalone Spark clusters. Also determines the region of the S3 bucket and SDB "
                             "domain for Toil's job store. The availability zone for spot instances may be chosen "
                             "independently from all zones in the region containing the specified zone.")
    for sp in cluster_sp, metric_sp:
        sp.add_argument('-t', '--instance-type', default='r3.8xlarge' if not mock_mode() else 'c3.large',
                        help='Worker instance type, e.g. m4.large or c3.8xlarge. Defaults to r3.8xlarge in production '
                             'mode. Will always use c3.large in mock mode, regardless of input value.')
        sp.add_argument('--spot-price', default=None, required=False,
                        help='Instance spot price if desired.')
    for sp in metric_sp, cluster_sp:
        sp.add_argument('-etc', '--add-to-etc-hosts', default=None, required=False,
                        help='Deprecated. Optional entry to add to /etc/hosts on Toil workers. This should *not* be '
                             'used to communicate the address of a standalone Spark master to driver jobs running on '
                             'Toil nodes. Use --master_ip=auto instead.')

    params = parser.parse_args()
    
    if params.command == 'launch-pipeline' and mock_mode() and params.master_ip:
        params.spark_sample_slots = 1

    if params.command == 'launch-cluster':
        launch_cluster(params)
        place_boto_on_leader(params)
    elif params.command == 'launch-pipeline':
        launch_pipeline(params)
    elif params.command == 'launch-metrics':
        manage_metrics_and_cluster_scaling(params)


def region_of_zone(availability_zone):
    # FIXME: cgcloud-lib's Context has something more robust for this
    return availability_zone[:-1]

if __name__ == '__main__':
    main()
