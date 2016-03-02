#!/usr/bin/env python2.7
"""
Author: John Vivian

Library of convenience functions for AWS

Assumes a ~/.boto credentials file:

[Credentials]
aws_access_key_id = AWS_ACCESS_KEY_ID
aws_secret_access_key = AWS_SECRET_ACCESS_KEY

"""
import base64
import hashlib
import logging
import os
import subprocess
import boto
import boto.ec2
import time
import datetime

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def dl_encrypted_file_from_s3(url, file_path):
    """
    url: str            S3 AWS string
    file_path: str      path where file will be downloaded
    """
    def generate_unique_key(master_key_path, url):
        """
        master_key_path: str    Path to the BD2K Master Key (for S3 Encryption)
        url: str                S3 URL (e.g. https://s3-us-west-2.amazonaws.com/bucket/file.txt)

        Returns: str            32-byte unique key generated for that URL
        """
        with open(master_key_path, 'r') as f:
            master_key = f.read()
        assert len(master_key) == 32, 'Invalid Key! Must be 32 characters. Key: {}, Length: {}'.format(master_key, len(master_key))
        new_key = hashlib.sha256(master_key + url).digest()
        assert len(new_key) == 32, 'New key is invalid and is not 32 characters: {}'.format(new_key)
        return new_key
    key = generate_unique_key('master.key', url)
    encoded_key = base64.b64encode(key)
    encoded_key_md5 = base64.b64encode(hashlib.md5(key).digest())
    h1 = 'x-amz-server-side-encryption-customer-algorithm:AES256'
    h2 = 'x-amz-server-side-encryption-customer-key:{}'.format(encoded_key)
    h3 = 'x-amz-server-side-encryption-customer-key-md5:{}'.format(encoded_key_md5)
    subprocess.check_call(['curl', '-fs', '--retry', '5', '-H', h1, '-H', h2, '-H', h3, url, '-o', file_path])


def get_instance_ips(filter_name=None, filter_cluster=None):
    import boto.ec2
    ips = []
    filters = {}
    if filter_name:
        filters['tag:Name'] = filter_name
    if filter_cluster:
        filters['tag:cluster_name'] = filter_cluster
    conn = boto.ec2.connect_to_region("us-west-2")
    reservations = conn.get_all_reservations(filters=filters)
    for r in reservations:
        for i in r.instances:
            if i.ip_address is not None:
                ips.append(i.ip_address)
                print str(i.ip_address)
    return ips


def list_s3_urls(bucket_name, directory=''):
    conn = boto.connect_s3()
    bucket = conn.get_bucket(bucket_name)
    prefix = 'https://s3-us-west-2.amazonaws.com/'
    for key in bucket.list(directory):
        print os.path.join(prefix, bucket_name, directory, key.name)


def calculate_ec2_spot_instance_cost(instance_id, instance_type='c3.8xlarge', avail_zone='us-west-2a'):
    """
    Computes the spot market cost for a stopped or terminated instance given 3 values:
    instance_type, instance_id, avail_zone

    Example:
    python calculate_ec2_spot_instance.py -t c3.8xlarge -i i-b3a1cd6a -a us-west-2a
    """

    def get_start_and_stop(instance_id, region='us-west-2'):
        """
        calculates start and stop time of an instance

        :returns: strs      startTime, endTime
        """
        logging.info('Acquiring start and stop time of instance...')
        start, stop = 0, 0
        conn = boto.ec2.connect_to_region(region)
        reservations = conn.get_all_instances(instance_id)
        for r in reservations:
            for i in r.instances:
                start = i.launch_time
                if i.state == 'terminated' or i.state == 'stopped':
                    # Convert stop to proper format
                    stop = i.reason.split('(')[1].split(')')[0]
                    stop = stop.split()
                    stop = stop[0] + 'T' + stop[1] + '.000Z'
                else:
                    logging.info('Instance not stopped or terminated yet. Using current UTC time.')
                    t = datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
                    stop = t[:4] + '-' + t[4:6] + '-' + t[6:8] + 'T' + \
                           t[8:10] + ':' + t[10:12] + ':' + t[12:14] + '.000Z'
        if start == 0:
            raise RuntimeError('Spot Instance {} not found'.format(instance_id))
        return start, stop

    def calculate_cost(instance_type, start_time, end_time, avail_zone, region='us-west-2'):
        # Some values
        logging.info('Calculating costs...')
        total, n = 0.0, 0
        # Connect to EC2 -- requires ~/.boto
        conn = boto.ec2.connect_to_region(region)
        # Get prices for instance, AZ and time range
        prices = conn.get_spot_price_history(instance_type=instance_type, start_time=start_time,
                                             end_time=end_time, availability_zone=avail_zone)
        # Output the prices
        for price in prices:
            total += price.price
            n += 1
        # Difference b/w first and last returned times
        stop = time.mktime(datetime.datetime.strptime(end_time, "%Y-%m-%dT%H:%M:%S.000Z").timetuple())
        start = time.mktime(datetime.datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S.000Z").timetuple())
        time_diff = (stop - start) / 3600
        # Output aggregate, average and max results
        print "For one: {} in Zone: {}".format(instance_type, avail_zone)
        print "From: {} to {}".format(start_time, end_time)
        print "\tTotal cost = $" + str(time_diff * (total/n))
        print "\tAvg hourly cost = $" + str(total / n)

    start_time, end_time = get_start_and_stop(instance_id)
    calculate_cost(instance_type, start_time, end_time, avail_zone)


def count_items_in_bucket(bucket_name, directory=''):
    conn = boto.connect_s3()
    bucket = conn.get_bucket(bucket_name)
    n = 0
    for key in bucket.list(directory):
        n += 1
    return n


def get_instance_ids(filter_name=None, filter_cluster=None):
    """
    returns a list of instance IDs given some filters

    filter_name: str        Name by which to filter
    filter_cluster : str    Cluster name by which to filter
    """
    import boto.ec2
    ids = []
    filters = {}
    if filter_name:
        filters['tag:Name'] = filter_name
    if filter_cluster:
        filters['tag:cluster_name'] = filter_cluster
    conn = boto.ec2.connect_to_region("us-west-2")
    reservations = conn.get_all_reservations(filters=filters)
    for r in reservations:
        for i in r.instances:
            if i.state == 'running':
                ids.append(str(i.id))
    return ids


def get_avail_zone(filter_name=None, filter_cluster=None):
    import boto.ec2
    ids = []
    filters = {}
    if filter_name:
        filters['tag:Name'] = filter_name
    if filter_cluster:
        filters['tag:cluster_name'] = filter_cluster
    conn = boto.ec2.connect_to_region("us-west-2")
    reservations = conn.get_all_reservations(filters=filters)
    for r in reservations:
        for i in r.instances:
            if i.placement is not None:
                ids.append(str(i.placement))
    return ids


def apply_alarm_to_instance(instance_id, namespace='AWS/EC2', metric='CPUUtilization', statistic='Average',
                            comparison='<', threshold=0.5, period=300, evaluation_periods=1, region='us-west-2'):
    """
    Applies an alarm to a given instance that terminates after a consecutive period of 1 hour at 0.5 CPU usage.

    instance_id: str        ID of the EC2 Instance
    region: str             AWS region
    """
    import boto.ec2.cloudwatch
    logging.info('Applying Cloudwatch alarm to: {}'.format(instance_id))
    cw = boto.ec2.cloudwatch.connect_to_region(region)
    alarm = boto.ec2.cloudwatch.MetricAlarm(name='CPUAlarm_{}'.format(instance_id),
                                            description='Terminate instance after low CPU.', namespace=namespace,
                                            metric=metric, statistic=statistic, comparison=comparison,
                                            threshold=threshold, period=period, evaluation_periods=evaluation_periods,
                                            dimensions={'InstanceId': [instance_id]},
                                            alarm_actions=['arn:aws:automate:{}:ec2:terminate'.format(region)])
    cw.put_metric_alarm(alarm)