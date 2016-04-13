import base64
import glob
import hashlib
import os
import subprocess
from urlparse import urlparse
import shutil
from toil_scripts.lib.programs import docker_call


def download_url(url, work_dir='.', name=None, s3_key_path=None, cghub_key_path=None):
    """
    Downloads URL, can pass in file://, http://, s3://, or ftp://, gnos://cghub/analysisID, or gnos:///analysisID

    :param str url: URL to download from
    :param str work_dir: Directory to download file to
    :param str name: Name of output file, if None, basename of URL is used
    :param str s3_key_path: Path to 32-byte encryption key if url points to S3 file that uses SSE-C
    :param str cghub_key_path: Path to cghub key used to download from CGHub.
    :return: Path to the downloaded file
    :rtype: str
    """
    file_path = os.path.join(work_dir, name) if name else os.path.join(work_dir, os.path.basename(url))
    if s3_key_path:
        _download_encrypted_file(url, file_path, s3_key_path)
    elif cghub_key_path:
        _download_from_genetorrent(url, file_path, cghub_key_path)
    elif urlparse(url).scheme == 's3':
        _download_s3_url(file_path, url)
    elif urlparse(url).scheme == 'file':
        shutil.copy(urlparse(url).path, file_path)
    else:
        subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
    assert os.path.exists(file_path)
    return file_path


def download_url_job(job, url, name=None, s3_key_path=None, cghub_key_path=None):
    """Job version of `download_url`"""
    work_dir = job.fileStore.getLocalTempDir()
    fpath = download_url(url, work_dir=work_dir, name=name,
                         s3_key_path=s3_key_path, cghub_key_path=cghub_key_path)
    return job.fileStore.writeGlobalFile(fpath)


def s3am_upload(fpath, s3_dir, num_cores=1, s3_key_path=None):
    """
    Uploads a file to s3 via S3AM
    For SSE-C encryption: provide a path to a 32-byte file

    :param str fpath: Path to file to upload
    :param str s3_dir: Ouptut S3 path. Format: s3://bucket/[directory]
    :param int num_cores: Number of cores to use for up/download with S3AM
    :param str s3_key_path: (OPTIONAL) Path to 32-byte key to be used for SSE-C encryption
    """
    if not s3_dir.startswith('s3://'):
        raise ValueError('Format of s3_dir (s3://) is incorrect: {}'.format(s3_dir))
    s3_dir = os.path.join(s3_dir, os.path.basename(fpath))
    if s3_key_path:
        _s3am_with_retry(num_cores, '--sse-key-is-master', '--sse-key-file', s3_key_path,
                         'file://{}'.format(fpath), s3_dir)
    else:
        _s3am_with_retry(num_cores, 'file://{}'.format(fpath), s3_dir)


def s3am_upload_job(job, file_id, file_name, s3_dir, num_cores, s3_key_path=None):
    """Job version of `s3am_upload`"""
    work_dir = job.fileStore.getLocalTempDir()
    fpath = job.fileStore.readGlobalFile(file_id, os.path.join(work_dir, file_name))
    s3am_upload(fpath=fpath, s3_dir=s3_dir, num_cores=num_cores, s3_key_path=s3_key_path)


def _download_s3_url(file_path, url, ceph=None):
    """
    Downloads from S3 URL via Boto

    :param str file_path: Path to file
    :param str url: S3 URL
    """
    from boto.s3.connection import S3Connection
    s3 = S3Connection()
    try:
        parsed_url = urlparse(url)
        if not parsed_url.netloc or not parsed_url.path.startswith('/'):
            raise ValueError("An S3 URL must be of the form s3:/BUCKET/ or "
                             "s3://BUCKET/KEY. '%s' is not." % url)
        bucket = s3.get_bucket(parsed_url.netloc)
        key = bucket.get_key(parsed_url.path[1:])
        key.get_contents_to_filename(file_path)
    finally:
        s3.close()


def _download_encrypted_file(url, file_path, key_path):
    """
    Downloads encrypted files from S3

    :param str url: URL to be downloaded
    :param str file_path: Output path to file
    :param str key_path: Path to the 32-byte key file
    """
    # Grab master key
    with open(key_path, 'r') as f:
        key = f.read()
    if len(key) != 32:
        raise ValueError('Bad key in {}. Must be 32 bytes'.format(key_path))

    key = _generate_unique_key(key_path, url)
    # Create necessary headers for SSE-C encryption and download
    encoded_key = base64.b64encode(key)
    encoded_key_md5 = base64.b64encode(hashlib.md5(key).digest())
    h1 = 'x-amz-server-side-encryption-customer-algorithm:AES256'
    h2 = 'x-amz-server-side-encryption-customer-key:{}'.format(encoded_key)
    h3 = 'x-amz-server-side-encryption-customer-key-md5:{}'.format(encoded_key_md5)
    subprocess.check_call(['curl', '-fs', '--retry', '5', '-H', h1, '-H', h2, '-H', h3, url, '-o', file_path])
    assert os.path.exists(file_path)


def _download_from_genetorrent(url, file_path, cghub_key_path):
    url = urlparse(url)
    analysis_id = url.path[1:]
    assert url.scheme == 'gnos', 'Improper format. gnos://cghub/ID. User supplied: {}'.format(url)
    work_dir = os.path.dirname(file_path)
    folder_path = os.path.join(work_dir, os.path.basename(analysis_id))
    parameters = ['-vv', '-c', cghub_key_path, '-d', analysis_id]
    docker_call(tool='quay.io/ucsc_cgl/genetorrent:3.8.7--9911761265b6f08bc3ef09f53af05f56848d805b',
                work_dir=work_dir, parameters=parameters)
    sample = glob.glob(os.path.join(folder_path, '*tar*'))
    assert len(sample) == 1, 'More than one sample tar in CGHub download: {}'.format(analysis_id)


def _s3am_with_retry(num_cores, *args):
    """
    Calls S3AM upload with retries

    :param int num_cores: Number of cores to pass to upload/download slots
    :param list[str] args: Additional arguments to append to s3am
    """
    retry_count = 3
    for i in xrange(retry_count):
        s3am_command = ['s3am', 'upload', '--force', '--part-size=50M', '--exists=overwrite',
                        '--upload-slots={}'.format(num_cores),
                        '--download-slots={}'.format(num_cores)] + list(args)
        ret_code = subprocess.call(s3am_command)
        if ret_code == 0:
            return
        else:
            print 'S3AM failed with status code: {}'.format(ret_code)
    raise RuntimeError('S3AM failed to upload after {} retries.'.format(retry_count))


def _generate_unique_key(master_key_path, url):
    """
    Generate a unique encryption key given a URL and a path to another "master" encrypion key

    :param str master_key_path: Path to the master key
    :param str url: URL used to generate unique encryption key
    :return: The new 32-byte key
    :rtype: str
    """
    with open(master_key_path, 'r') as f:
        master_key = f.read()
    if len(master_key) != 32:
        raise ValueError('Bad key in {}. Must be 32 bytes'.format(master_key_path))
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is not 32 bytes and is invalid! Check code'
    return new_key
