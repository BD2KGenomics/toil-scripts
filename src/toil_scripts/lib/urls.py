import glob
import os
import shutil
import subprocess
from urlparse import urlparse

from toil_scripts.lib import require
from toil_scripts.lib.programs import docker_call


def download_url(url, work_dir='.', name=None, num_cores=1, s3_key_path=None, cghub_key_path=None):
    """
    Downloads URL, can pass in file://, http://, s3://, or ftp://, gnos://cghub/analysisID, or gnos:///analysisID

    :param str url: URL to download from
    :param str work_dir: Directory to download file to
    :param str name: Name of output file, if None, basename of URL is used
    :param int num_cores: Number of cores to use if downloading with s3am
    :param str s3_key_path: Path to 32-byte encryption key if url points to S3 file that uses SSE-C
    :param str cghub_key_path: Path to cghub key used to download from CGHub.
    :return: Path to the downloaded file
    :rtype: str
    """
    file_path = os.path.join(work_dir, name) if name else os.path.join(work_dir, os.path.basename(url))
    if cghub_key_path:
        _download_from_genetorrent(url, file_path, cghub_key_path)
    elif urlparse(url).scheme == 's3':
        _s3am_download_with_retry(file_path, url, num_cores, s3_key_path)
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


def s3am_upload(fpath, s3_dir, num_cores=1, s3_key_path=None):
    """
    Uploads a file to s3 via S3AM
    For SSE-C encryption: provide a path to a 32-byte file

    :param str fpath: Path to file to upload
    :param str s3_dir: Ouptut S3 path. Format: s3://bucket/[directory]
    :param int num_cores: Number of cores to use for up/download with S3AM
    :param str s3_key_path: (OPTIONAL) Path to 32-byte key to be used for SSE-C encryption
    """
    require(s3_dir.startswith('s3://'), 'Format of s3_dir (s3://) is incorrect: {}'.format(s3_dir))
    s3_dir = os.path.join(s3_dir, os.path.basename(fpath))
    args = ['upload', '--force', '--upload-slots={}'.format(num_cores), '--exists=overwrite']
    if s3_key_path:
        args.extend(['--sse-key-is-master', '--sse-key-file', s3_key_path])
    args.extend(['file://{}'.format(fpath), s3_dir])
    _s3am_with_retry(num_cores, s3am_args=args)


def s3am_upload_job(job, file_id, file_name, s3_dir, num_cores, s3_key_path=None):
    """Job version of `s3am_upload`"""
    work_dir = job.fileStore.getLocalTempDir()
    fpath = job.fileStore.readGlobalFile(file_id, os.path.join(work_dir, file_name))
    s3am_upload(fpath=fpath, s3_dir=s3_dir, num_cores=num_cores, s3_key_path=s3_key_path)


def _s3am_download_with_retry(file_path, s3_url, num_cores=1, s3_key_path=None):
    """
    Calls s3am for downloading

    :param str file_path: Location to download file
    :param str s3_url: S3 URL to download
    :param int num_cores: Number of cores to pass
    :param str s3_key_path: (OPTIONAL) Path to 32-byte key to be used for SSE-C encryption
    """
    require(s3_url.startswith('s3://'), 'Format of s3_dir (s3://) is incorrect: {}'.format(s3_url))
    args = ['download', '--file-exists=overwrite', '--download-exists=discard']
    if s3_key_path:
        args.extend(['--sse-key-is-master', '--sse-key-file', s3_key_path])
    args.extend([s3_url, 'file://' + file_path])
    _s3am_with_retry(num_cores, s3am_args=args)


def _s3am_with_retry(num_cores, s3am_args):
    """
    Calls S3AM upload with retries

    :param int num_cores: Number of cores to pass to upload/download slots
    :param list[str] s3am_args: Additional arguments for s3am
    """
    retry_count = 3
    for i in xrange(retry_count):
        s3am_command = ['s3am'] + s3am_args + ['--part-size=50M','--download-slots={}'.format(num_cores)]
        ret_code = subprocess.call(s3am_command)
        if ret_code == 0:
            return
        else:
            print 'S3AM failed with status code: {}'.format(ret_code)
    raise RuntimeError('S3AM failed to upload after {} retries.'.format(retry_count))
