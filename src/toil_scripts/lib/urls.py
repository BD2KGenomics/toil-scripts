import glob
import os
import shutil
import subprocess
from urlparse import urlparse

from toil_scripts.lib import require
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
    if cghub_key_path:
        _download_with_genetorrent(url, file_path, cghub_key_path)
    elif urlparse(url).scheme == 's3':
        _s3am_with_retry(num_cores=1, file_path=file_path, s3_url=url, mode='download', s3_key_path=s3_key_path)
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


def _download_with_genetorrent(url, file_path, cghub_key_path):
    parsed_url = urlparse(url)
    analysis_id = parsed_url.path[1:]
    assert parsed_url.scheme == 'gnos', 'Improper format. gnos://cghub/ID. User supplied: {}'.format(parsed_url)
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
    _s3am_with_retry(num_cores, file_path=fpath, s3_url=s3_dir, mode='upload', s3_key_path=s3_key_path)


def s3am_upload_job(job, file_id, file_name, s3_dir, s3_key_path=None):
    """Job version of s3am_upload"""
    work_dir = job.fileStore.getLocalTempDir()
    fpath = job.fileStore.readGlobalFile(file_id, os.path.join(work_dir, file_name))
    s3am_upload(fpath=fpath, s3_dir=s3_dir, num_cores=job.cores, s3_key_path=s3_key_path)


def _s3am_with_retry(num_cores, file_path, s3_url, mode='upload', s3_key_path=None):
    """
    Run s3am with 3 retries

    :param int num_cores: Number of cores to pass to upload/download slots
    :param str file_path: Full path to the file
    :param str s3_url: S3 URL
    :param str mode: Mode to run s3am in. Either "upload" or "download"
    :param str s3_key_path: Path to the SSE-C key if using encryption
    """
    command = ['s3am']
    if mode == 'upload':
        command.extend(['upload', '--force', '--upload-slots={}'.format(num_cores),
                        '--exists=overwrite', 'file://' + file_path, s3_url])
    elif mode == 'download':
        command.extend(['download', '--file-exists=overwrite',
                        '--download-exists=discard', s3_url, 'file://' + file_path])
    else:
        raise ValueError('Improper mode specified. mode must be equal to "upload" or "download".')
    if s3_key_path:
        for arg in [s3_key_path, '--sse-key-file', '--sse-key-is-master']:
            command.insert(2, arg)
    for arg in ['--part-size=50M', '--download-slots={}'.format(num_cores)]:
        command.insert(2, arg)
    # Run s3am with retries
    retry_count = 3
    for i in xrange(retry_count):
        ret_code = subprocess.call(command)
        if ret_code == 0:
            return
        else:
            print 'S3AM failed with status code: {}'.format(ret_code)
    raise RuntimeError('S3AM failed to {} after {} retries.'.format(mode, retry_count))
