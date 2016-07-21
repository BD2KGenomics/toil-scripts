import os
import subprocess
import filecmp
from contextlib import closing
from uuid import uuid4

from toil.job import Job


def test_download_url_job(tmpdir):
    from toil_scripts.lib.urls import download_url_job
    options = Job.Runner.getDefaultOptions(os.path.join(str(tmpdir), 'test_store'))
    j = Job.wrapJobFn(download_url_job, 'www.google.com')
    Job.Runner.startToil(j, options)


def test_download_url(tmpdir):
    from toil_scripts.lib.urls import download_url
    work_dir = str(tmpdir)
    download_url(work_dir=work_dir, url='www.google.com', name='testy')
    assert os.path.exists(os.path.join(work_dir, 'testy'))


def test_upload_and_download_with_encryption(tmpdir):
    from toil_scripts.lib.urls import s3am_upload
    from toil_scripts.lib.urls import download_url
    from boto.s3.connection import S3Connection, Bucket, Key
    work_dir = str(tmpdir)
    # Create temporary encryption key
    key_path = os.path.join(work_dir, 'foo.key')
    subprocess.check_call(['dd', 'if=/dev/urandom', 'bs=1', 'count=32',
                           'of={}'.format(key_path)])
    # Create test file
    upload_fpath = os.path.join(work_dir, 'upload_file')
    with open(upload_fpath, 'wb') as fout:
        fout.write(os.urandom(1024))
    # Upload file
    random_key = os.path.join('test/', str(uuid4()), 'upload_file')
    s3_url = os.path.join('s3://cgl-driver-projects/', random_key)
    try:
        s3_dir = os.path.split(s3_url)[0]
        s3am_upload(fpath=upload_fpath, s3_dir=s3_dir, s3_key_path=key_path)
        # Download the file
        download_url(url=s3_url, name='download_file', work_dir=work_dir, s3_key_path=key_path)
        download_fpath = os.path.join(work_dir, 'download_file')
        assert os.path.exists(download_fpath)
        assert filecmp.cmp(upload_fpath, download_fpath)
    finally:
        # Delete the Key. Key deletion never fails so we don't need to catch any exceptions
        with closing(S3Connection()) as conn:
            b = Bucket(conn, 'cgl-driver-projects')
            k = Key(b)
            k.key = random_key
            k.delete()
