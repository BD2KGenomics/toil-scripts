#!/usr/bin/env python2.7
import os
from urlparse import urlparse

from bd2k.util.files import mkdir_p
from toil_lib.files import copy_files
from toil_lib.urls import s3am_upload


def output_file_job(job, filename, file_id, output_dir, s3_key_path=None):
    """
    Uploads a file from the FileStore to an output directory on the local filesystem or S3.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str filename: basename for file
    :param str file_id: FileStoreID
    :param str output_dir: Amazon S3 URL or local path
    :param str s3_key_path: (OPTIONAL) Path to 32-byte key to be used for SSE-C encryption
    :return:
    """
    job.fileStore.logToMaster('Writing {} to {}'.format(filename, output_dir))
    work_dir = job.fileStore.getLocalTempDir()
    filepath = job.fileStore.readGlobalFile(file_id, os.path.join(work_dir, filename))
    if urlparse(output_dir).scheme == 's3':
        s3am_upload(fpath=os.path.join(work_dir, filepath),
                    s3_dir=output_dir,
                    s3_key_path=s3_key_path)
    elif os.path.exists(os.path.join(output_dir, filename)):
        job.fileStore.logToMaster("File already exists: {}".format(filename))
    else:
        mkdir_p(output_dir)
        copy_files([filepath], output_dir)
