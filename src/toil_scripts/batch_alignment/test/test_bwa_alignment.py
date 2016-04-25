import os
import subprocess
import shutil
from boto.s3.connection import S3Connection, Bucket, Key


def test_rnaseq_cgl(tmpdir):
    work_dir = str(tmpdir)
    create_config(work_dir)
    subdir = '/mnt/ephemeral/toil-scripts/bwa'
    os.makedirs(os.path.join(subdir, 'workDir'))
    # URLs for
    ref = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa'
    amb = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.amb'
    ann = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.ann'
    bwt = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.bwt'
    fai = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.fai'
    pac = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.pac'
    sa = 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.sa'
    # Call Pipeline
    try:
        subprocess.check_call(['python', '-m', 'toil_scripts.batch_alignment.bwa_alignment',
                               os.path.join(subdir, 'jstore'),
                               '--config', os.path.join(work_dir, 'config.txt'),
                               '--retryCount', '1',
                               '--s3-dir', 's3://cgl-driver-projects/test/ci',
                               '--workDir', os.path.join(subdir, 'workDir'),
                               '--ref', ref,
                               '--amb', amb,
                               '--ann', ann,
                               '--bwt', bwt,
                               '--fai', fai,
                               '--pac', pac,
                               '--sa', sa,
                               '--library', 'test',
                               '--file-size', '1G'])
    finally:
        shutil.rmtree(subdir)
        conn = S3Connection()
        b = Bucket(conn, 'cgl-driver-projects')
        k = Key(b)
        k.key = 'test/ci/ci_test.bam'
        k.delete()


def create_config(path):
    """Creates config file for test at path"""
    fpath = os.path.join(path, 'config.txt')
    with open(fpath, 'w') as f:
        f.write('ci_test,'
                'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/r1.fq.gz,'
                'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/ci/r2.fq.gz')
