import os
import subprocess
import shutil
from boto.s3.connection import S3Connection, Bucket, Key


def test_bwa(tmpdir):
    work_dir = str(tmpdir)
    create_config(work_dir)
    create_manifest(work_dir)
    subdir = '/mnt/ephemeral/toil-scripts/bwa'
    os.makedirs(os.path.join(subdir, 'workDir'))
    # Call Pipeline
    try:
        subprocess.check_call(['python', '-m', 'toil_scripts.batch_alignment.bwa_alignment', 'run',
                               os.path.join(subdir, 'jstore'),
                               '--manifest', os.path.join(work_dir, 'manifest.txt'),
                               '--config', os.path.join(work_dir, 'config.txt'),
                               '--retryCount', '1',
                               '--output-dir', 's3://cgl-driver-projects/test/ci',
                               '--workDir', os.path.join(subdir, 'workDir')])
    finally:
        shutil.rmtree(subdir)
        conn = S3Connection()
        b = Bucket(conn, 'cgl-driver-projects')
        k = Key(b)
        k.key = 'test/ci/ci_test.bam'
        k.delete()

def create_config(path):
    """Creates manifest file for test at path"""
    fpath = os.path.join(path, 'config.txt')
    with open(fpath, 'w') as f:
        f.write('ssec:\n'
                'library: test\n'
                'program_unit: 12345\n'
                'platform: ILLUMINA\n'
                'rg_line:\n'
                'ref: s3://cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa\n'
                'amb: s3://cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.amb\n'
                'ann: s3://cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.ann\n'
                'bwt: s3://cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.bwt\n'
                'fai: s3://cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.fai\n'
                'pac: s3://cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.pac\n'
                'sa: s3://cgl-pipeline-inputs/alignment/ci/hg38_chr6.fa.sa\n'
                'alt:\n'
                'file_size: 1G\n'
                'sort: True\n'
                'trim_adapters: False\n')


def create_manifest(path):
    """Creates config file for test at path"""
    fpath = os.path.join(path, 'manifest.txt')
    with open(fpath, 'w') as f:
        f.write('ci_test\t'
                's3://cgl-pipeline-inputs/alignment/ci/r1_trunc.fq.gz\t'
                's3://cgl-pipeline-inputs/alignment/ci/r2_trunc.fq.gz\n')
