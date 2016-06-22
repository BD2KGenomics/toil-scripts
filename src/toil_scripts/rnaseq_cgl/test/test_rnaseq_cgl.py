import os
import shutil
import subprocess
import textwrap
from uuid import uuid4

from bd2k.util.files import mkdir_p
from boto.s3.connection import S3Connection, Bucket, Key


def test_rnaseq_cgl(tmpdir):
    subdir = '/mnt/ephemeral/toil-scripts/rnaseq/'
    work_dir = os.path.join(subdir, str(uuid4()))
    mkdir_p(work_dir)
    create_config_and_manifest(work_dir)
    sample = 's3://cgl-pipeline-inputs/rnaseq_cgl/ci/chr6_sample.tar.gz'
    # Call Pipeline
    try:
        base_command = ['toil-rnaseq', 'run',
                        '--config', os.path.join(work_dir, 'config-toil-rnaseq.yaml'),
                        os.path.join(work_dir, 'jstore'),
                        '--retryCount', '1',
                        '--workDir', work_dir]
        # Run with --samples
        subprocess.check_call(base_command + ['--samples', sample])
        # Run with manifest
        subprocess.check_call(base_command + ['--manifest', os.path.join(work_dir, 'manifest-toil-rnaseq.tsv')])
    finally:
        shutil.rmtree(subdir)
        conn = S3Connection()
        b = Bucket(conn, 'cgl-driver-projects')
        k = Key(b)
        k.key = 'test/ci/chr6_sample.tar.gz'
        k.delete()


def create_config_and_manifest(path):
    """Creates config file for test at path"""
    config_path = os.path.join(path, 'config-toil-rnaseq.yaml')
    manifest_path = os.path.join(path, 'manifest-toil-rnaseq.tsv')
    with open(config_path, 'w') as f:
        f.write(generate_config())
    with open(manifest_path, 'w') as f:
        f.write(generate_manifest())


def generate_config():
    return textwrap.dedent("""
        star-index: s3://cgl-pipeline-inputs/rnaseq_cgl/ci/starIndex_chr6.tar.gz
        kallisto-index: s3://cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx
        rsem-ref: s3://cgl-pipeline-inputs/rnaseq_cgl/ci/rsem_ref_chr6.tar.gz
        output-dir:
        s3-output-dir: s3://cgl-driver-projects/test/ci
        fastqc: true
        cutadapt:
        ssec:
        gtkey:
        wiggle:
        save-bam:
        fwd-3pr-adapter: AGATCGGAAGAG
        rev-3pr-adapter: AGATCGGAAGAG
        ci-test: true
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        tar\tpaired\tchr6_sample\ts3://cgl-pipeline-inputs/rnaseq_cgl/ci/chr6_sample.tar.gz
        """[1:])
