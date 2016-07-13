import logging
import shutil
import subprocess
import textwrap
from contextlib import closing
from urlparse import urlparse

import os
import posixpath
from boto.s3.connection import S3Connection, Bucket
from toil_scripts.lib import get_work_directory

log = logging.getLogger(__name__)

input_dir = urlparse('s3://cgl-pipeline-inputs/rnaseq_cgl/ci')
output_dir = urlparse('s3://cgl-driver-projects/test/ci')
sample = urlparse(input_dir.geturl() + '/chr6_sample.tar.gz')


def test_rnaseq_cgl():
    logging.basicConfig(level=logging.INFO)  # FIXME: I wish we didn't have to do this here
    workdir = get_work_directory()
    try:
        # Call Pipeline
        base_command = ['toil-rnaseq', 'run',
                        '--config', generate_config(workdir),
                        os.path.join(workdir, 'jstore'),
                        '--retryCount', '1',
                        '--workDir', workdir]
        # Run with --samples
        try:
            subprocess.check_call(base_command + ['--samples', sample.geturl()])
        finally:
            cleanup_and_validate()
        # Run with manifest
        try:
            subprocess.check_call(base_command + ['--manifest', generate_manifest(workdir)])
        finally:
            cleanup_and_validate()
    finally:
        shutil.rmtree(workdir)


def cleanup_and_validate():
    valid_output = False
    expected_output = posixpath.basename(sample.path)
    with closing(S3Connection()) as s3:
        b = Bucket(s3, output_dir.netloc)
        path = output_dir.path[1:]
        for k in b.list(prefix=path):
            assert k.name.startswith(path)
            if k.name[len(path):] == '/' + expected_output:
                # FIXME: We may want to validate the output a bit more
                valid_output = True
            else:
                log.warn('Unexpected output file %s/%s', output_dir.geturl(), k.name)
            k.delete()
    assert valid_output, 'Did not find expected output file'


def generate_config(workdir):
    path = os.path.join(workdir, 'config-toil-rnaseq.yaml')
    with open(path, 'w') as f:
        f.write(textwrap.dedent("""
            star-index: {input_dir}/starIndex_chr6.tar.gz
            kallisto-index: s3://cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx
            rsem-ref: {input_dir}/rsem_ref_chr6.tar.gz
            output-dir:
            s3-output-dir: {output_dir}
            fastqc: true
            cutadapt:
            ssec:
            gtkey:
            wiggle:
            save-bam:
            fwd-3pr-adapter: AGATCGGAAGAG
            rev-3pr-adapter: AGATCGGAAGAG
            ci-test: true
            """[1:]).format(output_dir=output_dir.geturl(),
                            input_dir=input_dir.geturl()))
    return path


def generate_manifest(workdir):
    path = os.path.join(workdir, 'manifest-toil-rnaseq.tsv')
    with open(path, 'w') as f:
        f.write('\t'.join(['tar', 'paired', 'chr6_sample', sample.geturl()]))
    return path
