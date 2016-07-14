import logging
import shutil
import subprocess
import textwrap
from contextlib import closing
from unittest import TestCase
from urlparse import urlparse
from uuid import uuid4

import os
import posixpath
from boto.s3.connection import S3Connection, Bucket
from toil_scripts.lib import get_work_directory

log = logging.getLogger(__name__)


class RNASeqCGLTest(TestCase):
    @classmethod
    def setUpClass(cls):
        super(RNASeqCGLTest, cls).setUpClass()
        # FIXME: pull up into common base class
        logging.basicConfig(level=logging.INFO)

    def setUp(self):
        self.input_dir = urlparse('s3://cgl-pipeline-inputs/rnaseq_cgl/ci')
        self.output_dir = urlparse('s3://cgl-driver-projects/test/ci/%s' % uuid4())
        self.sample = urlparse(self.input_dir.geturl() + '/chr6_sample.tar.gz')
        self.workdir = get_work_directory()
        self.base_command = ['toil-rnaseq', 'run',
                             '--config', self._generate_config(),
                             os.path.join(self.workdir, 'jobstore'),
                             '--retryCount', '1',
                             '--workDir', self.workdir]

    def test_with_samples_option(self):
        subprocess.check_call(self.base_command + ['--samples', self.sample.geturl()])
        self._assertOutput()

    def test_with_manifest(self):
        subprocess.check_call(self.base_command + ['--manifest', self._generate_manifest()])
        self._assertOutput()

    def _assertOutput(self):
        with closing(S3Connection()) as s3:
            bucket = Bucket(s3, self.output_dir.netloc)
            prefix = self.output_dir.path[1:]
            output_file = posixpath.basename(self.sample.path)
            bucket.get_key(posixpath.join(prefix, output_file), validate=True)
            # FIXME: We may want to validate the output a bit more

    def tearDown(self):
        shutil.rmtree(self.workdir)
        with closing(S3Connection()) as s3:
            bucket = Bucket(s3, self.output_dir.netloc)
            prefix = self.output_dir.path[1:]
            for key in bucket.list(prefix=prefix):
                assert key.name.startswith(prefix)
                key.delete()

    def _generate_config(self):
        path = os.path.join(self.workdir, 'config-toil-rnaseq.yaml')
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
                    """[1:]).format(output_dir=self.output_dir.geturl(),
                                    input_dir=self.input_dir.geturl()))
        return path

    def _generate_manifest(self):
        path = os.path.join(self.workdir, 'manifest-toil-rnaseq.tsv')
        with open(path, 'w') as f:
            f.write('\t'.join(['tar', 'paired', 'chr6_sample', self.sample.geturl()]))
        return path
