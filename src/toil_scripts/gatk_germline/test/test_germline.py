import logging
import shlex
import shutil
import subprocess
import tempfile
import textwrap
from unittest import TestCase
from uuid import uuid4

import os
from bd2k.util.iterables import concat

log = logging.getLogger(__name__)


class GermlineTest(TestCase):
    """
    These tests *can* be parameterized with the following optional environment variables:

    TOIL_SCRIPTS_TEST_TOIL_OPTIONS - a space-separated list of additional command line arguments to pass to Toil via
    the script entry point. Default is the empty string.

    TOIL_SCRIPTS_TEST_JOBSTORE - the job store locator to use for the tests. The default is a file: locator pointing
    at a local temporary directory.

    TOIL_SCRIPTS_TEST_NUM_SAMPLES - the number of sample lines to generate in the input manifest
    """

    @classmethod
    def setUpClass(cls):
        super(GermlineTest, cls).setUpClass()
        # FIXME: pull up into common base class
        logging.basicConfig(level=logging.INFO)

    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.fastq_url = 's3://cgl-pipeline-inputs/germline/ci/NIST7035_NIST7086.aln21.ci.1.fq'
        self.bam_sample = ['bam_test', 's3://cgl-pipeline-inputs/germline/ci/NIST7035_NIST7086.aln21.ci.bam']
        jobStore = os.getenv('TOIL_SCRIPTS_TEST_JOBSTORE', os.path.join(self.workdir, 'jobstore-%s' % uuid4()))
        toilOptions = shlex.split(os.environ.get('TOIL_SCRIPTS_TEST_TOIL_OPTIONS', ''))
        self.base_command = concat('toil-germline', 'run',
                                   toilOptions,
                                   jobStore)
        self.expected_files = set()

    def test_samples_option(self):
        """
        Tests passing a BAM sample using the command line interface
        """
        self.expected_files |= {'bam_test.raw.ci_test.gvcf',
                                'bam_test.hard_filter.ci_test.vcf',
                                'config-toil-germline-bam.yaml'}
        self._run(self.base_command,
                  '--sample', self.bam_sample,
                  '--config', self._generate_bam_config())
        self._assertOutput()

    def _test_manifest_option(self):
        """
        Tests using manifest containing FASTQ files
        """
        num_samples = int(os.environ.get('TOIL_SCRIPTS_TEST_NUM_SAMPLES', '3'))
        for i in range(1, num_samples+1):
            self.expected_files |= {'fastq_test_%s.raw.ci_test.gvcf' % i,
                                    'fastq_test_%s.hard_filter.ci_test.vcf' % i,
                                    'config-toil-germline-fastq.yaml',
                                    'manifest-toil-germline.tsv'}

        self._run(self.base_command,
                  '--config', self._generate_fastq_config(),
                  '--manifest', self._generate_manifest(num_samples))
        self._assertOutput()

    def _run(self, *args):
        args = list(concat(*args))
        log.info('Running %r', args)
        subprocess.check_call(args)

    def tearDown(self):
        shutil.rmtree(self.workdir)

    def _generate_bam_config(self):
        path = os.path.join(self.workdir, 'config-toil-germline-bam.yaml')
        with open(path, 'w') as f:
            f.write(textwrap.dedent("""
                    genome-fasta: s3://cgl-pipeline-inputs/germline/ci/b37_21.fa
                    assembly: b37
                    phase: s3://cgl-pipeline-inputs/germline/ci/1000G_phase1.indels.b37.21.recode.vcf
                    mills: s3://cgl-pipeline-inputs/germline/ci/Mills_and_1000G_gold_standard.indels.b37.21.recode.vcf
                    dbsnp: s3://cgl-pipeline-inputs/germline/ci/dbsnp_138.b37.21.recode.vcf
                    hapmap: s3://cgl-pipeline-inputs/germline/ci/hapmap_3.3.b37.21.recode.vcf
                    omni: s3://cgl-pipeline-inputs/germline/ci/1000G_omni2.5.b37.21.recode.vcf
                    run-bwa: True
                    trim: False
                    amb: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.amb
                    ann: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.ann
                    bwt: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.bwt
                    pac: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.pac
                    sa: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.sa
                    alt:
                    preprocess: True
                    ssec:
                    file-size: 1G
                    cores: 2
                    xmx: 5G
                    suffix: .ci_test
                    output-dir: {output_dir}
                    unsafe-mode:
                    run-vqsr:
                    joint: False
                    preprocess-only:
                    run-oncotator: False
                    synapse-name:
                    synapse-pwd:
                    """[1:]).format(output_dir=self.workdir))
        return path


    def _generate_fastq_config(self):
        path = os.path.join(self.workdir, 'config-toil-germline-fastq.yaml')
        with open(path, 'w') as f:
            f.write(textwrap.dedent("""
                    genome-fasta: s3://cgl-pipeline-inputs/germline/ci/b37_21.fa
                    assembly: b37
                    phase: s3://cgl-pipeline-inputs/germline/ci/1000G_phase1.indels.b37.21.recode.vcf
                    mills: s3://cgl-pipeline-inputs/germline/ci/Mills_and_1000G_gold_standard.indels.b37.21.recode.vcf
                    dbsnp: s3://cgl-pipeline-inputs/germline/ci/dbsnp_138.b37.21.recode.vcf
                    hapmap: s3://cgl-pipeline-inputs/germline/ci/hapmap_3.3.b37.21.recode.vcf
                    omni: s3://cgl-pipeline-inputs/germline/ci/1000G_omni2.5.b37.21.recode.vcf
                    run-bwa: True
                    trim: False
                    preprocess: True
                    amb: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.amb
                    ann: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.ann
                    bwt: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.bwt
                    pac: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.pac
                    sa: s3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.sa
                    alt:
                    ssec:
                    file-size: 1G
                    cores: 2
                    xmx: 5G
                    suffix: .ci_test
                    output-dir: {output_dir}
                    unsafe-mode: False
                    run-vqsr: False
                    joint: False
                    preprocess-only: False
                    run-oncotator: False
                    synapse-name:
                    synapse-pwd:
                    """[1:]).format(output_dir=self.workdir))
        return path

    def _generate_manifest(self, num_samples):
        path = os.path.join(self.workdir, 'manifest-toil-germline.tsv')
        with open(path, 'w') as f:
            f.write('\n'.join('\t'.join(['fastq_test_%s' % i,
                                         self.fastq_url,
                                         '@RG\\tID:foo\\tSM:bar\\tPL:ILLUMINA'])
                          for i in range(1, num_samples+1)))
        return path

    def _assertOutput(self):
        """
        Checks that all output files are expected and non-zero in size
        """
        for root, dirs, files in os.walk(self.workdir, topdown=False):
            for name in files:
                self.assertTrue(name in self.expected_files)
                self.assertTrue(os.stat(os.path.join(root, name)).st_size > 0)
