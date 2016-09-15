from __future__ import print_function

import argparse
import logging
import os
import shlex
import shutil
import subprocess
import tempfile
import textwrap
from unittest import TestCase
from uuid import uuid4

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
        self.jobStore = os.getenv('TOIL_SCRIPTS_TEST_JOBSTORE', os.path.join(self.workdir, 'jobstore-%s' % uuid4()))
        self.toilOptions = shlex.split(os.environ.get('TOIL_SCRIPTS_TEST_TOIL_OPTIONS', ''))
        self.base_command = concat('toil-germline', 'run',
                                   self.toilOptions,
                                   self.jobStore)

    def test_pipeline_with_vqsr(self):
        """
        Tests entire pipeline with VQSR using BAM input.

        Skips HaplotypeCaller step by swapping in a pre-cooked GVCF file.
        """
        expected_files = {'bam_test.processed.ci_test.bam',
                          'bam_test.ci_test.g.vcf',
                          'bam_test.genotyped.ci_test.vcf',
                          'bam_test.vqsr.ci_test.vcf',
                          'config-toil-germline.yaml'}

        inputs = self._get_default_inputs()
        inputs.run_bwa = True
        inputs.preprocess = True
        inputs.run_vqsr = True
        inputs.hc_output = 's3://cgl-pipeline-inputs/germline/ci/NA12878_21_Variants_CI_25.sep1VQSR.g.vcf'

        bam_sample = ['bam_test', 's3://cgl-pipeline-inputs/germline/ci/NIST7035_NIST7086.aln21.ci.bam']

        self._run(self.base_command,
                  '--sample', bam_sample,
                  '--config', self._generate_config(inputs))
        self._assertOutput(expected_files)

    def test_joint_genotype(self):
        """
        Aligns paired FASTQ files, joint genotypes, and hard filters.
        """
        num_samples = int(os.environ.get('TOIL_SCRIPTS_TEST_NUM_SAMPLES', '3'))
        expected_files = {'joint_genotyped.genotyped.ci_test.vcf',
                          'joint_genotyped.hard_filter.ci_test.vcf',
                          'config-toil-germline.yaml',
                          'manifest-toil-germline.tsv'}

        for i in range(1, num_samples+1):
            expected_files |= {'fastq_test_%s.processed.ci_test.bam' % i,
                               'fastq_test_%s.ci_test.g.vcf' % i}

        inputs = self._get_default_inputs()
        inputs.run_bwa = True
        inputs.preprocess = True
        inputs.joint_genotype = True

        self._run(self.base_command,
                  '--config', self._generate_config(inputs),
                  '--manifest', self._generate_manifest(num_samples))
        self._assertOutput(expected_files)

    def test_preprocess_only(self):
        """
        Tests preprocess only option
        """
        num_samples = int(os.environ.get('TOIL_SCRIPTS_TEST_NUM_SAMPLES', '3'))
        expected_files = {'config-toil-germline.yaml', 'manifest-toil-germline.tsv'}
        for i in range(1, num_samples+1):
            expected_files.add('fastq_test_%s.processed.ci_test.bam' % i)

        inputs = self._get_default_inputs()
        inputs.run_bwa = True
        inputs.preprocess = True
        inputs.preprocess_only = True

        self._run(self.base_command,
                  '--config', self._generate_config(inputs),
                  '--manifest', self._generate_manifest(num_samples),
                  '--preprocess-only')
        self._assertOutput(expected_files)

    def _run(self, *args):
        args = list(concat(*args))
        log.info('Running %r', args)
        subprocess.check_call(args)

    def tearDown(self):
        shutil.rmtree(self.workdir)

    def _generate_config(self, inputs=None):
        if inputs is None:
            inputs = self._get_default_inputs()

        path = os.path.join(self.workdir, 'config-toil-germline.yaml')
        with open(path, 'w') as f:
            f.write(textwrap.dedent("""
                    genome-fasta: {genome_fasta}
                    g1k-snp: {g1k_snp}
                    g1k-indel: {g1k_indel}
                    mills: {mills}
                    dbsnp: {dbsnp}
                    hapmap: {hapmap}
                    omni: {omni}
                    run-bwa: {run_bwa}
                    trim: False
                    preprocess: {preprocess}
                    amb: {amb}
                    ann: {ann}
                    bwt: {bwt}
                    pac: {pac}
                    sa: {sa}
                    alt:
                    snp-annotations: {snp_annotations}
                    indel-annotations: {indel_annotations}
                    ssec:
                    file-size: 1G
                    cores: 2
                    xmx: 5G
                    suffix: .ci_test
                    sorted:
                    output-dir: {output_dir}
                    unsafe-mode: False
                    run-vqsr: {run_vqsr}
                    joint-genotype: {joint_genotype}
                    preprocess-only:
                    run-oncotator: {run_oncotator}

                    # Special parameters for testing
                    hc-output: {hc_output}
                    """[1:]).format(**vars(inputs)))
        return path

    def _generate_manifest(self, num_samples):
        """
        Makes a Toil Germline manifest containing FASTQ sample data

        :param int num_samples: Number of replicate samples
        :return: Path to manifest
        :rtype: str
        """
        path = os.path.join(self.workdir, 'manifest-toil-germline.tsv')
        with open(path, 'w') as f:
            f.write('\n'.join('\t'.join(['fastq_test_%s' % i,
                                         self.fastq_url,
                                         '@RG\\tID:foo\\tSM:bar\\tPL:ILLUMINA'])
                              for i in range(1, num_samples + 1)))
        return path

    def _assertOutput(self, expected_files):
        """
        Checks that all output files are expected and non-zero in size
        """
        for root, dirs, files in os.walk(self.workdir, topdown=False):
            for name in files:
                self.assertTrue(name in expected_files)
                self.assertTrue(os.stat(os.path.join(root, name)).st_size > 0)

    def _get_default_inputs(self):
        """
        Creates a Namespace object with default parameters
        :return: Namespace object
        """
        inputs = argparse.Namespace()
        inputs.run_bwa = False
        inputs.preprocess = False
        inputs.run_vqsr = False
        inputs.run_oncotator = False
        inputs.joint_genotype = False
        inputs.ssec = None
        inputs.sorted = False
        inputs.cores = 4
        inputs.xmx = '8G'
        inputs.output_dir = self.workdir
        inputs.suffix = ''
        inputs.unsafe_mode = False
        inputs.genome_fasta = 's3://cgl-pipeline-inputs/germline/ci/b37_21.fa'
        inputs.g1k_snp = 's3://cgl-pipeline-inputs/germline/ci/1000G_phase1.snps.high_confidence.b37.21.recode.vcf'
        inputs.g1k_indel = 's3://cgl-pipeline-inputs/germline/ci/1000G_phase1.indels.b37.21.recode.vcf'
        inputs.mills = 's3://cgl-pipeline-inputs/germline/ci/Mills_and_1000G_gold_standard.indels.b37.21.recode.vcf'
        inputs.dbsnp = 's3://cgl-pipeline-inputs/germline/ci/dbsnp_138.b37.21.recode.vcf'
        inputs.hapmap = 's3://cgl-pipeline-inputs/germline/ci/hapmap_3.3.b37.21.recode.vcf'
        inputs.omni = 's3://cgl-pipeline-inputs/germline/ci/1000G_omni2.5.b37.21.recode.vcf'
        inputs.amb = 's3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.amb'
        inputs.ann = 's3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.ann'
        inputs.bwt = 's3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.bwt'
        inputs.pac = 's3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.pac'
        inputs.sa = 's3://cgl-pipeline-inputs/germline/ci/bwa_index_b37_21.sa'
        inputs.snp_annotations = ['QualByDepth', 'FisherStrand', 'StrandOddsRatio',
                                  'ReadPosRankSumTest', 'MappingQualityRankSumTest',
                                  'RMSMappingQuality']
        inputs.indel_annotations = ['QualByDepth', 'FisherStrand', 'StrandOddsRatio',
                                    'ReadPosRankSumTest', 'MappingQualityRankSumTest']
        # Special attributes for testing
        inputs.hc_output = ''
        return inputs

