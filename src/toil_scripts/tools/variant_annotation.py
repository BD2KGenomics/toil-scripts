#!/usr/bin/env python2.7
from __future__ import print_function
import os

from toil_scripts.lib.files import get_files_from_filestore
from toil_scripts.lib.programs import docker_call


def run_oncotator(job, vcf_id, config):
    """
    Runs Oncotator on VCF file. This tool only uses HG19 annotations.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str vcf_id: VCF FileStoreID
    :param Namespace config: Pipeline configuration and shared FileStoreIDs
    :return: Annotated VCF FileStoreID
    :rtype: str
    """
    job.fileStore.logToMaster('Running Oncotator')
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'oncotator_db': config.oncotator_db,
              'input.vcf': vcf_id}
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-i', 'VCF',
               '-o', 'VCF',
               '--db-dir', os.path.basename(inputs['oncotator_db']),
               'input.vcf',
               'annotated.vcf',
               'hg19']

    outputs={'annotated.vcf': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters = command,
                tool = 'jpfeil/oncotator:1.9--8fffc356981862d50cfacd711b753700b886b605',
                inputs=inputs.keys(),
                outputs=outputs)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'annotated.vcf'))


def gatk_variant_annotator(job, bam_id, bai_id, vcf_id, annotations, config):
    """
    Annotates a VCF file using GATK Annotations

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str bam: Sample BAM FileStoreID
    :param str bai: Bam Index FileStoreID
    :param str gvcf_id: GVCF FileStoreID
    :param list annotations: GATK Annotations
    :param Namespace config: Pipeline configuration and shared FileStoreIDs
    :return: Annotated VCF FileStoreID
    :rtype: str
    """
    job.fileStore.logToMaster('Running GATK VariantAnnotator: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    inputs['input.bam'] = bam_id
    inputs['input.bam.bai'] = bai_id
    inputs['input.vcf'] = vcf_id
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-R', 'genome.fa',
               '-T', 'VariantAnnotator',
               '-I', 'input.bam',
               '-V', 'input.vcf',
               '-o', 'output.vcf']

    for annotation in annotations:
        command.extend(['-A', annotation])

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs={'output.vcf': None}
    docker_call(work_dir = work_dir,
                env={'JAVA_OPTS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.vcf'))
