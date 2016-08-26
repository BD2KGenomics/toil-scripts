#!/usr/bin/env python2.7
from __future__ import print_function
import os

from toil.job import PromisedRequirement
from bd2k.util.humanize import human2bytes

from toil_scripts.tools.variant_filters import gatk_variant_recalibrator, \
    gatk_apply_variant_recalibration
from toil_scripts.lib.files import upload_or_move_job


def vqsr_pipeline(job, uuid, vcf_id, config):
    """
    Runs GATK Variant Quality Score Recalibration. Writes VQSR VCF to an output directory defined
    in the config dictionary. Multiple GVCFs are joined and recalibrated into a single VCF.

    0: Start                        0 --> 1 --> 3 --> 4 --> 5
    1: Recalibrate SNPs                   |      |
    2: Recalibrate INDELS                 +-> 2 -+
    3: Apply SNP Recalibration
    4: Apply INDEL Recalibration
    5: Write VCF to output directory

    :param JobFunctionWrappingJob job: Toil Job instance
    :param str vcf_id: VCF FileStoreID
    :param Namespace config: Input parameters
    :return: SNP and INDEL VQSR VCF FileStoreID
    :rtype: str
    """
    # Estimate disk resource as twice the input VCF plus 10G for genome data and SNP databases
    snp_recal_disk = PromisedRequirement(lambda x: 2 * x.size + human2bytes('10G'), vcf_id)
    snp_recal = job.wrapJobFn(gatk_variant_recalibrator,
                              'SNP',
                              vcf_id,
                              config,
                              disk=snp_recal_disk, cores=config.cores)

    # Estimate disk resource as twice the input VCF plus 10G for genome data and INDEL databases
    indel_recal_disk = PromisedRequirement(lambda x: 2 * x.size + human2bytes('10G'), vcf_id)
    indel_recal = job.wrapJobFn(gatk_variant_recalibrator,
                                'INDEL',
                                vcf_id,
                                config,
                                disk=indel_recal_disk, cores=config.cores, memory=config.xmx)

    # Estimate disk resource as the twice the sum of the input VCF, the recalibration table
    # and the tranche file plus 5G for reference genome data.
    apply_snp_recal_disk = PromisedRequirement(
        lambda x, y, z: 2 * (x.size + y.size + z.size) + human2bytes('5G'),
        vcf_id,
        snp_recal.rv(0),
        snp_recal.rv(1))
    apply_snp_recal = job.wrapJobFn(gatk_apply_variant_recalibration,
                                    'SNP',
                                    vcf_id,
                                    snp_recal.rv(0), snp_recal.rv(1),
                                    config,
                                    disk=apply_snp_recal_disk,
                                    cores=config.cores,
                                    memory=config.xmx)

    apply_indel_recal_disk = PromisedRequirement(
        lambda x, y, z: 2 * (x.size + y.size + z.size) + human2bytes('5G'),
        apply_snp_recal.rv(), indel_recal.rv(0), indel_recal.rv(1))
    apply_indel_recal = job.wrapJobFn(gatk_apply_variant_recalibration,
                                      'INDEL',
                                      apply_snp_recal.rv(),
                                      indel_recal.rv(0), indel_recal.rv(1),
                                      config,
                                      disk=apply_indel_recal_disk,
                                      cores=config.cores,
                                      memory=config.xmx)

    job.addChild(snp_recal)
    job.addChild(indel_recal)
    job.addFollowOn(apply_snp_recal)
    apply_snp_recal.addFollowOn(apply_indel_recal)

    # Output input VCF and recalibrated VCF
    output_dir = config.output_dir
    output_dir = os.path.join(output_dir, uuid)
    vqsr_name = '%s.vqsr%s.vcf' % (uuid, config.suffix)
    output_vqsr = job.wrapJobFn(upload_or_move_job,
                                vqsr_name,
                                apply_indel_recal.rv(),
                                output_dir,
                                s3_key_path=config.ssec)
    apply_indel_recal.addChild(output_vqsr)
    return apply_indel_recal.rv()
