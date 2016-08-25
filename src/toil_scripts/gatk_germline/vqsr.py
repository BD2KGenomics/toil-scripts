#!/usr/bin/env python2.7
from __future__ import print_function
import os

from toil.job import PromisedRequirement
from bd2k.util.humanize import human2bytes

from toil_scripts.tools.variant_annotation import gatk_genotype_gvcfs
from toil_scripts.tools.variant_filters import gatk_variant_recalibrator, \
    gatk_apply_variant_recalibration
from toil_scripts.lib.files import upload_or_move_job


def vqsr_pipeline(job, gvcfs, config):
    """
    Runs GATK Variant Quality Score Recalibration. Writes VQSR VCF to an output directory defined
    in the config dictionary. Multiple GVCFs are joined and recalibrated into a single VCF.

    0: GenotypeGVCFs                0 --> 1 --> 3 --> 4 --> 5
    1: Recalibrate SNPs                   |      |
    2: Recalibrate INDELS                 +-> 2 -+
    3: Apply SNP Recalibration
    4: Apply INDEL Recalibration
    5: Write VCF to output directory

    :param JobFunctionWrappingJob job: Toil Job instance
    :param dict gvcfs: Dictionary of UUIDs and GVCF FileStoreIDs
    :param Namespace config: Input parameters
    :return: SNP and INDEL VQSR VCF FileStoreID
    :rtype: str
    """
    # Estimate disk resource as twice the sum of all input GVCFs plus 5G for reference genome data
    genotype_gvcf_disk = PromisedRequirement(
        lambda lst: 2 * sum(x.size for x in lst) + human2bytes('5G'),
        gvcfs.values())
    genotype_gvcf = job.wrapJobFn(gatk_genotype_gvcfs,
                                  gvcfs,
                                  config,
                                  disk=genotype_gvcf_disk, cores=config.cores)

    # Estimate disk resource as twice the input VCF plus 10G for genome data and SNP databases
    snp_recal_disk = PromisedRequirement(lambda x: 2 * x.size + human2bytes('10G'),
                                         genotype_gvcf.rv())
    snp_recal = job.wrapJobFn(gatk_variant_recalibrator,
                              'SNP',
                              genotype_gvcf.rv(),
                              config,
                              disk=snp_recal_disk, cores=config.cores)

    # Estimate disk resource as twice the input VCF plus 10G for genome data and INDEL databases
    indel_recal_disk = PromisedRequirement(lambda x: 2 * x.size + human2bytes('10G'),
                                           genotype_gvcf.rv())
    indel_recal = job.wrapJobFn(gatk_variant_recalibrator,
                                'INDEL',
                                genotype_gvcf.rv(),
                                config,
                                disk=indel_recal_disk, cores=config.cores, memory=config.xmx)

    # Estimate disk resource as the twice the sum of the input VCF, the recalibration table
    # and the tranche file plus 5G for reference genome data.
    apply_snp_recal_disk = PromisedRequirement(
        lambda x, y, z: 2 * (x.size + y.size + z.size) + human2bytes('5G'),
        genotype_gvcf.rv(),
        snp_recal.rv(0),
        snp_recal.rv(1))
    apply_snp_recal = job.wrapJobFn(gatk_apply_variant_recalibration,
                                    'SNP',
                                    genotype_gvcf.rv(),
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

    job.addFollowOn(genotype_gvcf)
    genotype_gvcf.addChild(snp_recal)
    genotype_gvcf.addChild(indel_recal)
    job.addFollowOn(apply_snp_recal)
    apply_snp_recal.addFollowOn(apply_indel_recal)

    # Output input VCF and recalibrated VCF
    output_dir = config.output_dir
    if len(gvcfs) == 1:
        uuid = gvcfs.keys()[0]
        output_dir = os.path.join(output_dir, uuid)
        input_filename = '%s.genotyped%s.vcf' % (uuid, config.suffix)
        vqsr_name = '%s.vqsr%s.vcf' % (uuid, config.suffix)
    else:
        input_filename = 'joint.genotyped%s.vcf' % config.suffix
        vqsr_name = 'joint.vqsr%s.vcf' % config.suffix

    output_input = job.wrapJobFn(upload_or_move_job,
                                 input_filename,
                                 genotype_gvcf.rv(),
                                 output_dir,
                                 s3_key_path=config.ssec)
    output_vqsr = job.wrapJobFn(upload_or_move_job,
                                vqsr_name,
                                apply_indel_recal.rv(),
                                output_dir,
                                s3_key_path=config.ssec)
    apply_indel_recal.addChild(output_input)
    apply_indel_recal.addChild(output_vqsr)
    return apply_indel_recal.rv()
