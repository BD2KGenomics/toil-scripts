#!/usr/bin/env python2.7
from __future__ import print_function
import os

from toil.job import PromisedRequirement
from toil_lib.tools.variant_manipulation import gatk_variant_recalibrator, \
    gatk_apply_variant_recalibration

from toil_scripts.gatk_germline.common import output_file_job


def vqsr_pipeline(job, uuid, vcf_id, config):
    """
    Runs GATK Variant Quality Score Recalibration.

    0: Start                        0 --> 1 --> 3 --> 4 --> 5
    1: Recalibrate SNPs                   |      |
    2: Recalibrate INDELS                 +-> 2 -+
    3: Apply SNP Recalibration
    4: Apply INDEL Recalibration
    5: Write VCF to output directory

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str uuid: unique sample identifier
    :param str vcf_id: VCF FileStoreID
    :param Namespace config: Pipeline configuration options and shared files
        Requires the following config attributes:
        config.genome_fasta             FilesStoreID for reference genome fasta file
        config.genome_fai               FilesStoreID for reference genome fasta index file
        config.genome_dict              FilesStoreID for reference genome sequence dictionary file
        config.cores                    Number of cores for each job
        config.xmx                      Java heap size in bytes
        config.suffix                   Suffix for output filename
        config.output_dir               URL or local path to output directory
        config.ssec                     Path to key file for SSE-C encryption

        SNP VQSR attributes:
        config.snp_filter_annotations   List of GATK variant annotations
        config.hapmap                   FileStoreID for HapMap resource file
        config.omni                     FileStoreID for Omni resource file
        config.dbsnp                    FileStoreID for dbSNP resource file
        config.g1k_snp                  FileStoreID for 1000G SNP resource file

        INDEL VQSR attributes:
        config.indel_filter_annotations List of GATK variant annotations
        config.dbsnp                    FileStoreID for dbSNP resource file
        config.mills                    FileStoreID for Mills resource file

    :return: SNP and INDEL VQSR VCF FileStoreID
    :rtype: str
    """
    # Get the total size of the genome reference
    genome_ref_size = config.genome_fasta.size + config.genome_fai.size + config.genome_dict.size

    # The VariantRecalibator disk requirement depends on the input VCF, the resource files,
    # the genome reference files, and the output recalibration table, tranche file, and plots.
    # The sum of these output files are less than the input VCF.
    snp_resources = ['hapmap', 'omni', 'dbsnp', 'g1k_snp']
    snp_resource_size = sum(getattr(config, resource).size for resource in snp_resources)
    snp_recal_disk = PromisedRequirement(lambda in_vcf, ref_size, resource_size:
                                         2 * in_vcf.size + ref_size + resource_size,
                                         vcf_id,
                                         genome_ref_size,
                                         snp_resource_size)

    snp_recal = job.wrapJobFn(gatk_variant_recalibrator,
                              'SNP',
                              vcf_id,
                              config.genome_fasta,
                              config.genome_fai,
                              config.genome_dict,
                              get_short_annotations(config.snp_filter_annotations),
                              hapmap=config.hapmap,
                              omni=config.omni,
                              phase=config.g1k_snp,
                              dbsnp=config.dbsnp,
                              unsafe_mode=config.unsafe_mode,
                              disk=snp_recal_disk,
                              cores=config.cores,
                              memory=config.xmx)

    indel_resource_size = config.mills.size + config.dbsnp.size
    indel_recal_disk = PromisedRequirement(lambda in_vcf, ref_size, resource_size:
                                           2 * in_vcf.size + ref_size + resource_size,
                                           vcf_id,
                                           genome_ref_size,
                                           indel_resource_size)

    indel_recal = job.wrapJobFn(gatk_variant_recalibrator,
                                'INDEL',
                                vcf_id,
                                config.genome_fasta,
                                config.genome_fai,
                                config.genome_dict,
                                get_short_annotations(config.indel_filter_annotations),
                                dbsnp=config.dbsnp,
                                mills=config.mills,
                                unsafe_mode=config.unsafe_mode,
                                disk=indel_recal_disk,
                                cores=config.cores,
                                memory=config.xmx)

    # The ApplyRecalibration disk requirement depends on the input VCF size, the variant
    # recalibration table, the tranche file, the genome reference file, and the output VCF.
    # This step labels variants as filtered, so the output VCF file should be slightly larger
    # than the input file. Estimate a 10% increase in the VCF file size.
    apply_snp_recal_disk = PromisedRequirement(lambda in_vcf, recal, tranche, ref_size:
                                               int(2.1 * in_vcf.size + recal.size + tranche.size + ref_size),
                                               vcf_id,
                                               snp_recal.rv(0),
                                               snp_recal.rv(1),
                                               genome_ref_size)

    apply_snp_recal = job.wrapJobFn(gatk_apply_variant_recalibration,
                                    'SNP',
                                    vcf_id,
                                    snp_recal.rv(0), snp_recal.rv(1),
                                    config.genome_fasta,
                                    config.genome_fai,
                                    config.genome_dict,
                                    unsafe_mode=config.unsafe_mode,
                                    disk=apply_snp_recal_disk,
                                    cores=config.cores,
                                    memory=config.xmx)

    apply_indel_recal_disk = PromisedRequirement(lambda in_vcf, recal, tranche, ref_size:
                                                 int(2.1 * in_vcf.size + recal.size + tranche.size + ref_size),
                                                 vcf_id,
                                                 indel_recal.rv(0),
                                                 indel_recal.rv(1),
                                                 genome_ref_size)

    apply_indel_recal = job.wrapJobFn(gatk_apply_variant_recalibration,
                                      'INDEL',
                                      apply_snp_recal.rv(),
                                      indel_recal.rv(0), indel_recal.rv(1),
                                      config.genome_fasta,
                                      config.genome_fai,
                                      config.genome_dict,
                                      unsafe_mode=config.unsafe_mode,
                                      disk=apply_indel_recal_disk,
                                      cores=config.cores,
                                      memory=config.xmx)

    job.addChild(snp_recal)
    job.addChild(indel_recal)
    snp_recal.addChild(apply_snp_recal)
    indel_recal.addChild(apply_indel_recal)
    apply_snp_recal.addChild(apply_indel_recal)

    # Output recalibrated VCF
    output_dir = config.output_dir
    output_dir = os.path.join(output_dir, uuid)
    vqsr_name = '%s.vqsr%s.vcf' % (uuid, config.suffix)
    output_vqsr = job.wrapJobFn(output_file_job,
                                vqsr_name,
                                apply_indel_recal.rv(),
                                output_dir,
                                s3_key_path=config.ssec,
                                disk=PromisedRequirement(lambda x: x.size, apply_indel_recal.rv()))
    apply_indel_recal.addChild(output_vqsr)
    return apply_indel_recal.rv()


def get_short_annotations(annotations):
    """
    Converts full GATK annotation name to the shortened version
    :param annotations:
    :return:
    """
    # Annotations need to match VCF header
    short_name = {'QualByDepth': 'QD',
                  'FisherStrand': 'FS',
                  'StrandOddsRatio': 'SOR',
                  'ReadPosRankSumTest': 'ReadPosRankSum',
                  'MappingQualityRankSumTest': 'MQRankSum',
                  'RMSMappingQuality': 'MQ',
                  'InbreedingCoeff': 'ID'}

    short_annotations = []
    for annotation in annotations:
        if annotation in short_name:
            annotation = short_name[annotation]
        short_annotations.append(annotation)
    return short_annotations
