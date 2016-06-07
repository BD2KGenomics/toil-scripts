#!/usr/bin/env python2.7
from __future__ import print_function
import os

from bd2k.util.humanize import human2bytes
from toil.job import PromisedRequirement

from toil_scripts.lib.files import upload_or_move_job
from toil_scripts.tools.variant_filters import gatk_genotype_gvcf, gatk_select_variants, \
    gatk_variant_filtration


def hard_filter_pipeline(job, uuid, gvcf_id, config):
    """
    Runs GATK Hard Filtering on HaplotypeCaller output. Writes separate SNP and INDEL
    VCF files to output directory.

    0: GenotypeGVCFs            0 --> 1 --> 3 --> 5
    1: Select SNPs                |
    2: Select INDELs              +-> 2 --> 4 --> 5
    3: Apply SNP Filter
    4: Apply INDEL Filter
    5: Write filtered VCFs to output directory

    :param job: Toil Job instance
    :param str uuid: Unique sample identifier
    :param str gvcf_id: GVCF FileStoreID
    :param Namespace config: Pipeline configuration options and shared files
    :return: SNP and INDEL FileStoreIDs
    :rtype: tuple
    """
    job.fileStore.logToMaster('Running Hard Filter on {}'.format(uuid))

    genotype_gvcf_disk = PromisedRequirement(lambda x: 2*x.size + human2bytes('5G'), gvcf_id)
    genotype_gvcf = job.wrapJobFn(gatk_genotype_gvcf,
                                  dict(uuid=gvcf_id),
                                  config,
                                  memory=config.xmx, disk=genotype_gvcf_disk)

    select_snps_disk = PromisedRequirement(lambda x: 2*x.size + human2bytes('5G'),
                                           genotype_gvcf.rv())
    select_snps = job.wrapJobFn(gatk_select_variants,
                                'SNP',
                                genotype_gvcf.rv(),
                                config,
                                disk=select_snps_disk)

    snp_filter_disk = PromisedRequirement(lambda x: 2*x.size + human2bytes('5G'),
                                  select_snps.rv())
    snp_filter = job.wrapJobFn(gatk_variant_filtration,
                               'SNP',
                               select_snps.rv(),
                               config,
                               disk=snp_filter_disk)

    select_indels_disk = PromisedRequirement(lambda x: 2*x.size + human2bytes('5G'),
                                  genotype_gvcf.rv())
    select_indels = job.wrapJobFn(gatk_select_variants,
                                  'INDEL',
                                  genotype_gvcf.rv(),
                                  config,
                                  disk=select_indels_disk)

    indel_filter_disk = PromisedRequirement(lambda x: 2*x.size + human2bytes('5G'),
                                  select_snps.rv())
    indel_filter = job.wrapJobFn(gatk_variant_filtration,
                                 'INDEL',
                                 select_indels.rv(),
                                 config,
                                 disk=indel_filter_disk)

    job.addChild(genotype_gvcf)
    genotype_gvcf.addChild(select_snps)
    genotype_gvcf.addChild(select_indels)
    select_snps.addChild(snp_filter)
    select_indels.addChild(indel_filter)

    output_dir = os.path.join(config.output_dir, uuid)
    snp_filename = '%s_filtered_snps%s.vcf' % (uuid, config.suffix)
    indel_filename = '%s_filtered_indels%s.vcf' % (uuid, config.suffix)
    output_snps = job.wrapJobFn(upload_or_move_job, snp_filename, snp_filter.rv(),
                               output_dir, ssec=config.ssec)
    output_indels = job.wrapJobFn(upload_or_move_job, indel_filename, indel_filter.rv(),
                                  output_dir, ssec=config.ssec)
    snp_filter.addChild(output_snps)
    indel_filter.addChild(output_indels)
    return snp_filter.rv(), indel_filter.rv()
