#!/usr/bin/env python2.7
import os

from bd2k.util.humanize import human2bytes
from toil.job import PromisedRequirement

from toil_scripts.gatk_germline.common import output_file_job
from toil_scripts.tools.variant_filters import gatk_select_variants, \
    gatk_variant_filtration, gatk_combine_variants


def hard_filter_pipeline(job, uuid, vcf_id, config):
    """
    Runs GATK Hard Filtering on a Genomic VCF file and uploads the results.

    0: Start                0 --> 1 --> 3 --> 5 --> 6
    1: Select SNPs                |           |
    2: Select INDELs              +-> 2 --> 4 +
    3: Apply SNP Filter
    4: Apply INDEL Filter
    5: Merge SNP and INDEL VCFs
    6: Write filtered VCF to output directory

    :param job: Toil Job instance
    :param str uuid: Unique sample identifier
    :param str vcf_id: VCF FileStoreID
    :param Namespace config: Pipeline configuration options and shared files
    :return: SNP and INDEL FileStoreIDs
    :rtype: tuple
    """
    job.fileStore.logToMaster('Running Hard Filter on {}'.format(uuid))

    # Estimate all disk resources as twice the input VCF size plus 5G for reference genome data.
    select_snps_disk = PromisedRequirement(lambda vcf: 2 * vcf.size + human2bytes('5G'),
                                           vcf_id)
    select_snps = job.wrapJobFn(gatk_select_variants,
                                'SNP',
                                vcf_id,
                                config,
                                memory=config.xmx, disk=select_snps_disk)

    snp_filter_disk = PromisedRequirement(lambda vcf: 2 * vcf.size + human2bytes('5G'),
                                          select_snps.rv())
    snp_filter = job.wrapJobFn(gatk_variant_filtration,
                               'SNP',
                               select_snps.rv(),
                               config,
                               memory=config.xmx, disk=snp_filter_disk)

    select_indels_disk = PromisedRequirement(lambda vcf: 2 * vcf.size + human2bytes('5G'),
                                             vcf_id)
    select_indels = job.wrapJobFn(gatk_select_variants,
                                  'INDEL',
                                  vcf_id,
                                  config,
                                  memory=config.xmx, disk=select_indels_disk)

    indel_filter_disk = PromisedRequirement(lambda vcf: 2 * vcf.size + human2bytes('5G'),
                                            select_indels.rv())
    indel_filter = job.wrapJobFn(gatk_variant_filtration,
                                 'INDEL',
                                 select_indels.rv(),
                                 config,
                                 memory=config.xmx, disk=indel_filter_disk)

    combine_variants_disk = PromisedRequirement(
        lambda vcf1, vcf2: 2 * (vcf1.size + vcf2.size) + human2bytes('5G'),
        indel_filter.rv(), snp_filter.rv())
    merge_variants = job.wrapJobFn(gatk_combine_variants,
                                   {'SNPs': snp_filter.rv(), 'INDELs': indel_filter.rv()},
                                   config,
                                   memory=config.xmx, disk=combine_variants_disk)

    job.addChild(select_snps)
    job.addChild(select_indels)

    select_snps.addChild(snp_filter)
    snp_filter.addChild(merge_variants)

    select_indels.addChild(indel_filter)
    indel_filter.addChild(merge_variants)

    # Output the hard filtered VCF
    output_dir = os.path.join(config.output_dir, uuid)
    output_filename = '%s.hard_filter%s.vcf' % (uuid, config.suffix)
    output_vcf = job.wrapJobFn(output_file_job,
                               output_filename,
                               merge_variants.rv(),
                               output_dir, s3_key_path=config.ssec)
    merge_variants.addChild(output_vcf)
    return merge_variants.rv()

def main():
    """
    Simple command line interface to run hard filtering on a single sample.
    """
    import argparse

    from toil.job import Job

    from toil_scripts.gatk_germline.germline import download_shared_files
    from toil_scripts.lib.urls import download_url_job

    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--sample',
                        default=None,
                        nargs=2,
                        type=str,
                        help='Space delimited sample UUID and GVCF file in the format: uuid url')

    parser.add_argument('--genome-fasta',
                        required=True,
                        type=str,
                        help='Path or URL to genome fasta file')

    parser.add_argument('--output-dir',
                        default=None,
                        help='Path/URL to output directory')

    Job.Runner.addToilOptions(parser)
    options = parser.parse_args()
    options.cores = 8
    options.xmx = '10G'
    options.run_bwa = False
    options.preprocess = False
    options.run_vqsr = False
    options.run_oncotator = False
    options.ssec = None
    options.suffix = ''
    options.unsafe_mode = False
    options.annotations = ['QualByDepth',
                           'FisherStrand',
                           'StrandOddsRatio',
                           'ReadPosRankSumTest',
                           'MappingQualityRankSumTest',
                           'RMSMappingQuality',
                           'InbreedingCoeff']

    uuid, url = options.sample

    shared_files = Job.wrapJobFn(download_shared_files, options).encapsulate()
    download_sample = shared_files.addFollowOnJobFn(download_url_job,
                                                    url,
                                                    name='toil.g.vcf',
                                                    synapse_login=None,
                                                    s3_key_path=None,
                                                    disk='25G')

    download_sample.addFollowOnJobFn(hard_filter_pipeline,
                                     uuid,
                                     download_sample.rv(),
                                     shared_files.rv())

    Job.Runner.startToil(shared_files, options)

if __name__ == '__main__':
    main()

