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

    # Get the total size of the genome reference
    genome_ref_size = config.genome_fasta.size + config.genome_fai.size + config.genome_dict.size

    # The SelectVariants disk requirement depends on the input VCF, the genome reference files,
    # and the output VCF. The output VCF is smaller than the input VCF. The disk requirement
    # is identical for SNPs and INDELs.
    select_variants_disk = PromisedRequirement(lambda vcf, ref_size: 2 * vcf.size + ref_size,
                                               vcf_id,
                                               genome_ref_size)
    select_snps = job.wrapJobFn(gatk_select_variants,
                                'SNP',
                                vcf_id,
                                config,
                                memory=config.xmx,
                                disk=select_variants_disk)

    # The VariantFiltration disk requirement depends on the input VCF, the genome reference files,
    # and the output VCF. The filtered VCF is smaller than the input VCF.
    snp_filter_disk = PromisedRequirement(lambda vcf, ref_size: 2 * vcf.size + ref_size,
                                          select_snps.rv(),
                                          genome_ref_size)
    snp_filter = job.wrapJobFn(gatk_variant_filtration,
                               'SNP',
                               select_snps.rv(),
                               config,
                               memory=config.xmx,
                               disk=snp_filter_disk)

    select_indels = job.wrapJobFn(gatk_select_variants,
                                  'INDEL',
                                  vcf_id,
                                  config,
                                  memory=config.xmx,
                                  disk=select_variants_disk)

    indel_filter_disk = PromisedRequirement(lambda vcf, ref_size: 2 * vcf.size + ref_size,
                                            select_indels.rv(),
                                            genome_ref_size)
    indel_filter = job.wrapJobFn(gatk_variant_filtration,
                                 'INDEL',
                                 select_indels.rv(),
                                 config,
                                 memory=config.xmx,
                                 disk=indel_filter_disk)

    # The CombineVariants disk requirement depends on the SNP and INDEL input VCFs and the
    # genome reference files. The combined VCF is approximately the same size as the input files.
    combine_vcfs_disk = PromisedRequirement(lambda vcf1, vcf2, ref_size:
                                            2 * (vcf1.size + vcf2.size) + ref_size,
                                            indel_filter.rv(),
                                            snp_filter.rv(),
                                            genome_ref_size)
    combine_vcfs = job.wrapJobFn(gatk_combine_variants,
                                 {'SNPs': snp_filter.rv(), 'INDELs': indel_filter.rv()},
                                 config,
                                 merge_option='UNSORTED',  # Merges variants from a single sample
                                 memory=config.xmx,
                                 disk=combine_vcfs_disk)

    job.addChild(select_snps)
    job.addChild(select_indels)

    select_snps.addChild(snp_filter)
    snp_filter.addChild(combine_vcfs)

    select_indels.addChild(indel_filter)
    indel_filter.addChild(combine_vcfs)

    # Output the hard filtered VCF
    output_dir = os.path.join(config.output_dir, uuid)
    output_filename = '%s.hard_filter%s.vcf' % (uuid, config.suffix)
    output_vcf = job.wrapJobFn(output_file_job,
                               output_filename,
                               combine_vcfs.rv(),
                               output_dir,
                               s3_key_path=config.ssec,
                               disk=PromisedRequirement(lambda x: x.size, combine_vcfs.rv()))
    combine_vcfs.addChild(output_vcf)
    return combine_vcfs.rv()

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

