#!/usr/bin/env python2.7
from __future__ import print_function
import os

from toil.job import PromisedRequirement as PrR
from bd2k.util.humanize import human2bytes

from toil_scripts.gatk_germline.common import output_file_job
from toil_scripts.tools.variant_filters import gatk_variant_recalibrator, \
    gatk_apply_variant_recalibration


def vqsr_pipeline(job, uuid, vcf_id, config):
    """
    Runs GATK Variant Quality Score Recalibration.

    0: Start                        0 --> 1 --> 3 --> 4 --> 5
    1: Recalibrate SNPs                   |      |
    2: Recalibrate INDELS                 +-> 2 -+
    3: Apply SNP Recalibration
    4: Apply INDEL Recalibration
    5: Write VCF to output directory

    :param JobFunctionWrappingJob job: Toil Job instance
    :param str uuid: unique sample identifier
    :param str vcf_id: VCF FileStoreID
    :param Namespace config: Input parameters
    :return: SNP and INDEL VQSR VCF FileStoreID
    :rtype: str
    """
    # Get the total size of the genome reference
    genome_ref_size = config.genome_fasta.size + config.genome_fai.size + config.genome_dict.size

    # The VariantRecalibator disk requirement depends on the input VCF, the resource files,
    # the genome reference files, and the output recalibration table, tranche file, and plots.
    # The sum of these output files are less than the input VCF.
    snp_resources = ['hapmap', 'omni', 'dbsnp', 'phase']
    snp_resource_size = sum(getattr(config, resource).size for resource in snp_resources)
    snp_recal_disk = PrR(lambda in_vcf, ref_size, resource_size:
                         2 * in_vcf.size + ref_size + resource_size,
                         vcf_id,
                         genome_ref_size,
                         snp_resource_size)
    snp_recal = job.wrapJobFn(gatk_variant_recalibrator,
                              'SNP',
                              vcf_id,
                              config,
                              disk=snp_recal_disk,
                              cores=config.cores,
                              memory=config.xmx)

    indel_resource_size = config.mills.size + config.dbsnp.size
    indel_recal_disk = PrR(lambda in_vcf, ref_size, resource_size:
                           2 * in_vcf.size + ref_size + resource_size,
                           vcf_id,
                           genome_ref_size,
                           indel_resource_size)
    indel_recal = job.wrapJobFn(gatk_variant_recalibrator,
                                'INDEL',
                                vcf_id,
                                config,
                                disk=indel_recal_disk,
                                cores=config.cores,
                                memory=config.xmx)

    # The ApplyRecalibration disk requirement depends on the input VCF size, the variant
    # recalibration table, the tranche file, the genome reference file, and the output VCF.
    # This step labels variants as filtered, so the output VCF file should be slightly larger
    # than the input file. Estimate a 10% increase in the VCF file size.
    apply_snp_recal_disk = PrR(lambda in_vcf, recal, tranche, ref_size:
                               int(2.1 * in_vcf.size + recal.size + tranche.size + ref_size),
                               vcf_id,
                               snp_recal.rv(0),
                               snp_recal.rv(1),
                               genome_ref_size)
    apply_snp_recal = job.wrapJobFn(gatk_apply_variant_recalibration,
                                    'SNP',
                                    vcf_id,
                                    snp_recal.rv(0), snp_recal.rv(1),
                                    config,
                                    disk=apply_snp_recal_disk,
                                    cores=config.cores,
                                    memory=config.xmx)

    apply_indel_recal_disk = PrR(lambda in_vcf, recal, tranche, ref_size:
                                 int(2.1 * in_vcf.size + recal.size + tranche.size + ref_size),
                                 vcf_id,
                                 indel_recal.rv(0),
                                 indel_recal.rv(1),
                                 genome_ref_size)
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
    snp_recal.addChild(apply_snp_recal)
    indel_recal.addChild(apply_indel_recal)
    apply_snp_recal.addChild(apply_indel_recal)

    # Output input VCF and recalibrated VCF
    output_dir = config.output_dir
    output_dir = os.path.join(output_dir, uuid)
    vqsr_name = '%s.vqsr%s.vcf' % (uuid, config.suffix)
    output_vqsr = job.wrapJobFn(output_file_job,
                                vqsr_name,
                                apply_indel_recal.rv(),
                                output_dir,
                                s3_key_path=config.ssec,
                                disk=PrR(lambda x: x.size, apply_indel_recal.rv()))
    apply_indel_recal.addChild(output_vqsr)
    return apply_indel_recal.rv()


def main():
    """
    Runs GATK VQSR for a cohort of GVCFs
    """
    import argparse
    from multiprocessing import cpu_count

    from toil.job import Job
    import yaml

    from toil_scripts.gatk_germline.germline import download_shared_files
    from toil_scripts.lib.urls import download_url_job

    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--sample',
                        default=None,
                        nargs=2,
                        type=str,
                        help='Space delimited sample UUID and VCF file in the format: uuid url')

    parser.add_argument('--config',
                        required=True,
                        type=str,
                        help='Path or URL to Toil germline config file')

    parser.add_argument('--output-dir',
                        default=None,
                        help='Path/URL to output directory')

    Job.Runner.addToilOptions(parser)
    options = parser.parse_args()

    # Parse inputs
    inputs = {x.replace('-', '_'): y for x, y in
              yaml.load(open(options.config).read()).iteritems()}

    inputs = argparse.Namespace(**inputs)

    inputs.run_bwa = False
    inputs.preprocess = False
    inputs.run_oncotator = False
    inputs.annotations = ['QualByDepth',
                          'FisherStrand',
                          'StrandOddsRatio',
                          'ReadPosRankSumTest',
                          'MappingQualityRankSumTest',
                          'RMSMappingQuality',
                          'InbreedingCoeff']

    shared_files = Job.wrapJobFn(download_shared_files, inputs).encapsulate()

    uuid, url = options.sample
    gvcf = shared_files.addChildJobFn(download_url_job,
                                      url,
                                      name='toil.g.vcf',
                                      s3_key_path=None)

    shared_files.addFollowOnJobFn(vqsr_pipeline,
                                  uuid,
                                  gvcf.rv(),
                                  shared_files.rv(),
                                  cores=cpu_count())

    Job.Runner.startToil(shared_files, options)

if __name__ == '__main__':
    main()
