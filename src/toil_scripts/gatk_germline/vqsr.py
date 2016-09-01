#!/usr/bin/env python2.7
from __future__ import print_function
import os

from toil.job import PromisedRequirement
from bd2k.util.humanize import human2bytes

from toil_scripts.gatk_germline.common import output_file_job
from toil_scripts.tools.variant_filters import gatk_variant_recalibrator, \
    gatk_apply_variant_recalibration


def vqsr_pipeline(job, uuid, vcf_id, config):
    """
    Runs GATK Variant Quality Score Recalibration. Writes VQSR VCF to an output directory defined
    in the config dictionary. Multiple GVCFs are joint genotyped and recalibrated into a single VCF.

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
                              disk=snp_recal_disk, cores=config.cores, memory=config.xmx)

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
    output_vqsr = job.wrapJobFn(output_file_job,
                                vqsr_name,
                                apply_indel_recal.rv(),
                                output_dir,
                                s3_key_path=config.ssec)
    apply_indel_recal.addChild(output_vqsr)
    return apply_indel_recal.rv()


def main():
    """
    Simple command line interface to run VQSR pipeline
    """
    import argparse

    from toil.job import Job
    import yaml

    from toil_scripts.gatk_germline.germline import download_shared_files, joint_genotype_and_filter
    from toil_scripts.lib.urls import download_url_job

    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--manifest',
                        required=True,
                        type=str,
                        help='Path to UUID and GVCF file in the format: uuid url')

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

    inputs.cores = 8
    inputs.xmx = '10G'
    inputs.run_bwa = False
    inputs.preprocess = False
    inputs.joint = True
    inputs.run_vqsr = True
    inputs.run_oncotator = False
    inputs.ssec = None
    inputs.suffix = ''
    inputs.unsafe_mode = False
    inputs.annotations = ['QualByDepth',
                           'FisherStrand',
                           'StrandOddsRatio',
                           'ReadPosRankSumTest',
                           'MappingQualityRankSumTest',
                           'RMSMappingQuality',
                           'InbreedingCoeff']

    shared_files = Job.wrapJobFn(download_shared_files, inputs).encapsulate()

    gvcfs = {}

    with open(options.manifest, 'r') as f:
        for line in f:
            uuid, url = line.strip().split()
            gvcfs[uuid] = shared_files.addChildJobFn(download_url_job,
                                                     url,
                                                     name='toil.g.vcf',
                                                     s3_key_path=None).rv()

    shared_files.addFollowOnJobFn(joint_genotype_and_filter,
                                  gvcfs.items(),
                                  shared_files.rv(),
                                  cores=8)

    Job.Runner.startToil(shared_files, options)

if __name__ == '__main__':
    main()
