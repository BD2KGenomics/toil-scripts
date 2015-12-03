#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Precautionary step: Create location where jobStore and tmp files will exist
mkdir -p ${HOME}/toil_mnt
# Execution of pipeline
python exome_variant_pipeline.py \
${HOME}/toil_mnt/jobStore \
--retryCount 1 \
--config 'exome_variant_config.csv' \
--reference 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/Homo_sapiens_assembly19.fasta' \
--phase 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/1000G_phase1.indels.hg19.sites.fixed.vcf' \
--mills 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf' \
--dbsnp 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/dbsnp_132_b37.leftAligned.vcf' \
--cosmic 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/b37_cosmic_v54_120711.vcf' \
--output_dir ${HOME} \
--ssec ${HOME}/master.key \
--s3_dir 'cgl-driver-projects/test/variants/' \
--workDir ${HOME}/toil_mnt
--sudo \
#--restart
