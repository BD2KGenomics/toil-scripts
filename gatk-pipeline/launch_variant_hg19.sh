#!/usr/bin/env bash
# John Vivian
#
# This script provides all configuration necessary to run the Toil pipeline locally.
# It assumes there is a local file: exome_variant_config.csv.  One sample per line: uuid,url_normal,url_tumor.
#
# If --ssec is used, the program assumes input files are encrypted in S3 when retrieving them.
# If --sudo flag is used, 'sudo' will be prepended to the Docker subprocess call
# If --output_dir is used, the final VCF will be placed in the directory specified
# If --s3_dir is used, the final VCF will be uploaded to S3 using S3AM (pip install --pre s3am, need ~/.boto)
#
# Modify TMPDIR parameter to change location of tmp files.
# Modify first argument to change location of the local fileStore
# Uncomment the final line to resume your Toil job in the event of job failure.
mkdir -p ${HOME}/toil_mnt
export TMPDIR=${HOME}/toil_mnt
python gatk_pipeline.py \
${HOME}/toil_mnt/jstore \
--retryCount 3 \
--config "exome_variant_config.csv" \
--reference "https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa" \
--phase "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/1000G_phase1.indels.hg19.sites.vcf" \
--mills "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf" \
--dbsnp "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/dbsnp_138.hg19.vcf" \
--cosmic "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/Cosmic.hg19" \
--output_dir '${HOME}/' \
--ssec '${HOME}/master.key' \
--s3_dir 'cgl-driver-projects/wcdt/variants/' \
--sudo \
#--restart