#!/usr/bin/env bash
export TMPDIR=/home/ubuntu/fake_mnt
python gatk_pipeline.py \
/home/ubuntu/fake_mnt/jstore \
--retryCount 3 \
--config "config.txt" \
--reference "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/hg19.fa" \
--phase "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/1000G_phase1.indels.hg19.sites.vcf" \
--mills "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf" \
--dbsnp "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/dbsnp_138.hg19.vcf" \
--cosmic "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/Cosmic.hg19" \
--output_dir '/home/ubuntu/' \
--ssec '/home/ubuntu/master.key' \
--s3_dir 'cgl-driver-projects/wcdt/variants/'
