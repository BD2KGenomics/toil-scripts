#!/usr/bin/env bash
export TMPDIR=/home/ubuntu/fake_mnt
python gatk_pipeline.py \
/home/ubuntu/fake_mnt/jstore \
--retryCount 3 \
--config "config.txt" \
--reference "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/Homo_sapiens_assembly19.fasta" \
--phase "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/1000G_phase1.indels.hg19.sites.fixed.vcf" \
--mills "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf" \
--dbsnp "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/dbsnp_132_b37.leftAligned.vcf" \
--cosmic "https://s3-us-west-2.amazonaws.com/cgl-variant-inputs/b37_cosmic_v54_120711.vcf" \
--output_dir '/home/ubuntu/' \
--s3_dir 'cgl-driver-projects/wcdt/variants/'
