#!/usr/bin/env bash
export TMPDIR=/mnt
mkdir -p /mnt
python gpm.py \
--jobTree="aws:us-west-2:gatk-test-aba0325f-0a10-43b4-930b-0f4ea6d89259" \
--reference "https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/Homo_sapiens_assembly19.fasta" \
--n6 "https://s3-us-west-2.amazonaws.com/bd2k-test-data/pair6.normal.bam" \
--t6 "https://s3-us-west-2.amazonaws.com/bd2k-test-data/pair6.tumor.bam" \
--n7 "https://s3-us-west-2.amazonaws.com/bd2k-test-data/pair7.normal.bam" \
--t7 "https://s3-us-west-2.amazonaws.com/bd2k-test-data/pair7.tumor.bam" \
--n8 "https://s3-us-west-2.amazonaws.com/bd2k-test-data/pair8.normal.bam" \
--t8 "https://s3-us-west-2.amazonaws.com/bd2k-test-data/pair8.tumor.bam" \
--n9 "https://s3-us-west-2.amazonaws.com/bd2k-test-data/pair9.normal.bam" \
--t9 "https://s3-us-west-2.amazonaws.com/bd2k-test-data/pair9.tumor.bam" \
--phase "https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/1000G_phase1.indels.hg19.sites.fixed.vcf" \
--mills "https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf" \
--dbsnp "https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/dbsnp_132_b37.leftAligned.vcf" \
--cosmic "https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/b37_cosmic_v54_120711.vcf" \
--output_dir '/home/ubuntu/'
