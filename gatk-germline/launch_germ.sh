#!/usr/bin/env bash
export TMPDIR=/home/ubuntu/fake_mnt
python germ.py \
$TMPDIR/jstore \
--config "config.txt" \
--reference "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/hs37d5.fa" \
--phase "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/1000G_phase1.indels.b37.vcf" \
--mills "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/Mills_and_1000G_gold_standard.indels.b37.vcf" \
--dbsnp "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/dbsnp_138.b37.vcf" \
--hapmap "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/hapmap_3.3.b37.vcf" \
--omni "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/1000G_omni2.5.b37.vcf" \
-o '/home/ubuntu/' \
