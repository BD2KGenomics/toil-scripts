#!/usr/bin/env bash
data=/mnt/ephemeral/toil
mkdir -p $data/tmp
mkdir -p $data/output
export TMPDIR=$data/tmp
python -m toil_scripts.gatk_processing.gatk_preprocessing \
$data/store \
--config "config.txt" \
--genome "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/hs37d5.fa" \
--phase "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/1000G_phase1.indels.b37.vcf" \
--mills "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/Mills_and_1000G_gold_standard.indels.b37.vcf" \
--dbsnp "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/germline/dbsnp_138.b37.vcf" \
--output_dir $data/output \
--workDir $TMPDIR \
#--restart
