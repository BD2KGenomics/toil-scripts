#!/usr/bin/env bash
export TMPDIR=/home/ubuntu/fake_mnt
mkdir -p "/home/ubuntu/final_output"
python rna-seq_pipeline_single_sample.py \
--retryCount 3 \
-r1 "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/LNCapCtrl_S16_L005_R1_001.fastq.gz" \
-r2 "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/LNCapCtrl_S16_L005_R2_001.fastq.gz" \
-u "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/unc_hg19.bed" \
-f "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/hg19_M_rCRS_ref.transcripts.fa" \
-c "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/composite_exons.bed" \
-n "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/normalizeBedToolsExonQuant.pl" \
-o "/home/ubuntu/final_output" \
-i "test_UUID"
