#!/usr/bin/env bash
export TMPDIR=/home/ubuntu/fake_mnt
mkdir -p "/home/ubuntu/final_output"
python rna-seq_pipeline_multi_sample.py \
--retryCount 3 \
-c config.txt \
--unc "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/unc_hg19.bed" \
--fasta "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/hg19_M_rCRS_ref.transcripts.fa" \
--composite_exons "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/composite_exons.bed" \
--normalize "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/normalizeBedToolsExonQuant.pl" \
--rsem_ref "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/rsem_ref.zip" \
--chromosomes "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/chromosomes.zip" \
--ebwt "https://s3-us-west-2.amazonaws.com/rna-seq-pipeline/ebwt.zip" \
-o "/home/ubuntu/final_output" \
