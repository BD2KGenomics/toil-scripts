#!/usr/bin/env bash
export TMPDIR=/mnt/
mkdir -p "/home/ubuntu/final_output"
python rna-seq_pipeline.py \
--retryCount 3 \
-c config.txt \
--unc "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/unc_hg19.bed" \
--fasta "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/hg19_M_rCRS_ref.transcripts.fa" \
--composite_exons "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/composite_exons.bed" \
--normalize "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/normalizeBedToolsExonQuant.pl" \
--rsem_ref "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/rsem_ref.zip" \
--chromosomes "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/chromosomes.zip" \
--ebwt "https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/ebwt.zip" \
--ssec "/mnt/rna-seq/shared/foo.key" \
-o "/home/ubuntu/final_output" \
--s3_dir "cgl-driver-projects/wcdt/rna-seq-outputs/"
