#!/usr/bin/env bash
python rna-seq_pipeline.py \
aws:us-west-2:rna-seq-toil-1 \
--retryCount 3 \
-c "/home/mesosbox/shared/config.txt" \
--unc "https://s3-us-west-2.amazonaws.com/cgl-rna-seq-inputs/unc_hg19.bed" \
--fasta "https://s3-us-west-2.amazonaws.com/cgl-rna-seq-inputs/hg19_M_rCRS_ref.transcripts.fa" \
--composite_exons "https://s3-us-west-2.amazonaws.com/cgl-rna-seq-inputs/composite_exons.bed" \
--normalize "https://s3-us-west-2.amazonaws.com/cgl-rna-seq-inputs/normalizeBedToolsExonQuant.pl" \
--rsem_ref "https://s3-us-west-2.amazonaws.com/cgl-rna-seq-inputs/rsem_ref.zip" \
--chromosomes "https://s3-us-west-2.amazonaws.com/cgl-rna-seq-inputs/chromosomes.zip" \
--ebwt "https://s3-us-west-2.amazonaws.com/cgl-rna-seq-inputs/ebwt.zip" \
--ssec "/home/mesosbox/shared/master.key" \
--sseKey="/home/mesosbox/shared/master.key" \
--batchSystem="mesos" \
--masterIP=mesos-master:5050 \
--workDir=/var/lib/toil \
--s3_dir "cgl-driver-projects/wcdt/rna-seq-outputs/"