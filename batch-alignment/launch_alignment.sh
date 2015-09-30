#!/usr/bin/env bash
export TMPDIR=/mnt
mkdir -p "/mnt/final_output"
python batch_align.py \
/mnt/jstore \
--retryCount 3 \
-c config.txt \
--ref https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa \
--amb https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.amb \
--ann https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.ann \
--bwt https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.bwt \
--pac https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.pac \
--sa https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.sa \
--fai https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.fai \
--ssec "/home/ubuntu/master.key" \
-o "/mnt/final_output" \
--s3_dir "cgl-driver-projects-encrypted/wcdt/exome_bams/"
