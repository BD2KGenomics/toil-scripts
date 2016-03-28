#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Precautionary step: Create location where jobStore and tmp files will exist
mkdir -p ${HOME}/toil_mnt
# Execution of pipeline
python bwa_alignment.py \
${HOME}/toil_mnt/jobStore \
--retryCount 3 \
--config bwa_config.csv \
--lb KapaHyper \
--ref https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg19.fa \
--amb https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg19.fa.amb \
--ann https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg19.fa.ann \
--bwt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg19.fa.bwt \
--pac https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg19.fa.pac \
--sa https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg19.fa.sa \
--fai https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg19.fa.fai \
--workDir ${HOME}/toil_mnt \
--ssec ${HOME}/master.key \
--output_dir ${HOME} \
--s3_dir cgl-driver-projects/test/alignment \
--sudo \
#--restart
