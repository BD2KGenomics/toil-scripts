#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Precautionary step: Create location where jobStore and tmp files will exist
mkdir -p ${HOME}/toil_mnt
# Execution of pipeline
python rnaseq_unc_pipeline.py \
${HOME}/toil_mnt/jstore \
--config unc_config.csv \
--retryCount 2 \
--ssec ${HOME}/master.key \
--output_dir ${HOME}/rnaseq_output \
--s3_dir cgl-driver-projects/test/rna-test/ \
--workDir ${HOME}/toil_mnt \
--sudo \
#--restart
