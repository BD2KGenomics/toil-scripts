#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Execution of pipeline
python rnaseq_cgl_pipeline.py \
/mnt/jstore \
--genetorrent ${HOME}/config.txt \
--genetorrent_key ${HOME}/haussl_new.key \
--retryCount 2 \
--ssec ${HOME}/master.key \
--s3_dir tcga-test-output/edge-tests \
--workDir /mnt/ \
--sudo \
#--restart
