#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Precautionary step: Create location where jobStore and tmp files will exist
mkdir -p ${HOME}/toil_mnt
# Execution of pipeline
python rnaseq_cgl_pipeline.py \
${HOME}/toil_mnt/jstore \
--config rnaseq_cgl_config.csv \
--retryCount 2 \
--ssec /home/mesosbox/shared/master.key \
--s3_dir cgl-driver-projects/test/rna_cgl-test/ \
--sseKey=/home/mesosbox/shared/master.key \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
--workDir=/var/lib/toil \
#--restart
