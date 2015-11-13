#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
python rnaseq_unc_pipeline.py \
aws:us-west-2:unc_pipeline-run-1 \
--retryCount 3 \
--config toil_rnaseq_config.csv \
--ssec "/home/mesosbox/shared/master.key" \
--output_dir "/home/mesosbox/rnaseq_output" \
--s3_dir "cgl-driver-projects/test/rna-test/" \
--sudo \
--sseKey=/home/mesosbox/shared/master.key \
--batchSystem="mesos" \
--masterIP=mesos-master:5050 \
--workDir=/var/lib/toil \
#--restart
