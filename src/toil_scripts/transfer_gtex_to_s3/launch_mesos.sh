#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
# Execution of pipeline
python transfer_gtex_to_s3.py \
aws:us-west-2:gtex-transfer-run-1 \
--sra /home/mesosbox/shared/config.txt \
--dbgap_key /home/mesosbox/shared/prj_8908.ngc \
--ssec /home/mesosbox/shared/master.key \
--s3 gtex-data \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
--workDir /var/lib/toil \
--retryCount 2 \
#--restart