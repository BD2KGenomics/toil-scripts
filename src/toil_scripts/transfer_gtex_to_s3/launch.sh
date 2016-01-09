#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
# Execution of pipeline
python transfer_gtex_to_s3.py \
/var/lib/toil \
--sra /home/mesosbox/shared/config. \
--dbgap_key /home/mesosbox/shared/prj_8908.ngc \
--ssec /home/mesosbox/shared/master.key \
--sudo \
--s3 gtex-data \
--workDir /var/lib/toil \
--retryCount 2 \
#--restart