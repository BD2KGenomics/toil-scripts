#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
# Execution of pipeline
python transfer_tcga_to_s3.py \
/var/lib/toil/jstore \
--genetorrent /home/mesosbox/shared/config.txt \
--genetorrent_key /home/mesosbox/shared/haussl_new.key \
--ssec /home/mesosbox/shared/master.key \
--s3 tcga-data-cgl-recompute \
--workDir /var/lib/toil \
--retryCount 2 \
#--restart