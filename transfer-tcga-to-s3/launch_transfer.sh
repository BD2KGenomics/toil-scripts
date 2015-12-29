#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
# Execution of pipeline
python transfer_tcga_to_s3.py \
/home/ubuntu/jstore \
--genetorrent /home/ubuntu/config.txt \
--genetorrent_key /home/ubuntu/haussl_new.key \
--ssec /home/ubuntu/master.key \
--s3 tcga-data \
--sudo
#--restart