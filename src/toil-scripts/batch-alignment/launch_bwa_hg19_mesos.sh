#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
python bwa_alignment.py \
aws:us-west-2:alignment-wcdt-run-1 \
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
--ssec /home/mesosbox/shared/master.key \
--output_dir /home/mesosbox/shared \
--s3_dir cgl-driver-projects/test/alignment \
--sseKey=/home/mesosbox/shared/master.key \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
--workDir=/var/lib/toil \
#--restart