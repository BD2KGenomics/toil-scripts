#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
python bwa_alignment.py \
/home/mesosbox/shared/toil_mnt/jobStore \
--retryCount 3 \
--config bwa_config.csv \
--lb KapaHyper \
--ref https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa \
--amb https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.amb \
--ann https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.ann \
--bwt https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.bwt \
--pac https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.pac \
--sa https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.sa \
--fai https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.fai \
--ssec /home/mesosbox/shared/master.key \
--output_dir /home/mesosbox/shared \
--s3_dir cgl-driver-projects/test/alignment \
--sudo \
--sseKey=/home/mesosbox/shared/master.key \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
--workDir=/var/lib/toil \
#--restart