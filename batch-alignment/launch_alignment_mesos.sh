#!/usr/bin/env bash
export TMPDIR=/var/lib/toil/ # unnecessary
python toil_alignment.py \
aws:us-west-2:wcdt-toil-alignment-run-1 \
--retryCount 1 \
-c /home/mesosbox/shared/config.txt \
--ref https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa \
--amb https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.amb \
--ann https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.ann \
--bwt https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.bwt \
--pac https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.pac \
--sa https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.sa \
--fai https://s3-us-west-2.amazonaws.com/cgl-alignment-inputs/genome.fa.fai \
--ssec /home/mesosbox/shared/master.key \
-o /var/lib/toil/final_output \
--s3_dir "cgl-driver-projects-encrypted/wcdt/exome_bams/" \
--sseKey=/home/mesosbox/shared/master.key \ 
--batchSystem="mesos" \ 
--masterIP=mesos-master:5050 \
--workDir=/var/lib/toil
