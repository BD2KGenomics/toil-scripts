#!/usr/bin/env bash
# John Vivian
#
# This script provides all configuration necessary to run the Toil pipeline on a Mesos cluster.
# It assumes there is a local file: bwa_config.csv.  One sample per line: uuid,url_R1,url_R2.
#
#   --ssec          the program assumes input files are encrypted in S3 when retrieving them.
#   --output_dir    the final BAM will be placed in the directory specified
#   --s3_dir        the final BAM will be uploaded to S3 using S3AM (pip install --pre s3am, need ~/.boto)
#   --sudo          'sudo' will be prepended to the Docker subprocess call
#
# Modify TMPDIR parameter to change location of tmp files.
# Modify first argument to change location of the local fileStore
# Uncomment the final line to resume your Toil job in the event of job failure.
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
--masterIP=mesos-master:5050 \
--workDir=/var/lib/toil \
#--restart