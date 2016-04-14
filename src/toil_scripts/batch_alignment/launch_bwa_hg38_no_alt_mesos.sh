#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.batch_alignment.bwa_alignment \
aws:us-west-2:alignment-run-1 \
--retryCount 2 \
--config /home/mesosbox/shared/bwa-config.csv \
--lb KapaHyper \
--ref https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa \
--amb https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.amb \
--ann https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.ann \
--bwt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.bwt \
--pac https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.pac \
--sa https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.sa \
--fai https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.fai \
--workDir /var/lib/toil \
--ssec /home/mesosbox/shared/master.key \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
--s3_dir cgl-driver-projects/test/alignment \
--use_bwakit
