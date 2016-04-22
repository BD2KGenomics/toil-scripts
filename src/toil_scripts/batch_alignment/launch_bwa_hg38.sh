#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.batch_alignment.bwa_alignment \
/data/alignment-jobStore \
--retryCount 2 \
--config bwa-config.csv \
--library KapaHyper \
--ref https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38.fa \
--amb https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38.fa.amb \
--ann https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38.fa.ann \
--bwt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38.fa.bwt \
--pac https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38.fa.pac \
--sa https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38.fa.sa \
--fai https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38.fa.fai \
--workDir /data \
--ssec /data/master.key \
--output-dir /data
