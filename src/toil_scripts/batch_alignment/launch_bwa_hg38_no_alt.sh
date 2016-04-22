#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.batch_alignment.bwa_alignment \
/data/alignment-jstore \
--retryCount 2 \
--config /data/bwa-config.cvs \
--library KapaHyper \
--ref https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa \
--amb https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.amb \
--ann https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.ann \
--bwt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.bwt \
--pac https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.pac \
--sa https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.sa \
--fai https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/alignment/hg38_no_alt.fa.fai \
--workDir /data \
--ssec /data/master.key \
--output-dir /data
