#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline \
/data/rnseq-jstore \
--config /data/rnaseq-config.csv \
--retryCount 1 \
--ssec /data/master.key \
--output_dir /data \
--workDir /data
