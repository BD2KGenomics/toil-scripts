#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Precautionary step: Create location where jobStore and tmp files will exist
# Execution of pipeline
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline \
/data/rnseq-jstore \
--config /data/config.txt \
--retryCount 1 \
--ssec /data/master.key \
--s3_dir cgl-driver-projects/test/releases \
--workDir /data
#--restart
