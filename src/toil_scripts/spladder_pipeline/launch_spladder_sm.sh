#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Execution of pipeline
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.spladder_pipeline.spladder_pipeline \
/data/jstore \
--config /data/config.txt \
--retryCount 1 \
--ssec /data/master.key \
--output-dir /data \
--sseKey=/data/master.key \
--workDir=/data \
#--restart \
