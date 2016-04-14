#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.rnaseq_unc.rnaseq_unc_pipeline \
/data/unc-jobStore \
--config unc-config.csv \
--retryCount 2 \
--ssec /data/master.key \
--output_dir /data/rnaseq_output \
--workDir /data
