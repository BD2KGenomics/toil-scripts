#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Execution of pipeline
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.spladder_pipeline.spladder_pipeline \
aws:us-west-2:spladder-run-1 \
--config /home/mesosbox/shared/config.txt \
--retryCount 1 \
--ssec /home/mesosbox/shared/master.key \
--output-s3-dir s3://bucket/dir \
--sseKey=/home/mesosbox/shared/master.key \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
--workDir=/var/lib/toil \
# --restart
