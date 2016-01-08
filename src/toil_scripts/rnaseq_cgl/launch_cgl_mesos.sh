#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# Execution of pipeline
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline \
aws:us-west-2:cgl-pipeline-run-1 \
--config /home/mesosbox/shared/rnaseq_cgl_config.csv \
--retryCount 1 \
--ssec /home/mesosbox/shared/master.key \
--s3_dir cgl-driver-projects/test/spot-test/ \
--sseKey=/home/mesosbox/shared/master.key \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
--workDir=/var/lib/toil \
#--restart
