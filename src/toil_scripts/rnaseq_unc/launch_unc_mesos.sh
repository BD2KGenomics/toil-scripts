#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.rnaseq_unc.rnaseq_unc_pipeline \
aws:us-west-2:unc-pipeline-run-1 \
--retryCount 1 \
--config toil_rnaseq_config.csv \
--ssec "/home/mesosbox/shared/master.key" \
--output_dir "/home/mesosbox/rnaseq_output" \
--s3_dir "cgl-driver-projects/test/rna-test/" \
--sudo \
--sseKey=/home/mesosbox/shared/master.key \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
--workDir=/var/lib/toil \
#--restart
