#!/usr/bin/env bash
# Frank Austin Nothaft, fnothaft@berkeley.edu
#
# Pipeline for alt-aware alignment against the GRCh38 build used in the 1000G vs. b38 recompute.
# Precautionary step: Create location where jobStore and tmp files will exist and set TOIL_HOME.
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.batch_alignment.bwa_alignment \
/data/alignment-jobStore \
--retryCount 2 \
--config bwa-config.csv \
--library LIB \
--ref https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--amb https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb \
--ann https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann \
--bwt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt \
--pac https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.pac \
--sa https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.sa \
--fai https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
--alt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt \
--workDir /data \
--output-dir /data
