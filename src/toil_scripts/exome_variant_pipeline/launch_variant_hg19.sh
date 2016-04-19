#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python toil_scripts.exome_variant_pipeline.exome_variant_pipeline \
/data/exome-jobStore\
--retryCount 2 \
--config "exome-config.csv" \
--reference 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/hg19.fa' \
--phase 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/1000G_phase1.indels.hg19.sites.vcf' \
--mills 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf' \
--dbsnp 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/dbsnp_138.hg19.vcf' \
--cosmic 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/cosmic.hg19.vcf' \
--output_dir /data \
--ssec /data/master.key \
--workDir /data
