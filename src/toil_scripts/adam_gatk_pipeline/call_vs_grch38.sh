#!/usr/bin/env bash

set -x -v

# Frank Austin Nothaft, fnothaft@berkeley.edu
#
# Pipeline for alt-aware alignment against the GRCh38 build used in the 1000G vs. b38 recompute,
# followed by preprocessing using ADAM, and variant calling using the HaplotypeCaller.
#
# Precautionary step: Create location where jobStore and tmp files will exist and set TOIL_HOME.

# Execution of pipeline
python -m toil_scripts.adam_gatk_pipeline.align_and_call \
    aws:us-west-2:fnothaft-toil-jobstore \
    --retryCount 1 \
    --uuid SRR062640 \
    --s3_bucket fnothaft-fc-test-west-2 \
    --bucket_region us-west-2 \
    --aws_access_key ${FC_AWS_ACCESS_KEY_ID} \
    --aws_secret_key ${FC_AWS_SECRET_ACCESS_KEY} \
    --ref https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --amb https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb \
    --ann https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann \
    --bwt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt \
    --pac https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.pac \
    --sa https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.sa \
    --fai https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
    --alt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt \
    --use_bwakit \
    --num_nodes 1 \
    --driver_memory 50g \
    --executor_memory 50g \
    --phase https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz \
    --mills https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz \
    --dbsnp https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz \
    --omni https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz \
    --hapmap https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz \
    --batchSystem=mesos \
    --mesosMaster $(hostname -i):5050 \
    --workDir /var/lib/toil \
    --file_size 1G \
    --logInfo
