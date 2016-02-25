#!/usr/bin/env bash

set -x -v

# Frank Austin Nothaft, fnothaft@berkeley.edu
#
# Pipeline for alt-aware alignment against the GRCh38 build used in the 1000G vs. b38 recompute,
# preprocessing using ADAM followed by variant calling using GATK HaplotypeCaller, and preprocessing
# using GATK followed by variant calling using GATK HaplotypeCaller.

python -m toil_scripts.adam_gatk_pipeline.align_and_call \
    ${JOBSTORE} \
    --retryCount 1 \
    --uuid_manifest my_manifest_file \
    --s3_bucket ${BUCKET} \
    --bucket_region ${REGION} \
    --ref https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/GRCh38_full_analysis_set_plus_decoy_hla.reordered.fa \
    --amb https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/GRCh38_full_analysis_set_plus_decoy_hla.reordered.fa.amb \
    --ann https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/GRCh38_full_analysis_set_plus_decoy_hla.reordered.fa.ann \
    --bwt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/GRCh38_full_analysis_set_plus_decoy_hla.reordered.fa.bwt \
    --pac https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/GRCh38_full_analysis_set_plus_decoy_hla.reordered.fa.pac \
    --sa https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/GRCh38_full_analysis_set_plus_decoy_hla.reordered.fa.sa \
    --fai https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/GRCh38_full_analysis_set_plus_decoy_hla.reordered.fa.fai \
    --alt https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/GRCh38_full_analysis_set_plus_decoy_hla.reordered.fa.alt \
    --use_bwakit \
    --num_nodes 2 \
    --driver_memory 200 \
    --executor_memory 200 \
    --phase https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/1000G_phase1.snps.high_confidence.grch38.reordered.vcf \
    --mills https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/Mills_and_1000G_gold_standard.indels.grch38.reordered.vcf \
    --dbsnp https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf \
    --omni https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/1000G_omni2.5.grch38.reordered.vcf \
    --hapmap https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_grch38_reordered/hapmap_3.3.grch38.reordered.vcf \
    --batchSystem=mesos \
    --mesosMaster $(hostname -i):5050 \
    --workDir /var/lib/toil \
    --file_size 1G \
    --logInfo \
    --clean never
