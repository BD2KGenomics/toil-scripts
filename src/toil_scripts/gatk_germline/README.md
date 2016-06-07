## University of California, Santa Cruz Genomics Institute
### Guide: Running GATK Best Practices Variant Pipeline using Toil

## Overview
The GATK germline pipeline includes GATK preprocessing, variant calling, and filtering. HaplotypeCaller is run per sample and then individual VCFs are joined to increase statistical power for variant quality score recalibration. If there are fewer than thirty exome samples, it is recommended that hard filters are used instead. The pipeline also supports clinical variant annotation using Oncotator.

#### General Dependencies
    1. Python 2.7
    2. Curl         apt-get install curl
    3. Docker       http://docs.docker.com/engine/installation/

#### Python Dependencies
    1. pip          apt-get install python-pip
    2. virtualenv   pip install virtualenv
    3. Toil         pip install toil

#### Installation
Toil-scripts is now pip installable! `pip install toil-scripts` for a toil-stable version 
or `pip install --pre toil-scripts` for cutting edge development version.

Type: `toil-germline` to get basic help menu and instructions

To decrease the chance of versioning conflicts, install toil-scripts into a virtualenv: 

- `virtualenv ~/toil-scripts` 
- `source ~/toil-scripts/bin/activate`
- `pip install toil`
- `pip install toil-scripts`

If Toil is already installed globally (true for CGCloud users), or there are global dependencies (like Mesos),
use virtualenv's `--system-site-packages` flag.

## General Usage
    1. Type `toil-germline generate` to create an editable manifest and config in the current working directory.
    2. Parameterize the pipeline by editing the config.
    3. Fill in the manifest with information pertaining to your samples.
    4. Type `toil-germline run [jobStore]` to execute the pipeline.

## Example Commands
Run sample(s) locally using the manifest
    1. `toil-germline generate`
    2. Fill in config and manifest
    3. `toil-germline run --config config-toil-germline.yaml \
        --manifest manifest-toil-germline.tsv ./example-jobstore`

Toil options can be appended to `toil-germline run`, for example:
`toil-germline run ./example-jobstore --retryCount=1 --workDir=/data`

For a complete list of Toil options, just type `toil-germline run -h`

Run a single sample locally
    1. `toil-germline generate-config`
    2. Fill in config
    3. `toil-germline run ./example-jobstore --workDir /data --samples \
        UUID https://sample-depot.com/sample.bam`
        
## GATK Variant Annotations:
Currently, this pipeline uses the following variant annotations:

* QualByDepth
* FisherStrand
* StrandOddsRatio
* ReadPosRankSumTest
* MappingQualityRankSumTest
* RMSMappingQuality
* InbreedingCoeff
    
## Example Config
```
genome-fasta: s3://cgl-pipeline-inputs/hg19/ucsc.hg19.fasta
genome-fai: s3://cgl-pipeline-inputs/hg19/ucsc.hg19.fasta.fai
genome-dict: s3://cgl-pipeline-inputs/hg19/ucsc.hg19.dict
phase: s3://cgl-pipeline-inputs/hg19/1000G_phase1.indels.hg19.sites.vcf
mills: s3://cgl-pipeline-inputs/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
dbsnp: s3://cgl-pipeline-inputs/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf
hapmap: s3://cgl-pipeline-inputs/hg19/hapmap_3.3.hg19.sites.vcf
omni: s3://cgl-pipeline-inputs/hg19/1000G_omni2.5.hg19.sites.vcf
run-bwa: True
trim: True
amb: s3://cgl-pipeline-inputs/alignment/hg19.fa.amb
ann: s3://cgl-pipeline-inputs/alignment/hg19.fa.ann
bwt: s3://cgl-pipeline-inputs/alignment/hg19.fa.bwt
pac: s3://cgl-pipeline-inputs/alignment/hg19.fa.pac
sa: s3://cgl-pipeline-inputs/alignment/hg19.fa.sa
run-vqsr: True
joint: True
file-size: 10G
xmx: 30G
suffix: .toil
output_dir: /data/my-toil-run
ssec:
unsafe_mode: False
```
