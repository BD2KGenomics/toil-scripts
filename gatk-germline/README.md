## University of California, Santa Cruz Genomics Institute
### Guide: Running GATK Best Practices Variant Pipeline using Toil

## Overview
This pipeline takes an exhaustive list of bam files and uses GATK tools to identify germline variants.

#### General Dependencies
    1. Python 2.7
    2. Curl         apt-get install curl
    3. Docker       apt-get install docker.io

#### Python Dependencies
    1. Toil         pip install toil==3.1.0b1.dev56

## Getting Started
#### Running a single sample locally
From the BD2KGenomics toil-scripts Github repository, download the following files to the same directory.

    1. toil-scripts/gatk-germline/germline.py
    2. toil-scripts/gatk-germline/launch_germ.py
    
The bash script `launch_germ.sh` contains all of the default parameters required to run the pipeline.

| Parameter                 | Function                                                                                                                              |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| 1st argument (unlabelled) | Path to where the jobStore will exist. The jobStore hosts intermediate files during runtime.                                          |
| `--config`		    | Path to the config csv file. Each line of config file has UUID, URL for sample bam						    |
| `--reference'  	    | Reference genome in FASTA format                                                                                                     |
| `--phase'		    | URL 1000G_phase1.indels.b37.vcf                                                                                                           |
| `--mills`		    | URL Mills_and_1000G_gold_standard.indels.b37.sites.vcf										    |
| `--dbsnp`		    | URL dbsnp_138.b37.vcf															    |
| `--hapmap`	            | URL hapmap_3.3.b37.vcf														    |	
| `--omni`		    | URL 1000G_omni2.5.b37.vcf														    |					
| `--retryCount`            | OPTIONAL: Number of times a failed job will retried. Useful for non-systemic failures (HTTP requests, etc)                            |
| `--output_dir`            | OPTIONAL: Directory where final output of pipeline will be placed                                                                     |
| `--workDir`               | OPTIONAL: Location where tmp files will be placed during pipeline run. If not used, defaults to TMPDIR environment variable.          |
| `--restart`               | OPTIONAL: Restarts pipeline after failure, requires presence of an existing jobStore.                                                 |
