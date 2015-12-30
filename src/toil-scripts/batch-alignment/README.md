## University of California, Santa Cruz Genomics Institute
### Guide: Running GATK-compatible Alignment 

This guide attempts to walk the user through running this pipeline from start to finish. If there are any questions
please contact John Vivian (jtvivian@gmail.com). If you find any errors or corrections please feel free to make a 
pull request.  Feedback of any kind is appreciated.

## Overview
This pipeline accepts two fastq files (by URL) to be aligned into a BAMFILE, which is the final output of the pipeline.
A launch script is provided for 3 different references (b37, hg19, hg38).

## Dependencies
This pipeline has been tested on Ubuntu 14.04, but should also run on other unix based systems.  `apt-get` and `pip`
often require `sudo` privilege, so if the below commands fail, try prepending `sudo`.  If you do not have sudo 
privileges you will need to build these tools from source, or bug a sysadmin (they don't mind). 

#### General Dependencies
    1. Python 2.7
    2. Curl     apt-get install curl
    3. Docker   http://docs.docker.com/engine/installation/
    
#### Python Dependencies
    1. Toil     pip install toil
    2. S3AM     pip install --pre s3am (optional, for upload of BAMFILE to S3)
 
## Getting Started
#### Aligning one sample locally
From the BD2KGenomics toil-scripts Github repository, download the following files to the same directory.

    1. toil-scripts/batch-alignment/bwa_alignment.py
    2. toil-scripts/batch-alignment/launch_bwa_hg19.sh
    3. toil-scripts/batch-alignment/bwa_config.csv
    
The bash script `launch_bwa_hg19.sh` contains all of the parameters required to run this pipeline, although you 
will likely want to modify a couple lines as it assumes everything will be staged from your home directory.

| Parameter                 | Function                                                                                                                              |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| 1st argument (unlabelled) | Path to where the jobStore will exist.  The jobStore hosts intermediate files during runtime                                          |
| `--config`                | Path to the config file. With the format:  UUID,url_to_fastq1,url_to_fastq2                                                           |
| `--retryCount`            | OPTIONAL: Number of times a failed job will retried. Useful for non-systemic failures (HTTP requests, etc)                            |
| `--ssec`                  | OPTIONAL: Path to a master key if input files are encrypted in S3                                                                     |
| `--output_dir`            | OPTIONAL: Directory where final output of pipeline will be placed                                                                     |
| `--s3_dir`                | OPTIONAL: S3 "Directory" (bucket + directories)                                                                                       |
| `--workDir`               | OPTIONAL: Location where tmp files will be placed during pipeline run.,If not used, defaults to TMPDIR environment variable.          |
| `--sudo`                  | OPTIONAL: Prepends "sudo" to all docker commands. Necessary if user is not a member of a docker group or does not have root privilege |
| `--restart`               | OPTIONAL: Restarts pipeline after failure, requires presence of an existing jobStore.                                                 |


The first argument (location of the jobStore) and the directory set in `--workDir`, need *plenty* of space to store 
intermediate files during pipeline execution.  Change those parameters to point to the appropriate scratch space or
wherever there exists sufficient storage. The servers I have tested on have 700GB of disk space, which is plenty,
but ultimately this is contingent upon sample size.

#### Running a sample on a batch system (gridEngine, Parasol, etc).
To run your pipeline using the gridEngine batch system, simply add the argument `--batchSystem=gridEngine` to the launch
script.  We currently support Grid Engine, Parasol, and Mesos. 
 
#### Running batches of samples
`--config` accepts a path to a CSV file that **must** follow the format of one sample 
per line: UUID,url_to_fastq1,url_to_fastq2

# Advanced: Running the Pipeline on a Distributed Cloud Cluster (using Mesos)
From the BD2KGenomics toil-scripts Github repository, download the following files which will run on the head node.

    1. toil-scripts/batch-alignment/bwa_alignment.py
    2. toil-scripts/batch-alignment/launch_bwa_hg19_mesos.sh
    3. toil-scripts/batch-alignment/bwa_config.csv
    
It is outside the scope of this guide to explain how to setup a distributed cloud cluster.  I recommend taking a 
look at the BD2KGenomics tool: [CGCloud](https://github.com/BD2KGenomics/cgcloud), which can setup a distributed 
cloud cluster using the Mesos batch system in AWS.  Please do not direct questions related to CGCloud or 
setting up a distributed cluster to the author of this pipeline. 

A launch script (`launch_bwa_hg19_mesos.sh`) has been prepared that will run on the head node of the Mesos cluster, scheduling jobs to the worker
nodes that exist within the cluster.

Explanation of additional parameters

| Parameter     | Function                                                                                                                 |
|---------------|--------------------------------------------------------------------------------------------------------------------------|
| 1st argument  | This now points to an AWS jobStore                                                                                       |
| `--batchSystem` | Path to the config csv file OR the sample.tar.  UUID for the sample is based off the filename before the .tar extension. |
| `--masterIP`    | A boilerplate argument that indicates what port to use                                                                   |
| `--sseKey`      | OPTIONAL: Encrypts intermediate files when using cloud jobStore.
 
## Additional Information
Launch scripts are provided for alignment to hg19, hg38, and hg38 with no alternative loci. 