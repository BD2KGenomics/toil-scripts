## University of California, Santa Cruz Genomics Institute
### Guide: Running UNC Best Practices RNA-seq Pipeline using Toil

This guide attempts to walk the user through running this pipeline from start to finish. If there are any questions
please contact John Vivian (jtvivian@gmail.com). If you find any errors or corrections please feel free to make a 
pull request.  Feedback of any kind is appreciated.

## Overview
This pipeline accepts a sample.tar file (either by URL or locally) that contains RNA-seq fq.gz data files.  It assumes
that every fastq file has either **R1** or **R2** in the file name.  If there are several fastq files in the tarfile, 
(typically in the context of sequencing lanes), they will be concatenated together into one **R1.fq.gz** and 
one **R2.fq.gz** file.

This pipeline produces a tarball (tar.gz) file for a given sample that contains:

    Gene & isoform: counts, TPM, FPKM, and raw counts.
    Exon quantification
    Directory of QC metrics produced by RseqQC
 
The output tarball is *stamped* with the UUID for the sample (e.g. UUID.tar.gz). If a config file is being used, the
UUID is specified on each line per sample, otherwise the UUID is derived from the input file (UUID.tar). 

## Dependencies
This pipeline has been tested on Ubuntu 14.04, but should also run on other unix based systems.  `apt-get` and `pip`
often require `sudo` privilege, so if the below commands fail, try prepending `sudo`.  If you do not have sudo 
privileges you will need to build these tools from source, or bug a sysadmin (they don't mind). 

#### General Dependencies
    1. Python 2.7
    2. Curl         apt-get install curl
    3. Docker       http://docs.docker.com/engine/installation/
    4. Samtools     apt-get install samtools (attempting to remove)
    5. Unzip        apt-get install unzip

#### Python Dependencies
    1. Toil         pip install toil
    2. Boto         pip install boto (optional, only needed if uploading results to S3)


## Getting Started
#### Running a single sample locally
From the BD2KGenomics toil-scripts Github repository, download the following files to the same directory.

    1. toil-scripts/rnaseq-unc/rnaseq_unc_pipeline.py
    2. toil-scripts/rnaseq-unc/launch_unc.sh
    
The bash script `launch_unc.sh` contains all of the parameters required to run this pipeline, although you 
will likely want to modify a couple lines as it assumes everything will be staged from your home directory.

| Parameter                 | Function                                                                                                                              |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| 1st argument (unlabelled) | Path to where the jobStore will exist. The jobStore hosts intermediate files during runtime.                                          |
| `--config` OR `--input`   | Path to the config csv file OR the sample.tar.  UUID for the sample is based off the filename before the .tar extension               |
| `--single_end_reads`      | OPTIONAL: Set this flag if the data is non-paired (single end reads)                                                                  |
| `--retryCount`            | OPTIONAL: Number of times a failed job will retried. Useful for non-systemic failures (HTTP requests, etc)                            |
| `--ssec`                  | OPTIONAL: Path to a master key if input files are encrypted in S3                                                                     |
| `--output_dir`            | OPTIONAL: Directory where final output of pipeline will be placed                                                                     |
| `--s3_dir`                | OPTIONAL: S3 "Directory" (bucket + directories)                                                                                       |
| `--workDir`               | OPTIONAL: Location where tmp files will be placed during pipeline run.,If not used, defaults to TMPDIR environment variable.          |
| `--sudo`                  | OPTIONAL: Prepends "sudo" to all docker commands. Necessary if user is not a member of a docker group or does not have root privilege |
| `--restart`               | OPTIONAL: Restarts pipeline after failure, requires presence of an existing jobStore.                                                 |

For users *outside* of the BD2K group at UC Santa Cruz, here is an example of a modified launch script that assumes the 
RNA-seq sample is local, the user has sudo privilege, and wants the output of the rna-seq pipeline locally.

``` shell
#!/usr/bin/env bash
# Ensure directory where jobStore and temp files will go exists.
mkdir -p ${HOME}/toil_mnt
# Execution of pipeline
python rnaseq_unc_pipeline.py \
${HOME}/toil_mnt/jstore \
--retryCount 2 \
--input /path/to/sample.tar \
--output_dir ${HOME}/rnaseq_output \
--workDir ${HOME}/toil_mnt \
--sudo 
```

The first argument (location of the jobStore) and the directory set in `--workDir`, need *plenty* of space to store 
intermediate files during pipeline execution.  Change those parameters to point to the appropriate scratch space or
wherever there exists sufficient storage. The servers I have tested on have 700GB of disk space, which is plenty,
but ultimately this is contingent upon sample size.

#### Running a sample on a batch system (gridEngine, Parasol, etc).
To run your pipeline using the gridEngine batch system, simply add the argument `--batchSystem=gridEngine` to the launch
script.  We currently support Grid Engine, Parasol, and Mesos. 
 
#### Running batches of samples
Instead of using the `--input` argument, use mutually exclusive argument `--config`. This accepts a path to a 
CSV file that **must** follow the format of one sample per line:  UUID,url_to_sample.tar.  Support does not
currently exist for running batches of local samples, as we are making an effort to move most of the data we 
process to cloud storage.

## Advanced: Running the Pipeline on a Distributed Cloud Cluster (using Mesos)
From the BD2KGenomics toil-scripts Github repository, download the following files which will run on the head node.

    1. toil-scripts/rnaseq-unc/rnaseq_unc_pipeline.py
    2. toil-scripts/rnaseq-unc/launch_unc_mesos.sh
    
It is outside the scope of this guide to explain how to setup a distributed cloud cluster.  I recommend taking a 
look at the BD2KGenomics tool: [CGCloud](https://github.com/BD2KGenomics/cgcloud), which can setup a distributed 
cloud cluster using the Mesos batch system in AWS.  Please do not direct questions related to CGCloud or 
setting up a distributed cluster to the author of this pipeline. 

A launch script has been prepared that will run on the head node of the Mesos cluster, scheduling jobs to the worker
nodes that exist within the cluster.

``` shell
#!/usr/bin/env bash
python rnaseq_unc_pipeline.py \
aws:us-west-2:unc_pipeline-run-1 \
--retryCount 3 \
--config toil_rnaseq_config.csv \
--ssec "/home/mesosbox/shared/master.key" \
--output_dir "/home/mesosbox/rnaseq_output" \
--s3_dir "cgl-driver-projects/test/rna-test/" \
--workDir=/var/lib/toil \
--batchSystem="mesos" \
--masterIP=mesos-master:5050 \
--sseKey=/home/mesosbox/shared/master.key \
```

Explanation of additional parameters

| Parameter     | Function                                                                                                                 |
|---------------|--------------------------------------------------------------------------------------------------------------------------|
| 1st argument  | This now points to an AWS jobStore                                                                                       |
| `--batchSystem` | Path to the config csv file OR the sample.tar.  UUID for the sample is based off the filename before the .tar extension. |
| `--masterIP`    | A boilerplate argument that indicates what port to use                                                                   |
| `--sseKey`      | OPTIONAL: Encrypts intermediate files when using cloud jobStore.   

**NOTE:** Every worker node must have all of the required dependencies, as well as the inputs for the arguments 
`--sseKey`, `--ssec`, and `--config` must be on *every* worker node, otherwise as the pipeline runs, jobs wil fail 
as those files will not be found.  If `--s3_dir` is used, a ~/.boto config file with credentials must also be on every
worker.

## Troubleshooting
This section is a work in progress that will be amended as issues arise.
#### Encryption Key Errors
If encryption keys of any type will be used, they must be **exactly** 32-byte keys.

## Additional Information
This pipeline aligns data to the HG19 human reference genome before performing RNA-seq analysis. If you want to align
to a different genome then you will need to look at the help menu `python rnaseq_unc_pipeline.py -h` to see all
of the different inputs that you will need to change.  Support is not offered if you would like to use a different
reference genome. Our group plans on developing a new RNA-seq pipeline for HG38 that uses STAR+RSEM and Kallisto.