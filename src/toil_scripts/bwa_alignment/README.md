## University of California, Santa Cruz Genomics Institute
### Guide: Running the BWA Pipeline using Toil

This guide attempts to walk the user through running this pipeline from start to finish. 

If you find any errors or corrections please feel free to make a pull request.  Feedback of any kind is appreciated.


## Overview

Fastqs are aligned to create a BAM that is compatible with GATK.

## Dependencies

This pipeline has been tested on Ubuntu 14.04, but should also run on other unix based systems.  `apt-get` and `pip`
often require `sudo` privilege, so if the below commands fail, try prepending `sudo`.  If you do not have `sudo` 
privileges you will need to build these tools from source, or bug a sysadmin about how to get them. 

#### General Dependencies

    1. Python 2.7
    2. Curl         apt-get install curl
    3. Docker       http://docs.docker.com/engine/installation/

#### Python Dependencies

    1. Toil         pip install toil
    2. S3AM         pip install --pre s3am (optional, needed for uploading output to S3)

## Installation

Toil-scripts is now pip installable! `pip install toil-scripts` for a toil-stable version. 

If there is an existing, system-wide installation of Toil, as is the case when using CGCloud, 
the `pip install toil` step should be skipped and virtualenv should be invoked with `--system-site-packages`. 
This way the existing Toil installation will be available inside the virtualenv.

To decrease the chance of versioning conflicts, install toil-scripts into a virtualenv: 

- `virtualenv ~/toil-scripts` 
- `source ~/toil-scripts/bin/activate`
- `pip install toil`, but see next paragraph
- `pip install toil-scripts`

The reason that Toil isn't installed automatically as a dependency of toil-scripts is to
 give users the opportunity to customize the Toil installation by adding optional extras, 
 e.g. via `pip install toil[aws,mesos]`.

## Inputs

The BWA pipeline requires input files in order to run. The only required input, aside from the sample(s), is a 
reference genome.  The pipeline can be sped up by specifying URLs for the reference index files, which are generated 
with `bwa index` and `samtools faidx`.

## General Usage

Type `toil-bwa` to get basic help menu and instructions
 
1. Type `toil-bwa generate` to create an editable manifest and config in the current working directory.
2. Parameterize the pipeline by editing the config.
3. Fill in the manifest with information pertaining to your samples.
4. Type `toil-bwa run [jobStore]` to execute the pipeline.

## Example Commands

Run sample(s) locally using the manifest
1. `toil-bwa generate`
2. Fill in config and manifest
3. `toil-bwa run ./example-jobstore`

Toil options can be appended to `toil-bwa run`, for example:
`toil-bwa run ./example-jobstore --retryCount=1 --workDir=/data`

For a complete list of Toil options, just type `toil-bwa run -h`

Run a variety of samples locally
1. `toil-bwa generate-config`
2. Fill in config
3. `toil-bwa run ./example-jobstore --retryCount=1 --workDir=/data --sample \
    test-uuid file:///full/path/to/read1.fq.gz file:///full/path/to/read2.fq.gz`

## Example Config

   ``` 
    # BWA Alignment Pipeline configuration file
    # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
    # Edit the values in this configuration file and then rerun the pipeline: "toil-bwa run"
    # URLs can take the form: http://, ftp://, file://, s3://, gnos://.
    # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank
    ##############################################################################################################
    # Required: Reference fasta file
    ref: s3://cgl-pipeline-inputs/alignment/hg19.fa
    
    # Required: Output location of sample. Can be full path to a directory or an s3:// URL
    output-dir: /data/
    
    # Required: The library entry to go in the BAM read group.
    library: Illumina
    
    # Required: Platform to put in the read group
    platform: Illumina
    
    # Required: Program Unit for BAM header. Required for use with GATK.
    program_unit: 12345
    
    # Required: Approximate input file size. Provided as a number followed by (base-10) [TGMK]. E.g. 10M, 150G
    file-size: 50G
    
    # Optional: If true, sorts bam
    sort: True
    
    # Optional. If true, trims adapters
    trim: false
    
    # Optional: Reference fasta file (amb) -- if not present will be generated
    amb: s3://cgl-pipeline-inputs/alignment/hg19.fa.amb
    
    # Optional: Reference fasta file (ann) -- If not present will be generated
    ann: s3://cgl-pipeline-inputs/alignment/hg19.fa.ann
    
    # Optional: Reference fasta file (bwt) -- If not present will be generated
    bwt: s3://cgl-pipeline-inputs/alignment/hg19.fa.bwt
    
    # Optional: Reference fasta file (pac) -- If not present will be generated
    pac: s3://cgl-pipeline-inputs/alignment/hg19.fa.pac
    
    # Optional: Reference fasta file (sa) -- If not present will be generated
    sa: s3://cgl-pipeline-inputs/alignment/hg19.fa.sa
    
    # Optional: Reference fasta file (fai) -- If not present will be generated
    fai: s3://cgl-pipeline-inputs/alignment/hg19.fa.fai
    
    # Optional: (string) Path to Key File for SSE-C Encryption
    ssec:
    
    # Optional: Use instead of library, program_unit, and platform.
    rg-line:
    
    # Optional: Alternate file for reference build (alt). Necessary for alt aware alignment
    alt:
    
    # Optional: If true, runs the pipeline in mock mode, generating a fake output bam
    mock-mode:
```

## Distributed Run

To run on a distributed AWS cluster, see [CGCloud](https://github.com/BD2KGenomics/cgcloud) for instance provisioning, 
then run `toil-bwa run aws:us-west-2:example-jobstore-bucket --batchSystem=mesos --mesosMaster mesos-master:5050`
to use the AWS job store and mesos batch system. 

## Dockerized Pipeline
To run the dockerized bwa alignment pipeline, please see [this link](https://github.com/BD2KGenomics/cgl-docker-lib/tree/alex-dockerized-pipelines/bwa-alignment-cgl-pipeline) in cgl-docker-lib.

Note: this link points to a branch off of master. Once the branch has been merged, this link will be changed to a URL to master.
