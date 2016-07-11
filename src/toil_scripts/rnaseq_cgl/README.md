## University of California, Santa Cruz Genomics Institute
### Guide: Running the CGL HG38 RNA-seq Pipeline using Toil

This guide attempts to walk the user through running this pipeline from start to finish. If there are any questions
please contact John Vivian (jtvivian@gmail.com). If you find any errors or corrections please feel free to make a 
pull request.  Feedback of any kind is appreciated.


## Overview

RNA-seq fastqs are combined, aligned, and quantified with 2 different methods (RSEM and Kallisto)

This pipeline produces a tarball (tar.gz) file for a given sample that contains:

    RSEM: TPM, FPKM, counts and raw counts (parsed from RSEM output)
    Kallisto: abundance.tsv, abundance.h5, and a JSON of run information

The output tarball is *stamped* with the UUID for the sample (e.g. UUID.tar.gz). 

## Installation

Toil-scripts is now pip installable! `pip install toil-scripts` for a toil-stable version 
or `pip install --pre toil-scripts` for cutting edge development version.

Type: `toil-rnaseq` to get basic help menu and instructions

To decrease the chance of versioning conflicts, install toil-scripts into a virtualenv: 

- `virtualenv ~/toil-scripts` 
- `source ~/toil-scripts/bin/activate`
- `pip install toil`
- `pip install toil-scripts`

If Toil is already installed globally (true for CGCloud users), or there are global dependencies (like Mesos),
use virtualenv's `--system-site-packages` flag.
 

## Dependencies

This pipeline has been tested on Ubuntu 14.04, but should also run on other unix based systems.  `apt-get` and `pip`
often require `sudo` privilege, so if the below commands fail, try prepending `sudo`.  If you do not have `sudo` 
privileges you will need to build these tools from source, or bug a sysadmin about how to get them (they don't mind). 

#### General Dependencies

    1. Python 2.7
    2. Curl         apt-get install curl
    3. Docker       http://docs.docker.com/engine/installation/

#### Python Dependencies

    1. Toil         pip install toil
    2. S3AM         pip install --pre s3am (optional, needed for uploading output to S3)

## Inputs

The CGL RNA-seq pipeline requires input files in order to run. These files are hosted on Synapse and can 
be downloaded after creating an account which takes about 1 minute and is free. 

* Register for a [Synapse account](https://www.synapse.org/#!RegisterAccount:0)
* Either download the samples from the [website GUI](https://www.synapse.org/#!Synapse:syn5886029) or use the Python API
* `pip install synapseclient`
* `python`
    * `import synapseclient`
    * `syn = synapseclient.Synapse()`
    * `syn.login('foo@bar.com', 'password')`
    * Get the RSEM reference (1 G)
        * `syn.get('syn5889216', downloadLocation='.')`
    * Get the Kallisto index (2 G)
        * `syn.get('syn5886142', downloadLocation='.')`
    * Get the STAR index (25 G)
        * `syn.get('syn5886182', downloadLocation='.')`


## General Usage
 
1. Type `toil-rnaseq generate` to create an editable manifest and config in the current working directory.
2. Parameterize the pipeline by editing the config.
3. Fill in the manifest with information pertaining to your samples.
4. Type `toil-rnaseq run [jobStore]` to execute the pipeline.

## Example Commands

Run sample(s) locally using the manifest

1. `toil-rnaseq generate`
2. Fill in config and manifest
3. `toil-rnaseq run ./example-jobstore`

Toil options can be appended to `toil-rnaseq run`, for example:
`toil-rnaseq run ./example-jobstore --retryCount=1 --workDir=/data`

For a complete list of Toil options, just type `toil-rnaseq run -h`

Run a variety of samples locally

1. `toil-rnaseq generate-config`
2. Fill in config
3. `toil-rnaseq run ./example-jobstore --retryCount=1 --workDir=/data --samples \
    s3://example-bucket/sample_1.tar file:///full/path/to/sample_2.tar https://sample-depot.com/sample_3.tar`

## Example Config

```
star-index: s3://cgl-pipeline-inputs/rnaseq_cgl/ci/starIndex_chr6.tar.gz
kallisto-index: s3://cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx
rsem-ref: s3://cgl-pipeline-inputs/rnaseq_cgl/ci/rsem_ref_chr6.tar.gz
output-dir: /data/my-toil-run
s3-dir: s3://my-bucket/test/rnaseq
ssec: 
gt-key: 
wiggle: true
save-bam: true
ci-test:
fwd-3pr-adapter: AGATCGGAAGAG
rev-3pr-adapter: AGATCGGAAGAG
```

Example with local input files

```
star-index: file://data/starIndex_chr6.tar.gz
kallisto-index: file://data/kallisto_hg38.idx
rsem-ref: file://data/rsem_ref_chr6.tar.gz
output-dir: /data/my-toil-run
s3-dir: s3://my-bucket/test/rnaseq
ssec: 
gt-key: 
wiggle: true
save-bam: true
ci-test:
fwd-3pr-adapter: AGATCGGAAGAG
rev-3pr-adapter: AGATCGGAAGAG
```

## Distributed Run

To run on a distributed AWS cluster, see [CGCloud](https://github.com/BD2KGenomics/cgcloud) for instance provisioning, 
then run `toil-rnaseq run aws:us-west-2:example-jobstore-bucket --batchSystem=mesos --mesosMaster mesos-master:5050`
to use the AWS job store and mesos batch system. 
