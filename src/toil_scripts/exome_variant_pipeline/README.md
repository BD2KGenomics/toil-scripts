## University of California, Santa Cruz Genomics Institute
### Guide: Running the CGL Exome Pipeline using Toil

This guide attempts to walk the user through running this pipeline from start to finish. If there are any questions
please contact John Vivian (jtvivian@gmail.com). If you find any errors or corrections please feel free to make a 
pull request.  Feedback of any kind is appreciated.

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Inputs](#inputs)
- [Usage](#general-usage)
- [Methods](#methods)

## Overview

A pair of Tumor/Normal exome BAMs are preprocessed (GATK), indels are found (Pindel), and variants are
called with two mutation callers (MuTect and MuSe).  This pipeline is modular — any part of the 
pipeline can be run on it's own. If preprocessing is selected, it will always occur before any of the other tools.

This pipeline produces a tarball (tar.gz) file for a given sample that contains:

- MuTect: Mutect.vcf, Mutect.cov, Mutect.out
- Pindel: pindel-config.txt pindel_BP pindel_CloseEndMapped pindel_D pindel_INT pindel_INT_final pindel_INV pindel_LI pindel_RP pindel_SI pindel_TD
- MuSe: Muse.vcf

The output tarball is *stamped* with the UUID for the sample (e.g. UUID.tar.gz). 

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
    
#### System Dependencies

This pipeline needs approximately 15G of RAM in order to run various GATK steps.  

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

The CGL exome pipeline requires input files in order to run. These files are hosted on Synapse and can 
be downloaded after creating an account which takes about 1 minute and is free. 

* Register for a [Synapse account](https://www.synapse.org/#!RegisterAccount:0)
* Either download the samples from the [website GUI](https://www.synapse.org/#!Synapse:syn5886029) or use the Python API
* `pip install synapseclient`
* `python`
    * `import synapseclient`
    * `syn = synapseclient.Synapse()`
    * `syn.login('foo@bar.com', 'password')`
    * Get the Reference Genome (3 G)
        * `syn.get('syn6128232', downloadLocation='.')`
    * Get the Phase VCF (0.3 G)
        * `syn.get('syn6128233', downloadLocation='.')`
    * Get the Mills VCF (0.1 G)
        * `syn.get('syn6128236', downloadLocation='.')`
    * Get the DBSNP VCF (10 G)
        * `syn.get('syn6128237', downloadLocation='.')`
    * Get the Cosmic VCF (0.01 G)
        * `syn.get('syn6128235', downloadLocation='.')`
        
A sample, which consists of a tumor and normal BAM file, can be passed via the command line options
`--normal`, `--tumor`, and `--uuid`. If wanting to run more than one sample, then the use the `toil-exome --generate-manifest`
command and fill in the manifest as instructed. 
All samples and inputs must be submitted as URLs with support for the following schemas: 
`http://`, `file://`, `s3://`, `ftp://`.


## General Usage
 
Type `toil-exome` to get basic help menu and instructions
 
1. Type `toil-exome generate` to create an editable manifest and config in the current working directory.
2. Parameterize the pipeline by editing the config.
3. Fill in the manifest with information pertaining to your samples.
4. Type `toil-exome run [jobStore]` to execute the pipeline.

## Example Commands

Run sample(s) locally using the manifest
1. `toil-exome generate`
2. Fill in config and manifest
3. `toil-exome run ./example-jobstore`

Toil options can be appended to `toil-exome run`, for example:
`toil-exome run ./example-jobstore --retryCount=1 --workDir=/data`

For a complete list of Toil options, just type `toil-exome run -h`

Run a variety of samples locally
1. `toil-exome generate-config`
2. Fill in config
3. `toil-exome run \
        ./example-jobstore \
        --retryCount 1 \
        --workDir /data \
        --normal s3://example-bucket/normal.bam \
        --tumor s3://example-bucket/tumor.bam \ 
        --uuid test-sample`

## Example Config

HG19
```
reference: s3://cgl-pipeline-inputs/variant_hg19/hg19.fa                     
phase: s3://cgl-pipeline-inputs/variant_hg19/1000G_phase1.indels.hg19.sites.vcf                  
mills: s3://cgl-pipeline-inputs/variant_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
dbsnp: s3://cgl-pipeline-inputs/variant_hg19/dbsnp_138.hg19.vcf
cosmic: s3://cgl-pipeline-inputs/variant_hg19/cosmic.hg19.vcf                 
run-mutect: true        
run-pindel: true        
run-muse: true          
preprocessing: true     
output-dir: /data/my-toil-run          
s3-dir: s3://my-bucket/test/exome
ssec:                   
gtkey:                  
ci-test:
```

B37
```
reference: https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/Homo_sapiens_assembly19.fasta
phase: https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/1000G_phase1.indels.hg19.sites.fixed.vcf
mills: https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf
dbsnp: https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/dbsnp_132_b37.leftAligned.vcf
cosmic: https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_b37/b37_cosmic_v54_120711.vcf
run-mutect: true        
run-pindel: true        
run-muse: true          
preprocessing: true     
output-dir:          
s3-output-dir:                 
ssec:                   
gtkey:                  
ci-test:
```

## Distributed Run

To run on a distributed AWS cluster, see [CGCloud](https://github.com/BD2KGenomics/cgcloud) for instance provisioning, 
then run `toil-exome run aws:us-west-2:example-jobstore-bucket --batchSystem=mesos --mesosMaster mesos-master:5050`
to use the AWS job store and mesos batch system. 

# Methods

## Tools

| Tool     | Version | Description                                                                                       |
|----------|---------|---------------------------------------------------------------------------------------------------|
| Samtools | 0.1.19  | Indexes bam files and creates reference index.                                                    |
| Picard   | 1.95    | Creates sequence dictionary for reference.                                                        |
| GATK     | 3.5     | Used for GATK preprocessing. Indel realignment (IR) and base quality score recalibration (BQSR).  |
| Pindel   | 0.2.5b6 | Computes insertions and deletions for each BAM file.                                              |
| MuSe     | 1.0     | Finds somatic variants from a pair of normal and tumor BAM files.                                 |
| MuTect   | 1.1.7   | Finds somatic variants from a pair of normal and tumor BAM files.                                 |

All tool containers can be found on our [quay.io account](quay.io/organization/ucsc_cgl).

## Reference Data

This pipeline is designed to work with HG19 and GRCh37. HG38 / GRCh38 support is in the works. All other input files
to this pipeline are part of the [GATK bundle](http://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it).

## Tool options

- Samtools index/faidx are run with default options
- Picard dictionary is run with default options
- GATK BQSR is run with default options
- GATK PR is run with `--emit_original_quals`

#### GATK - RTC

```
'-T', 'RealignerTargetCreator',
'-nt', cores,
'-R', ref_fasta,
'-I', sample_bam,
'-known', phase_vcf,
'-known', mills_vcf,
'--downsampling_type', 'NONE',
'-o', sample_intervals
```

#### GATK - IR

```
'-T', 'IndelRealigner',
'-R', ref_fasta,
'-I', sample_bam
'-known', phase_vcf,
'-known', mills_vcf,
'-targetIntervals', sample_intervals,
'--downsampling_type', 'NONE',
'-maxReads', '720000',
'-maxInMemory', '5400000', 
'-o', sample_indel_bam
```

#### MuTect

```
'--analysis_type', 'MuTect',
'--reference_sequence', ref_fasta,
'--cosmic', cosmic_vcf,
'--dbsnp', dbsnp_vcf,
'--input_file:normal', normal_bam,
'--input_file:tumor', tumor_bam,
'--tumor_lod', '10', 
'--initial_tumor_lod', '4.0',
'--out', 'mutect.out',
'--coverage_file', 'mutect.cov',
'--vcf', 'mutect.vcf'
```

#### MuSe

```
'--mode', 'wxs',
'--dbsnp', dbsnp_vcf,
'--fafile', ref_fasta,
'--tumor-bam', tumor_bam,
'--tumor-bam-index', tumor_bai,
'--normal-bam', normal_bam,
'--normal-bam-index', normal_bai,
'--outfile', muse_vcf,
'--cpus', cores
```

#### Pindel

pindel-config.txt contains mean insert size for tumor and normal bam
```
'-f', ref_fasta,
'-i', pindel_config
'--number_of_threads', cores,
'--minimum_support_for_event', '3',
'--report_long_insertions', 'true',
'--report_breakpoints', 'true',
'-o', 'pindel'
```

## Dockerized Pipeline
To run the dockerized exome var pipeline, please see [this link](https://github.com/BD2KGenomics/cgl-docker-lib/tree/alex-dockerized-pipelines/exome-cgl-pipeline) in cgl-docker-lib.

Note: this link points to a branch off of master. Once the branch has been merged, this link will be changed to a URL to master.
