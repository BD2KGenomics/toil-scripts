## University of California, Santa Cruz Genomics Institute
### Repository of Toil-based Genomic Pipelines
 
This repository contains a collection of genomic pipelines developed by the UCSC Computational Genomics Lab.

To install these pipelines, just type: `pip install toil-scripts`

Entrypoints are supplied to pipelines that work with the current stable release of Toil:
    
- Fastq to BAM alignment: `toil-bwa`
- CGL RNA-seq Pipeline (hg38): `toil-rnaseq`
- UNC best practices RNA-seq pipeline (hg19): `toil-rnaseq-unc` 
- SplAdder alternative splicing pipeline: `toil-spladder`
- GATK exome variant pipeline: `toil-exome`

Each pipeline has its own README that provides instructions on how to get up and running. 
The general dependencies for these pipelines are:

1. [Toil](https://github.com/BD2KGenomics/toil)
2. [Docker](https://www.docker.com/)

Our group utilizes genomics tools encapsulated within Docker containers for portability.  Each of these
pipelines can be run locally on your laptop, on a baremetal cluster, or on a cloud provider. 

If there are any questions please contact the Toil team at: Toil@googlegroups.com 
If you find any errors or corrections please feel free to make a pull request.  Feedback of any kind is appreciated.
