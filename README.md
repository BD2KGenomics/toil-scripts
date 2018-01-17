## [DEPRECATED]

## University of California, Santa Cruz Genomics Institute
### Repository of Toil-based Genomic Pipelines

Toil-Scripts is now deprecated and no longer being maintained.  It may be run and used at the user's discretion.  Dependencies are locked at: toil-lib==1.2.0a1.dev126, and toil==3.5.0a1.dev277 .  Other versions may work, but are not gauranteed to run.

To install these pipelines, type: `pip install toil-scripts`

Entrypoints are supplied to pipelines that work with the current stable release of Toil:
    
- Fastq to BAM alignment: `toil-bwa`
- GATK exome variant pipeline: `toil-exome`

Source code and installation instructions for the **CGL RNA-seq pipeline** have moved to: 
https://github.com/BD2KGenomics/toil-rnaseq

Each pipeline has its own README that provides instructions on how to get up and running. 
The general dependencies for these pipelines are:

1. [Toil](https://github.com/BD2KGenomics/toil)
2. [Docker](https://www.docker.com/)

Our group utilizes genomics tools encapsulated within Docker containers for portability.  Each of these
pipelines can be run locally on your laptop, on a baremetal cluster, or on a cloud provider. 

If there are any questions please contact the Toil team at: Toil@googlegroups.com 
If you find any errors or corrections please feel free to make a pull request.  Feedback of any kind is appreciated.
