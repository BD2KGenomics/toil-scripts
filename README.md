## University of California, Santa Cruz Genomics Institute
### Repository of Toil-based Genomic Pipelines
 
 This repository contains a collection of genomic pipelines developed by the UCSC Computational Genomics Lab.
 These pipelines are used on a regular basis in-house and are freely available to anyone who would like to use
 them.

Pipelines currently supported for the latest version of Toil:
    
- Fastq to BAM alignment (GATK compatible format)
- Exome variant pipeline
- UNC best practices RNA-seq pipeline 
- Our own best practices RNA-seq pipeline (STAR + RSEM, Kallisto)
- Germline variant pipeline

Each pipeline has its own README that provides explicit instructions on how to get up and running. 
The general dependencies for these pipelines are:

1. [Toil](https://github.com/BD2KGenomics/toil)
2. [Docker](https://www.docker.com/)

Our group utilizes genomics tools encapsulated within Docker containers for portability.  Each of these
pipelines can be run locally on your laptop, on a baremetal cluster, or on a cloud provider. 

If there are any questions please contact John Vivian (jtvivian@gmail.com).
If you find any errors or corrections please feel free to make a pull request.  Feedback of any kind is appreciated.
