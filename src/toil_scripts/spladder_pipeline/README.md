## University of California, Santa Cruz Genomics Institute
### SplAdder Pipeline (HG38)

This README provides an overview of the SplAdder pipeline. If there are any questions
please contact John Vivian (jtvivian@gmail.com). If you find any errors or corrections please feel free to make a 
pull request.  Feedback of any kind is appreciated.

## Overview

This pipeline accepts URLS (`http://`, `https://`, `ftp://`, `file://`, `gnos://` prefixes) of sample **tarfiles**
that contain RNA-seq fastq data files. It assumes that every fastq file follows the standard of having either 
 **R1** or **R2** in the file name.  If there are several fastq files in the tarfile, (typically in the context 
 of sequencing lanes), they will be concatenated together into one **R1.fq.gz** and one **R2.fq.gz** file.  Note,
 local files follow the syntax of `file://` followed by the full path to that file's location, e.g. 
 `file:///home/ubuntu/sample.tar`.
 
## Output

This pipeline produces a tar.gz file for a given sample that contains two subdirs:

- splAdder
    - alignment.filt.hdf5  
    - alignment.hdf5  
    - genes_graph_conf3.alignment.pickle
- variants_and_qc
    - output.vcf.gz
    - qccounts.tsv
        
The output tarball is *stamped* with the UUID for the sample (e.g. UUID.tar.gz). The UUID of each sample is 
specified in the config file. 
   
   
## Dependencies

This pipeline has been tested on Ubuntu 14.04, but should also run on other unix based systems.  `apt-get` and `pip`
often require `sudo` privilege, so if the below commands fail, try prepending `sudo`.  If you do not have `sudo` 
privileges you will need to build these tools from source, or bug a sysadmin about how to get them (they don't mind).  
 
#### General Dependencies

1. Python 2.7
2. Curl — apt-get install curl
3. Docker — http://docs.docker.com/engine/installation/

#### Python Dependencies

1. Toil — pip install toil
2. S3AM — pip install --pre s3am (optional, needed for uploading output to S3)
    
    
## Running / Help

It is recommended to use the associated launch scripts (launch_spladder_sm.sh for single machine mode and 
launch_spladder_mesos.sh for distributed clusters) which take care of filling in the default arguments
needed to run a pipeline. To run a pipeline after dependencies have been installed, simply:

- `git clone https://github.com/BD2KGenomics/toil-scripts`
- `/toil-scripts/src/toil_scripts/spladder_pipeline/launch_spladder_sm.sh`

Due to PYTHONPATH issues, help can be found by:

- `cd toil-scripts/src`
- `python -m stoil_scripts.spladder_pipeline.spladder_pipeline --help`
