## University of California, Santa Cruz Genomics Institute
### GATK-compatible Alignment

If there are any questions please contact John Vivian (jtvivian@gmail.com). 
If you find any errors or corrections please feel free to make a pull request.  
Feedback of any kind is appreciated.

## Overview

This pipeline accepts two fastq files (by URL) to be aligned into a BAMFILE, which is the final output of the pipeline.
A launch script is provided for 4 different references (b37, hg19, hg38, and hg38 no alternative loci).

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
 
## Output

This pipeline produces a BAMFILE for a given sample.

## Running / Help

It is recommended to use the associated launch scripts which provide default arguments needed to run the pipeline. 
It is likely that the job store positional argument, `--workDir`, and `--output-dir` arguments will need to be modified.
To run a pipeline after dependencies have been installed, simply:

* `git clone https://github.com/BD2KGenomics/toil-scripts`
* `/toil-scripts/src/toil_scripts/batch_alignment/launch_bwa_hg38_no_alt.sh`

Due to PYTHONPATH issues, help can be found by typing:

* `cd toil-scripts/src`
* `python -m toil_scripts.batch_alignment.bwa_alignment --help`
 