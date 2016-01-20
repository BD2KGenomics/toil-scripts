## University of California, Santa Cruz Genomics Institute
### Guide: Running TCGA Data on PODK (single machine mode)

This guide attempts to walk the user through running this pipeline from start to finish. If there are any questions
please contact John Vivian (jtvivian@gmail.com). If you find any errors or corrections please feel free to make a 
pull request.

## Dependencies

1. PODK Access and Docker permissions
    1. For docker permissions contact Erich Weiler (weiler@soe.ucsc.edu)
2. Python 2.7.3 (install from source into your home directory: /pod/home/USER)
3. virtualenv ([install from source](http://stackoverflow.com/questions/9348869/how-to-install-virtualenv-without-using-sudo))
4. Acquire a [CGHub permissions key](https://cghub.ucsc.edu/keyfile/keyfile.html) to access TCGA and place in your home directory

## Setup Steps

Go to your home directory (/pod/home/USER) and execute the following commands:

1. Create an instance of virtualenv: `virtualenv venv`
2. Activate venv:  `source venv/bin/activate`
3. Install toil:  `pip install toil`
4. Export LD_LIBRARY_PATH. Replace USER with your username: `export LD_LIBRARY_PATH=/pod/home/USER/python-2.7.3/lib`


## Run the Pipeline

1. Create a folder for your run: `mkdir test_run`
2. Download the CGL pipeline to your **test_run** folder.
    1. rnaseq_cgl_pipeline.py

Create a configuration file for your run.  `nano config.txt`, then paste one TCGA analysis_id per line.  
An example would be:

```
5557a728-1827-4aff-b28b-f004d835f9d6
2826301c-5d33-465a-99fa-401aea553a7f
bb9ecd73-ded7-4c4e-9674-d05647be7a22
f23c04eb-6f22-4b19-a6e3-90e7887f535e
```

Now make the launch script executable: `chmod u+x launch.sh`.

Create a launch script for your run.  `nano launch.sh` then paste the below code and change the appropriate lines.

``` bash
#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
#
# CHANGE THE BELOW LINES
USER=your_name_here
KEYNAME=CGHub_key_here
# Make sure pstore folder exists
mkdir -p /pod/pstore/${USER}
# Execution of pipeline
python rnaseq_cgl_pipeline.py \
/pod/pstore/${USER}/jstore \
--genetorrent /pod/home/${USER}/test_run/config.txt \
--genetorrent_key /pod/home/${USER}/${KEYNAME} \
--output_dir /pod/pstore/${USER}
--workDir=/scratch \
--logLevel=debug  \
#--restart
```

Since Docker isn't installed on the head node, you'll need to SSH into a single node to run your pipeline. 
For example: `ssh USER@podk-1-2`.  I generally run my pipelines in a screen, so type `screen`.
Now, make sure your **venv** is active and **LD_LIBRARY_PATH** is set.  If everything has been setup properly
you should be able to type `time ./launch.sh` to kick off the pipeline.
 
Save any stack traces to send to the author if you encounter issues.  