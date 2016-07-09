## ADAM Überscript

The ADAM überscript runs the `adam_gatk_pipeline` on a cluster whose
size is automatically managed. The überscript itself manages the cluster
by collecting metrics on job execution, and looking at the cluster occupancy.
The process of running the script is somewhat complex and is described in this
document.

## Setup

This document assumes that all machines are run in a separate
[cgcloud](https://github.com/BD2KGenomics/cgcloud) namespace, given as `<CGCLOUD_NAMESPACE>`.

1. Create the head node:
    - `cgcloud create -t m3.medium -n /<CGCLOUD_NAMESPACE>/ generic-ubuntu-trusty-box`
2. On the head node (ssh via `cgcloud ssh -n /<CGCLOUD_NAMESPACE>/ generic-ubuntu-trusty-box`):
    - `sudo apt-get install git make libxml2-dev libxslt1-dev python-dev`
    - `virtualenv cgcloud-master`
    - `source cgcloud-master/bin/activate`
    - `pip install boto tqdm`
    - Add to bash profile:
        - `export CGCLOUD_ZONE=us-west-2a`
        - `export CGCLOUD_PLUGINS=cgcloud.mesos:cgcloud.toil`
        - `export CGCLOUD_KEYPAIRS="__me__ @@developers"`
        - `export CGCLOUD_NAMESPACE=/<CGCLOUD_NAMESPACE>/`
        - `export CGCLOUD_ME="ubuntu@$(hostname)"`
    - `source ~/.profile`, or log out and log in
    - Generate a throw-away ssh key via `ssh-keygen -C "$CGCLOUD_ME" -f ~/.ssh/id_rsa -N ""`
    - create a `~/.boto`
    - `git clone https://github.com/BD2KGenomics/cgcloud.git`
    - from the cgcloud repository:
        - `make develop sdist`
    - `cgcloud register-key ~/.ssh/id_rsa.pub`
    - `cgcloud create -IT toil-box`
    - `git clone https://github.com/BD2KGenomics/toil-scripts.git`
3. Running the pipelines:
    These should be run in succession, each in a separate screen, from the head node.     
    - Launch Cluster:
        `PYTHONPATH=$PYTHONPATH:~/toil-scripts/src python -m toil_scripts.adam_uberscript.adam_uberscript launch-cluster -c <cluster name> -S ~/toil-scripts/ -t r3.8xlarge -b ~/.boto -M <manifest file> -C <path to config file>`
    - Launch Pipeline:
        `PYTHONPATH=$PYTHONPATH:~/toil-scripts/src python -m toil_scripts.adam_uberscript.adam_uberscript launch-pipeline -c <cluster name> -j <jobstore> -C <path to config file>`
    - Launch Metrics:
        `PYTHONPATH=$PYTHONPATH:~/toil-scripts/src python -m toil_scripts.adam_uberscript.adam_uberscript launch-metrics -c "<cluster name>" -j <jobstore> --namespace $CGCLOUD_NAMESPACE -t r3.8xlarge`