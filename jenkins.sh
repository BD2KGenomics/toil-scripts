#!/usr/bin/env bash
# Create s3am venv
virtualenv s3am
. s3am/bin/activate
pip install --pre s3am
deactivate
# Create Toil venv
virtualenv toil
. toil/bin/activate
pip install pytest toil boto tqdm
# Expose binaries to the PATH
mkdir bin
ln -snf ${PWD}/s3am/bin/s3am bin/
export PATH=$PATH:${PWD}/bin
# Set PYTHONPATH
export PYTHONPATH=$(python -c 'from os.path import abspath as a;import sys;print a("src")' $0)
py.test src --doctest-modules --junitxml=test-report.xml
