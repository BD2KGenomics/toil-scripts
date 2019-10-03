# Copyright (C) 2015 UCSC Computational Genomics Lab
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function

import sys

from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand
from version import version


kwargs = dict(
    name='toil-scripts',
    version=version,
    description='A repository of genomic workflows developed by the UCSC Computational Genomics lab ',
    author='UCSC Computational Genomics Lab',
    author_email='cgl-toil@googlegroups.com',
    url="https://github.com/BD2KGenomics/toil-scripts",
    install_requires=[
        'toil-lib==1.2.0a1.dev126',
        'pyyaml==5.1'],
    tests_require=[
        'pytest==2.8.3'],
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={
        'console_scripts': [
            'toil-bwa = toil_scripts.bwa_alignment.bwa_alignment:main',
            'toil-exome = toil_scripts.exome_variant_pipeline.exome_variant_pipeline:main',
            'toil-germline = toil_scripts.gatk_germline.germline:main']})


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        # Sanitize command line arguments to avoid confusing Toil code attempting to parse them
        sys.argv[1:] = []
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

kwargs['cmdclass'] = {'test': PyTest}

setup(**kwargs)

toil_version = '3.5.0a1.dev277'

print("\n\n"
      "Thank you for installing toil-scripts! If you want to run Toil on a cluster in a cloud, please reinstall it "
      "with the appropriate extras. To install AWS/EC2 support for example, run "
      "\n\n"
      "pip install toil[aws,mesos]==%s"
      "\n\n"
      "on every EC2 instance. Refer to Toil's documentation at http://toil.readthedocs.io/en/latest/installation.html "
      "for more information."
      % toil_version)
