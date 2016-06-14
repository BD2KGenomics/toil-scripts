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
from pkg_resources import parse_version

# Toil version check -- Raise warning instead of using intall_requires to avoid virtualenv conflicts
toil_min_version = '3.1.6'
toil_max_version = '3.2.0'

try:
    from toil.version import version as toil_version
except ImportError:
    raise RuntimeError('Cannot find Toil version {}.x for Toil. Read about installing Toil at: '
                       'http://toil.readthedocs.io/en/latest/installation.html'.format(toil_min_version))

if not parse_version(str(toil_min_version)) <= parse_version(toil_version) < parse_version(toil_max_version):
    raise RuntimeError('Need Toil version within range [{},{}). Read about installing Toil at: '
                       'http://toil.readthedocs.io/en/latest/installation.html'
                       ''.format(toil_min_version, toil_max_version))

kwargs = dict(
    name='toil-scripts',
    version=version,
    description='A repository of genomic workflows developed by the UCSC Computational Genomics lab ',
    author='UCSC Computational Genomics Lab',
    author_email='cgl-toil@googlegroups.com',
    url="https://github.com/BD2KGenomics/toil-scripts",
    install_requires=[
        'boto==2.38.0', # FIXME: Make an extra
        'tqdm==3.8.0', # FIXME: Remove once ADAM stops using it (superfluous import)
        'pyyaml==3.11'],
    tests_require=[
        'pytest==2.8.3'],
    package_dir={'': 'src'},
    packages=find_packages('src', exclude=['*.test']),
    entry_points={
        'console_scripts': [
            'toil-bwa = toil_scripts.batch_alignment.bwa_alignment:main',
            'toil-rnaseq = toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline:main',
            'toil-rnaseq-unc = toil_scripts.rnaseq_unc.rnaseq_unc_pipeline:main',
            'toil-spladder = toil_scripts.spladder_pipeline.spladder_pipeline:main',
            'toil-exome = toil_scripts.exome_variant_pipeline.exome_variant_pipeline:main']})


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

print("\n\nThank you for installing toil-scripts! If you want to run Toil in a cloud environment or on a distributed "
      "system, please install Toil with the desired extras, for example:\n\n"
      "'pip install toil[aws,mesos,azure,encryption]==%s'\n\n"
      "Please take a look at Toil's documentation for more information: http://toil.readthedocs.io/en/releases-3.1.x/"
      % toil_version)
