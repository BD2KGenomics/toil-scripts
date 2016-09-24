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
from pkg_resources import parse_version, require, DistributionNotFound

def check_provided(distribution, min_version, max_version=None, optional=False):
    min_version = parse_version(min_version)
    if max_version is not None:
        max_version = parse_version(max_version)

    messages = []

    toil_missing = 'Cannot find a valid installation of Toil.'
    dist_missing = 'Cannot find an installed copy of the %s distribution, typically provided by Toil.' % distribution
    version_too_low = 'The installed copy of %s is out of date. It is typically provided by Toil.' % distribution
    version_too_high = 'The installed copy of %s is too new. It is typically provided by Toil.' % distribution
    required_version = 'Setup requires version %s or higher' % (min_version,)
    required_version += '.' if max_version is None else ', up to but not including %s.' % (max_version,)
    install_toil = 'Installing Toil should fix this problem.'
    upgrade_toil = 'Upgrading Toil should fix this problem.'
    reinstall_dist = 'Uninstalling %s and reinstalling Toil should fix this problem.' % distribution
    reinstall_toil = 'Uninstalling Toil and reinstalling it should fix this problem.'
    footer = ("Setup doesn't install Toil automatically to give you a chance to choose any of the optional extras "
              "that Toil provides. More on installing Toil at http://toil.readthedocs.io/en/latest/installation.html.")
    try:
        # This check will fail if the distribution or any of its dependencies are missing.
        version = require(distribution)[0].version
    except DistributionNotFound:
        version = None
        if not optional:
            messages.extend([toil_missing if distribution == 'toil' else dist_missing, install_toil])
    else:
        if parse_version(version) < min_version:
            messages.extend([version_too_low, required_version,
                             upgrade_toil if distribution == 'toil' else reinstall_dist])
        elif max_version is not None and max_version < parse_version(version):
            messages.extend([version_too_high, required_version,
                             reinstall_toil if distribution == 'toil' else reinstall_dist])
    if messages:
        messages.append(footer)
        raise RuntimeError(' '.join(messages))
    else:
        return version


toil_version = check_provided('toil', min_version='3.3.0', max_version='3.5.0')
check_provided('bd2k-python-lib', min_version='1.14a1.dev29' )
check_provided('boto', min_version='2.38.0', optional=True)

kwargs = dict(
    name='toil-scripts',
    version=version,
    description='A repository of genomic workflows developed by the UCSC Computational Genomics lab ',
    author='UCSC Computational Genomics Lab',
    author_email='cgl-toil@googlegroups.com',
    url="https://github.com/BD2KGenomics/toil-scripts",
    install_requires=[
        'toil-lib==1.1.0a1.dev65',
        'pyyaml==3.11'],
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

print("\n\n"
      "Thank you for installing toil-scripts! If you want to run Toil on a cluster in a cloud, please reinstall it "
      "with the appropriate extras. To install AWS/EC2 support for example, run "
      "\n\n"
      "pip install toil[aws,mesos]==%s"
      "\n\n"
      "on every EC2 instance. Refer to Toil's documentation at http://toil.readthedocs.io/en/latest/installation.html "
      "for more information."
      % toil_version)
