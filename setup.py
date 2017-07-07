# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

import os
import sys
from setuptools import find_packages, setup
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

rootpath = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    return open(os.path.join(rootpath, *parts), 'r').read()


LICENSE = read('LICENSE.txt')
long_description = read('README.md')

# Dependencies.
with open('requirements.txt') as f:
    tests_require = f.readlines()
install_requires = [t.strip() for t in tests_require]

package_data = {
    '': ['constituents.npz',
         'hawaii_coast/*',
         'roms/cdl/*.cdl']
}

config = Configuration('')
config.add_extension('oalib', sources='src/oalib.F')
config.add_extension('hindices', sources='src/hindices.F')

config = dict(
    name='seapy',
    version='0.1.0',
    description='State Estimation and Analysis in PYthon',
    long_description=long_description,
    author='Brian Powell',
    author_email='powellb@hawaii.edu',
    url='https://github.com/powellb/seapy',
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: MIT License',
        ],
    packages=find_packages(),
    package_data=package_data,
    ext_package='seapy.external',
    scripts=['bin/convert_clim.py', 'bin/convert_frc.py'],
    license=LICENSE,
    install_requires=install_requires,
    zip_safe=False,
    **config.todict(),
)


setup(**config)
