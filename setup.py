# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

import os
import sys
from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

# Only way to get --noopt defaulted into the build environment
# since f2py_options is not working below
sys.argv[:] = sys.argv[:1] + ['config_fc', '--noopt'] + sys.argv[1:]

rootpath = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    return open(os.path.join(rootpath, *parts), 'r').read()


with open("README.md", "r") as f:
    long_description = f.read()

# Dependencies.
with open('requirements.txt') as f:
    tests_require = f.readlines()
install_requires = [t.strip() for t in tests_require]

package_data = {
    '': ['constituents.npz',
         'hawaii_coast/*',
         'roms/cdl/*.cdl',
         'roms/cobalt/*.cdl']
}

config = Configuration('')
flags = [] if os.name == 'nt' else ['-fPIC']
# ifort generated libraries produce invalid results in interpolation (NOT
# OBVIOUS)
config.add_extension('oalib', sources='src/oalib.f',
                     # f2py_options=["noopt"],
                     extra_f77_compile_args=flags)
config.add_extension('hindices', sources='src/hindices.f',
                     # f2py_options=["noopt"],
                     extra_f77_compile_args=flags)

config = dict(
    name=os.getenv('PACKAGE_NAME', 'seapy'),
    version='0.5.2',
    license='MIT',
    description='State Estimation and Analysis in PYthon',
    long_description=long_description,
    long_description_content_type='text/markdown',
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
    install_requires=install_requires,
    zip_safe=False,
    **config.todict()
)


setup(**config)
