#!/usr/bin/env python
"""
  Module to store all environment variables in an easy dictionary

  Written by Brian Powell on 04/26/13
  Copyright (c)2017 University of Hawaii under the BSD-License.
"""


import os


class options(object):
    pass

opt = options()
for key in os.environ.keys():
    setattr(opt, key, os.environ[key])

