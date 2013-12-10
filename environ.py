#!/usr/bin/env python
"""
  Class to abstract the netcdf information for a ROMS netcdf file/url.

  This does not load any actual data information until accessed. The accessor
  method will attempt to extract the variable of interest.

  Written by Brian Powell on 04/26/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import os

class options(object):
    pass
    
opt = options()
for key in os.environ.keys():
    setattr(opt,key,os.environ[key])

