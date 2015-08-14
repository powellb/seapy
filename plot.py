#!/usr/bin/env python
"""
  plot.py
  
  State Estimation and Analysis for PYthon

  Module with plotting utilities

  Written by Brian Powell on 10/18/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np
from scipy import ndimage
import os
import re
from matplotlib import pyplot as plt

def stackbar(x, y, colors=None, **kwargs):
    """
    Given an array of vectors in y, draw a bar chart for each one stacked on
    the prior.
    """
    s=y[0,:]
    if colors is None:
        colors = [ "" for i in range(0,y.shape[0]) ]
    plt.bar(x, y[0,:], color=colors[0], **kwargs)
    for i in range(1,y.shape[0]):
        plt.bar(x, y[i,:], color=colors[i], bottom=s, **kwargs)
        s=s+y[i,:]

    
