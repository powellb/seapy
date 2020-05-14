#!/usr/bin/env python
"""
  plot.py

  State Estimation and Analysis for PYthon

  Module with plotting utilities

  Written by Brian Powell on 10/18/13
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


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
    x = np.asarray(x)
    if colors is None:
        colors = ["" for i in range(0, y.shape[0])]
    # Stack positive and negative separately
    for op in ("__ge__", "__lt__"):
        d = getattr(y, op)(0)
        s = y[0, :] * 0
        l = np.where(d[0, :])[0]
        if np.any(l):
            plt.bar(x[l], y[0, l], color=colors[0], **kwargs)
        s[l] = y[0, l]
        for i in range(1, y.shape[0]):
            l = np.where(d[i, :])[0]
            if np.any(l):
                plt.bar(x[l], y[i, l], color=colors[i], bottom=s[l], **kwargs)
                s[l] += y[i, l]
