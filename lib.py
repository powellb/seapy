#!/usr/bin/env python
"""
  lib.py
  
  State Estimation and Analysis for PYthon

  Library of utilities for general seapy module, imported into the namespace
  when importing the seapy module

  Written by Brian Powell on 10/18/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np

def adddim(fld,size=1):
    """
    Given a 2D field, fld, add a new first dimension of given size
    """
    return np.transpose( \
       np.kron(fld,np.ones(size)).reshape([fld.shape[0],fld.shape[1],size]), \
       [2,0,1])


def convolve_mask(fld, ksize=3, kernel=None):
    """
    Given a masked array, convolve data over the masking using the
    specified kernel size, or the provided kernel
    """
    if kernel == None:
        center=np.round(ksize/2)
        kernel=np.ones([ksize,ksize])
        kernel[center,center]=0.0
    # Convolve the mask
    msk=np.ma.getmaskarray(fld)
    count=np.convolve( (~msk).view(np.int8).ravel(), kernel.ravel(), \
                      'same').reshape(msk.shape)
    nfld=np.convolve(fld.ravel(),kernel.ravel(),'same').reshape(msk.shape)
    lst=np.nonzero(msk is True and count>0)
    msk[lst] = False
    fld[lst] = nfld[lst] / count[lst]

def earth_distance(lon1, lat1, lon2, lat2):
    """
    earth_distance(lon1, lat1, lon2, lat2)

    Compute the distance between lat/lon points
    """
    epsilon = 0.99664718940443;  # This is Sqrt(1-epsilon^2)
    radius = 6378137; # Radius in meters
    d2r = np.pi/180.0
    
    # Using trig identities of tan(atan(b)), cos(atan(b)), sin(atan(b)) for 
    # working with geocentric where lat_gc = atan(epsilon * tan(lat))
    tan_lat = epsilon * np.tan(d2r*lat1)
    cos_lat = 1 / np.sqrt(1.0+tan_lat**2)
    sin_lat = tan_lat / np.sqrt(1.0+tan_lat**2)
    tan_lat = epsilon * np.tan(d2r*lat2)
    cos_latj = 1 / np.sqrt(1.0+tan_lat**2)
    sin_latj = tan_lat / np.sqrt(1.0+tan_lat**2)

    return radius*np.sqrt(2.0*(1.0-cos_lat*cos_latj* \
                          np.cos(d2r*(lon1-lon2))-sin_lat*sin_latj))

pass