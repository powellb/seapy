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
from scipy import ndimage

secs2day = 1.0/86400.0

def adddim(fld,size=1):
    """
    Given a field, replicate with a new first dimension of given size
    """
    fld=np.asanyarray(fld)
    return np.transpose( \
       np.kron(fld,np.ones(size)).reshape(np.append(np.asarray(fld.shape),size)),
       np.append(fld.ndim,np.arange(0,fld.ndim))) 


def convolve_mask(fld, ksize=3, kernel=None):
    """
    Given a masked array, convolve data over the masking using the
    specified kernel size, or the provided kernel
    """
    fld = np.ma.array(fld, copy=False)
    if fld.ndim > 3 or fld.ndim < 2:
        raise AttributeError("Can only convolve 2- or 3-D fields")
        
    if kernel == None:
        center=np.round(ksize/2)
        kernel=np.ones([ksize,ksize])
        kernel[center,center]=0.0

    # Convolve the mask
    msk=np.ma.getmaskarray(fld)
    if fld.ndim == 2:
        count=ndimage.convolve((~msk).view(np.int8), kernel)
        nfld=ndimage.convolve(fld.data*(~msk).view(np.int8), kernel)
    else:
        kernel=np.expand_dims(kernel, axis=3)
        count=np.transpose( ndimage.convolve( 
                (~msk).view(np.int8).transpose(1,2,0), kernel),(2,0,1))
        nfld=np.transpose( ndimage.convolve( 
            (fld.data*(~msk).view(np.int8)).transpose(1,2,0), kernel),(2,0,1))
        
    lst=np.nonzero(np.logical_and(msk, count>0))
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
    
    lon1 = np.asanyarray(lon1)
    lat1 = np.asanyarray(lat1)
    lon2 = np.asanyarray(lon2)
    lat2 = np.asanyarray(lat2)
    
    # Using trig identities of tan(atan(b)), cos(atan(b)), sin(atan(b)) for 
    # working with geocentric where lat_gc = atan(epsilon * tan(lat))
    tan_lat = epsilon * np.tan(d2r*lat1.astype(np.float64))
    cos_lat = 1 / np.sqrt(1.0+tan_lat**2)
    sin_lat = tan_lat / np.sqrt(1.0+tan_lat**2)
    tan_lat = epsilon * np.tan(d2r*lat2.astype(np.float64))
    cos_latj = 1 / np.sqrt(1.0+tan_lat**2)
    sin_latj = tan_lat / np.sqrt(1.0+tan_lat**2)

    return radius*np.sqrt(2.0*(1.0-cos_lat*cos_latj* \
                  np.cos(d2r*(lon1-lon2))-sin_lat*sin_latj))

def rotate(u, v, angle):
    """
    rotate(u,v,angle)
    
    Rotate a vector field, given by u and v, by the angle given.
    """
    u=np.asanyarray(u)
    v=np.asanyarray(v)
    sa=np.sin(np.asanyarray(angle))
    ca=np.cos(np.asanyarray(angle))

    return u*ca - v*sa, u*sa + v*ca
    
    
pass