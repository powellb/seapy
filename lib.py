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
import os
import re
import datetime

secs2day = 1.0/86400.0

def adddim(fld,size=1):
    """
    replicate a field and add a new first dimension with given size

    Parameters
    ----------
    fld : array_like
        Input field.
    size : int, optional
        Size of additional first dimension

    Returns
    -------
    fld : array
    """
    fld=np.atleast_1d(fld)
    s=np.ones(fld.ndim+1)
    s[0]=size
    return np.tile(fld,s)


def convolve_mask(data, ksize=3, kernel=None, copy=True):
    """
    Convolve data over the missing regions of a mask
    
    Parameters
    ----------
    data : masked array_like
        Input field.
    ksize : int, optional
        Size of square kernel
    kernel : ndarray, optional
        Define a convolution kernel. Default is averaging
    copy : bool, optional
        If true, a copy of input array is made

    Returns
    -------
    fld : masked array
    """
    fld = np.ma.array(data, copy=copy)
    # Make sure ksize is odd
    ksize = int(ksize+1) if int(ksize)%2==0 else int(ksize)
    if fld.ndim > 3 or fld.ndim < 2:
        raise AttributeError("Can only convolve 2- or 3-D fields")
    if ksize < 3:
        raise ValueError("ksize must be greater than or equal to 3")
        
    if kernel is None:
        center=np.round(ksize/2)
        kernel=np.ones([ksize,ksize])
        kernel[center,center]=0.0

    # Convolve the mask
    msk=np.ma.getmaskarray(fld)
    if fld.ndim == 2:
        count=ndimage.convolve((~msk).view(np.int8), kernel, 
                               mode="constant", cval=0.0)
        nfld=ndimage.convolve(fld.data*(~msk).view(np.int8), kernel,
                               mode="constant", cval=0.0)
    else:
        kernel=np.expand_dims(kernel, axis=3)
        count=np.transpose( ndimage.convolve( 
                (~msk).view(np.int8).transpose(1,2,0), kernel,
                mode="constant", cval=0.0),(2,0,1))
        nfld=np.transpose( ndimage.convolve( 
            (fld.data*(~msk).view(np.int8)).transpose(1,2,0), kernel,
            mode="constant", cval=0.0),(2,0,1))
        
    lst=np.nonzero(np.logical_and(msk, count>0))
    msk[lst] = False
    fld[lst] = nfld[lst] / count[lst]
    return fld

def earth_distance(lon1, lat1, lon2, lat2):
    """
    Compute the distance between lat/lon points
    
    Parameters
    ----------
    lon1 : array_like or scalar
        Input array of source longitude(s)
    lat1 : array_like or scalar
        Input array of source latitude(s)
    lon2 : array_like or scalar
        Input array of destination longitude(s)
    lon2 : array_like or scalar
        Input array of destination longitude(s)

    Returns
    -------
    distance : array or scalar of distance in meters
    
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
    cos_lat = 1.0 / np.sqrt(1.0+tan_lat**2)
    sin_lat = tan_lat / np.sqrt(1.0+tan_lat**2)
    tan_lat = epsilon * np.tan(d2r*lat2.astype(np.float64))
    cos_latj = 1.0 / np.sqrt(1.0+tan_lat**2)
    sin_latj = tan_lat / np.sqrt(1.0+tan_lat**2)

    return radius*np.sqrt(2.0*(1.0-cos_lat*cos_latj* \
                  np.cos(d2r*(lon1-lon2))-sin_lat*sin_latj))

def rotate(u, v, angle):
    """
    Rotate a vector field by the given angle

    Parameters
    ----------
    u : array like
        Input u component
    v : array like
        Input v component
    angle : array like
        Input angle of rotation

    Returns
    -------
    rotated_u, rotated_v : array
    """
    u=np.asanyarray(u)
    v=np.asanyarray(v)
    angle=np.asanyarray(angle)
    sa=np.sin(angle)
    ca=np.cos(angle)

    return u*ca - v*sa, u*sa + v*ca
    
def vecfind(a, b, tolerance=0):
    """
    Find all occurences of b in a within the given tolerance and return
    the indices into a and b that correspond.

    Parameters
    ----------
    a : array
        Input vector
    b : array
        Input vector
    tolerance : float, optional
        Input tolerance for how close a==b

    Returns
    -------
    index_a, index_b : arrays of indices for each vector where values are equal,
            such that a[index_a] == b[index_b]
    
    """
    a=np.asanyarray(a)
    b=np.asanyarray(b)
    index_a=[]
    index_b=[]
    for i,c in enumerate(b):
        d=np.abs(a-c)
        x=np.nonzero(d<tolerance)[0]
        if x:
            index_a.append(np.where(d==np.min(d))[0][0])
            index_b.append(i)
    return index_a, index_b
    
def list_files(path=".", regex=".*"):
    """
    list all file names in the given path that conform to the regular
    expression pattern.

    Parameters
    ----------
    path : string, optional
        Input directory path
    regex : string, optional
        Input regular expression string to filter filenames

    Returns
    -------
    files : array

    """
    files=[]
    prog=re.compile(regex)
    for file in os.listdir(path):
        if prog.search(file) != None:
            files.append(file)
    return files
    
def day2date(day=0,epoch=datetime.datetime(2000,1,1)):
    """
    Return a datetime object from the number of days since the epoch

    Parameters
    ----------
    day : scalar
        Input day number
    epoch : datetime
        Date of epoch

    Returns
    -------
    date : datetime
    """
    return epoch + datetime.timedelta(days=day)

def date2day(date=datetime.datetime(2000,1,1),epoch=datetime.datetime(2000,1,1)):
    """
    Compute the fractional number of days elapsed since the epoch to the date
    given.

    Parameters
    ----------
    date : datetime
        Input date
    epoch : datetime
        Date of epoch

    Returns
    -------
    numdays : scalar
    """
    return (date-epoch).total_seconds() * secs2day

def today2day(epoch=datetime.datetime(2000,1,1)):
    """
    Return the day number of today (UTC time) since the epoch.

    Parameters
    ----------
    epoch : datetime
        Date of epoch

    Returns
    -------
    numdays : scalar
    """
    return date2day(datetime.datetime.utcnow(),epoch)

pass