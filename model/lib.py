#!/usr/bin/env python
"""
  lib.py
  
  Library of utilities for ocean models, imported into the namespace
  when importing the model module

  Written by Brian Powell on 10/18/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np
import pudb

def convolve_land( field, mask=None, kernel=3 ):
    """
    Given a field, mask [optional], and kernel size, use a convolution kernel
    to replace unknown (or masked) values.
    """
    if field.ndim == 2:
        t=1
        m,n = field.shape
    else:
        t,m,n=field.shape
    
    # Apply the mask if given one
    if mask != None:
        mm,mn=mask.shape
        if mm != m or mn != n:
            raise ValueError("mask is different size from field")
        fld = fld * adddim(mask,t)
    
def _cgrid_rho_vel( field, dim ):
    shp = np.array(field.shape)
    fore = np.product([shp[i] for i in np.arange(0,dim)])
    aft =  np.product([shp[i] for i in np.arange(dim+1,field.ndim)])
    nfld = 0.5 * (field.reshape([fore,shp[dim],aft])[:,0:-1,:] + 
                  field.reshape([fore,shp[dim],aft])[:,1:,:])
    shp[dim] = shp[dim]-1
    return nfld.reshape(shp)
    
def rho2u( field ):
    return _cgrid_rho_vel(field, field.ndim-1)

def rho2v( field ):
    return _cgrid_rho_vel(field, field.ndim-2)
    
def u2rho( field ):
    shp = np.array(field.shape)
    nshp = shp.copy()
    nshp[-1]=nshp[-1]+1
    fore = np.product([shp[i] for i in np.arange(0,field.ndim-1)])
    nfld = np.ones([fore,nshp[-1]])
    nfld[:,1:-1] = 0.5 * \
                 (field.reshape([fore,shp[-1]])[:,0:-1] + 
                  field.reshape([fore,shp[-1]])[:,1:])
    nfld[:,0] = nfld[:,1] + (nfld[:,2]-nfld[:,3])
    nfld[:,-1] = nfld[:,-2] + (nfld[:,-2]-nfld[:,-3])
    return nfld.reshape(nshp)

def v2rho( field ):
    shp = np.array(field.shape)
    nshp = shp.copy()
    nshp[-2]=nshp[-2]+1
    fore = np.product([shp[i] for i in np.arange(0,field.ndim-2)])
    nfld = np.ones([fore,nshp[-2],nshp[-1]])
    nfld[:,1:-1,:] = 0.5 * \
                 (field.reshape([fore,shp[-2],shp[-1]])[:,0:-1,:] + 
                  field.reshape([fore,shp[-2],shp[-1]])[:,1:,:])
    nfld[:,0,:] = nfld[:,1,:] + (nfld[:,2,:]-nfld[:,3,:])
    nfld[:,-1,:] = nfld[:,-2,:] + (nfld[:,-2,:]-nfld[:,-3,:])
    return nfld.reshape(nshp)
    
    