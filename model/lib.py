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
    
def _cgrid_rho_vel( field, dim ):
    """
    Compute the u- or v-grid velocity from a rho-field for a c-grid
    """
    shp = np.array(field.shape)
    fore = np.product([shp[i] for i in np.arange(0,dim)])
    aft =  np.product([shp[i] for i in np.arange(dim+1,field.ndim)])
    nfld = 0.5 * (field.reshape([fore,shp[dim],aft])[:,0:-1,:] + 
                  field.reshape([fore,shp[dim],aft])[:,1:,:])
    if np.ma.isMA(field):
        nmsk = np.logical_or(field.mask.reshape([fore,shp[dim],aft])[:,0:-1,:], 
                             field.mask.reshape([fore,shp[dim],aft])[:,1:,:])
        shp[dim] = shp[dim]-1
        return np.ma.array(nfld.reshape(shp), mask=nmsk.reshape(shp), copy=False)
    else:
        shp[dim] = shp[dim]-1
        return nfld.reshape(shp)
    
def rho2u( field ):
    """
    Put the rho field onto the u field for the c-grid
    """
    return _cgrid_rho_vel(field, field.ndim-1)

def rho2v( field ):
    """
    Put the rho field onto the v field for the c-grid
    """
    return _cgrid_rho_vel(field, field.ndim-2)
    
def u2rho( field ):
    """
    Put the u field onto the rho field for the c-grid
    """
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
    if np.ma.isMA(field):
        nmsk = np.zeros([fore,nshp[-1]])
        nmsk[:,1:-1] = np.logical_or(field.mask.reshape([fore,shp[-1]])[:,0:-1], 
                  field.mask.reshape([fore,shp[-1]])[:,1:])
                 
        nmsk[:,0] = nmsk[:,1]
        nmsk[:,-1] = nmsk[:,-2]
        return np.ma.array(nfld.reshape(nshp), mask=nmsk.reshape(nshp), copy=False)
    else:
        return nfld.reshape(nshp)

def v2rho( field ):
    """
    Put the v field onto the rho field for the c-grid
    """
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
    if np.ma.isMA(field):
        nmsk = np.zeros([fore,nshp[-2],nshp[-1]])
        nmsk[:,1:-1,:] = np.logical_or( \
                 field.mask.reshape([fore,shp[-2],shp[-1]])[:,0:-1,:], \
                 field.mask.reshape([fore,shp[-2],shp[-1]])[:,1:,:])
        nmsk[:,0,:] = nmsk[:,1,:]
        nmsk[:,-1,:] = nmsk[:,-2,:]
        return np.ma.array(nfld.reshape(nshp), mask=nmsk.reshape(nshp), copy=False)
    else:
        return nfld.reshape(nshp)
    
    