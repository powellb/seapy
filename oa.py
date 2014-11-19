#!/usr/bin/env python
"""
  oa
  
  Objective analysis.  This function will interpolate data using the 
  fortran routines written by Emanuelle Di Lorenzo and Bruce Cornuelle

  Written by Brian Powell on 10/08/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function
import numpy as np
import seapy.oalib

def oasurf(x,y,d,xx,yy,pmap=None,weight=10,nx=2,ny=2):
    # Do some error checking
    nx = ny if nx==0 else nx
    ny = nx if ny==0 else ny
    d = np.ma.fix_invalid(d, copy=False, fill_value=-999999.0)
    
    # Generate a mapping weight matrix if not passed
    if pmap is None:
        pmap=np.zeros([xx.size,weight],order="F")
    
    # Call FORTRAN library to objectively map
    vv, err = seapy.oalib.oa2d(x.ravel(),y.ravel(),d.ravel(),
                                 xx.ravel(), yy.ravel(), nx, ny, pmap)
    
    # Reshape the results and return
    return np.ma.fix_invalid(vv.reshape(xx.shape), copy=False, 
               fill_value=-999999.0), pmap
    
def oavol(x,y,z,v,xx,yy,zz,pmap=None,weight=10,nx=2,ny=2):
    # Do some error checking
    nx = ny if nx==0 else nx
    ny = nx if ny==0 else ny
    z = np.ma.fix_invalid(z, copy=False, fill_value=-999999.0)

    # Generate a mapping weight matrix if not passed
    if pmap is None:
        pmap=np.zeros([xx.size,weight],order="F")
        # Build the map
        seapy.oalib.oa2d(x.ravel(),y.ravel(),ones(x.shape),
                           xx.ravel(), yy.ravel(), nx, ny, pmap)
        
    # Call FORTRAN library to objectively map
    vv, err = seapy.oalib.oa3d(x.ravel(),y.ravel(),
                                 z.data.reshape(z.shape[0],-1).transpose(),
                                 v.reshape(v.shape[0],-1).transpose(),
                                 xx.ravel(), yy.ravel(), 
                                 zz.reshape(zz.shape[0],-1).transpose(), 
                                 nx, ny, pmap)
    
    # Reshape the results and return
    return np.ma.fix_invalid(vv.transpose().reshape(zz.shape), copy=False, 
               fill_value=-999999.0), pmap

    
