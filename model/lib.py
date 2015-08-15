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
from seapy import convolve_mask

def _cgrid_rho_vel( rho, dim, fill ):
    """
    Private Method: Compute the u- or v-grid velocity from a rho-field for a c-grid
    """
    rho = np.ma.array(rho, copy=False)
    if fill:
        rho = convolve_mask(rho, copy=True)
    shp = np.array(rho.shape)
    fore = np.product([shp[i] for i in np.arange(0,dim)])
    aft =  np.product([shp[i] for i in np.arange(dim+1,rho.ndim)])
    nfld = 0.5 * (rho.reshape([fore,shp[dim],aft])[:,0:-1,:].filled(np.nan) +
                  rho.reshape([fore,shp[dim],aft])[:,1:,:].filled(np.nan))
    shp[dim] = shp[dim]-1
    return np.ma.fix_invalid(nfld.reshape(shp), copy=False, fill_value=1e+37)

def rho2u( rho, fill=False ):
    """
    Put the rho field onto the u field for the c-grid

    Parameters
    ----------
    rho : masked array like
        Input rho field
    fill : bool, optional
        Fill the masked data before moving grids

    Returns
    -------
    u : masked array
    """
    return _cgrid_rho_vel(rho, rho.ndim-1, fill)

def rho2v( rho, fill=False ):
    """
    Put the rho field onto the v field for the c-grid

    Parameters
    ----------
    rho : masked array like
        Input rho field
    fill : bool, optional
        Fill the masked data before moving grids

    Returns
    -------
    v : masked array
    """
    return _cgrid_rho_vel(rho, rho.ndim-2, fill)

def u2rho( u, fill=False ):
    """
    Put the u field onto the rho field for the c-grid

    Parameters
    ----------
    u : masked array like
        Input u field
    fill : bool, optional
        Fill the masked data before moving grids

    Returns
    -------
    rho : masked array
    """
    u = np.ma.array(u, copy=False)
    if fill:
        u = convolve_mask(u, copy=True)
    shp = np.array(u.shape)
    nshp = shp.copy()
    nshp[-1]=nshp[-1]+1
    fore = np.product([shp[i] for i in np.arange(0,u.ndim-1)])
    nfld = np.ones([fore,nshp[-1]])
    nfld[:,1:-1] = 0.5 * \
                 (u.reshape([fore,shp[-1]])[:,0:-1].filled(np.nan) +
                  u.reshape([fore,shp[-1]])[:,1:].filled(np.nan))
    nfld[:,0] = nfld[:,1] + (nfld[:,2]-nfld[:,3])
    nfld[:,-1] = nfld[:,-2] + (nfld[:,-2]-nfld[:,-3])
    return np.ma.fix_invalid(nfld.reshape(nshp), copy=False, fill_value=1e+37)

def v2rho( v, fill=False ):
    """
    Put the v field onto the rho field for the c-grid

    Parameters
    ----------
    v : masked array like
        Input v field
    fill : bool, optional
        Fill the masked data before moving grids

    Returns
    -------
    rho : masked array
    """
    v = np.ma.array(v, copy=False)
    if fill:
        v = convolve_mask(v, copy=True)
    shp = np.array(v.shape)
    nshp = shp.copy()
    nshp[-2]=nshp[-2]+1
    fore = np.product([shp[i] for i in np.arange(0,v.ndim-2)])
    nfld = np.ones([fore,nshp[-2],nshp[-1]])
    nfld[:,1:-1,:] = 0.5 * \
                 (v.reshape([fore,shp[-2],shp[-1]])[:,0:-1,:].filled(np.nan) +
                  v.reshape([fore,shp[-2],shp[-1]])[:,1:,:].filled(np.nan))
    nfld[:,0,:] = nfld[:,1,:] + (nfld[:,2,:]-nfld[:,3,:])
    nfld[:,-1,:] = nfld[:,-2,:] + (nfld[:,-2,:]-nfld[:,-3,:])
    return np.ma.fix_invalid(nfld.reshape(nshp), copy=False, fill_value=1e+37)


