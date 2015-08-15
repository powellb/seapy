#!/usr/bin/env python
"""
  analysis.py

  Methods to assist in the analysis of ROMS fields

  Written by Brian Powell on 05/24/15
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function
import numpy as np
from joblib import Parallel, delayed
from seapy import oavol, oasurf, adddim
from seapy.model import asgrid

# Create a function to do the multiprocessing vertical interpolation.
# This function is used because we don't want to construct an array of
# tuples from the oa routine.
def __dinterp(x,y,z,dat,fz,pmap):
    ndat, pm = oavol(x, y, z, dat, x, y, fz, pmap, 5, 1, 1)
    ndat = np.squeeze(ndat)
    return np.ma.masked_where(np.abs(ndat)>9e10, ndat, copy=False)

def constant_depth(field, grid, depth, zeta=None, threads=-2):
    """
    Find the values of a 3-D field at a constant depth for all times given.

    Parameters
    ----------
    field : ndarray,
        ROMS 3-D field to interpolate onto a constant depth level. Can be
        two- or three-dimensional array (first dimension assumed to be time).
    grid : seapy.model.grid or string,
        Grid that defines the depths and stretching for the field given
    depth : float,
        Depth (in meters) to find all values
    zeta : ndarray, optional,
        ROMS zeta field corresponding to field if you wish to apply the SSH
        correction to the depth calculations.
    threads : int, optional,
        Number of threads to use for processing

    Returns
    -------
    nfield : ndarray,
        Values from ROMS field on the given constant depth
    """

    # Make sure our inputs are all valid
    grid = asgrid(grid)
    if np.ndim(field)==3:
        field = adddim(field)
    if zeta is not None and np.ndim(zeta==2):
        zeta = adddim(zeta)
    depth = depth if depth < 0 else -depth

    # Set up some arrays
    x, y = np.meshgrid(np.arange(field.shape[-1]),np.arange(field.shape[-2]))
    fz, pmap = oasurf(x,y,x,x,y,None,5,1,1)
    # pu.db
    fz = adddim(np.ones(x.shape))*depth
    # Loop over all times, generate new field at depth
    nfield = np.ma.array(Parallel(n_jobs=threads,verbose=2)\
                  (delayed(__dinterp)(x,y,grid.depth_rho,
                                      np.squeeze(field[i,:,:,:]),fz,pmap)
               for i in range(field.shape[0])), copy=False)

    return nfield

