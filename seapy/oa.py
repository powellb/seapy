#!/usr/bin/env python
"""
  oa

  Objective analysis.  This function will interpolate data using the
  fortran routines written by Emanuelle Di Lorenzo and Bruce Cornuelle

  Written by Brian Powell on 10/08/13
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""

import numpy as np
from seapy.external import oalib

__bad_val = -999999.0


def oasurf(x, y, d, xx, yy, pmap=None, weight=10, nx=2, ny=2, verbose=False):
    """
    Objective analysis interpolation for 2D fields

    Parameters
    ----------
    x: array [2-D]
        x-values of source data
    y: array [2-D]
        y-values of source data
    d: array [2-D]
        data values of source
    xx: array [2-D]
        x-values of destination
    yy: array [2-D]
        y-values of destination
    pmap: array, optional
        weighting array to map between source and destination.
        NOTE: best to save this after using to prevent recomputing
        weights for every interpolate
    weight: int, optional
        number of neighbor points to consider for every destination point
    nx: int, optional
        decorrelation lengthscale in x [same units as x]
    ny: int, optional
        decorrelation lengthscale in y [same units as y]
    verbose : bool, optional
        display information within the OA routine

    Returns
    -------
    new_data: ndarray
        data interpolated onto the new grid
    pmap: ndarray
        weighting map used in the interpolation

    """
    # Do some error checking
    nx = ny if nx == 0 else nx
    ny = nx if ny == 0 else ny
    d = np.ma.masked_invalid(d, copy=False)

    # Generate a mapping weight matrix if not passed
    if pmap is None:
        pmap = np.zeros([xx.size, weight], order='F')

    # Call FORTRAN library to objectively map
    vv, err = oalib.oa2d(x.ravel(), y.ravel(),
                         d.filled(__bad_val).ravel(),
                         xx.ravel(), yy.ravel(), nx, ny,
                         pmap, verbose)

    # Reshape the results and return
    return np.ma.masked_equal(vv.reshape(xx.shape), __bad_val, copy=False), \
        pmap


def oavol(x, y, z, v, xx, yy, zz, pmap=None, weight=10, nx=2, ny=2,
          verbose=False):
    """
    Objective analysis interpolation for 3D fields

    Parameters
    ----------
    x: array [2-D]
        x-values of source data
    y: array [2-D]
        y-values of source data
    z: array [3-D]
        z-values of source data
    v: array [3-D]
        data values of source
    xx: array [2-D]
        x-values of destination
    yy: array [2-D]
        y-values of destination
    zz: array [3-D]
        z-values of destination
    pmap: array, optional
        weighting array to map between source and destination.
        NOTE: best to save this after using to prevent recomputing
        weights for every interpolate
    weight: int, optional
        number of neighbor points to consider for every destination point
    nx: int, optional
        decorrelation lengthscale in x [same units as x]
    ny: int, optional
        decorrelation lengthscale in y [same units as y]
    verbose : bool, optional
        display information within the OA routine

    Returns
    -------
    new_data: ndarray
        data interpolated onto the new grid
    pmap: ndarray
        weighting map used in the interpolation

    """
    # Do some error checking
    nx = ny if nx == 0 else nx
    ny = nx if ny == 0 else ny
    v = np.ma.masked_invalid(v, copy=False)

    # Generate a mapping weight matrix if not passed
    if pmap is None:
        # Build the map
        tmp = np.ones(x.ravel().shape, order='F')
        pmap = np.zeros([xx.size, weight], order='F')
        oalib.oa2d(x.ravel(), y.ravel(), tmp,
                   xx.ravel(), yy.ravel(), nx, ny, pmap, verbose)

    # Call FORTRAN library to objectively map
    vv, err = oalib.oa3d(x.ravel(), y.ravel(),
                         z.reshape(z.shape[0], -1).transpose(),
                         v.filled(__bad_val).reshape(
                             v.shape[0], -1).transpose(),
                         xx.ravel(), yy.ravel(),
                         zz.reshape(zz.shape[0], -1).transpose(),
                         nx, ny, pmap, verbose)

    # Reshape the results and return
    return np.ma.masked_equal(vv.transpose().reshape(zz.shape), __bad_val,
                              copy=False), pmap
