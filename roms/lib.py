#!/usr/bin/env python
"""
  lib.py

  General ROMS utils

  Written by Brian Powell on 05/24/13
  Copyright (c)2016 University of Hawaii under the BSD-License.
"""

import numpy as np
from seapy.lib import default_epoch
import netCDF4

fields = {"zeta": {"grid": "rho", "dims": 2},
          "ubar": {"grid": "u", "dims": 2, "rotate": "vbar"},
          "vbar": {"grid": "v", "dims": 2, "rotate": "ubar"},
          "u": {"grid": "u", "dims": 3, "rotate": "v"},
          "v": {"grid": "v", "dims": 3, "rotate": "u"},
          "temp": {"grid": "rho", "dims": 3},
          "salt": {"grid": "rho", "dims": 3}}


def stretching(vstretching=2, theta_s=2, theta_b=0.1, hc=100, s_rho=10,
               w_grid=False):
    """
    Compute the stretching function for ROMS

    Parameters
    ----------
    vstretching : int, optional
        stretching algorithm type
    theta_s: float, optional
        value of surface theta
    theta_b: float, optional
        value of bottom theta
    hc: int, optional
        critical depth
    s_rho: int, optional
        number of s-levels
    w_grid: bool, optional
        solve stretching on the w-grid

    Returns
    -------
    s, cs: array

    """
    ds = 1.0 / s_rho
    if w_grid:
        lev = np.arange(1, s_rho + 1)
    else:
        lev = np.arange(1, s_rho + 1) - 0.5
    s = (lev - s_rho) * ds

    if vstretching == 1:
        if theta_s > 0:
            ptheta = np.sinh(theta_s * s) / np.sinh(theta_s)
            rtheta = np.tanh(theta_s * (s + 0.5)) / \
                (2.0 * np.tanh(0.5 * theta_s)) - 0.5
            cs = (1.0 - theta_b) * ptheta + theta_b * rtheta
        else:
            cs = s
        pass
    elif vstretching == 2:
        if theta_s > 0:
            csur = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1.0)
            if theta_b > 0:
                cbot = -1.0 + np.sinh(theta_b * (s + 1.0)) / np.sinh(theta_b)
                weight = (s + 1.0) * (1.0 - s)
                cs = weight * csur + (1.0 - weight) * cbot
            else:
                cs = csur
        else:
            cs = s
        pass
    elif vstretching == 3:
        if theta_s > 0:
            exp_s = theta_s
            exp_b = theta_b
            alpha = 3
            cbot = np.log(np.cosh(alpha * (s + 1.0)**exp_b)) / \
                np.log(np.cosh(alpha)) - 1.0
            csur = -np.log(np.cosh(alpha * np.abs(s)**exp_s)) / \
                np.log(np.cosh(alpha))
            weight = (1.0 - np.tanh(alpha * (s + .5))) / 2.0
            cs = weight * cbot + (1.0 - weight) * csur
        else:
            cs = s
        pass
    elif vstretching == 4:
        if theta_s > 0:
            csur = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1.0)
        else:
            csur = -(s * s)
        pass
        if theta_b > 0:
            cs = (np.exp(theta_b * csur) - 1.0) / (1.0 - np.exp(-theta_b))
        else:
            cs = s
        pass
    elif vstretching == 5:
        s = -(lev * lev - 2 * lev * s_rho + lev + s_rho * s_rho - s_rho) / \
            (1.0 * s_rho * s_rho - s_rho) - \
            0.01 * (lev * lev - lev * s_rho) / (1.0 - s_rho)
        if theta_s > 0:
            csur = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1)
        else:
            csur = -(s * s)
        if theta_b > 0:
            cs = (np.exp(theta_b * (csur + 1.0)) - 1.0) / \
                (np.exp(theta_b) - 1.0) - 1.0
        else:
            cs = csur
        pass
    else:
        raise ValueError("stretching value must be between 1 and 5")

    return s, cs


def depth(vtransform=1, h=None, hc=100, scoord=None,
          stretching=None, zeta=0, w_grid=False):
    """
    Solve the depth of the given bathymetry in s-levels.

    Parameters
    ----------
    vtransform : int, optional
        transform algorithm type
    h: array, optional
        value of bottom depths
    hc: int, optional
        critical depth
    scoord: array
        s coordinates from stretching method
    stretching: array
        stretching values from stretching method
    zeta: array
        sea surface height to add to bottom
    w_grid: bool, optional
        solve stretching on the w-grid

    Returns
    -------
    s, cs: array

    """
    if h is None or scoord is None or stretching is None:
        raise AttributeError("you must supply h, scoord, and stretching")
    if scoord.size != stretching.size:
        raise ValueError(
            "the stretching and scoord arrays must be the same size")
    N = scoord.size
    hinv = 1 / h
    h = np.asanyarray(h)
    r = np.arange(0, N)
    if w_grid:
        N = N + 1
    z = np.zeros(np.hstack((N, h.shape)))

    if vtransform == 1:
        cff = hc * (scoord - stretching)
        for k in r:
            z0 = cff[k] + stretching[k] * h
            z[k, :] = z0 + zeta * (1.0 + z0 * hinv)
    elif vtransform == 2:
        cff = 1 / (hc + h)
        for k in r:
            cff1 = hc * scoord[k] + h * stretching[k]
            z[k, :] = zeta + (zeta + h) * cff * cff1
    else:
        raise ValueError("transform value must be between 1 and 2")
    if w_grid:
        z[1:-1, :] = z[0:-2, :]
        z[0, :] = -h

    return z


def thickness(vtransform=1, h=None, hc=100, scoord=None,
              stretching=None, zeta=0):
    """
     Given the transform method, the bathyemtry, the sea surface height
     [optional], the critical depth, the s-coordinates, and stretching
     function (returned by the roms.utils.stretching function), determine
     the thickness of the grid cells
    """
    # We need the w-coordinate depths
    z_w = depth(vtransform, h, hc, scoord, stretching, zeta, True)
    return z_w[1:, :, :] - z_w[0:-1, :, :]


def get_time(nc, tvar=None, epoch=None):
    """
    Load the time vector from a netCDF file as a datetime array.

    Parameters
    ----------
    nc : netCDF4.Dataset,
      netcdf input file
    tvar : string, optional
      time variable to load. If not specified, it will find the
      time variable from predefined

    Returns
    -------
    ndarray,
       Array of datetimes if no epoch is supplied. If epoch, array
       is in days since epoch
    """
    tvar = tvar if tvar else get_timevar(nc)
    times = netCDF4.num2date(nc.variables[tvar][:],
                             nc.variables[tvar].units)
    if not epoch:
        return times
    else:
        return np.asarray([ (t - epoch).total_seconds() * secs2day for t in times ])


def get_timevar(nc):
    """
    Find the appropriate time variable (bry_time, ocean_time, etc.) from a
    given netcdf file

    Parameters
    ----------
    nc : netCDF4.Dataset netcdf input file

    Returns
    -------
    time: string

    """
    for time in ("ocean_time", "time", "bry_time", "wind_time",
                 "clim_time", "zeta_time"):
        if time in nc.variables:
            return time
    return None


def get_reftime(nc, epoch=default_epoch):
    """
    Given a ROMS netCDF4 file, return the reference time for the file. This
    is the timebase of the record dimension in the format:
    "<units> since <reftime>"

    Parameters
    ----------
    nc : netCDF4 dataset
        Input ROMS file
    epoch_str : string, optional
        If lacking units, use this string as the units

    Returns
    -------
    timebase : datetime
        datetime of the origin for the file
    time : string
        name of variable used to generate the base (None if default)
    """
    try:
        time = get_timevar(nc)
        return netCDF4.num2date(0, nc.variables[time].units), time
    except AttributeError:
        return epoch, None

pass

