#!/usr/bin/env python
"""
  lib.py

  General ROMS utils

  Written by Brian Powell on 05/24/13
  Copyright (c)2017 University of Hawaii under the MIT-License.
"""

import numpy as np
import seapy
from seapy.lib import default_epoch, secs2day
import netCDF4

fields = {"zeta": {"grid": "rho", "dims": 2},
          "ubar": {"grid": "u", "dims": 2, "rotate": "vbar"},
          "vbar": {"grid": "v", "dims": 2, "rotate": "ubar"},
          "u": {"grid": "u", "dims": 3, "rotate": "v"},
          "v": {"grid": "v", "dims": 3, "rotate": "u"},
          "temp": {"grid": "rho", "dims": 3},
          "salt": {"grid": "rho", "dims": 3}}
ids = {1: "zeta", 2: "ubar", 3: "vbar", 4: "u", 5: "v", 6: "temp", 7: "salt"}


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
    z: ndarray,
      depth of grid cells

    """
    if h is None or scoord is None or stretching is None:
        raise AttributeError("you must supply h, scoord, and stretching")
    if scoord.size != stretching.size:
        raise ValueError(
            "the stretching and scoord arrays must be the same size")
    N = scoord.size
    hinv = 1 / h
    h = np.asanyarray(h)
    wk = 0
    r = range(N)
    if w_grid:
        N = N + 1
        wk = 1
    z = np.zeros(np.hstack((N, h.shape)))

    if vtransform == 1:
        cff = hc * (scoord - stretching)
        for k in r:
            z0 = cff[k] + stretching[k] * h
            z[k + wk, :] = z0 + zeta * (1.0 + z0 * hinv)
    elif vtransform == 2:
        cff = 1 / (hc + h)
        for k in r:
            cff1 = hc * scoord[k] + h * stretching[k]
            z[k + wk, :] = zeta + (zeta + h) * cff * cff1
    else:
        raise ValueError("transform value must be between 1 and 2")
    if w_grid:
        z[0, :] = -h

    return z


def thickness(vtransform=1, h=None, hc=100, scoord=None,
              stretching=None, zeta=0):
    """
    Get the thickness of the grid cells for the given sigma-parameters.

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
    hz : array,
      thickness
    """
    # Get the w-coordinate depths and return the differenc
    z_w = depth(vtransform, h, hc, scoord, stretching, zeta, True)
    return z_w[1:, :, :] - z_w[0:-1, :, :]


def gen_boundary_region(shp, north=None, east=None, west=None, south=None,
                        kind='linear'):
    """
    Generate a masked field varying from 1 at the boundary to 0 in the
    middle along each of the specified boundaries. This is used to create
    nudging and sponge fields to save into the respective ROMS files.

    Parameters
    ----------
    shp : tuple,
      The shape of the grid to use
    north : int, optional,
      The size of the region in the north boundary
    south : int, optional,
      The size of the region in the south boundary
    east : int, optional,
      The size of the region in the east boundary
    west : int, optional,
      The size of the region in the west boundary
    kind : string, optional,
      The type of transition:
         'linear' (default)
         'cosine'

    Returns
    -------
    fld : np.ma.array,
      array containing boundary values ranging from 0 to 1. masked values
      were not set by the routine, but the fill_value is set to 0.
    """
    fld = np.ma.zeros(shp, fill_value=0)
    fld[:] = np.ma.masked

    # Set up a dictionary to define how to deal with each boundary.
    # The tuple is (dimension, array_end, rotate)
    dirs = {"north": (shp[1], True, True),
            "south": (shp[1], False, True),
            "east": (shp[0], True, False),
            "west": (shp[0], False, False)}
    ref = locals()
    for d in dirs:
        # Set the factor to generate the values
        nx = ref[d]
        if nx is None or nx == 0:
            continue
        x = np.arange(nx)
        if kind == "cosine":
            x = np.cos(np.pi / (2.0 * nx) * x)[::-1]
        else:
            x = 1.0 / nx * x
        x = np.tile(x[::-1], [dirs[d][0], 1])
        # If the boundary is the end, flip it
        sl = np.array([slice(None, None, None), slice(None, nx, None)])
        if dirs[d][1]:
            x = np.fliplr(x)
            sl[1] = slice(-nx, None, None)
        # If the dimensions are rotated, transpose
        if dirs[d][2]:
            x = np.transpose(x)
            sl = sl[::-1]
        sl = (sl[0], sl[1])
        fld[sl] = np.maximum(fld.filled()[sl], x)

    return fld


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
        return np.asarray([(t - epoch).total_seconds() * secs2day for t in times])


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


def omega(grid, u, v, zeta=0, scale=False):
    """
    Compute the vertical velocity on s-grid.

    Parameters
    ----------
    grid : seapy.model.grid,
      The grid to use for the calculations
    u : ndarray,
      The u-field at a given time
    v : ndarray,
      The v-field at a given time
    zeta : ndarray, optional,
      The zeta-field at a given time
    scale : bool, optional,
      If True, return omega in [m s**-1];
      If [False], return omega in [m**3 s**-1]

    Returns
    -------
    omega : ndarray,
      Vertical Velocity on s-grid
    """
    u = np.ma.array(u, fill_value=0, copy=False).filled()
    v = np.ma.array(v, fill_value=0, copy=False).filled()
    zeta = np.ma.array(zeta, fill_value=0, copy=False).filled()

    z_w = depth(grid.vtransform, grid.h, grid.hc, grid.s_rho,
                grid.cs_r, zeta=zeta, w_grid=True)
    thick = z_w[1:, :, :] - z_w[0:-1, :, :]
    Huon = u * seapy.model.rho2u(thick / grid.pn)
    Hvon = v * seapy.model.rho2v(thick / grid.pm)
    W = np.zeros((grid.n + 1, grid.eta_rho, grid.xi_rho))
    W[:-1, 1:-1, 1:-1] = np.cumsum((Huon[:, 1:-1, :-1] - Huon[:, 1:-1, 1:]) +
                                   (Hvon[:, :-1, 1:-1] - Hvon[:, 1:, 1:-1]), axis=0)
    W -= ((z_w - z_w[0, :, :]) / (z_w[-1, :, :] - z_w[0, :, :])) * W[-1, :, :]
    if scale:
        return W * grid.pn * grid.pm
    return W


def wvelocity(grid, u, v, zeta=0, scale=False):
    """
    Compute "true" vertical velocity

    Parameters
    ----------
    grid : seapy.model.grid,
      The grid to use for the calculations
    u : ndarray,
      The u-field at a given time
    v : ndarray,
      The v-field at a given time
    zeta : ndarray, optional,
      The zeta-field at a given time
    scale : bool, optional,
      If True, return omega in [m s**-1];
      If [False], return omega in [m**3 s**-1]

    Returns
    -------
    w : ndarray,
      Vertical Velocity
    """
    u = np.ma.array(u, fill_value=0, copy=False).filled()
    v = np.ma.array(v, fill_value=0, copy=False).filled()
    zeta = np.ma.array(zeta, fill_value=0, copy=False).filled()

    z_r = depth(grid.vtransform, grid.h, grid.hc, grid.s_rho,
                grid.cs_r, zeta=zeta)
    z_w = depth(grid.vtransform, grid.h, grid.hc, grid.s_rho,
                grid.cs_r, zeta=zeta, w_grid=True)
    thick = z_w[1:, :, :] - z_w[0:-1, :, :]
    totdep = z_w[-1, :, :] - z_w[0, :, :]

    # Get omega
    sw = omega(grid, u, v, zeta)

    # Compute quasi-horizontal motions (Ui + Vj)*GRAD s(z)
    udel = u * seapy.model.rho2u(z_r * grid.pm)
    vdel = v * seapy.model.rho2v(z_r * grid.pn)
    vert = 0.25 * (udel[:, 1:-1, 1:] + udel[:, 1:-1, :-1] +
                   vdel[:, 1:, 1:-1] + vdel[:, :-1, 1:-1])

    # Compute barotropic velocity
    ubar = np.sum(u * seapy.model.rho2u(thick), axis=0) / \
        seapy.model.rho2u(totdep)
    vbar = np.sum(v * seapy.model.rho2v(thick), axis=0) / \
        seapy.model.rho2v(totdep)
    wrk = (ubar[1:-1, :-1] - ubar[1:-1, 1:] + vbar[:-1, 1:-1] -
           vbar[1:, 1:-1]) / totdep[1:-1, 1:-1]

    # Cubic interpolation
    cff1 = 0.375
    cff2 = 0.75
    cff3 = 0.125
    cff4 = 0.5625
    cff5 = 0.0625

    w = np.zeros((grid.n + 1, grid.eta_rho, grid.xi_rho))

    # Do the bottom two layers
    slope = (z_r[0, 1:-1, 1:-1] - z_w[0, 1:-1, 1:-1]) / \
        (z_r[1, 1:-1, 1:-1] - z_r[0, 1:-1, 1:-1])
    w[0, 1:-1, 1:-1] = cff1 * (vert[0, :, :] -
                               slope * (vert[1, :, :] - vert[0, :, :])) + \
        cff2 * vert[0, :, :] + cff3 * vert[1, :, :]
    w[1, 1:-1, 1:-1] = grid.pm[1:-1, 1:-1] * grid.pn[1:-1, 1:-1] * \
        (sw[1, 1:-1, 1:-1] + wrk[:, :] *
         (z_w[1, 1:-1, 1:-1] - z_w[0, 1:-1, 1:-1])) + \
        cff1 * vert[0, :, :] + cff2 * vert[1, :, :] + cff3 * vert[2, :, :]

    # Do the middle bit
    w[2:-2, 1:-1, 1:-1] = grid.pm[1:-1, 1:-1] * grid.pn[1:-1, 1:-1] * \
        (sw[2:-2, 1:-1, 1:-1] + wrk[:, :] *
         (z_w[2:-2, 1:-1, 1:-1] - z_w[0, 1:-1, 1:-1])) + \
        cff4 * (vert[1:-2, :, :] + vert[2:-1, :, :]) + \
        cff5 * (vert[:-3, :, :] + vert[3:, :, :])

    # Do the upper two layers
    slope = (z_w[-1, 1:-1, 1:-1] - z_r[0, 1:-1, 1:-1]) / \
        (z_r[-1, 1:-1, 1:-1] - z_r[-2, 1:-1, 1:-1])
    w[-2, 1:-1, 1:-1] = grid.pm[1:-1, 1:-1] * grid.pn[1:-1, 1:-1] * \
        (sw[-2, 1:-1, 1:-1] + wrk[:, :] *
         (z_w[-1, 1:-1, 1:-1] - z_w[0, 1:-1, 1:-1])) + \
        cff1 * vert[-1, :, :] + cff2 * vert[-2, :, :] + cff3 * vert[-3, :, :]
    w[-1, 1:-1, 1:-1] = grid.pm[1:-1, 1:-1] * grid.pn[1:-1, 1:-1] * \
        wrk[:, :] * (z_w[-1, 1:-1, 1:-1] - z_w[0, 1:-1, 1:-1]) + \
        cff1 * (vert[-1, :, :] + slope * (vert[-1, :, :] - vert[-2, :, :])) + \
        cff2 * vert[-1, :, :] - cff3 * vert[-2, :, :]

    # Done
    return w


pass

