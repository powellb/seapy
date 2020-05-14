#!/usr/bin/env python
"""
  lib.py

  General ROMS utils

  Written by Brian Powell on 05/24/13
  Copyright (c)2020 University of Hawaii under the MIT-License.
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


def _get_calendar(var):
    """
    Get the proper calendar string from a netcdf file

    Parameters
    ----------
    var : netCDF4.variable

    Returns
    -------
    calendar type: string,
      The type of calendar system used
    convert : bool
      True if the calendar needs to be converted to datetime
    """
    # Set up the mapping for calendars
    default = 1
    calendar_conv = [False, False, False,
                     True, True, True, False, False, False]
    calendar_types = ['standard', 'gregorian', 'proleptic_gregorian', 'noleap',
                      'julian', 'all_leap', '365_day', '366_day', '360_day']
    cals = {v: v for v in calendar_types}
    cals['gregorian_proleptic'] = 'proleptic_gregorian'

    # Load the calendar type. If it is incorrectly specified (*cough* ROMS), change it
    for cal in ('calendar', 'calendar_type'):
        if hasattr(var, cal):
            cal = cals.get(str(getattr(var, cal)).lower(),
                           calendar_types[default])
            return cal, calendar_conv[calendar_types == cal]
    return calendar_types[default], calendar_conv[default]


def date2num(dates, nc, tvar=None):
    """
    Convert the datetime vector to number for the given netcdf files considering
    the units and the calendar type used. This is a wrapper to the netCDF4.date2num
    function to account for calendar strangeness in ROMS

    Parameters
    ----------
    dates : array of datetime.datetime
      Values to convert
    nc : netCDF4.Dataset,
      netcdf input file
    tvar : string, optional
      time variable to load. If not specified, it will find the
      time variable from predefined

    Returns
    -------
    ndarray,
       Array of values in the correct units/calendar of the netCDF file
    """
    tvar = tvar if tvar else get_timevar(nc)

    calendar, _ = _get_calendar(nc.variables[tvar])
    # Convert the times
    return netCDF4.date2num(dates,
                            nc.variables[tvar].units,
                            calendar=calendar)


def num2date(nc, tvar=None, records=None, epoch=None):
    """
    Load the time vector from a netCDF file as a datetime array, accounting
    for units and the calendar type used. This is a wrapper to the netCDF4.num2date
    function to account for calendar strangeness in ROMS

    Parameters
    ----------
    nc : netCDF4.Dataset,
      netcdf input file
    tvar : string, optional
      time variable to load. If not specified, it will find the
      time variable from predefined
    records : array or slice, optional
      the indices of records to load
    epoch : datetime.datetime, optional
      if you would like the values relative to an epoch, then
      specify the epoch to remove.

    Returns
    -------
    ndarray,
       Array of datetimes if no epoch is supplied. If epoch, array
       is in days since epoch
    """
    import datetime
    records = records if records is not None else np.s_[:]
    tvar = tvar if tvar else get_timevar(nc)
    calendar, convert = _get_calendar(nc.variables[tvar])

    # Load the times
    times = netCDF4.num2date(nc.variables[tvar][records],
                             nc.variables[tvar].units,
                             calendar=calendar)
    if convert:
        times = [datetime.datetime.strptime(t.strftime(), '%Y-%m-%d-%H:%M:%S')
                 for t in times]

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
                 "clim_time", "frc_time", "zeta_time"):
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
        tvar = get_timevar(nc)
        calendar, _ = _get_calendar(nc.variables[tvar])

        return netCDF4.num2date(0, nc.variables[tvar].units,
                                calendar=calendar), tvar
    except AttributeError:
        return epoch, None


def omega(grid, u, v, zeta=0, scale=True, work=False):
    """
    Compute the vertical velocity on s-grid.

    Parameters
    ----------
    grid : seapy.model.grid,
      The grid to use for the calculations
    u : ndarray,
      The u-field in time
    v : ndarray,
      The v-field in time
    zeta : ndarray, optional,
      The zeta-field in time
    scale : bool, optional,
      If [True], return omega in [m s**-1];
      If False, return omega in [m**3 s**-1]
    work : bool, optional,
      If True, return the work arrays:
        z_r : ndarray,
          Depth on rho-grid (time-varying if zeta != 0)
        z_w : ndarray,
          Depth on w-grid (time-varying if zeta != 0)
        thick_u : ndarray
          Thickness of the u-grid
        thick_v : ndarray
          Thickness of the v-grid
      If False, return only omega

    Returns
    -------
    omega : ndarray,
      Vertical Velocity on s-grid
    """
    grid = seapy.model.asgrid(grid)
    u = np.ma.array(u)
    v = np.ma.array(v)
    zeta = np.ma.array(zeta)

    # Check the sizes
    while u.ndim < 4:
        u = u[np.newaxis, ...]
    while v.ndim < 4:
        v = v[np.newaxis, ...]
    while zeta.ndim < 3:
        zeta = zeta[np.newaxis, ...]

    # Get the model grid parameters for the given thickness
    thick_u = u * 0
    thick_v = v * 0
    z_r = np.ma.zeros((u.shape[0], u.shape[1], zeta.shape[1], zeta.shape[2]))
    z_w = np.ma.zeros((u.shape[0], u.shape[1] + 1,
                       zeta.shape[1], zeta.shape[2]))
    for i in range(zeta.shape[0]):
        s_w, cs_w = seapy.roms.stretching(
            grid.vstretching, grid.theta_s, grid.theta_b, grid.hc,
            grid.n, w_grid=True)
        z_r[i, ...] = seapy.roms.depth(grid.vtransform, grid.h, grid.hc,
                                       s_w, cs_w, zeta=zeta[i, ...], w_grid=False)
        z_w[i, ...] = seapy.roms.depth(grid.vtransform, grid.h, grid.hc,
                                       s_w, cs_w, zeta=zeta[i, ...], w_grid=True)
        thick_rho = np.squeeze(z_w[i, 1:, :, :] - z_w[i, :-1, :, :])
        thick_u[i, ...] = seapy.model.rho2u(thick_rho)
        thick_v[i, ...] = seapy.model.rho2v(thick_rho)
    z_r[z_r > 50000] = np.ma.masked
    z_w[z_w > 50000] = np.ma.masked

    # Compute W (omega)
    Huon = u * thick_u * seapy.model.rho2u(grid.dn)
    Hvom = v * thick_v * seapy.model.rho2v(grid.dm)
    W = z_w * 0
    for k in range(grid.n):
        W[:, k + 1, :-2, :-2] = W[:, k, :-2, :-2] - \
            (Huon[:, k, 1:-1, 1:] - Huon[:, k, 1:-1, :-1]
             + Hvom[:, k, 1:, 1:-1] - Hvom[:, k, :-1, 1:-1])
    wrk = W[:, -1:, :, :] / (z_w[:, -1:, :, :] - z_w[:, 0:1, :, :])
    W[:, :-1, :, :] = W[:, :-1, :, :] - wrk * \
        (z_w[:, :-1, :, :] - z_w[:, 0:1, :, :])
    W[:, -1, :, :] = 0

    if scale:
        W *= grid.pn * grid.pm
    if work:
        return W, z_r, z_w, thick_u, thick_v
    else:
        return W


def wvelocity(grid, u, v, zeta=0):
    """
    Compute "true" vertical velocity

    Parameters
    ----------
    grid : seapy.model.grid,
      The grid to use for the calculations
    u : ndarray,
      The u-field in time
    v : ndarray,
      The v-field in time
    zeta : ndarray, optional,
      The zeta-field in time

    Returns
    -------
    w : ndarray,
      Vertical Velocity
    """
    grid = seapy.model.asgrid(grid)
    u = np.ma.array(u)
    v = np.ma.array(v)
    zeta = np.ma.array(zeta)

    # Check the sizes
    while u.ndim < 4:
        u = u[np.newaxis, ...]
    while v.ndim < 4:
        v = v[np.newaxis, ...]
    while zeta.ndim < 3:
        zeta = zeta[np.newaxis, ...]

    # Get omega
    W, z_r, z_w, thick_u, thick_v = omega(grid, u, v, zeta, scale=True,
                                          work=True)

    # Compute quasi-horizontal motions (Ui + Vj)*GRAD s(z)
    vert = z_r * 0
    # U-contribution
    wrk = u * (z_r[:, :, :, 1:] - z_r[:, :, :, :-1]) * \
        (grid.pm[:, 1:] - grid.pm[:, :-1])
    vert[:, :, :, 1:-1] = 0.25 * (wrk[:, :, :, :-1] + wrk[:, :, :, 1:])
    # V-contribution
    wrk = v * (z_r[:, :, 1:, :] - z_r[:, :, :-1, :]) * \
        (grid.pn[1:, :] - grid.pn[:-1, :])
    vert[:, :, 1:-1, :] += 0.25 * (wrk[:, :, :-1, :] + wrk[:, :, 1:, :])

    # Compute barotropic velocity [ERROR IN FORMULATION RIGHT NOW]
    wrk = np.zeros((vert.shape[0], vert.shape[2], vert.shape[3]))
    ubar = np.sum(u * thick_u, axis=1) / np.sum(thick_u, axis=1)
    vbar = np.sum(v * thick_v, axis=1) / np.sum(thick_v, axis=1)
    # wrk[:, 1:-1, 1:-1] = (ubar[:, 1:-1, :-1] - ubar[:, 1:-1, 1:] +
    #                       vbar[:, :-1, 1:-1] - vbar[:, 1:, 1:-1])

    # Shift vert from rho to w
    wvel = z_w * 0
    # First two layers
    slope = (z_r[:, 0, :, :] - z_w[:, 0, :, :]) / \
        (z_r[:, 1, :, :] - z_r[:, 0, :, :])
    wvel[:, 0, :, :] = 0.375 * (vert[:, 0, :, :] - slope *
                                (vert[:, 1, :, :] - vert[:, 0, :, :])) + \
        0.75 * vert[:, 0, :, :] - \
        0.125 * vert[:, 1, :, :]
    wvel[:, 1, :, :] = W[:, 1, :, :] + wrk + \
        0.375 * vert[:, 0, :, :] + \
        0.75 * vert[:, 1, :, :] - 0.125 * vert[:, 2, :, :]

    # Middle of the grid
    wvel[:, 2:-2, :, :] = W[:, 2:-2, :, :] + \
        wrk[:, np.newaxis, :, :] + \
        0.5625 * (vert[:, 1:-2, :, :] + vert[:, 2:-1, :, :]) - \
        0.0625 * (vert[:, :-3, :, :] + vert[:, 3:, :, :])

    # Upper two layers
    slope = (z_w[:, -1, :, :] - z_r[:, -1, :, :]) / \
        (z_r[:, -1, :, :] - z_r[:, -2, :, :])
    wvel[:, -1, :, :] = wrk + 0.375 * (vert[:, -1, :, :] + slope *
                                       (vert[:, -1, :, :] - vert[:, -2, :, :])) + \
        0.75 * vert[:, -1, :, :] - \
        0.0625 * vert[:, -2, :, :]
    wvel[:, -2, :, :] = W[:, -2, :, :] + 0.375 * vert[:, -1, :, :] + \
        wrk + 0.75 * vert[:, -2, :, :] - \
        0.125 * vert[:, -3, :, :]

    # No gradient at the boundaries
    wvel[:, :, 0, :] = wvel[:, :, 1, :]
    wvel[:, :, -2:, :] = wvel[:, :, -3:-2, :]
    wvel[:, :, :, 0] = wvel[:, :, :, 1]
    wvel[:, :, :, -2:] = wvel[:, :, :, -3:-2]

    return wvel


pass
