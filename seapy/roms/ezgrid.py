#!/usr/bin/env python
"""
  ezgrid.py

  Functions for generating ROMS grid files

  Copyright (c)2020 University of Hawaii under the MIT-License.
"""
import numpy as np
from datetime import datetime
import netCDF4
import seapy
from collections import namedtuple


def create_grid(grid_file, lat, lon):
    """
    Create a new, basic grid. This will construct the ROMS grid netcdf file
    given a set of latitude and longitude coordinates for the rho-grid.
    The coordinates can be spaced however specified, and the grid will be
    created; however, this does not guarantee the grid will be valid.
    After the grid is created, the bathymetry and mask will need to be
    generated independently.

    Parameters
    ----------
    grid_file : string,
      Name of the grid file to create. NOTE: this file will be
      overwritten.
    lat : ndarray,
      latitude of the grid cells. The array must be the size
      of the desired grid.
    lon : ndarray,
      longitude of the grid cells. The array must be the size
      of the desired grid.

    Returns
    -------
    netCDF4 :
      The netcdf object of the new grid

    Examples
    --------
    To create a basic, evenly spaced grid:

    >>> lat = np.linspace(10,20,0.25)
    >>> lon = np.linspace(-100,-80,0.25)
    >>> lon, lat = np.meshgrid(lon, lat)
    >>> create_grid('mygrid.nc', lat, lon)

    To create more advanced grids, simply generate the
    2D arrays of lat and lon in the manner you want your cells
    and call create_grid:

    >>> create_grid('mygrid.nc', lat, lon)

    """

    # Put lat/lon into proper arrays
    lat = np.atleast_2d(lat)
    lon = np.atleast_2d(lon)

    if lat.shape != lon.shape:
        raise AttributeError("lat and lon shapes are not equal")

    # Calculate the angle between the points
    angle = np.zeros(lat.shape)
    angle[:, :-1] = seapy.earth_angle(lon[:, :-1],
                                      lat[:, :-1], lon[:, 1:], lat[:, 1:])
    angle[:, -1] = angle[:, -2]

    # Calculate distances/parameters
    f = 2.0 * 7.2921150e-5 * np.sin(lat * np.pi / 180.0)
    dx = np.zeros(f.shape)
    dy = np.zeros(f.shape)
    dx[:, 1:] = seapy.earth_distance(
        lon[:, 1:], lat[:, 1:], lon[:, :-1], lat[:, :-1])
    dy[1:, :] = seapy.earth_distance(
        lon[1:, :], lat[1:, :], lon[:-1, :], lat[:-1, :])
    dx[:, 0] = dx[:, 1]
    dy[0, :] = dy[1, :]
    pm = 1.0 / dx
    pn = 1.0 / dy
    dndx = np.zeros(dx.shape)
    dmde = np.zeros(dx.shape)
    dndx[:, 1:-1] = 0.5 * (dy[:, 2:] - dy[:, :-2])
    dmde[1:-1, :] = 0.5 * (dx[2:, :] - dx[:-2, :])
    xl = seapy.earth_distance(
        np.min(lon), np.mean(lat), np.max(lon), np.mean(lat))
    el = seapy.earth_distance(
        np.mean(lon), np.min(lat), np.mean(lon), np.max(lat))

    # Generate rho-grid coordinates
    x_rho = np.zeros(lat.shape)
    y_rho = np.zeros(lat.shape)
    x_rho[:, 1:] = seapy.earth_distance(
        lon[:, :-1], lat[:, :-1], lon[:, 1:], lat[:, 1:])
    x_rho = np.cumsum(x_rho, axis=1)
    y_rho[1:, :] = seapy.earth_distance(
        lon[:-1, :], lat[:-1, :], lon[1:, :], lat[1:, :])
    y_rho = np.cumsum(y_rho, axis=0)

    # Create u-grid
    lat_u = 0.5 * (lat[:, 1:] + lat[:, :-1])
    lon_u = 0.5 * (lon[:, 1:] + lon[:, 0:-1])
    x_u = np.zeros(lat_u.shape)
    y_u = np.zeros(lat_u.shape)
    x_u[:, 1:] = seapy.earth_distance(
        lon_u[:, :-1], lat_u[:, :-1], lon_u[:, 1:], lat_u[:, 1:])
    x_u = np.cumsum(x_u, axis=1)
    y_u[1:, :] = seapy.earth_distance(
        lon_u[:-1, :], lat_u[:-1, :], lon_u[1:, :], lat_u[1:, :])
    y_u = np.cumsum(y_u, axis=0)

    # Create v-grid
    lat_v = 0.5 * (lat[1:, :] + lat[0:-1, :])
    lon_v = 0.5 * (lon[1:, :] + lon[0:-1, :])
    x_v = np.zeros(lat_v.shape)
    y_v = np.zeros(lat_v.shape)
    x_v[:, 1:] = seapy.earth_distance(
        lon_v[:, :-1], lat_v[:, :-1], lon_v[:, 1:], lat_v[:, 1:])
    x_v = np.cumsum(x_v, axis=1)
    y_v[1:, :] = seapy.earth_distance(
        lon_v[:-1, :], lat_v[:-1, :], lon_v[1:, :], lat_v[1:, :])
    y_v = np.cumsum(y_v, axis=0)

    # Create psi-grid
    lat_psi = lat_v[:, :-1]
    lon_psi = lon_u[:-1, :]
    x_psi = np.zeros(lat_psi.shape)
    y_psi = np.zeros(lat_psi.shape)
    x_psi[:, 1:] = seapy.earth_distance(
        lon_psi[:, :-1], lat_psi[:, :-1], lon_psi[:, 1:], lat_psi[:, 1:])
    x_psi = np.cumsum(x_psi, axis=1)
    y_psi[1:, :] = seapy.earth_distance(
        lon_psi[:-1, :], lat_psi[:-1, :], lon_psi[1:, :], lat_psi[1:, :])
    y_psi = np.cumsum(y_psi, axis=0)

    # Create the new grid
    nc = seapy.roms.ncgen.create_grid(
        grid_file, lat.shape[0], lat.shape[1], clobber=True)

    nc.variables["xl"][:] = xl
    nc.variables["el"][:] = el
    nc.variables["spherical"][:] = 1
    nc.variables["f"][:] = f
    nc.variables["pm"][:] = pm
    nc.variables["pn"][:] = pn
    nc.variables["dndx"][:] = dndx
    nc.variables["dmde"][:] = dmde
    nc.variables["x_rho"][:] = x_rho
    nc.variables["y_rho"][:] = y_rho
    nc.variables["x_psi"][:] = x_psi
    nc.variables["y_psi"][:] = y_psi
    nc.variables["x_u"][:] = x_u
    nc.variables["y_u"][:] = y_u
    nc.variables["x_v"][:] = x_v
    nc.variables["y_v"][:] = y_v
    nc.variables["lat_rho"][:] = lat
    nc.variables["lon_rho"][:] = lon
    nc.variables["lat_psi"][:] = lat_psi
    nc.variables["lon_psi"][:] = lon_psi
    nc.variables["lat_u"][:] = lat_u
    nc.variables["lon_u"][:] = lon_u
    nc.variables["lat_v"][:] = lat_v
    nc.variables["lon_v"][:] = lon_v
    nc.variables["N"][:] = 1
    nc.variables["mask_rho"][:] = np.ones(lat.shape)
    nc.variables["mask_u"][:] = np.ones(lat_u.shape)
    nc.variables["mask_v"][:] = np.ones(lat_v.shape)
    nc.variables["mask_psi"][:] = np.ones(lat_psi.shape)
    nc.variables["angle"][:] = angle
    nc.variables["rdrag"][:] = np.ones(lon.shape) * 0.0003
    nc.variables["rdrag2"][:] = np.ones(lon.shape) * 0.003
    nc.variables["visc_factor"][:] = np.ones(lon.shape)
    nc.variables["diff_factor"][:] = np.ones(lon.shape)
    nc.variables["ZoBot"][:] = np.ones(lon.shape) * 0.02
    nc.sync()

    return nc


def calc_latlon(llcrnrlat, llcrnrlon, reseta, resxi=None, rotate=0):
    """
    Generate arrays for latitude and longitude for use in
    creating simple grids.

    NOTE: You can specify variational resolutions and rotations;
    however, this uses a simple algorithm to step along from the
    origin, and discrepencies are averaged together. THEREFORE,
    if you are not specifying inconsistent resolutions or
    rotations within single row or columns, you may get results
    slightly different than specified.

    Parameters
    ----------
    llcrnrlat : float,
      Latitude for the lower, left of the grid
    llcrnrlon : float,
      Longitude for the lower, left of the grid
    reseta : ndarray,
      Resolution in meters in the eta-direction of each grid cell.
      A 2D array that is the horizontal size of the grid (e.g.,
      resolution.shape = (100,70)) and the values stored are
      the desired resolution of the grid cell in m. Hence,
      a 100 x 70 array of ones would create a grid that is
      100 in the eta-direction, 70 in the xi-direction, and
      a constant resolution of 1km. If the lat and lon is
      specified as arrays, this is optional.
    resxi : ndarray, optional,
      Resolution in meters in the xi-direction of each grid cell.
      This is the same as reseta; however, it is optional, and
      is set to the same as reseta if not specified.
    rotate : float or ndarray,
      Amount to rotate the grid in degrees. If given as a scalar,
      the entire grid is rotated at the same angle; otherwise, the
      grid my have curvilinear shape. The angle is geometric
      (counter-clockwise).

    Returns
    -------
    lat, lon : ndarray
       Two arrays of size resolution containing the computed lat
       and lons

    Examples
    --------
    Create a grid of 1km resolution in both eta and xi,
    rotated toward the SE by 33 degrees, with the lower left
    point at 20N, 150E:

    >>> res = np.ones((100,70)) * 1000
    >>> lat, lon = calc_latlon(20, 150, res, rotate=-33)
    """

    # Set up the resolutions
    if resxi is None:
        resxi = reseta
    else:
        if resxi.shape != reseta.shape:
            raise AttributeError(
                "xi and eta resolutions must be same size array")
    # Since we are taking steps, average the values together as we
    # can't take the final step
    reseta[:-1, :] = 0.5 * (reseta[:-1, :] + reseta[1:, :])
    resxi[:, :-1] = 0.5 * (resxi[:, :-1] + resxi[:, 1:])

    # Set up the rotations
    if np.isscalar(rotate):
        rotate = np.ones(reseta.shape) * rotate
    rotate = np.radians(rotate)

    # Set up the lat and lon arrays
    lat = np.ones(reseta.shape) * np.nan
    lon = np.ones(reseta.shape) * np.nan
    lat[0, 0] = llcrnrlat
    lon[0, 0] = llcrnrlon
    pi2 = np.pi / 2.0

    # Loop over each row, and within the row, step along in the
    # xi-direction before moving to the next row
    for j in seapy.progressbar.progress(range(lat.shape[0] - 1)):
        for i in range(lat.shape[1] - 1):
            # Compute the local deltas
            dlat = seapy.earth_distance(
                lon[j, i], lat[j, i] - 0.5, lon[j, i], lat[j, i] + 0.5)
            dlon = seapy.earth_distance(
                lon[j, i] - 0.5, lat[j, i], lon[j, i] + 0.5, lat[j, i])

            # Compute how far to step in xi
            xdx = resxi[j, i + 1] * np.cos(rotate[j, i]) / dlon
            xdy = resxi[j, i + 1] * np.sin(rotate[j, i]) / dlat

            # Take the step
            lat[j, i + 1] = lat[j, i] + xdy
            lon[j, i + 1] = lon[j, i] + xdx

            # Compute how far to step in eta
            edx = reseta[j + 1, i] * np.cos(rotate[j, i] + pi2) / dlon
            edy = reseta[j + 1, i] * np.sin(rotate[j, i] + pi2) / dlat

            # Take the step
            lat[j + 1, i] = lat[j, i] + edy
            lon[j + 1, i] = lon[j, i] + edx

    lon[-1, -1] = lon[-1, -2] + xdx
    lat[-1, -1] = lat[-2, -1] + edy

    return lat, lon
