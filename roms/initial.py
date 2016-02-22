#!/usr/bin/env python
"""
  initial.py

  ROMS initial conditions utilities

  Written by Brian Powell on 01/15/14
  Copyright (c)2016 University of Hawaii under the BSD-License.
"""


import seapy
import numpy as np
import netCDF4


def from_roms(roms_file, ini_file, record=0, time=None, grid=None):
    """
    Given a ROMS history, average, or climatology file, generate
    initial conditions on the same grid.

    Parameters
    ----------
    roms_file: string or list
        Input ROMS source (history, average, climatology file)
    ini_file: string
        Input name for output initial condition file
    record: int
        Input index to use as initial condition
    time: datetime optional
        Input datetime to use for the initial condition (default to record time)
    grid: seapy.model.grid or string, optional
        Input ROMS grid: specify the grid if loaded or filename to load

    Returns
    -------
    None

    """
    # Load the grid
    if grid is None:
        grid = seapy.model.asgrid(roms_file)
    else:
        grid = seapy.model.asgrid(grid)
    ncroms = seapy.netcdf4(roms_file)
    src_ref, romstime = seapy.roms.get_reftime(ncroms)

    # Create the initial file and fill up the descriptive data
    ncini = seapy.roms.ncgen.create_ini(ini_file,
                                        eta_rho=grid.eta_rho, xi_rho=grid.xi_rho, s_rho=grid.n,
                                        reftime=src_ref, title="generated from " + roms_file)
    grid.to_netcdf(ncini)
    if time is None:
        time = netCDF4.num2date(ncroms.variables[romstime][record],
                                ncroms.variables[romstime].units)
    ncini.variables["ocean_time"][:] = netCDF4.date2num(time,
                                                        ncini.variables["ocean_time"].units)

    # Fill up the initial state with the roms file data
    for var in seapy.roms.fields:
        ncini.variables[var][0, :] = ncroms.variables[var][record, :]

    # Close up
    ncini.close()
    ncroms.close()
    pass


