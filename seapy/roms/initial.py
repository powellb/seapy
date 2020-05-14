#!/usr/bin/env python
"""
  initial.py

  ROMS initial conditions utilities

  Written by Brian Powell on 01/15/14
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import seapy
import numpy as np
import netCDF4


def from_roms(roms_file, ini_file, record=0, time=None, grid=None,
              clobber=False, cdl=None):
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
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.

    Returns
    -------
    None

    """
    # Load the grid
    if grid is None:
        grid = seapy.model.asgrid(roms_file)
    else:
        grid = seapy.model.asgrid(grid)
    ncroms = seapy.netcdf(roms_file)
    src_ref, romstime = seapy.roms.get_reftime(ncroms)

    # Create the initial file and fill up the descriptive data
    ncini = seapy.roms.ncgen.create_ini(ini_file,
                                        eta_rho=grid.eta_rho, xi_rho=grid.xi_rho, s_rho=grid.n,
                                        reftime=src_ref, clobber=clobber, cdl=cdl,
                                        title="generated from " + roms_file)
    grid.to_netcdf(ncini)
    if time is None:
        time = seapy.roms.num2date(ncroms, romstime, record)
    ncini.variables["ocean_time"][:] = seapy.roms.date2num(
        time, ncini, "ocean_time")

    # Fill up the initial state with the roms file data
    for var in seapy.roms.fields:
        if var in ncini.variables and var in ncroms.variables:
            ncini.variables[var][0, :] = ncroms.variables[var][record, :]

    # Close up
    ncini.close()
    ncroms.close()
    pass
