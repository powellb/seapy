#!/usr/bin/env python
"""
  initial.py
  
  ROMS initial conditions utilities

  Written by Brian Powell on 01/15/14
  Copyright (c)2014 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import seapy
import numpy as np
import netCDF4
import netcdftime

def from_roms(roms_file, ini_file, record=0, time=None, grid=None):
    """
    Given a ROMS history, average, or climatology file, generate 
    initial conditions on the same grid.
    
    Parameters
    ----------
    roms_file: string
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
    grid=seapy.model.asgrid(grid)
    ncroms = netCDF4.Dataset(roms_file)
    romstime = seapy.roms.get_timevar(ncroms)
    try:
        src_time=netcdftime.utime(ncroms.variables[romstime].units)
    except AttributeError:
        src_time=netcdftime.utime(seapy.roms.default_epoch)

    # Create the initial file and fill up the descriptive data
    ncini=seapy.roms.ncgen.create_ini(ini_file, 
             eta_rho=grid.eta_rho,xi_rho=grid.xi_rho,s_rho=grid.n,
             timebase=src_time.origin,title="generated from "+roms_file)
    ini_time=netcdftime.utime(ncini.variables[seapy.roms.get_timevar(ncini)].units)
    grid.to_netcdf(ncini)
    if time is None:
        time=src_time.num2date(ncroms.variables[romstime][record])
    ncini.variables["ocean_time"][:]=ini_time.date2num(time)

    # Fill up the initial state with the roms file data
    for var in seapy.roms.fields.keys():
        ncini.variables[var][0,:]=ncroms.variables[var][record,:]

    # Close up
    ncini.close()
    ncroms.close()

    pass


