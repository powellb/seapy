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
    from_roms(roms_file, ini_file, grid=None, recors=None)
    
    Given a ROMS history, average, or climatology file, generate 
    initial conditions on the same grid.
    
    Parameters
    ----------
    roms_file : string ROMS source (history, average, climatology file)
    ini_file : string name for output initial condition file
    record :  record index to use as initial condition
    [time] : datetime to use for the initial condition (default to record time)
    [grid] : seapy.model.grid or string for ROMS grid
    
    Returns
    -------
    None
    
    """
    # Load the grid
    if grid != None:
        if isinstance(grid,basestring):
            grid = seapy.model.grid(grid)
    else:
        # If we weren't given a grid, try to construct from the climate file
        grid = seapy.model.grid(roms_file)
        
    ncroms = netCDF4.Dataset(roms_file)
    romstime = seapy.roms.get_timevar(ncroms)
    src_time=netcdftime.utime(ncroms.variables[romstime].units)

    # Create the initial file and fill up the descriptive data
    ncini=seapy.roms.ncgen.create_ini(ini_file, 
             eta_rho=grid.ln,xi_rho=grid.lm,N=grid.n,
             timebase=src_time.origin,title="generated from "+roms_file)
    ini_time=netcdftime.utime(ncini.variables[seapy.roms.get_timevar(ncini)].units)
    ncini.variables["lat_rho"][:]=grid.lat_rho
    ncini.variables["lon_rho"][:]=grid.lon_rho
    ncini.variables["lat_u"][:]=grid.lat_u
    ncini.variables["lon_u"][:]=grid.lon_u
    ncini.variables["lat_v"][:]=grid.lat_v
    ncini.variables["lon_v"][:]=grid.lon_v
    ncini.variables["Vtransform"][:]=grid.vtransform
    ncini.variables["Vstretching"][:]=grid.vstretching
    ncini.variables["theta_s"][:]=grid.theta_s
    ncini.variables["theta_b"][:]=grid.theta_b
    ncini.variables["hc"][:]=grid.hc
    ncini.variables["Tcline"][:]=grid.tcline
    ncini.variables["s_rho"][:]=grid.s_rho
    ncini.variables["Cs_r"][:]=grid.cs_r
    ncini.variables["h"][:]=grid.h
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


