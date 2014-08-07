#!/usr/bin/env python
"""
  boundary.py
  
  ROMS boundary utilities

  Written by Brian Powell on 01/15/14
  Copyright (c)2014 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import seapy
import os
import numpy as np
import netCDF4
import netcdftime

def from_roms(roms_file, bry_file, grid=None, records=None):
    """
    from_roms(roms_file, bry_file, grid=None, recors=None)
    
    Given a ROMS history, average, or climatology file, generate 
    boundary conditions on the same grid.
    
    Parameters
    ----------
    roms_file : string ROMS source (history, average, romsatology file)
    bry_file : string output boundary file
    [grid] : seapy.model.grid or string for ROMS grid
    [records] : array of the record indices to put into the boundary
    
    Returns
    -------
    None
    
    """
    if grid != None:
        if isinstance(grid,basestring):
            grid = seapy.model.grid(grid, z=True)
    else:
        # If we weren't given a grid, try to construct from the romsate file
        grid = seapy.model.grid(roms_file, z=True)
        
    ncroms = netCDF4.Dataset(roms_file)
    time = seapy.roms.get_timevar(ncroms)
    records = np.arange(0, len(ncroms.variables[time][:])) \
             if records == None else records
    src_time=netcdftime.utime(ncroms.variables[time].units)

    # Create the boundary file and fill up the descriptive data
    if os.path.isfile(bry_file):
        ncbry=netCDF4.Dataset(bry_file,"a")
    else:
        ncbry=seapy.roms.ncgen.create_bry(bry_file, 
             eta_rho=grid.ln,xi_rho=grid.lm,N=grid.n,
             timebase=src_time.origin,title="generated from "+roms_file)
    brytime = seapy.roms.get_timevar(ncbry)
    bry_time=netcdftime.utime(ncbry.variables[brytime].units)
    ncbry.variables["lat_rho"][:]=grid.lat_rho
    ncbry.variables["lon_rho"][:]=grid.lon_rho
    ncbry.variables["lat_u"][:]=grid.lat_u
    ncbry.variables["lon_u"][:]=grid.lon_u
    ncbry.variables["lat_v"][:]=grid.lat_v
    ncbry.variables["lon_v"][:]=grid.lon_v
    ncbry.variables["Vtransform"][:]=grid.vtransform
    ncbry.variables["Vstretching"][:]=grid.vstretching
    ncbry.variables["theta_s"][:]=grid.theta_s
    ncbry.variables["theta_b"][:]=grid.theta_b
    ncbry.variables["hc"][:]=grid.hc
    ncbry.variables["Tcline"][:]=grid.tcline
    ncbry.variables["s_rho"][:]=grid.s_rho
    ncbry.variables["Cs_r"][:]=grid.cs_r
    ncbry.variables["h"][:]=grid.h
    ncbry.variables["bry_time"][:]=bry_time.date2num(
        src_time.num2date(ncroms.variables[time][records]))

    # Go over the variables for each side
    sides={"north":[-2,"-1"], "south":[-2,"0"], "east":[-1,"-1"], "west":[-1,"0"]}
    index=["records",":",":",":"]
    for var in seapy.roms.fields.keys():
        if var in ncroms.variables:
            for side in sides: 
                outvar=var+"_"+side
                ndim=seapy.roms.fields[var]["dims"]
                if ndim==3:
                    myindex=list(index)
                elif ndim==2:
                    myindex=list(index[:-1])
                myindex[sides[side][0]]=sides[side][1]
                indexstr=",".join(myindex)
                ncbry.variables[outvar][:]=eval("ncroms.variables[var]["+indexstr+"]")
    ncbry.close()
    
    pass

