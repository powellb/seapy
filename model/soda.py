#!/usr/bin/env python
"""
  soda.py
  
  Functions for dealing with the SODA model for importation into ROMS

  Written by Brian Powell on 05/24/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np
import datetime
import netCDF4 as nc
from progressbar import ProgressBar
import seapy.Null as Null
from seapy.roms.ncgen import create_zlevel

_url="http://apdrc.soest.hawaii.edu:80/dods/public_data/SODA/soda_pop2.2.4";

"""
    Given a ROMS grid, output file name, and time period, load the data
    from SODA into a new file.
"""
def load_history(grid=None, file_name=None, 
                 time=[datetime.datetime(1,1,1),datetime.datetime(1,1,1)],
                 epoch=datetime.datetime(1,1,1), url=_url):
    if grid is None:
        raise Exception("Missing input grid size to define the boundaries")
    if file_name is None:
        raise Exception("Missing output file name")
    
    # Open the SODA data
    soda = nc.Dataset(url,"r")
    
    # Figure out the time records that are required
    stime = nc.num2date(soda.variables["time"][:], soda.variables["time"].units)
    # stime = np.datetime64("0001-01-01") + np.array(soda.time-2, "timedelta64")
    time[0] = np.datetime64(time[0])
    time[1] = np.datetime64(time[1])
    tlist = np.nonzero(np.logical_and(stime>=time[0],stime<=time[1]))[0]
    if len(tlist) == 0:
        raise Exception("Cannot find valid times")

    # Get the latitude and longitude ranges
    minlat = np.min(grid.lat_rho.flatten())-0.25
    maxlat = np.max(grid.lat_rho.flatten())+0.25
    minlon = np.min(grid.lon_rho.flatten())-0.25
    maxlon = np.max(grid.lon_rho.flatten())+0.25

    # Check if we are using negative values versus positive
    slon = soda.variables["lon"][:]
    if minlon < 0:
        l=slon > 180
        slon[l] = slon[l] - 360
        
    latlist = np.nonzero( np.logical_and( soda.variables["lat"][:]>=minlat, \
                                          soda.variables["lat"][:]<=maxlat))[0]
    lonlist = np.nonzero( np.logical_and( slon>=minlon, slon<=maxlon))[0]
    if len(latlist) == 0 or len(lonlist) == 0:
        raise Exception("Bounds not found")

    # Build the history file
    his = create_zlevel(file_name, len(latlist), len(lonlist), 
                len(soda.variables["lev"][:])), epoch, 
                "SODA history from "+url)

    # Write out the data
    his.variables["lat"] = soda.variables["lat"][latlist]
    his.variables["lon"] = soda.variables["lon"][lonlist]
    his.variables["depth"] = soda.variables["lev"]
    his.variables["time"] = soda.variables["time"][tlist]
    
    # Loop over the variables
    dims = [3,4,4,4,4]
    sodavars = ("ssh", "u", "v", "temp", "salt")
    hisvars = ("zeta", "u", "v", "temp", "salt")
    p = ProgressBar()
    for i in p(range(len(sodavars))):
        if dims[i]==3:
            his.put(hisvars[i], 
             soda.variables[sodavars[i]][tlist, latlist, lonlist].filled(fill_value=9.99E10))
        else:
            his.put(hisvars[i], 
             soda.variables[sodavars[i]][tlist, :, latlist, lonlist].filled(fill_value=9.99E10))
    pass


