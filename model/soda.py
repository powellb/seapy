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
from seapy.roms.ncgen import create_zlevel

_url="http://apdrc.soest.hawaii.edu:80/dods/public_data/SODA/soda_pop2.2.4";

def load_history(filename, 
                 time=[datetime.datetime(1,1,1),datetime.datetime(1,1,1)],
                 grid=None,
                 epoch=datetime.datetime(1,1,1), url=_url):
    """
    Download soda data and save into local file

    Parameters
    ----------
    filename: string
        name of output file
    time: array of datetime
        start and end times to load SODA data
    grid: seapy.model.grid, optional
        if specified, only load SODA data that covers the grid
    epoch: datetime, optional
        timebase for new file
    url: string, optional
        URL to load SODA data from
        
    Returns
    -------
    None
    """
    if grid is None:
        raise Exception("Missing input grid size to define the boundaries")
    if filename is None:
        raise Exception("Missing output file name")
    
    # Open the SODA data
    soda = nc.Dataset(url,"r")
    
    # Figure out the time records that are required
    stime = nc.num2date(soda.variables["time"][:], soda.variables["time"].units)
    # stime = np.datetime64("0001-01-01") + np.array(soda.time-2, "timedelta64")
    time[0] = np.datetime64(time[0])
    time[1] = np.datetime64(time[1])
    tlist = np.nonzero(np.logical_and(stime>=time[0],stime<=time[1]))[0]
    if not tlist:
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
    if not latlist or not lonlist:
        raise Exception("Bounds not found")

    # Build the history file
    his = create_zlevel(filename, len(latlist), len(lonlist), 
                len(soda.variables["lev"][:]), epoch, 
                "SODA history from "+url)

    # Write out the data
    his.variables["lat"] = soda.variables["lat"][latlist]
    his.variables["lon"] = soda.variables["lon"][lonlist]
    his.variables["depth"] = soda.variables["lev"]
    his.variables["time"] = soda.variables["time"][tlist]
    
    # Loop over the variables
    sodavars = {"ssh":3, "u":4, "v":4, "temp":4, "salt":4}
    hisvars = {"ssh":"zeta", "u":"u", "v":"v", "temp":"temp", "salt":"salt"}
    for var in seapy.progressbar.progress(sodavars):
        if sodavars[var]==3:
            his.put(hisvars[var], 
             soda.variables[var][tlist, latlist, lonlist].filled(fill_value=9.99E10))
        else:
            his.put(hisvars[var], 
             soda.variables[var][tlist, :, latlist, lonlist].filled(fill_value=9.99E10))
    pass


