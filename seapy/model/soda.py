#!/usr/bin/env python
"""
  soda.py

  Functions for dealing with the soda model for importation into ROMS

  Written by Brian Powell on 07/24/15
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import numpy as np
from datetime import datetime
import netCDF4
from seapy.lib import default_epoch, chunker
from seapy.model.grid import asgrid
from seapy.roms import ncgen, num2date, date2num

_url = "http://apdrc.soest.hawaii.edu:80/dods/public_data/SODA/soda_pop2.2.4"
_maxrecs = 5


def load_history(filename,
                 start_time=datetime(1, 1, 1),
                 end_time=datetime(1, 1, 1),
                 grid=None,
                 epoch=default_epoch, url=_url, load_data=False):
    """
    Download soda data and save into local file

    Parameters
    ----------
    filename: string
        name of output file
    start_time: datetime
        starting date to load soda data
    end_time: datetime
        ending date for loading soda data
    grid: seapy.model.grid, optional
        if specified, only load SODA data that covers the grid
    epoch: datetime, optional
        reference time for new file
    url: string, optional
        URL to load SODA data from
    load_data: bool, optional
        If true (default) actually load the data. If false, it
        displays the information needed to load the data using ncks

    Returns
    -------
    None
    """
    # Load the grid
    grid = asgrid(grid)

    # Open the soda data
    soda = netCDF4.Dataset(url)

    # Figure out the time records that are required
    soda_time = num2date(soda, "time")

    time_list = np.where(np.logical_and(soda_time >= start_time,
                                        soda_time <= end_time))
    if not any(time_list):
        raise Exception("Cannot find valid times")

    # Get the latitude and longitude ranges
    minlat = np.min(grid.lat_rho) - 0.5
    maxlat = np.max(grid.lat_rho) + 0.5
    minlon = np.min(grid.lon_rho) - 0.5
    maxlon = np.max(grid.lon_rho) + 0.5
    soda_lon = soda.variables["lon"][:]
    soda_lat = soda.variables["lat"][:]

    # Ensure same convention
    if not grid.east():
        soda_lon[soda_lon > 180] -= 360

    latlist = np.where(np.logical_and(soda_lat >= minlat,
                                      soda_lat <= maxlat))
    lonlist = np.where(np.logical_and(soda_lon >= minlon,
                                      soda_lon <= maxlon))
    if not np.any(latlist) or not np.any(lonlist):
        raise Exception("Bounds not found")

    # Build the history file
    if load_data:
        his = ncgen.create_zlevel(filename, len(latlist[0]),
                                  len(lonlist[0]),
                                  len(soda.variables["lev"][:]), epoch,
                                  "soda history from " + url, dims=1)

        # Write out the data
        his.variables["lat"][:] = soda_lat[latlist]
        his.variables["lon"][:] = soda_lon[lonlist]
        his.variables["depth"][:] = soda.variables["lev"]
        his.variables["time"][:] = date2num(
            soda_time[time_list], his, 'time')
    # Loop over the variables
    sodavars = {"ssh": 3, "u": 4, "v": 4, "temp": 4, "salt": 4}
    hisvars = {"ssh": "zeta", "u": "u", "v": "v", "temp": "temp",
               "salt": "salt"}

    if not load_data:
        print("ncks -v {:s} -d time,{:d},{:d} -d lat,{:d},{:d} -d lon,{:d},{:d} {:s} {:s}".format(
            ",".join(sodavars.keys()),
            time_list[0][0], time_list[0][-1], latlist[0][0],
            latlist[0][-1], lonlist[0][0], lonlist[0][-1], _url, filename))
    else:
        for rn, recs in enumerate(chunker(time_list[0], _maxrecs)):
            print("{:s}-{:s}: ".format(soda_time[recs[0]].strftime("%m/%d/%Y"),
                                       soda_time[recs[-1]].strftime("%m/%d/%Y")),
                  end='', flush=True)
            for var in sodavars:
                print("{:s} ".format(var), end='', flush=True)
                hisrange = np.arange(
                    rn * _maxrecs, (rn * _maxrecs) + len(recs))
                if sodavars[var] == 3:
                    his.variables[hisvars[var]][hisrange, :, :] = \
                        np.ma.array(
                        soda.variables[var][recs, latlist[0], lonlist[0]]). \
                        filled(fill_value=9.99E10)
                else:
                    his.variables[hisvars[var]][hisrange, :, :, :] = \
                        soda.variables[var][recs, :, latlist[0],
                                            lonlist[0]].filled(fill_value=9.99E10)
            his.sync()
            print("", flush=True)
    pass
