#!/usr/bin/env python
"""
  hycom.py

  Functions for dealing with the HYCOM model for importation into ROMS

  Written by Brian Powell on 07/24/15
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import numpy as np
from datetime import datetime
import netCDF4
from seapy.lib import default_epoch, chunker
from seapy.model.grid import asgrid
from seapy.roms import ncgen, num2date

# _url = "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/2010"
_url = "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1"
_maxrecs = 5


def load_history(filename,
                 start_time=datetime(1, 1, 1),
                 end_time=datetime(1, 1, 1),
                 grid=None,
                 epoch=default_epoch, url=_url, load_data=False):
    """
    Download HYCOM data and save into local file

    Parameters
    ----------
    filename: string
        name of output file
    start_time: datetime
        starting date to load HYCOM data
    end_time: datetime
        ending date for loading HYCOM data
    grid: seapy.model.grid, optional
        if specified, only load SODA data that covers the grid
    epoch: datetime, optional
        reference time for new file
    url: string, optional
        URL to load SODA data from
    load_data: bool, optional
        If true actually load the data. If false (default), it
        displays the information needed to load the data using ncks

    Returns
    -------
    None
    """
    # Load the grid
    grid = asgrid(grid)

    # Open the HYCOM data
    hycom = netCDF4.Dataset(url)

    # Figure out the time records that are required
    hycom_time = num2date(hycom, "time")

    time_list = np.where(np.logical_and(hycom_time >= start_time,
                                        hycom_time <= end_time))
    if not np.any(time_list):
        raise Exception("Cannot find valid times")

    # Get the latitude and longitude ranges
    minlat = np.min(grid.lat_rho) - 0.5
    maxlat = np.max(grid.lat_rho) + 0.5
    minlon = np.min(grid.lon_rho) - 0.5
    maxlon = np.max(grid.lon_rho) + 0.5
    hycom_lon = hycom.variables["lon"][:]
    hycom_lat = hycom.variables["lat"][:]

    # Ensure same convention
    if not grid.east():
        hycom_lon[hycom_lon > 180] -= 360

    latlist = np.where(np.logical_and(hycom_lat >= minlat,
                                      hycom_lat <= maxlat))
    lonlist = np.where(np.logical_and(hycom_lon >= minlon,
                                      hycom_lon <= maxlon))
    if not np.any(latlist) or not np.any(lonlist):
        raise Exception("Bounds not found")

    # Build the history file
    if load_data:
        his = ncgen.create_zlevel(filename, len(latlist[0]),
                                  len(lonlist[0]),
                                  len(hycom.variables["depth"][:]), epoch,
                                  "HYCOM history from " + url, dims=1)

        # Write out the data
        his.variables["lat"][:] = hycom_lat[latlist]
        his.variables["lon"][:] = hycom_lon[lonlist]
        his.variables["depth"][:] = hycom.variables["depth"]
        his.variables["time"][:] = seapy.roms.date2num(
            hycom_time[time_list], his, 'time')

    # Loop over the variables
    hycomvars = {"surf_el": 3, "water_u": 4, "water_v": 4, "water_temp": 4,
                 "salinity": 4}
    hisvars = {"surf_el": "zeta", "water_u": "u", "water_v": "v",
               "water_temp": "temp", "salinity": "salt"}

    if not load_data:
        print("ncks -v {:s} -d time,{:d},{:d} -d lat,{:d},{:d} -d lon,{:d},{:d} {:s} {:s}".format(
            ",".join(hycomvars.keys()),
            time_list[0][0], time_list[0][-1], latlist[0][0],
            latlist[0][-1], lonlist[0][0], lonlist[0][-1], url, filename))
    else:
        for rn, recs in enumerate(chunker(time_list[0], _maxrecs)):
            print("{:s}-{:s}: ".format(hycom_time[recs[0]].strftime("%m/%d/%Y"),
                                       hycom_time[recs[-1]].strftime("%m/%d/%Y")),
                  end='', flush=True)
            for var in hycomvars:
                print("{:s} ".format(var), end='', flush=True)
                hisrange = np.arange(
                    rn * _maxrecs, (rn * _maxrecs) + len(recs))
                if hycomvars[var] == 3:
                    his.variables[hisvars[var]][hisrange, :, :] = \
                        hycom.variables[var][recs, latlist[0], lonlist[0]].filled(
                        fill_value=9.99E10)
                else:
                    his.variables[hisvars[var]][hisrange, :, :, :] = \
                        hycom.variables[var][recs, :, latlist[0],
                                             lonlist[0]].filled(fill_value=9.99E10)
            his.sync()
            print("", flush=True)
    pass
