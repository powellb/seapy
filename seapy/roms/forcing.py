#!/usr/bin/env python
"""
  forcing.py

  Functions for generating ROMS forcing files (atmosphere, tides,
  rivers, etc.)

  Written by Brian Powell on 02/09/16
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""
import numpy as np
from datetime import datetime
import netCDF4
import seapy
from collections import namedtuple

# Define a named tuple to store raw data for the gridder
forcing_data = namedtuple('forcing_data', 'field ratio offset')
fields = ("Tair", "Qair", "Pair", "rain",
          "Uwind", "Vwind", "lwrad_down", "swrad")

gfs_url = "http://oos.soest.hawaii.edu/thredds/dodsC/hioos/model/atm/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd"

gfs_map = {
    "pad": 1.0,
    "frc_lat": "latitude",
    "frc_lon": "longitude",
    "frc_time": "time",
    "Tair": forcing_data("tmp2m", 1, -273.15),
    "Pair": forcing_data("prmslmsl", 0.01, 0),
    "Qair": forcing_data("rh2m", 1, 0),
    "rain": forcing_data("pratesfc", 1, 0),
    "Uwind": forcing_data("ugrd10m", 1, 0),
    "Vwind": forcing_data("vgrd10m", 1, 0),
    "lwrad_down": forcing_data("dlwrfsfc", 1, 0),
    "swrad": forcing_data("dswrfsfc", 1, 0)
}

ncep_map = {
    "pad": 2.0,
    "frc_lat": "lat",
    "frc_lon": "lon",
    "frc_time": "time",
    "Tair": forcing_data("air", 1, -273.15),
    "Pair": forcing_data("slp", 0.01, 0),
    "Qair": forcing_data("rhum", 1, 0),
    "rain": forcing_data("prate", 1, 0),
    "Uwind": forcing_data("uwnd", 1, 0),
    "Vwind": forcing_data("vwnd", 1, 0),
    "lwrad_down": forcing_data("dlwrf", 1, 0),
    "swrad": forcing_data("dswrf", 1, 0)
}


def gen_bulk_forcing(infile, fields, outfile, grid, start_time, end_time,
                     epoch=seapy.default_epoch, clobber=False, cdl=None):
    """
    Given a source file (or URL), a dictionary that defines the
    source fields mapped to the ROMS fields, then it will generate
    a new bulk forcing file for ROMS.

    Parameters
    ----------
    infile: string,
      The source file (or URL) to load the data from
    fields: dict,
        A dictionary of the fields to load and process. The dictionary
        is composed of:
         "frc_lat":STRING name of latitude field in forcing file
         "frc_lon":STRING name of longitude field in forcing file
         "frc_time":STRING name of time field in forcing file
         "frc_time_units":STRING optional, supply units of frc time
                                field in forcing file
         keys of ROMS bulk forcing field names (Tair, Pair, Qair,
         rain, Uwind, Vwind, lwrad_down, swrad) each with an
         array of values of a named tuple (forcing_data) with the
         following fields:
            field: STRING value of the forcing field to use
            ratio: FLOAT value to multiply with the source data
            offset: FLOAT value to add to the source data
    outfile: string,
        Name of the output file to create
    grid: seapy.model.grid or string,
        Grid to use for selecting spatial region
    start_time: datetime,
        Starting time of data to process
    end_time: datetime,
        Ending time of data to process
    epoch: datetime,
        Epoch to use for ROMS times
    clobber: bool optional,
        Delete existing file or not, default False
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.

    Returns
    -------
    None: Generates an output file of bulk forcings

    Examples
    --------
    To generate GFS forcing for the grid "mygrid.nc" for the year
    2014, then use the standard GFS map definitions (and even
    the built-in GFS archive url):

    >>> seapy.roms.forcing.gen_bulk_forcing(seapy.roms.forcing.gfs_url,
            seapy.roms.forcing.gfs_map, 'my_forcing.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2014,1,1))

    NCEP reanalysis is trickier as the files are broken up by
    variable type; hence, a separate file will be created for
    each variable. We can use the wildcards though to put together
    multiple time period (e.g., 2014 through 2015).

    >>> seapy.roms.forcing.gen_bulk_forcing("uwnd.10m.*nc",
            seapy.roms.forcing.ncep_map, 'ncep_frc.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2015,12,31)), clobber=True)
    >>> seapy.roms.forcing.gen_bulk_forcing("vwnd.10m.*nc",
            seapy.roms.forcing.ncep_map, 'ncep_frc.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2015,12,31)), clobber=False)
    >>> seapy.roms.forcing.gen_bulk_forcing("air.2m.*nc",
            seapy.roms.forcing.ncep_map, 'ncep_frc.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2015,12,31)), clobber=False)
    >>> seapy.roms.forcing.gen_bulk_forcing("dlwrf.sfc.*nc",
            seapy.roms.forcing.ncep_map, 'ncep_frc.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2015,12,31)), clobber=False)
    >>> seapy.roms.forcing.gen_bulk_forcing("dswrf.sfc.*nc",
            seapy.roms.forcing.ncep_map, 'ncep_frc.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2015,12,31)), clobber=False)
    >>> seapy.roms.forcing.gen_bulk_forcing("prate.sfc.*nc",
            seapy.roms.forcing.ncep_map, 'ncep_frc.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2015,12,31)), clobber=False)
    >>> seapy.roms.forcing.gen_bulk_forcing("rhum.sig995.*nc",
            seapy.roms.forcing.ncep_map, 'ncep_frc_rhum_slp.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2015,12,31)), clobber=True)
    >>> seapy.roms.forcing.gen_bulk_forcing("slp.*nc",
            seapy.roms.forcing.ncep_map, 'ncep_frc_rhum_slp.nc',
            'mygrid.nc', datetime.datetime(2014,1,1),
            datetime.datetime(2015,12,31)), clobber=False)

    Two forcing files, 'ncep_frc.nc' and 'ncep_frc_rhum_slp.nc', are
    generated for use with ROMS. NOTE: You will have to use 'ncks'
    to eliminate the empty forcing fields between the two files
    to prevent ROMS from loading them.
    """
    # Load the grid
    grid = seapy.model.asgrid(grid)

    # Open the Forcing data
    forcing = seapy.netcdf(infile)

    # Gather the information about the forcing
    if 'frc_time_units' in fields:
        frc_time = netCDF4.num2date(forcing.variables[fields['frc_time']][:],
                                    fields['frc_time_units'])
    else:
        frc_time = seapy.roms.num2date(forcing, fields['frc_time'])

    # Figure out the time records that are required
    time_list = np.where(np.logical_and(frc_time >= start_time,
                                        frc_time <= end_time))[0]
    if not np.any(time_list):
        raise Exception("Cannot find valid times")

    # Get the latitude and longitude ranges
    minlat = np.floor(np.min(grid.lat_rho)) - fields['pad']
    maxlat = np.ceil(np.max(grid.lat_rho)) + fields['pad']
    minlon = np.floor(np.min(grid.lon_rho)) - fields['pad']
    maxlon = np.ceil(np.max(grid.lon_rho)) + fields['pad']
    frc_lon = forcing.variables[fields['frc_lon']][:]
    frc_lat = forcing.variables[fields['frc_lat']][:]
    # Make the forcing lat/lon on 2D grid
    if frc_lon.ndim == 3:
        frc_lon = np.squeeze(frc_lon[0, :, :])
        frc_lat = np.squeeze(frc_lat[0, :, :])
    elif frc_lon.ndim == 1:
        frc_lon, frc_lat = np.meshgrid(frc_lon, frc_lat)

    # Find the values in our region
    if not grid.east():
        frc_lon[frc_lon > 180] -= 360
    region_list = np.where(np.logical_and.reduce((
        frc_lon <= maxlon,
        frc_lon >= minlon,
        frc_lat <= maxlat,
        frc_lat >= minlat)))
    if not np.any(region_list):
        raise Exception("Cannot find valid region")

    eta_list = np.s_[np.min(region_list[0]):np.max(region_list[0]) + 1]
    xi_list = np.s_[np.min(region_list[1]):np.max(region_list[1]) + 1]
    frc_lon = frc_lon[eta_list, xi_list]
    frc_lat = frc_lat[eta_list, xi_list]

    # Create the output file
    out = seapy.roms.ncgen.create_frc_bulk(outfile, lat=frc_lat.shape[0],
                                           lon=frc_lon.shape[1], reftime=epoch,
                                           clobber=clobber, cdl=cdl)
    out.variables['frc_time'][:] = seapy.roms.date2num(
        frc_time[time_list], out, 'frc_time')
    out.variables['lat'][:] = frc_lat
    out.variables['lon'][:] = frc_lon

    # Loop over the fields and fill out the output file
    for f in seapy.progressbar.progress(list(set(fields.keys()) & (out.variables.keys()))):
        if hasattr(fields[f], 'field') and fields[f].field in forcing.variables:
            out.variables[f][:] = \
                forcing.variables[fields[f].field][time_list, eta_list, xi_list] * \
                fields[f].ratio + fields[f].offset
            out.sync()

    out.close()


def gen_direct_forcing(his_file, frc_file, cdl=None):
    """
    Generate a direct forcing file from a history (or other ROMS output) file. It requires
    that sustr, svstr, shflux, and ssflux (or swflux) with salt be available. This will
    generate a forcing file that contains: sustr, svstr, swflux, and ssflux.

    Parameters
    ----------
    his_file: string,
      The ROMS history (or other) file(s) (can use wildcards) that contains the fields to
      make forcing from
    frc_file: string,
      The output forcing file
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.

    Returns
    -------
    None: Generates an output file of bulk forcings
    """
    import os

    infile = seapy.netcdf(his_file)
    ref, _ = seapy.roms.get_reftime(infile)

    # Create the output file
    nc = seapy.roms.ncgen.create_frc_direct(frc_file,
                                            eta_rho=infile.dimensions[
                                                'eta_rho'].size,
                                            xi_rho=infile.dimensions[
                                                'xi_rho'].size,
                                            reftime=ref,
                                            clobber=True,
                                            title="Forcing from " +
                                            os.path.basename(his_file),
                                            cdl=cdl)

    # Copy the data over
    time = seapy.roms.num2date(infile, 'ocean_time')
    nc.variables['frc_time'][:] = seapy.roms.date2num(time, nc, 'frc_time')
    for x in seapy.progressbar.progress(seapy.chunker(range(len(time)), 1000)):
        nc.variables['SSS'][x, :, :] = seapy.convolve_mask(
            infile.variables['salt'][x, -1, :, :], copy=False)
        if 'EminusP' in infile.variables:
            nc.variables['swflux'][x, :, :] = seapy.convolve_mask(
                infile.variables['EminusP'][x, :, :], copy=False) * 86400
        elif 'swflux' in infile.variables:
            nc.variables['swflux'][x, :, :] = seapy.convolve_mask(
                infile.variables['swflux'][x, :, :], copy=False)
        else:
            nc.variables['swflux'][x, :, :] = seapy.convolve_mask(
                infile.variables['ssflux'][x, :, :]
                / nc.variables['SSS'][x, :, :], copy=False)

        nc.sync()
        for f in ("sustr", "svstr", "shflux", "swrad"):
            if f in infile.variables:
                nc.variables[f][x, :, :] = seapy.convolve_mask(
                    infile.variables[f][x, :, :], copy=False)
                nc.sync()

    for f in ("lat_rho", "lat_u", "lat_v", "lon_rho", "lon_u", "lon_v"):
        if f in infile.variables:
            nc.variables[f][:] = infile.variables[f][:]
    nc.close()
