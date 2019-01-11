#!/usr/bin/env python
"""
   psource.py

   Routines for dealing with point source files in ROMS

   Author: Brian Powell <powellb@hawaii.edu>
   Copyright (c) 2019, Brian Powell, all rights reserved.
   Created: 10 January 2019

"""

import seapy.roms.ncgen
import numpy as np
import urllib.request as urllib
import json
import datetime

discharge_url = 'https://waterdata.usgs.gov/nwisweb/get_ratings?file_type=exsa&site_no='
river_url = 'https://waterservices.usgs.gov/nwis/iv/?format=json'


def make_psource(filename, rivers, cdl=None):
    """
    Construct a point source file with all of the point sources configured.

    Input
    -----
    filename : string
            Output filename
    s_rho : int,
            Number of s-levels in the point source file (should match the grid)
    rivers : array,
            array of dictionaries of rivers with the following information that
            defines a point source:
               The x and y values on the grid for where the point source is,
               the direction of the point source (0 for xi or 1 for eta),
               an identification number (any choice),
               flag (1=temp, 2=salt, 3=temp+salt, 4=temp+salt+sed, 5=temp+salt+sed+bio),
               and an array of values for the vertical shape (a value for each s-level)
               that sum to 1.
            { "x":grid_x, "y":grid_y, "direction":0 or 1, "flag":1,2,3,4,or 5, "id":value, "vshape":[vals] }
    cdl : string, optional,
            Name of CDL file to use

    Output
    ------
    nc : netCDF4 id
       If successful, returns the netcdf id for writing
    """
    river = np.asarray(river)
    s_rho = len(river[0]['vshape'])

    # Create an empty, new file
    nc = seapy.roms.ncgen.create_psource(
        fname, nriver=len(river), s_rho=s_rho, clobber=True, cdl=cdl)

    # Go through each river and set up the basics
    for i, r in enumerate(river):
        nc.variables['river'][i] = int(r['id'])
        nc.variables['river_Xposition'][i] = int(r['x'])
        nc.variables['river_Yposition'][i] = int(r['y'])
        nc.variables['river_direction'][i] = int(r['direction'])
        nc.variables['river_flag'][i] = int(r['flag'])
        try:
            vshape = np.asarray(r['vshape'][:, np.newaxis])
            nc.variables['river_Vshape'][:, i] = vshape
        except Exeption as err:
            print(f"Error in VSHAPE: {err}")
            print("Using default shape")
            vshape = np.ones((s_rho, 1)) / s_rho
            nc.variables['river_Vshape'][:, i] = vshape

    return nc


def stage2discharge(gage, usgs_id):
    '''
    Function to convert stage data to discharge using usgs reference table

    Input
    -------
    gage : masked array,
                array of stage data to be converted
    usgs_id : int,
                8 digit identifier for river on usgs

    Output
    ---------
    flux : masked array,
                discharge data in cubic feet
    '''
    import urllib.request as urllib

    # Load lookup table
    url = discharge_url + usgs_id
    stage = []
    discharge = []
    header = 0
    for line in urllib.urlopen(url):
        if not line.startswith(b'#'):
            if header < 2:
                header += 1
            else:
                a = line.split(b'\t')
                stage.append(float(line.split(b'\t')[0]))
                discharge.append(float(line.split(b'\t')[2]))
    stage = np.asarray(stage).astype(float)
    discharge = np.asarray(discharge).astype(float)

    # Convert data
    flux = np.ma.masked_array(np.zeros(len(gage)), mask=gage.mask)
    for i, f in enumerate(gage):
        if f:
            flux[i] = discharge[(np.abs(stage - f)).argmin()]
    return flux


def get_river_transport(usgs_id, times=1, source='discharge'):
    '''
    Function to get flux data from usgs and convert to cubic meters

    Input
    -------
    usgs_id : int,
                8 digit identifier for river on usgs
    times : int/list of datetimes, default = 1
                If int supplied, function will fetch last n days of data available
                If list of datetimes, function will fetch data between the
                start and end value
    source : string,
                  set to 'discharge' or 'stage' to access corresponding
                  parameter. Stage data is then converted to discharge
                  using stage2discharge function.
    Output
    ---------
    dates : list,
                list of datetimes associated with the flux data
    flux : array,
                flux data in cubic meters per second
    '''
    # Build url
    siteurl = '&sites=' + usgs_id
    if source == 'discharge':
        sourceurl = '&parameterCd=00060'
    elif source == 'stage':
        sourceurl = '&parameterCd=00065'
    else:
        print('Incorrect source type specified')
        return None
    if isinstance(times, int):
        timeurl = '&period=P%sD' % str(times)
    else:
        timeurl = '&startDT=%s&endDT=%s' % (times[0].strftime(
            '%Y-%m-%d'), times[1].strftime('%Y-%m-%d'))
    url = river_url + siteurl + sourceurl + timeurl

    # Access url
    print('Accessing ' + url)
    x = json.loads(urllib.urlopen(url).read().decode('utf-8'))
    dates = []
    flux = []
    if x['value']['timeSeries']:
        for l in x['value']['timeSeries'][0]['values'][0]['value']:
            dates.append(datetime.datetime.strptime(
                l['dateTime'][:16], '%Y-%m-%dT%H:%M'))
            flux.append(l['value'])
        flux = np.ma.masked_values(
            np.ma.masked_array(flux).astype(np.float), -999999)
        if flux.all() is np.ma.masked:
            print('No valid data found')
            return None
        else:
            if source == 'stage':
                flux = stage2discharge(flux, usgs_id)
            return np.array(dates)[~np.ma.getmaskarray(flux)], \
                flux.compressed() / 3.28**3
    else:
        print('Cannot process file from url')
        return None
