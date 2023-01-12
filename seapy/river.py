#!/usr/bin/env python
"""
  river.py

  Python functions to gather USGS stream data

  Written by Dale Partridge on 08/21/17
  Copyright (c)2017 University of Hawaii under the BSD-License.


    River Parameters
    ----------------------
    river : named tuple containing:
        usgs_id : int, 8 digit identifier for river
        direction : int, direction in grid- 0 for xi direction, 1 for eta
        flow_direction : int, Direction of flow, 1 for north/east, -1 for south/west
        flag : int, parameters to include in roms,
                1=temp, 2=salt, 3=temp+salt, 4=temp+salt+sed, 5=temp+salt+sed+bio
        x : int, xi rho-grid position of river
        y : int, eta rho-grid position of river
        vshape : array, vertical profile of river transport (should sum to 1)
        multiplier : float, value used to scale river transport to account for
                     other sources
        source : string, used to denote accessing 'stage' or 'discharge'
                 data from usgs site
        defaults : named tuple containing the following:
            transport : float, default river transport (m^3/s)
            temp : float, default temperature (Celsius)
            salt : float, default salinity (PSU)
            turb : float, default turbidity (NTU)

    For each domain build a dictionary of river tuples, e.g. example :-
    manoa_defaults = defaults(0.0236,20,1,0.8)
    manoa_vshape = [1/10]*10
    hiomsg = {'palolo-manoa' : river(16247100, 3, 1, -1, 112, 39, manoa_vshape, 1.3, 'stage' manoa_defaults)}

"""
import seapy
import numpy as np
import datetime
from collections import namedtuple
import urllib.request as urllib
import json

'''
'''
defaults = namedtuple('defaults', 'transport temp salt turb')
river_info = namedtuple(
    'river_info', 'usgs_id direction flow_direction flag x y vshape multiplier source defaults')
discharge_url = 'https://waterdata.usgs.gov/nwisweb/get_ratings?file_type=exsa&site_no='
river_url = 'https://waterservices.usgs.gov/nwis/iv/?format=json'


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
                flux data in cubic meters
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
            np.ma.masked_array(flux).astype(np.float32), -999999)
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


def get_turbidity(fname, rivers):
    '''
    Function to get turbidity from usgs river sensor

    Input
    -------
    fname : string
            River file name
    '''
    import netCDF4
    nc = netCDF4.Dataset(fname, 'a')
    time = list(seapy.day2date(nc.variables['river_time'][:]))
    for i, rd in enumerate(sorted(rivers)):
        r = rivers[rd]
        # Build url
        siteurl = '&sites=' + r.usgs_id
        sourceurl = '&parameterCd=63680'
        timeurl = f"&startDT={time[0].strftime('%Y-%m-%d')}&endDT={time[-1].strftime('%Y-%m-%d')}"
        url = river_url + siteurl + sourceurl + timeurl

        # Access url
        print(f'Adding turbidity to {rd} from {url}')
        x = json.loads(urllib.urlopen(url).read().decode('utf-8'))
        dates = []
        turb = []
        if x['value']['timeSeries']:
            for l in x['value']['timeSeries'][0]['values'][0]['value']:
                dates.append(datetime.datetime.strptime(
                    l['dateTime'][:16], '%Y-%m-%dT%H:%M'))
                turb.append(l['value'])
            turb = np.ma.masked_values(
                np.ma.masked_array(turb).astype(np.float32), -999999)

        for j, rt in enumerate(nc.variables['river'][:]):
            if int(rt) == r.usgs_id:
                turb_01 = np.squeeze(nc.variables['river_turb_01'][:, 0, j])
                for l in turb.nonzero()[0]:
                    p = time.index(min(time, key=lambda d: abs(dates[l] - d)))
                    turb_01[p] = turb[l]
                for k in range(nc.dimensions['s_rho'].size):
                    nc.variables['river_turb_01'][:, k, j] = turb_01
        nc.sync()
    nc.close()
    return
