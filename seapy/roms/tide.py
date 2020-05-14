#!/usr/bin/env python
"""
  tide.py

  Methods for working tidal forcing files in ROMS

  Written by Brian Powell on 04/05/16
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import numpy as np
import netCDF4
import seapy
import datetime
from warnings import warn


def create_forcing(filename, tide, title="Tidal Forcing", epoch=seapy.default_epoch):
    """
    Create a tidal forcing file from the given tidal values.

    Parameters
    ----------
    filename: string,
      File name of the tidal forcing to create
    tide: dict,
      Dictionary of the tidal forcing containing the following keys:
       Eamp : SSH amplitdue
       Ephase : SSH phase (radians)
       Cmajor : velocity major ellipse
       Cminor : velocity minor ellipse
       Cphase : velocity ellipse phase (radians)
       Cangle : velocity ellipse angle (radians)
       tide_start : datetime of the tide reference
       tides  : list of the tides
    title: string, optional,
      NetCDF title string to use
    epoch: datetime, optional,
      Epoch date for time reference

    Returns
    -------
    None
    """
    # Create the tide forcing file
    ntides, eta_rho, xi_rho = tide['Eamp'].shape
    if ntides != len(tide['tides']):
        raise ValueError(
            "The number of tidal data are different than the tides.")

    tideout = seapy.roms.ncgen.create_tide(filename, eta_rho=eta_rho,
                                           xi_rho=xi_rho,
                                           reftime=epoch,
                                           ntides=ntides,
                                           clobber=True, title=title)
    # Set the tide periods and attributes
    tideout.variables['tide_period'][:] = 1.0 / \
        seapy.tide.frequency(tide['tides'])
    tideout.setncattr("tidal_constituents", ", ".join(tide['tides']))
    tideout.setncattr("tide_start", "Day {:5.1f} ({:s})".format((tide['tide_start']
                                                                 - epoch).total_seconds() / 86400,
                                                                str(tide['tide_start'])))
    tideout.setncattr("base_date", "days since {:s}".format(
        str(tide['tide_start'])))
    tideout.variables['tide_Eamp'][:] = tide['Eamp']
    tideout.variables['tide_Ephase'][:] = np.degrees(tide['Ephase'])
    tideout.variables['tide_Cmax'][:] = tide['Cmajor']
    tideout.variables['tide_Cmin'][:] = tide['Cminor']
    tideout.variables['tide_Cphase'][:] = np.degrees(tide['Cphase'])
    tideout.variables['tide_Cangle'][:] = np.degrees(tide['Cangle'])
    tideout.close()


def load_forcing(filename):
    """
    Load a tidal forcing file into a dictionary

    Parameters
    ----------
    filename: string
      File name of the tidal forcing file to load

    Returns
    -------
    dict:
      Dictionary of the tidal forcing information with keys:
       Eamp : SSH amplitdue
       Ephase : SSH phase (radians)
       Cmajor : velocity major ellipse
       Cminor : velocity minor ellipse
       Cphase : velocity ellipse phase (radians)
       Cangle : velocity ellipse angle (radians)
       tide_start : datetime of the tide reference
       tides  : list of the tides
    """
    import re

    nc = seapy.netcdf(filename)
    frc = {}
    frc['Eamp'] = nc.variables['tide_Eamp'][:]
    frc['Ephase'] = np.radians(nc.variables['tide_Ephase'][:])
    frc['Cmajor'] = nc.variables['tide_Cmax'][:]
    frc['Cminor'] = nc.variables['tide_Cmin'][:]
    frc['Cphase'] = np.radians(nc.variables['tide_Cphase'][:])
    frc['Cangle'] = np.radians(nc.variables['tide_Cangle'][:])
    start_str = getattr(nc, 'tide_start', None) or \
        getattr(nc, 'base_date', None)
    tides = getattr(nc, 'tidal_constituents', None) or \
        getattr(nc, 'tides', None)
    frc['tides'] = tides.upper().split(", ")
    frc['tide_start'] = None
    nc.close()
    if start_str:
        try:
            frc['tide_start'] = datetime.datetime.strptime(
                re.sub('^.*since\s*', '', start_str),
                "%Y-%m-%d %H:%M:%S")
        except ValueError:
            pass

    return frc


def tide_error(his_file, tide_file, grid=None):
    """
    Calculates the tidal error for each point given a model history and the
    tidal file used

    Parameters
    ----------
    his_file : string,
      String of history file location. Can be multiple files using wildcard
    tide_file: string,
      String of tidal file location
    grid : string or grid, optional,
      If specified, use this grid. Default is to build a grid from the history
      file.

    Returns
    -------
      tide_error : masked_array,
        Array containing the tidal error at each point, with land points masked

    """
    if grid:
        grid = seapy.model.asgrid(grid)
    else:
        grid = seapy.model.asgrid(his_file)

    # Load tidal file data
    frc = load_forcing(tide_file)

    # Calculate tidal error for each point
    nc = seapy.netcdf(his_file)
    times = seapy.roms.num2date(nc)
    tide_error = np.ma.masked_where(
        grid.mask_rho == 0, np.zeros((grid.mask_rho.shape)))
    zeta = nc.variables['zeta'][:]
    nc.close()
    for i in seapy.progressbar.progress(range(grid.ln)):
        for j in range(grid.lm):
            if not tide_error.mask[i, j]:
                z = zeta[:, i, j]
                t_ap = seapy.tide.pack_amp_phase(frc['tides'],
                                                 frc['Eamp'][:, i, j], frc['Ephase'][:, i, j])
                mout = seapy.tide.fit(times, z, tides=frc['tides'],
                                      lat=grid.lat_rho[i, j], tide_start=frc['tide_start'])
                for c in t_ap:
                    m = mout['major'][c]
                    t = t_ap[c]
                    tide_error[i, j] += 0.5 * (m.amp**2 + t.amp**2) - \
                        m.amp * t.amp * np.cos(m.phase - t.phase)
                tide_error[i, j] = np.sqrt(tide_error[i, j])
    return tide_error
