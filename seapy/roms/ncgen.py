#!/usr/bin/env python
"""
  Functions to generate ROMS netcdf files

  Written by Brian Powell on 04/26/13
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import os
import re
import netCDF4
import numpy as np
from datetime import datetime
from seapy.lib import default_epoch
from seapy.cdl_parser import cdl_parser
from seapy.roms import lib
from warnings import warn

"""
    Module variables
"""
_cdl_dir = os.path.dirname(lib.__file__)
_cdl_dir = "/".join((('.' if not _cdl_dir else _cdl_dir), "cdl/"))
_format = "NETCDF4_CLASSIC"


def __number_or_string(val):
    """
    convert a string to a number if the string represents a number;
    otherwise, return the string.
    """
    try:
        val = float(val.strip())
    except ValueError:
        pass
    return val


def ncgen(filename, dims=None, vars=None, attr=None, title=None,
          clobber=False, format=_format):
    """
    Create a new netcdf file with the given definitions. Need to define
    the dimensions, the variables, and the attributes.

    Parameters
    ----------
    filename : string
        name and path of file to create
    dims : dict
        dictionary of dimensions with dimension name as keys, and the value
        as the length of the dimension. NOTE: 0 value means UNLIMITED.
    vars: list of dictionaries
        each variable to define is a dictionary that contains three keys:
            name: string name of variable
            type: string type (float, double, etc.)
            dims: comma separated string of dimensions ("ocean_time, eta_rho")
            attr: dictionary of variable attributes where the key is
                  the attribute name and the value is the attribute string
    attr: dict, optional
        optional dictionary of global attributes for the netcdf file:
        key is the attribute name and the value is the attribute string
    title: string, optional
        netcdf attribute title
    clobber: bool, optional
        If True, destroy existing file
    format: string, optional
        NetCDF format to use. Default is NETCDF4_CLASSIC

    Returns
    -------
    nc, netCDF4 object

    Examples
    --------
    >>> dims = {"ocean_time":0, "eta_rho":120, "xi_rho":100}
    >>> vars = [  {"name":"eta_slice", "type":"double",
                   "dims":"ocean_time, eta_rho",
                   "attr":{"units":"degrees Celcius"}},
                  {"name":"xi_slice", "type":"double",
                   "dims":"ocean_time, xi_rho",
                   "attr":{"units":"degrees Celcius"}} ]
    >>> seapy.roms.ncgen("test.nc", dims=dims, vars=vars, title="Test")
    """
    vars = np.atleast_1d(vars)
    if dims is None:
        dims = {}
    if attr is None:
        attr = {}
    # Create the file
    if not os.path.isfile(filename) or clobber:
        _nc = netCDF4.Dataset(filename, "w", format=format)
        # Loop over the dimensions and add them
        for dim in dims:
            _nc.createDimension(dim, dims[dim])
        # Loop over the variables and add them
        for var in vars:
            add_variable(_nc, var)

        # Add global attributes
        for a in attr:
            _nc.setncattr(a, attr[a])

        try:
            _nc.author = os.getenv('USER') or \
                os.getenv('LOGNAME') or \
                os.getenv('USERNAME') or \
                os.getlogin() or \
                'nobody'
        except (AttributeError, IOError, OSError, FileNotFoundError) as e:
            _nc.author = 'nobody'

        _nc.history = datetime.now().strftime(
            "Created on %a, %B %d, %Y at %H:%M")
        if title is not None:
            _nc.title = title
        _nc.close()
    else:
        warn(filename + " already exists. Using existing definition")
    return netCDF4.Dataset(filename, "a")
    pass


def _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho):
    """
        internal method: Set grid dimensions
    """
    if "xi_rho" in dims.keys():
        dims["xi_rho"] = xi_rho
    if "xi_u" in dims.keys():
        dims["xi_u"] = xi_rho - 1
    if "xi_v" in dims.keys():
        dims["xi_v"] = xi_rho
    if "xi_psi" in dims.keys():
        dims["xi_psi"] = xi_rho - 1
    if "eta_rho" in dims.keys():
        dims["eta_rho"] = eta_rho
    if "eta_u" in dims.keys():
        dims["eta_u"] = eta_rho
    if "eta_v" in dims.keys():
        dims["eta_v"] = eta_rho - 1
    if "eta_psi" in dims.keys():
        dims["eta_psi"] = eta_rho - 1
    if "s_rho" in dims.keys():
        dims["s_rho"] = s_rho
    if "s_w" in dims.keys():
        dims["s_w"] = s_rho + 1

    return dims


def _set_time_ref(vars, timevar, reftime, cycle=None):
    """
        internal method: Set time reference
    """
    if isinstance(timevar, str):
        timevar = [timevar]
    for tvar in timevar:
        for nvar in vars:
            if nvar["name"] == tvar:
                if "units" in nvar["attr"]:
                    t = re.findall('(\w+) since .*', nvar["attr"]["units"])
                    nvar["attr"]["units"] = \
                        "{:s} since {:s}".format(t[0], str(reftime))
                else:
                    nvar["attr"]["units"] = \
                        "days since {:s}".format(str(reftime))
                if cycle is not None:
                    nvar["attr"]["cycle_length"] = cycle
    return vars


def add_variable(nc, var):
    """
    Add a new variable with meta data to an existing netcdf file

    Parameters
    ----------
    filename : string or netCDF4 object
        name or file to add the variable to
    vars: dictionary
            name: string name of variable
            type: string type (float, double, etc.)
            dims: comma separated string of dimensions ("ocean_time, eta_rho")
            attr: dictionary of variable attributes where the key is
                  the attribute name and the value is the attribute string

    Returns
    -------
    nc, netCDF4 object

    Examples
    --------
    >>> var = {"name":"eta_slice", "type":"double",
               "dims":"ocean_time, eta_rho",
               "attr":{"units":"degrees Celcius"}}
    >>> nc = seapy.roms.ncgen.add_variable("test.nc", var)
    """
    if nc is None:
        raise AttributeError("No file was specified")
    if isinstance(nc, netCDF4._netCDF4.Dataset):
        pass
    else:
        nc = netCDF4.Dataset(nc, "a")

    # Handle the dimensions by enforcing a tuple list rather
    # than a list of strings, then add whatever we have
    try:
        dims = var['dims'].replace(" ", "").split(',')
    except:
        dims = var['dims']
    try:
        nvar = nc.createVariable(var["name"], var["type"], dims)
    except:
        nvar = nc.createVariable(var["name"], var["type"])

    try:
        for key in var["attr"]:
            # Check if it is a number and convert
            astr = __number_or_string(var["attr"][key])
            setattr(nvar, key, astr)
    except KeyError:
        pass

    return nc


def _create_generic_file(filename, cdl, eta_rho, xi_rho, s_rho,
                         reftime=None, clobber=False, title="ROMS"):
    """
        internal method: Generic file creator that uses ocean_time
    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    if reftime is not None:
        vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_psource(filename, nriver=1, s_rho=5,
                   reftime=default_epoch, clobber=False, cdl=None, title="My River"):
    """
    Create a new, blank point source file

    Parameters
    ----------
    filename : string
        name and path of file to create
    nriver : int, optional
        number of rivers to put in file
    s_rho: int, optional
        number of s-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "frc_rivers.cdl" if cdl is None else cdl)

    # Fill in the appropriate river values
    dims["river"] = nriver
    dims["s_rho"] = s_rho
    vars = _set_time_ref(vars, "river_time", reftime)

    # Create the river file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_grid(filename, eta_rho=10, xi_rho=10, s_rho=1, clobber=False,
                cdl=None, title="My Grid"):
    """
    Create a new, blank grid file


    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "roms_grid.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)

    print(dims)

    # Create the grid file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_adsen(filename, eta_rho=10, xi_rho=10, s_rho=1,
                 reftime=default_epoch, clobber=False, cdl=None, title="My Adsen"):
    """
    Create a new adjoint sensitivity file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Create the general file
    return _create_generic_file(filename, _cdl_dir + "adsen.cdl" if cdl is None else cdl,
                                eta_rho, xi_rho, s_rho, reftime, clobber, title)


def create_bry(filename, eta_rho=10, xi_rho=10, s_rho=1,
               reftime=default_epoch, clobber=False, cdl=None, title="My BRY"):
    """
    Create a bry forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "bry_unlimit.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "bry_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_clim(filename, eta_rho=10, xi_rho=10, s_rho=1,
                reftime=default_epoch, clobber=False, cdl=None, title="My CLIM"):
    """
    Create a climatology forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "clm_ts.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "clim_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_bulk(filename, lat=10, lon=10,
                    reftime=default_epoch, clobber=False, cdl=None,
                    title="My Forcing"):
    """
    Create a bulk flux forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "frc_bulk.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims["lat"] = lat
    dims["lon"] = lon
    vars = _set_time_ref(vars, "frc_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_direct(filename, eta_rho=10, xi_rho=10,
                      reftime=default_epoch, clobber=False, cdl=None,
                      title="My Forcing"):
    """
    Create a direct surface forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "frc_direct.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = {'y_rho': eta_rho,
            'y_u': eta_rho,
            'y_v': eta_rho - 1,
            'x_rho': xi_rho,
            'x_u': xi_rho - 1,
            'x_v': xi_rho,
            'frc_time': 0}
    vars = _set_time_ref(vars, 'frc_time', reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_flux(filename, eta_rho=10, xi_rho=10, ntimes=1,
                    cycle=None, reftime=default_epoch, clobber=False,
                    cdl=None, title="My Flux"):
    """
    Create a surface flux forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    ntimes: int, optional
        number of time records (climatology files do not have unlimited
        dimension)
    cycle: int or None, optional
        The number of days before cycling the forcing records
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "frc_fluxclm.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, 1)
    times = ("srf_time", "shf_time", "swf_time", "sss_time")
    for n in times:
        dims[n] = ntimes
    vars = _set_time_ref(vars, times, reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_srelax(filename, eta_rho=10, xi_rho=10, s_rho=1, cycle=None,
                      reftime=default_epoch, clobber=False, cdl=None,
                      title="My Srelaxation"):
    """
    Create a Salt Relaxation forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    cycle: int or None, optional
        The number of days before cycling the forcing records
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "frc_srelax.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "sss_time", reftime, cycle)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_qcorr(filename, eta_rho=10, xi_rho=10, s_rho=1, cycle=None,
                     reftime=default_epoch, clobber=False, cdl=None,
                     title="My Qcorrection"):
    """
    Create a Q Correction forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    cycle: int or None, optional
        The number of days before cycling the forcing records
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "frc_qcorr.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "sst_time", reftime, cycle)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_wind(filename, eta_rho=10, xi_rho=10, s_rho=1, cycle=None,
                    reftime=default_epoch, clobber=False, cdl=None,
                    title="My Winds"):
    """
    Create a surface wind stress forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    cycle: int or None, optional
        The number of days before cycling the forcing records
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "frc_windstress.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "sms_time", reftime, cycle)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_wave(filename, eta_rho=10, xi_rho=10, reftime=default_epoch,
                    clobber=False, cdl=None, title="My Waves"):
    """
    Create a surface wave forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "frc_wave.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho=1)
    vars = _set_time_ref(vars, "wave_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_tide(filename, eta_rho=10, xi_rho=10, s_rho=1, ntides=1,
                reftime=default_epoch, clobber=False,
                title="My Tides"):
    """
    Create a barotropic tide forcing file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    ntides: int, optional
        number of tidal frequencies to force with
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "frc_tides.cdl")

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    dims["tide_period"] = ntides

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_ini(filename, eta_rho=10, xi_rho=10, s_rho=1,
               reftime=default_epoch, clobber=False, cdl=None, title="My Ini"):
    """
    Create an initial condition file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "ini_hydro.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_nudge_coef(filename, eta_rho=10, xi_rho=10, s_rho=1, clobber=False,
                      cdl=None, title="My Nudging"):
    """
    Create a nudging coefficients file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "nudge_coef.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_da_obs(filename, state_variable=20, survey=1, provenance=None,
                  clobber=False, cdl=None, title="My Observations"):
    """
    Create an assimilation observations file

    Parameters
    ----------
    filename : string
        name and path of file to create
    survey: int, optional
        number of surveys in the file
    state_variable: int, optional
        number of state variables in the observations
    provenance: string, optional
        Description of the provenance values
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object


    """

    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "s4dvar_obs.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims["survey"] = survey
    dims["state_variable"] = state_variable

    # Set the provenance values in the global attributes
    if provenance is not None:
        attr["obs_provenance"] = str(provenance)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title, format="NETCDF3_64BIT")

    # Return the new file
    return _nc


def create_da_ray_obs(filename, ray_datum=1, provenance="None",
                      reftime=default_epoch, clobber=False,
                      cdl=None, title="My Observations"):
    """
    Create an acoustic ray assimilation observations file

    Parameters
    ----------
    filename : string
        name and path of file to create
    ray_datum: int, optional
        Number of rays to assimilate
    provenance: string, optional
        Description of the provenance values
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """

    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "s4dvar_obs_ray.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims["ray_datum"] = ray_datum
    vars = _set_time_ref(vars, "obs_time", reftime)

    # Set the provenance values in the global attributes
    attr["obs_provenance"] = provenance

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_da_bry_std(filename, eta_rho=10, xi_rho=10, s_rho=1, bry=4,
                      reftime=default_epoch, clobber=False, cdl=None,
                      title="My BRY STD"):
    """
    Create a boundaries standard deviation file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    bry: int, optional
        number of open boundaries to specify
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "s4dvar_std_b.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    dims["IorJ"] = max(eta_rho, xi_rho)
    dims["boundary"] = bry
    vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_da_frc_std(filename, eta_rho=10, xi_rho=10, s_rho=1,
                      reftime=default_epoch, clobber=False,
                      cdl=None, title="My FRC STD"):
    """
    Create a forcing standard deviation file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "s4dvar_std_f.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_da_ini_std(filename, eta_rho=10, xi_rho=10, s_rho=1,
                      reftime=default_epoch, clobber=False,
                      cdl=None, title="My INI STD"):
    """
    Create an initialization standard deviation file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "s4dvar_std_i.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_da_model_std(filename, eta_rho=10, xi_rho=10, s_rho=1,
                        reftime=default_epoch, clobber=False,
                        cdl=None, title="My Model STD"):
    """
    Create an time varying model standard deviation file

    Parameters
    ----------
    filename : string
        name and path of file to create
    eta_rho: int, optional
        number of rows in the eta direction
    xi_rho: int, optional
        number of columns in the xi direction
    s_rho: int, optional
        number of s-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(
        _cdl_dir + "s4dvar_std_m.cdl" if cdl is None else cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_zlevel_grid(filename, lat=10, lon=10, depth=1,
                       clobber=False, cdl=None,
                       title="Zlevel Grid", dims=2):
    """
    Create z-level grid file

    Parameters
    ----------
    filename : string
        name and path of file to create
    lat: int, optional
        number of latitudinal rows
    lon: int, optional
        number of longitudinal columns
    depth: int, optional
        number of z-levels
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title
    dims: int, optional
        number of dimensions to use for lat/lon

    Returns
    -------
    nc, netCDF4 object

    """
    if cdl == None:
        if dims == 1:
            cdlfile = _cdl_dir + "zlevel_1d_grid.cdl"
        else:
            cdlfile = _cdl_dir + "zlevel_2d_grid.cdl"
    else:
        cdlfile = cdl

    # Generate the Structure
    dims, vars, attr = cdl_parser(cdlfile)

    # Fill in the appropriate dimension values
    dims["lat"] = lat
    dims["lon"] = lon
    dims["depth"] = depth

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_zlevel(filename, lat=10, lon=10, depth=1,
                  reftime=default_epoch,
                  clobber=False, cdl=None,
                  title="Zlevel Model Data", dims=2):
    """
    Create an time varying model standard deviation file

    Parameters
    ----------
    filename : string
        name and path of file to create
    lat: int, optional
        number of latitudinal rows
    lon: int, optional
        number of longitudinal columns
    depth: int, optional
        number of z-levels
    reftime: datetime, optional
        date of epoch for time origin in netcdf
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.
    title: string, optional
        netcdf attribute title
    dims: int, optional
        number of dimensions to use for lat/lon

    Returns
    -------
    nc, netCDF4 object

    """
    if cdl == None:
        if dims == 1:
            cdlfile = _cdl_dir + "zlevel_1d.cdl"
        else:
            cdlfile = _cdl_dir + "zlevel_2d.cdl"
    else:
        cdlfile = cdl

    # Generate the Structure
    dims, vars, attr = cdl_parser(cdlfile)

    # Fill in the appropriate dimension values
    dims["lat"] = lat
    dims["lon"] = lon
    dims["depth"] = depth
    vars = _set_time_ref(vars, "time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


if __name__ == "__main__":
    grid = create_zlevel("test.nc")
