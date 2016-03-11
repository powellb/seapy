#!/usr/bin/env python
"""
  Functions to generate ROMS netcdf files

  Written by Brian Powell on 04/26/13
  Copyright (c)2016 University of Hawaii under the BSD-License.
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
            if var["dims"][0]:
                nvar = _nc.createVariable(var["name"], var["type"],
                                          var["dims"])
            else:
                nvar = _nc.createVariable(var["name"], var["type"])
            try:
                for key in var["attr"]:
                    setattr(nvar, key, var["attr"][key])
            except KeyError:
                pass
        # Add global attributes
        for a in attr:
            setattr(_nc, a, attr[a])
        _nc.author = os.getlogin()
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
    dims["xi_rho"] = xi_rho
    dims["xi_u"] = xi_rho - 1
    dims["xi_v"] = xi_rho
    dims["xi_psi"] = xi_rho - 1
    dims["eta_rho"] = eta_rho
    dims["eta_u"] = eta_rho
    dims["eta_v"] = eta_rho - 1
    dims["eta_psi"] = eta_rho - 1
    dims["N"] = s_rho

    # Fill in the appropriate dimension values
    try:
        dims["s_rho"] = s_rho
        dims["s_w"] = s_rho + 1
    except KeyError:
        pass

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


def _create_generic_file(filename, cdl, eta_rho, xi_rho, s_rho,
                         reftime=None, clobber=False, title="ROMS"):
    """
        internal method: Generic file creator that uses ocean_time
    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + cdl)

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    if reftime is not None:
        vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_river(filename, nriver=1, s_rho=5,
                 reftime=default_epoch, clobber=False, title="My River"):
    """
    Create a new, blank river file

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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "frc_rivers.cdl")

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
                title="My Grid"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "roms_grid.cdl")

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)

    # Create the grid file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_adsen(filename, eta_rho=10, xi_rho=10, s_rho=1,
                 reftime=default_epoch, clobber=False, title="My Adsen"):
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
    return _create_generic_file(filename, "adsen.cdl", eta_rho, xi_rho, s_rho,
                                reftime, clobber, title)


def create_bry(filename, eta_rho=10, xi_rho=10, s_rho=1,
               reftime=default_epoch, clobber=False, title="My BRY"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "bry_unlimit.cdl")

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "bry_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_clim(filename, eta_rho=10, xi_rho=10, s_rho=1,
                reftime=default_epoch, clobber=False, title="My CLIM"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "clm_ts.cdl")

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "clim_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_bulk(filename, lat=10, lon=10,
                    reftime=default_epoch, clobber=False,
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "frc_bulk.cdl")

    # Fill in the appropriate dimension values
    dims["lat"] = lat
    dims["lon"] = lon
    vars = _set_time_ref(vars, "time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_flux(filename, eta_rho=10, xi_rho=10, ntimes=1,
                    cycle=None, reftime=default_epoch, clobber=False,
                    title="My Flux"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "frc_fluxclm.cdl")

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


def create_frc_qcorr(filename, eta_rho=10, xi_rho=10, s_rho=1, cycle=None,
                     reftime=default_epoch, clobber=False,
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "frc_qcorr.cdl")

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "sst_time", reftime, cycle)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_frc_wind(filename, eta_rho=10, xi_rho=10, s_rho=1, cycle=None,
                    reftime=default_epoch, clobber=False,
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "frc_windstress.cdl")

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "sms_time", reftime, cycle)

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
               reftime=default_epoch, clobber=False, title="My Ini"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "ini_hydro.cdl")

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_da_obs(filename, state_variable=20, survey=1, provenance=None,
                  clobber=False, title="My Observations"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object


    """

    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "s4dvar_obs.cdl")

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
                      title="My Observations"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """

    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "s4dvar_obs_ray.cdl")

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
                      reftime=default_epoch, clobber=False,
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "s4dvar_std_b.cdl")

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
                      title="My FRC STD"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "s4dvar_std_f.cdl")

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
                      title="My INI STD"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "s4dvar_std_i.cdl")

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
                        title="My Model STD"):
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
    title: string, optional
        netcdf attribute title

    Returns
    -------
    nc, netCDF4 object

    """
    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + "s4dvar_std_m.cdl")

    # Fill in the appropriate dimension values
    dims = _set_grid_dimensions(dims, eta_rho, xi_rho, s_rho)
    vars = _set_time_ref(vars, "ocean_time", reftime)

    # Create the file
    _nc = ncgen(filename, dims=dims, vars=vars, attr=attr, clobber=clobber,
                title=title)

    # Return the new file
    return _nc


def create_zlevel_grid(filename, lat=10, lon=10, depth=1,
                       clobber=False,
                       title="Zlevel Grid", cdlfile=None, dims=2):
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
    title: string, optional
        netcdf attribute title
    cdlfile: string, optional
        name of CDL file to use for construction
    dims: int, optional
        number of dimensions to use for lat/lon

    Returns
    -------
    nc, netCDF4 object

    """
    if cdlfile == None:
        if dims == 1:
            cdlfile = "zlevel_1d_grid.cdl"
        else:
            cdlfile = "zlevel_2d_grid.cdl"

    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + cdlfile)

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
                  clobber=False,
                  title="Zlevel Model Data", cdlfile=None, dims=2):
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
    title: string, optional
        netcdf attribute title
    cdlfile: string, optional
        name of CDL file to use for construction
    dims: int, optional
        number of dimensions to use for lat/lon

    Returns
    -------
    nc, netCDF4 object

    """
    if cdlfile == None:
        if dims == 1:
            cdlfile = "zlevel_1d.cdl"
        else:
            cdlfile = "zlevel_2d.cdl"

    # Generate the Structure
    dims, vars, attr = cdl_parser(_cdl_dir + cdlfile)

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


