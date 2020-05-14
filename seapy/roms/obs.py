#!/usr/bin/env python
"""
  obs.py

  State Estimation and Analysis for PYthon

  Module to handle the observation structure within ROMS. The ROMS structure
  defines the obs_provenance, which is a numeric ID for tracking the source
  of the observations you use. This module defines a dictionary to translate
  between the numeric and string representations so that either can be used.
  Standard instruments are predefined; however, for your own applications,
  you will want to define additional provenances for the observations you use.
  This is accomplished via:

  >>> import seapy
  >>> seapy.roms.obs.obs_provenance.update({353:'MY_OBS1', 488:'MY_OBS2'})

  You can make your own module for importing seapy and adding your definitions
  easily.

  Written by Brian Powell on 08/05/14
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import numpy as np
import netCDF4
import seapy
from collections import namedtuple
from warnings import warn

# Define a named tuple to store raw data for the gridder
raw_data = namedtuple('raw_data', 'type provenance values error min_error')

# Define the observation type
obs_types = {
    1: "ZETA",
    2: "UBAR",
    3: "VBAR",
    4: "U",
    5: "V",
    6: "TEMP",
    7: "SALT",
    20: "RADIAL"
}
# Define the observation provenances used within my applications
obs_provenance = {
    0: "UNKNOWN",
    100: "GLIDER",
    102: "GLIDER_SG022",
    103: "GLIDER_SG023",
    114: "GLIDER_SG114",
    139: "GLIDER_SG139",
    146: "GLIDER_SG146",
    147: "GLIDER_SG147",
    148: "GLIDER_SG148",
    150: "GLIDER_SG500",
    151: "GLIDER_SG511",
    152: "GLIDER_SG512",
    153: "GLIDER_SG513",
    154: "GLIDER_SG523",
    155: "GLIDER_SG626",
    200: "CTD",
    210: "CTD_HOT",
    220: "CTD_ARGO",
    300: "SST",
    301: "SST_OSTIA",
    315: "SST_NAVO_MAP",
    317: "SST_AVHRR_17",
    318: "SST_AVHRR_18",
    330: "SST_MODIS_AQUA",
    331: "SST_MODIS_TERRA",
    332: "SST_VIIRS",
    340: "SST_REMSS",
    341: "SST_AMSRE",
    342: "SST_TMI",
    400: "SSH",
    410: "SSH_AVISO_MAP",
    411: "SSH_AVISO_TOPEX_POSEIDON",
    412: "SSH_AVISO_JASON1",
    413: "SSH_AVISO_JASON2",
    414: "SSH_AVISO_JASON3",
    420: "SSH_AVISO_GFO",
    421: "SSH_AVISO_ENVISAT",
    422: "SSH_AVISO_ERS1",
    423: "SSH_AVISO_ERS2",
    430: "SSH_AVISO_ALTIKA",
    431: "SSH_AVISO_CRYOSAT2",
    432: "SSH_AVISO_HAIYANG",
    433: "SSH_AVISO_SENTINEL3A",
    450: "SSH_HYCOM",
    460: "SSS_AQUARIUS",
    500: "DRIFTERS",
    600: "RADAR",
    610: "RADAR_KOK",
    620: "RADAR_KAK",
    630: "RADAR_KAL",
    640: "RADAR_KAP",
    650: "RADAR_KNA",
    660: "RADAR_KKH",
    670: "RADAR_PPK",
    700: "ADCP",
    800: "MOORING",
    810: "TAO_ARRAY"
}


def _type_from_string(s):
    """
    PRIVATE method: Search the type dictionary for the key of the
    given the value. If the key isn't a string or properly resolves, try to
    return the value as such
    """
    try:
        return list(obs_types.keys())[
            list(obs_types.values()).index(s.upper())]
    except AttributeError:
        return int(s)


def _provenance_from_string(s):
    """
    PRIVATE method: Search the provenance dictionary for the key of the
    given the value. If the key isn't a string or properly resolves, try to
    return the value as such
    """
    try:
        return list(obs_provenance.keys())[
            list(obs_provenance.values()).index(s.upper())]
    except AttributeError:
        return int(s)


def asobs(obs):
    """
    Return the input as an observation array if possible. If the parameter
    is already an observation, just return; otherwise, create a new class.

    Parameters
    ----------
    obs: obs class, string, or list
        what to cast as observation

    Output
    ------
    obs: seapy.roms.obs.obs
    """
    if obs is None:
        raise AttributeError("No obs were specified")
    if isinstance(obs, seapy.roms.obs.obs):
        return obs
    else:
        return seapy.roms.obs.obs(filename=obs)


def astype(otype):
    """
    Return the integer type of the given observation array.

    Input
    -----
    type: ndarray,
        List of types to put into integer form

    Output
    ------
    types: array,
        List of integer types
    """
    otype = np.atleast_1d(otype)
    if otype.dtype.type == np.str_:
        return np.array([_type_from_string(s) for s in otype])
    else:
        return otype


def asprovenance(prov):
    """
    Return the integer provenance of the given provenance array.

    Input
    -----
    prov: array,
        List of provenances to put into integer form

    Output
    ------
    provs: ndarray,
        List of integer provenances
    """
    prov = np.atleast_1d(prov)
    if prov.dtype.type == np.str_:
        return np.array([_provenance_from_string(s) for s in prov])
    else:
        return prov


class obs:

    def __init__(self, filename=None, time=None, x=None, y=None, z=None,
                 lat=None, lon=None, depth=None, value=None, error=None,
                 type=None, provenance=None, meta=None,
                 title="ROMS Observations"):
        """
        Class to deal with ROMS observations for data assimilation

        Parameters
        ----------
        filename : string or list, optional,
            if filename is given, the data are loaded from a netcdf file
        time : ndarray, optional,
          time of observation in days
        x : ndarray, optional,
          obs location on grid in x (eta)
        y : ndarray, optional,
          obs location on grid in y (xi)
        z : ndarray, optional,
          obs location on grid in z (positive layers or negative depth [m])
        lat : ndarray, optional,
          obs true latitude [deg]
        lon : ndarray, optional,
          obs true longitude [deg]
        depth : ndarray, optional,
          obs true depth [m]
        value : ndarray, optional,
          obs value [units]
        error : ndarray, optional,
          obs error [units**2]
        type : ndarray, optional,
          obs type [1-zeta, 2-ubar, 3-vbar, 4-u, 5-v, 6-temp, 7-salt]
        provenance : ndarray, optional,
          obs provenance
        meta : ndarray, optional,
          obs additional information
        """
        self.title = title
        if filename is not None:
            nc = seapy.netcdf(filename)
            # Construct an array from the data in the file. If obs_meta
            # exists in the file, then load it; otherwise, fill with zeros
            self.filename = filename
            self.time = nc.variables["obs_time"][:]
            self.x = nc.variables["obs_Xgrid"][:]
            self.y = nc.variables["obs_Ygrid"][:]
            self.z = nc.variables["obs_Zgrid"][:]
            self.lat = nc.variables["obs_lat"][:]
            self.lon = nc.variables["obs_lon"][:]
            self.depth = nc.variables["obs_depth"][:]
            self.value = nc.variables["obs_value"][:]
            self.error = nc.variables["obs_error"][:]
            self.type = nc.variables["obs_type"][:]
            self.provenance = nc.variables["obs_provenance"][:]
            # Update the provenance definitions
            try:
                obs_provenance.update(dict((int(k.strip()), v.strip())
                                           for v, k in
                                           (it.split(':') for it in
                                            nc.obs_provenance.split(','))))
            except (AttributeError, ValueError):
                pass
            try:
                self.meta = nc.variables["obs_meta"][:]
            except KeyError:
                self.meta = np.zeros(self.value.size)
            finally:
                nc.close()
        else:
            self.filename = None
            if time is not None:
                self.time = np.atleast_1d(time)
            if x is not None:
                self.x = np.atleast_1d(x)
            if y is not None:
                self.y = np.atleast_1d(y)
            if z is not None:
                self.z = np.atleast_1d(z)
            if lat is not None:
                self.lat = np.atleast_1d(lat)
            if lon is not None:
                self.lon = np.atleast_1d(lon)
            if depth is not None:
                self.depth = np.atleast_1d(depth)
            if value is not None:
                self.value = np.atleast_1d(value)
            if error is not None:
                self.error = np.atleast_1d(error)
            if type is not None:
                self.type = astype(type)
            if provenance is not None:
                self.provenance = asprovenance(provenance)
            else:
                self.provenance = 0
            if meta is not None:
                self.meta = np.atleast_1d(meta)
            self._consistent()

    def _consistent(self):
        """
        PRIVATE method: try to make the structure self-consistent. Throw
        an exception if not possible.
        """
        # Make sure required arrays are a 1d array
        self.time = self.time.ravel()
        self.x = self.x.ravel()
        self.y = self.y.ravel()
        self.value = self.value.ravel()
        self.error = self.error.ravel()
        self.type = astype(self.type.ravel())

        lt = self.time.size
        if not lt == self.x.size == self.y.size == \
                self.value.size == self.error.size == self.type.size:
            # If these lengths are not equal, then there is a serious issue
            raise ValueError(
                "Lengths of observation attributes are not equal.")
        else:
            # For the others, we can pad the information to ensure
            # consistency
            def _resizearr(key, n):
                arr = getattr(self, key, np.zeros(n))
                if arr.size == n:
                    return arr
                return np.resize(arr, n)

            self.z = _resizearr('z', lt)
            self.lat = _resizearr('lat', lt)
            self.lon = _resizearr('lon', lt)
            self.depth = _resizearr('depth', lt)
            self.provenance = asprovenance(_resizearr('provenance', lt))
            self.meta = _resizearr('meta', lt)

        # Eliminate bad values
        good_vals = np.logical_and.reduce((
            np.isfinite(self.value),
            np.isfinite(self.x),
            np.isfinite(self.y),
            np.isfinite(self.error),
            np.isfinite(self.time)))
        if np.any(~good_vals):
            self.delete(np.where(good_vals == False))

        # Set the shape parameter
        self.shape = self.value.shape

        # Ensure consistency in depth and z
        self.z[self.depth > 0] = self.depth[self.depth > 0]

    def __len__(self):
        self.shape = self.value.shape
        return self.value.size

    def __getitem__(self, l):
        return obs(time=self.time[l], x=self.x[l], y=self.y[l],
                   z=self.z[l], lon=self.lon[l], lat=self.lat[l],
                   depth=self.depth[l], value=self.value[l],
                   error=self.error[l], type=self.type[l],
                   provenance=self.provenance[l], meta=self.meta[l])

    def __setitem__(self, l, new_obs):
        if not isinstance(new_obs, seapy.roms.obs.obs):
            raise TypeError("Trying to assign obs to a non-obs type.")

        self.time[l] = new_obs.time
        self.x[l] = new_obs.x
        self.y[l] = new_obs.y
        self.z[l] = new_obs.z
        self.lon[l] = new_obs.lon
        self.lat[l] = new_obs.lat
        self.depth[l] = new_obs.depth
        self.value[l] = new_obs.value
        self.error[l] = new_obs.error
        self.type[l] = new_obs.type
        self.provenance[l] = new_obs.provenance
        self.meta[l] = new_obs.meta
        self._consistent()

    def __repr__(self):
        return "< {:d} obs: {:.1f} to {:.1f} >".format(self.value.size,
                                                       np.min(self.time), np.max(self.time))

    def __str__(self):
        return "\n".join([repr(self), "\n".join(
            "{:.3f}, [{:s}:{:s}] ({:.2f},{:.2f},{:.2f}) = {:.4f} +/- {:.4f}".format(
                t, obs_types[self.type[n]],
                obs_provenance.get(self.provenance[n], "UNKNOWN"),
                self.lon[n], self.lat[n], self.depth[n],
                self.value[n], self.error[n])
            for n, t in enumerate(self.time))])

    def add(self, new_obs):
        """
        Add another class of obs into this one

        Parameters
        ----------
        new_obs : obs,
            Class of obs to append to the existing

        Returns
        -------
        None

        Examples
        --------
        Load a list from netcdf, then append a new set of values

        >>> a=obs("test.nc")
        >>> b=obs(time=4,x=3.2,y=2.8,z=0,value=23.44,error=0.5,type="temp",
        >>>       provenance="glider")
        >>> a.add(b)

        """
        self._consistent()
        new_obs._consistent()
        self.time = np.append(self.time, new_obs.time)
        self.x = np.append(self.x, new_obs.x)
        self.y = np.append(self.y, new_obs.y)
        self.z = np.append(self.z, new_obs.z)
        self.lat = np.append(self.lat, new_obs.lat)
        self.lon = np.append(self.lon, new_obs.lon)
        self.depth = np.append(self.depth, new_obs.depth)
        self.value = np.append(self.value, new_obs.value)
        self.error = np.append(self.error, new_obs.error)
        self.type = np.append(self.type, new_obs.type)
        self.provenance = np.append(self.provenance, new_obs.provenance)
        self.meta = np.append(self.meta, new_obs.meta)

    def copy(self):
        """
        deep copy this class and return the new copy.

        Returns
        -------
        obs : obs,
            deep copy of the class
        """
        import copy
        return copy.deepcopy(self)

    def delete(self, obj):
        """
        delete observations from the record.

        Parameters
        ----------
        obj : slice, int or array of ints
            Indicate which sub-arrays to remove.

        Returns
        -------
        Nothing: updates the class arrays

        Examples
        --------
        Delete every other observation
        >>> myobs.delete(np.s_[::2])
        """
        self.time = np.delete(self.time, obj)
        self.x = np.delete(self.x, obj)
        self.y = np.delete(self.y, obj)
        self.z = np.delete(self.z, obj)
        self.lat = np.delete(self.lat, obj)
        self.lon = np.delete(self.lon, obj)
        self.depth = np.delete(self.depth, obj)
        self.value = np.delete(self.value, obj)
        self.error = np.delete(self.error, obj)
        self.type = np.delete(self.type, obj)
        self.provenance = np.delete(self.provenance, obj)
        self.meta = np.delete(self.meta, obj)

    def create_survey(self, dt=0):
        """
        Build the survey structure from the observations
        """
        # Generate the sort list
        self.sort = np.argsort(self.time, kind='mergesort')

        # Build the survey structure
        times, counts = np.unique(self.time[self.sort], return_counts=True)

        # Make sure everything is within dt
        if dt:
            delta = np.diff(times)
            while np.any(delta < dt):
                idx = np.argmin(delta)
                self.time[self.time == times[idx + 1]] = times[idx]
                times[idx + 1] = times[idx]
                times = np.unique(times)
                delta = np.diff(times)

            # Re-sort the surveys
            times, counts = np.unique(self.time[self.sort], return_counts=True)

        self.survey_time = times
        self.nobs = counts

    def to_netcdf(self, filename=None, dt=0, clobber=True):
        """
        Write out the observations into the specified netcdf file

        Parameters
        ----------
        filename : string, optional
            name of file to save. If obs were loaded from a file and filename
            is not specified, then write to the same.
        dt : float,
            ensure data are at least separated in time by dt; otherwise,
            make as part of same survey
        clobber : bool, optional
            if True, any existing file is overwritten
        """
        import os

        # Check filename
        if filename is None and self.filename is not None:
            filename = self.filename
        if filename is None:
            error("No filename given")

        # Save out the observations by survey
        self._consistent()
        self.create_survey(dt)
        if not self.value.size:
            warn(
                "No observations are available to be written to {:s}".format(filename))
            return None

        if not clobber and os.path.exists(filename):
            warn("{:s} exists with no clobber.".format(filename))
            return None

        state_vars = np.maximum(7, np.max(self.type))
        nc = seapy.roms.ncgen.create_da_obs(filename,
                                            survey=self.survey_time.size,
                                            state_variable=state_vars,
                                            provenance=','.join((':'.join(
                                                (obs_provenance.get(v, "UNKNOWN"), str(v)))
                                                for v in np.unique(self.provenance))),
                                            clobber=True, title=self.title)
        nc.variables["spherical"][:] = 1
        nc.variables["Nobs"][:] = self.nobs
        nc.variables["survey_time"][:] = self.survey_time
        nc.variables["obs_variance"][:] = np.ones(state_vars) * 0.1
        nc.variables["obs_time"][:] = self.time[self.sort]
        nc.variables["obs_Xgrid"][:] = self.x[self.sort]
        nc.variables["obs_Ygrid"][:] = self.y[self.sort]
        nc.variables["obs_Zgrid"][:] = self.z[self.sort]
        nc.variables["obs_lat"][:] = self.lat[self.sort]
        nc.variables["obs_lon"][:] = self.lon[self.sort]
        nc.variables["obs_depth"][:] = self.depth[self.sort]
        nc.variables["obs_value"][:] = self.value[self.sort]
        nc.variables["obs_error"][:] = self.error[self.sort]
        nc.variables["obs_type"][:] = self.type[self.sort]
        nc.variables["obs_provenance"][:] = self.provenance[self.sort]
        nc.variables["obs_meta"][:] = self.meta[self.sort]
        nc.close()


def gridder(grid, time, lon, lat, depth, data, dt, depth_adjust=False,
            title='ROMS Observations'):
    """
    Construct an observations set from raw observations by placing them
    onto a grid.

    Parameters
    ----------
    grid : seapy.model.grid or filename string,
        Grid to place the raw observations onto
    time : ndarray,
        Time of the observations. This can be a scalar and all values
        will be assigned to the single time; otherwise, there must be a
        corresponding time to each value in the data.
    lon : ndarray,
        longitude of the observations. This can be a scalar and all values
        will be assigned to the single location; otherwise, there must be a
        corresponding longitude to each value in the data.
    lat : ndarray,
        latitude of the observations. This can be a scalar and all values
        will be assigned to the single location; otherwise, there must be a
        corresponding latitude to each value in the data.
    depth : ndarray or None,
        depth of the observations. If None, then all values are placed on
        the surface; otherwise, must be a corresponding depth for each
        value in the data.
    data : list of named tuples of seapy.roms.obs.raw_data,
        This list is comprised of each set of observation data types that
        are to be gridded together. If there is only one type (e.g.,
        SSH observations, there is only one item). An Argo float would have
        two items in the list (temperature and salinity observations).
        The list is comprised of named tuples of the raw observations
        with the following fields:
            "type" : string (or integer) of the type from
                     seapy.roms.obs.obs_types
             "provenance"  : string (or integer) of the type from
                             seapy.roms.obs.obs_provenance
            "values" : ndarray of actual observed values in units
                       for type
            "error" : ndarray (or None) of individual observational
                      uncertainty (same units of values). If not known,
                      use None
            "min_error" : float of the minimum error that should be
                          prescribed to the observations (typically,
                          the instrument error) in the same units of
                          values.
    dt : float
        The bin size of time for observations to be considered at the
        same time. The units must be the same as the provided time.
    title : string, optional,
        Title to assign the observations structure for output

    Returns
    -------
    obs : seapy.obs class
        Resulting observations from the raw data as placed onto grid.

    Examples
    --------
    A profile of temp and salt observations at a given lat/lon:

    >>> obs = seapy.obs.gridder(grid, times, lon, lat,
            [ seapy.roms.obs.raw_data("TEMP", "CTD_ARGO", temp, None, 0.1),
              seapy.roms.obs.raw_data("SALT", "CTD_ARGO", salt, None, 0.05)],
            dt = 1/24, title="Argo")

    Satellite Data from a number of lat/lons at a single time

    >>> obs = seapy.obs.gridder(grid, time, lon, lat,
            seapy.roms.obs.raw_data("ZETA", "SSH_AVISO", sla, sla_err, 0.05),
            dt = 2/24, title="SSH")

    These will generate new observation structures from the raw data.
    """
    from numpy_groupies import aggregate

    # Make sure the input is of the proper form
    grid = seapy.model.asgrid(grid)
    time = np.atleast_1d(time)
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)

    # First, before relying on gridding, extract only the data that are
    # encompassed by the grid
    region_list = np.where(np.logical_and.reduce((
        lat >= np.min(grid.lat_rho), lat <= np.max(grid.lat_rho),
        lon >= np.min(grid.lon_rho), lon <= np.max(grid.lon_rho))))
    if not np.any(region_list):
        warn("No observations were located within grid region_list")
        return None
    lat = lat[region_list]
    lon = lon[region_list]

    # Get the appropriate k-dimension depending on whether the data
    # are 2-D or 3-D
    if depth is None:
        # Get the grid locations from the data locations
        subsurface_values = False
        (j, i) = grid.ij((lon, lat))
        depth = grid.n * np.ones(i.size)
        k = np.ma.array(np.resize(grid.n, i.size))
    else:
        # Get the grid locations from the data locations
        subsurface_values = True
        depth = np.atleast_1d(depth)[region_list]
        (k, j, i) = grid.ijk((lon, lat, depth), depth_adjust)

    # Sub-select only the points that lie on our grid
    valid_list = np.where((~i.mask * ~j.mask * ~k.mask) == True)
    i = i[valid_list].compressed()
    j = j[valid_list].compressed()
    k = k[valid_list].compressed()
    depth = depth[valid_list]

    # Make sure the times are consistent and in dt-space
    if time.size == 1:
        time = np.resize(time, valid_list[0].size)
    else:
        time = time[region_list][valid_list]
    dtime = np.floor(time / dt)

    # Loop over all time intervals putting everything together. NOTE: The
    # preference is to use aggregate over the time-dimension just as we do
    # in the spatial-dimension; however, this led to crashing.
    ot = list()
    ox = list()
    oy = list()
    oz = list()
    odep = list()
    olon = list()
    olat = list()
    oval = list()
    oerr = list()
    oprov = list()
    otype = list()
    for t in seapy.progressbar.progress(np.unique(dtime)):
        time_list = np.where(dtime == t)
        mtime = np.nanmean(time[time_list])

        for v in data:
            valid_data = np.s_[:]
            if isinstance(v.values, np.ma.core.MaskedArray):
                valid_data = \
                    (v.values[region_list][valid_list][time_list].nonzero())[0]
                if not valid_data.size:
                    continue

            # Put together the indices based on the type of data we have
            if subsurface_values:
                idx = (k[time_list][valid_data],
                       j[time_list][valid_data],
                       i[time_list][valid_data])
            else:
                idx = (j[time_list][valid_data],
                       i[time_list][valid_data])
            indices = np.floor(idx).astype(int)

            # Grid the data onto our grid and compute the mean and variance
            ii = aggregate(indices, i[time_list][valid_data], func='mean')
            jj = aggregate(indices, j[time_list][valid_data], func='mean')
            binned = np.where(ii * jj > 0)
            ii = ii[binned].ravel()
            jj = jj[binned].ravel()
            (latl, lonl) = grid.latlon((ii, jj))
            Nd = ii.size

            # Put the co-located values together
            nvalues = aggregate(indices,
                                v.values[region_list][valid_list][
                                    time_list][valid_data],
                                func='mean')

            # Get their variance
            vari = aggregate(indices,
                             v.values[region_list][valid_list][
                                 time_list][valid_data],
                             func='var')

            # Put together the known observation values
            if v.error is not None:
                errs = aggregate(indices,
                                 v.error[region_list][valid_list][
                                     time_list][valid_data]**2,
                                 func='mean')
                errs = errs[binned].flatten()
            else:
                errs = 0.0

            # Build the depth vectors
            if subsurface_values:
                dd = aggregate(indices, depth[time_list][valid_data],
                               func='mean')
                kk = aggregate(indices, k[time_list][valid_data],
                               func='mean')
                dd = dd[binned].ravel()
                # ROMS counts from 1 for depth layers
                kk = kk[binned].ravel() + 1
            else:
                kk = np.resize(grid.n, Nd)
                dd = kk

            # Put all of the data from this time into our lists
            ot.append(np.resize(mtime, Nd))
            ox.append(ii)
            oy.append(jj)
            oz.append(kk)
            odep.append(dd)
            olon.append(lonl)
            olat.append(latl)
            oval.append(nvalues[binned].flatten())
            otype.append(np.resize(seapy.roms.obs.astype(v.type), Nd))
            oprov.append(np.resize(
                seapy.roms.obs.asprovenance(v.provenance), Nd))
            oerr.append(np.maximum(v.min_error**2,
                                   np.maximum(vari[binned].flatten(),
                                              errs)))

    # Make sure that we have something relevant
    if not oval:
        return None

    # Put everything together and create an observation class
    return seapy.roms.obs.obs(time=np.hstack(ot).ravel(),
                              x=np.hstack(ox).ravel(),
                              y=np.hstack(oy).ravel(),
                              z=np.hstack(oz).ravel(),
                              lat=np.hstack(olat).ravel(),
                              lon=np.hstack(olon).ravel(),
                              depth=np.hstack(odep).ravel(),
                              value=np.hstack(oval).ravel(),
                              error=np.hstack(oerr).ravel(),
                              type=np.hstack(otype).ravel(),
                              provenance=np.hstack(oprov).ravel(),
                              title=title)


def merge_files(obs_files, out_files, days, dt, limits=None, clobber=True):
    """
    merge together a group of observation files into combined new files
    with observations that lie only within the corresponding dates

    Parameters
    ----------
    obs_files : list,
        List of files to merge together (a single file will work, it will
        just be filtered by the dates)
    out_files : list or string,
        list of the filenames to create for each of the output periods.
        If a single string is given, the character '#' will be replaced
        by the starting time of the observation (e.g. out_files="out_#.nc"
        will become out_03234.nc)
    days : list of tuples,
        List of starting and ending day numbers for each cycle to process.
        The first value is the start day, the second is the end day. The
        number of tuples is the number of files to output.
    dt : float,
        Time separation of observations. Observations that are less than
        dt apart in time will be set to the same time.
    limits : dict, optional
        Set the limits of the grid points that observations are allowed
        within, {'north':i, 'south':i, 'east':i, 'west':i }. As obs near
        the boundaries are not advisable, this allows you to specify the
        valid grid range to accept obs within.
    clobber: bool, optional
        If True, output files are overwritten. If False, they are skipped.

    Returns
    -------
    None

    Examples
    --------

    Put together three files into 5 separate files in two day intervals from
    day 10 through day 20:

    >>> merge_files(["obs_1.nc", "obs_2.nc", "obs_3.nc"], "new_#.nc",
                   [(i, i+2) for i in range(10, 20, 2)])

    Put together same three files into 3 overlapping separate files in five
    day intervals with one overlapping day:

    >>> merge_files(["obs_1.nc", "obs_2.nc", "obs_3.nc"], "new_#.nc",
                   [(i, i+5) for i in range(10, 20, 4)])

    """
    import re
    import os

    # Only unique files
    obs_files = set().union(seapy.flatten(obs_files))
    outtime = False
    if isinstance(out_files, str):
        outtime = True
        time = re.compile('\#')

    # Go through the files to determine which periods they cover
    myobs = list()
    sdays = list()
    edays = list()
    for file in obs_files:
        nc = seapy.netcdf(file)
        fdays = nc.variables['survey_time'][:]
        nc.close()
        l = np.where(np.logical_and(fdays >= np.min(days),
                                    fdays <= np.max(days)))[0]
        if not l.size:
            continue
        myobs.append(file)
        sdays.append(fdays[0])
        edays.append(fdays[-1])
    sdays = np.asarray(sdays)
    edays = np.asarray(edays)

    # Loop over the dates in pairs
    for n, t in enumerate(seapy.progressbar.progress(days)):
        # Set output file name
        if outtime:
            outfile = time.sub("{:05d}".format(t[0]), out_files)
        else:
            outfile = out_files[n]

        if os.path.exists(outfile) and not clobber:
            continue

        # Find the files that cover the current period
        fidx = np.where(np.logical_and(sdays <= t[1], edays >= t[0]))[0]
        if not fidx.size:
            continue

        # Create new observations for this time period
        nobs = obs(myobs[fidx[0]])
        l = np.where(np.logical_or(nobs.time < t[0], nobs.time > t[1]))
        nobs.delete(l)
        for idx in fidx[1:]:
            o = obs(myobs[idx])
            l = np.where(np.logical_and(o.time >= t[0], o.time <= t[1]))
            nobs.add(o[l])
        # Remove any limits
        if limits is not None:
            l = np.where(np.logical_or.reduce((
                nobs.x < limits['west'],
                nobs.x > limits['east'],
                nobs.y < limits['south'],
                nobs.y > limits['north'])))
            nobs.delete(l)

        # Save out the new observations
        nobs.to_netcdf(outfile, dt=dt)

        pass
