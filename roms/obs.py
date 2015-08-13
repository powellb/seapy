#!/usr/bin/env python
"""
  obs.py

  State Estimation and Analysis for PYthon

  Module to handle the observation structure within ROMS

  Written by Brian Powell on 08/05/14
  Copyright (c)2014 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np
import netCDF4
import seapy
import matplotlib.path
from collections import namedtuple
from warnings import warn

# Define a named tuple to store raw data for the gridder
raw_data = namedtuple('raw_data', 'type provenance values error min_error')

# Define the observation type
obs_types = {
    1:"ZETA",
    2:"UBAR",
    3:"VBAR",
    4:"U",
    5:"V",
    6:"TEMP",
    7:"SALT",
    20:"RADIAL"
}
# Define the observation provenances used within my applications
obs_provenance = {
    0:"UNKNOWN",
    100:"GLIDER",
    114:"GLIDER_SG022",
    114:"GLIDER_SG023",
    114:"GLIDER_SG114",
    139:"GLIDER_SG139",
    146:"GLIDER_SG146",
    147:"GLIDER_SG147",
    148:"GLIDER_SG148",
    150:"GLIDER_SG500",
    151:"GLIDER_SG511",
    152:"GLIDER_SG512",
    153:"GLIDER_SG513",
    154:"GLIDER_SG523",
    200:"CTD",
    210:"CTD_HOT",
    220:"CTD_ARGO",
    300:"SST",
    301:"SST_OSTIA",
    317:"SST_AVHRR_17",
    318:"SST_AVHRR_18",
    330:"SST_MODIS_AQUA",
    331:"SST_MODIS_TERRA",
    400:"SSH",
    410:"SSH_AVISO_MAP",
    411:"SSH_AVISO_ENVISAT",
    412:"SSH_AVISO_JASON1",
    413:"SSH_AVISO_JASON2",
    414:"SSH_AVISO_GFO",
    415:"SSH_HYCOM",
    460:"SSS_AQUARIUS",
    500:"DRIFTERS",
    600:"RADAR",
    610:"RADAR_KOKOHEAD",
    620:"RADAR_KAKAAKO",
    630:"RADAR_KALAELOA",
    700:"ADCP",
    711:"ADCP_KILOMOANA_os38nb",
    712:"ADCP_KILOMOANA_os38bb",
    713:"ADCP_KILOMOANA_os75nb",
    714:"ADCP_KILOMOANA_os75bb",
    715:"ADCP_KILOMOANA_wh300",
    721:"ADCP_WECOMA_os38nb",
    722:"ADCP_WECOMA_os38bb",
    723:"ADCP_WECOMA_os75nb",
    724:"ADCP_WECOMA_os75bb",
    725:"ADCP_WECOMA_wh300",
    731:"ADCP_KOK_os38nb",
    732:"ADCP_KOK_os38bb",
    733:"ADCP_KOK_os75nb",
    734:"ADCP_KOK_os75bb",
    735:"ADCP_KOK_wh300",
    741:"ADCP_HIIALAKAI_os38nb",
    742:"ADCP_HIIALAKAI_os38bb",
    743:"ADCP_HIIALAKAI_os75nb",
    744:"ADCP_HIIALAKAI_os75bb",
    745:"ADCP_HIIALAKAI_wh300",
    751:"ADCP_THOMPSON_os38nb",
    752:"ADCP_THOMPSON_os38bb",
    753:"ADCP_THOMPSON_os75nb",
    754:"ADCP_THOMPSON_os75bb",
    755:"ADCP_THOMPSON_wh300"
}

def _type_from_string(s):
    """
    PRIVATE method: Search the type dictionary for the key of the
    given the value. If the key isn't a string or properly resolves, try to
    return the value as such
    """
    try:
        return list(obs_types.keys())[ \
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
        return list(obs_provenance.keys())[ \
                list(obs_provenance.values()).index(s.upper())]
    except AttributeError:
        return int(s)

def astype(type):
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
    type = np.atleast_1d(type)
    return np.array([_type_from_string(s) for s in type])

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
    return np.array([_provenance_from_string(s) for s in prov])

class obs:
    def __init__(self, filename=None, time=None, x=None, y=None, z=None,
                 lat=None, lon=None, depth=None, value=None, error=None,
                 type=None, provenance=None, meta=None,
                 title="ROMS Observations"):
        """
        Class to deal with ROMS observations for data assimilation

        Parameters
        ----------
        filename : string, optional,
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
        self.title=title
        if filename is not None:
            nc=netCDF4.Dataset(filename)
            # Construct an array from the data in the file. If obs_meta
            # exists in the file, then load it; otherwise, fill with zeros
            self.time = nc.variables["obs_time"][:]
            self.x = nc.variables["obs_Xgrid"][:]
            self.y = nc.variables["obs_Ygrid"][:]
            self.z = nc.variables["obs_Zgrid"][:]
            self.lat = nc.variables["obs_lat"][:]
            self.lon = nc.variables["obs_lon"][:]
            self.depth = nc.variables["obs_depth"][:]
            self.value = nc.variables["obs_value"][:]
            self.error = nc.variables["obs_error"][:]
            self.type = nc.variables["obs_type"][:].astype(np.int)
            self.provenance = nc.variables["obs_provenance"][:].astype(np.int)
            if "obs_meta" in nc.variables:
                self.meta = nc.variables["obs_meta"][:]
            else:
                self.meta = self.value.copy() * 0
            nc.close()
        else:
            if time is not None: self.time = np.atleast_1d(time)
            if x is not None: self.x = np.atleast_1d(x)
            if y is not None: self.y = np.atleast_1d(y)
            if z is not None: self.z = np.atleast_1d(z)
            if lat is not None: self.lat = np.atleast_1d(lat)
            if lon is not None: self.lon = np.atleast_1d(lon)
            if depth is not None: self.depth = np.atleast_1d(depth)
            if value is not None: self.value = np.atleast_1d(value)
            if error is not None: self.error = np.atleast_1d(error)
            if type is not None: self.type = astype(type)
            if provenance is not None:
                self.provenance = asprovenance(provenance)
            else:
                self.provenance = 0
            if meta is not None: self.meta = np.atleast_1d(meta)

        # ensure everything is good
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

        lt=len(self.time)
        if not lt==len(self.x)==len(self.y)==\
               len(self.value)==len(self.error)==len(self.type):
            # If these lengths are not equal, then there is a serious issue
            raise ValueError("Lengths of observation attributes are not equal.")
        else:
            # For the others, we can pad the information to ensure
            # consistency
            def _resizearr(arr,n):
                try:
                    return np.resize(getattr(self,arr),n)
                except:
                    return np.zeros(n)

            self.z = _resizearr("z", lt)
            self.lat = _resizearr("lat", lt)
            self.lon = _resizearr("lon", lt)
            self.depth = _resizearr("depth", lt)
            self.provenance = asprovenance(_resizearr("provenance", lt))
            self.meta = _resizearr("meta", lt)

        # Eliminate bad values
        good_vals=np.logical_and(np.isfinite(self.value), np.isfinite(self.x))
        good_vals=np.logical_and(good_vals,np.isfinite(self.y))
        good_vals=np.logical_and(good_vals,np.isfinite(self.error))
        good_vals=np.logical_and(good_vals,np.isfinite(self.time))
        if np.any(~good_vals):
            self.delete(np.where(good_vals==False))

        # Set the shape parameter
        self.shape = self.value.shape

    def __len__(self):
        self.shape = self.value.shape
        return len(self.value)

    def __getitem__(self, l):
        return obs(time=self.time[l], x=self.x[l], y=self.y[l],
                            z=self.z[l],lon=self.lon[l],lat=self.lat[l],
                            depth=self.depth[l],value=self.value[l],
                            error=self.error[l],type=self.type[l],
                            provenance=self.provenance[l],meta=self.meta[l])

    def __setitem__(self, l, new_obs):
        if not isinstance(new_obs,seapy.roms.obs.obs):
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
        return self.__str__()

    def __str__(self):
        return "\n".join([ \
  "{:.2f}, [{:s}:{:s}] ({:.2f},{:.2f},{:.2f}) = {:.4f} +/- {:.4f}".format( \
                    t, obs_types[self.type[n]],
                    obs_provenance.get(self.provenance[n],"UNKNOWN"),
                    self.lon[n], self.lat[n], self.depth[n],
                    self.value[n], self.error[n] ) \
            for n,t in enumerate(self.time) ])

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

    def create_survey(self):
        """
        Build the survey structure from the observations
        """
        # Generate the sort list
        self.sort = np.argsort(self.time,kind='mergesort')

        # Build the survey structure
        times, counts = np.unique(self.time[self.sort], return_counts=True)
        self.survey_time=times
        self.nobs=counts

    def to_netcdf(self, filename):
        """
        Write out the observations into the specified netcdf file

        Parameters
        ----------
        filename : string,
            name of file to save
        """
        # Save out the observations by survey
        self._consistent()
        self.create_survey()
        state_vars = np.maximum(7,np.max(self.type))
        nc = seapy.roms.ncgen.create_da_obs(filename,
                survey=len(self.survey_time), state_variable=state_vars,
                provenance=','.join((':'.join((obs_provenance[v],str(v)))
                                     for v in np.unique(self.provenance))),
                clobber=True, title=self.title)
        nc.variables["spherical"][:] = 1
        nc.variables["Nobs"][:] = self.nobs
        nc.variables["survey_time"][:] = self.survey_time
        nc.variables["obs_variance"][:] = np.ones(state_vars)*0.1
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


def gridder(grid, time, lon, lat, depth, data, dt, title='ROMS Observations'):
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

    **Examples**

    A profile of temp and salt observations at a given lat/lon:
    >>> obs = seapy.obs.gridder(grid, times, lon, lat,
            [ seapy.roms.obs.raw_data("TEMP", "CTD_ARGO", temp, None, 0.1),
              seapy.roms.obs.raw_data("SALT", "CTD_ARGO", salt, None, 0.05)],
            dt = 1/24, title="Argo")

    Satellite Data from a number of lat/lons at a single time
    >>> obs = seapy.obs.gridder(grid, time, lon, lat,
            seapy.roms.obs.raw_data("ZETA", "SSH_AVISO", sla, sla_err, 0.05),
            dt = 2/24, title="SSH")

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
                lat>=np.min(grid.lat_rho), lat<=np.max(grid.lat_rho),
                lon>=np.min(grid.lon_rho), lon<=np.max(grid.lon_rho))))
    if len(region_list[0]) == 0:
        warn("No observations were located within grid region_list")
        return None
    lat = lat[region_list]
    lon = lon[region_list]

    # Get the appropriate k-dimension depending on whether the data
    # are 2-D or 3-D
    if depth is None:
        # Get the grid locations from the data locations
        subsurface_values = False
        (j,i) = grid.ij((lon,lat))
        depth = np.zeros(len(i))
        k = np.ma.array(np.resize(grid.n, len(i)))
    else:
        # Get the grid locations from the data locations
        subsurface_values = True
        depth=np.atleast_1d(depth)[region_list]
        (k,j,i) = grid.ijk((lon,lat,depth))

    # Sub-select only the points that lie on our grid
    valid_list = np.where((~i.mask * ~j.mask * ~k.mask) == True)
    i = i.compressed()
    j = j.compressed()
    k = k[valid_list]
    depth = depth[valid_list]

    # Make sure the times are consistent and in dt-space
    if len(time) == 1:
        time = np.resize(time, len(valid_list[0]))
    else:
        time = time[region_list][valid_list]
    dtime = np.floor(time/dt)

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
                valid_data = ~np.ma.getmaskarray(
                                v.values[region_list][valid_list][time_list])

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
            (latl, lonl) = grid.latlon((ii,jj))
            Nd = len(ii)

            # Put the co-located values together
            nvalues = aggregate(indices,
                    v.values[region_list][valid_list][time_list][valid_data],
                    func='mean')
            # Get their variance
            vari = aggregate(indices,
                    v.values[region_list][valid_list][time_list][valid_data],
                    func='var')
            # Put together the known observation values
            if v.error is not None:
                errs = aggregate(indices,
                    v.error[region_list][valid_list][time_list][valid_data]**2,
                    func='mean')
                errs = errs[binned].flatten()
            else:
                errs = 0.0

            # Build the depth vectors
            if subsurface_values:
                dd = aggregate(indices, depth[time_list], func='mean')
                kk = aggregate(indices, k[time_list], func='mean')
                dd = dd[binned].ravel()
                kk = kk[binned].ravel()
            else:
                kk = np.resize(grid.n, Nd)
                dd = np.zeros(ii.shape)

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

    # Put everything together and create an observation class
    return seapy.roms.obs.obs(time=np.hstack(ot).ravel(),
                             x=np.hstack(ox).ravel(),
                             y=np.hstack(oy).ravel(),
                             z=np.hstack(odep).ravel(),
                             lat=np.hstack(olat).ravel(),
                             lon=np.hstack(olon).ravel(),
                             depth=np.hstack(oz).ravel(),
                             value=np.hstack(oval).ravel(),
                             error=np.hstack(oerr).ravel(),
                             type=np.hstack(otype).ravel(),
                             provenance=np.hstack(oprov).ravel(),
                             title=title)
