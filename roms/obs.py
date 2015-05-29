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
from joblib import Parallel, delayed
import matplotlib.path

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
    114:"GLIDER_SG114",
    139:"GLIDER_SG139",
    146:"GLIDER_SG146",
    147:"GLIDER_SG147",
    148:"GLIDER_SG148",
    150:"GLIDER_SG500",
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
    411:"SSH_AVISO_ENVISAT",
    412:"SSH_AVISO_JASON1",
    413:"SSH_AVISO_JASON2",
    414:"SSH_AVISO_GFO",
    415:"SSH_HYCOM",
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
    # Search the dictionary for the key of the given the value
    return list(obs_types.keys())[ \
                list(obs_types.values()).index(s.upper())]

def _provenance_from_string(s):
    # Search the dictionary for the key of the given the value
    return list(obs_provenance.keys())[ \
                list(obs_provenance.values()).index(s.upper())]

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
            if time: self.time = np.atleast_1d(time)
            if x: self.x = np.atleast_1d(x)
            if y: self.y = np.atleast_1d(y)
            if z: self.z = np.atleast_1d(z)
            if lat: self.lat = np.atleast_1d(lat)
            if lon: self.lon = np.atleast_1d(lon)
            if depth: self.depth = np.atleast_1d(depth)
            if value: self.value = np.atleast_1d(value)
            if error: self.error = np.atleast_1d(error)
            if type: self._astype(type)
            if provenance: self._asprovenance(provenance)
            if meta: self.meta = np.atleast_1d(meta)

        # ensure everything is good
        self._consistent()

    def _astype(self, type):
        """
        PRIVATE method: build the type array from either string or values
        """
        if type is None:
            raise ValueError("you must provide a valid observation type")

        type=np.atleast_1d(type)
        # If the type given is a string, convert to ID
        if type.dtype.char is np.dtype(np.str).char:
            self.type = np.array([_type_from_string(s) for s in type])
        else:
            self.type = type.astype(np.int)

    def _asprovenance(self, prov):
        """
        PRIVATE method: build the provenance array from either string or values
        """
        if prov is None:
            prov=np.zeros(1)
        else:
            prov=np.atleast_1d(prov)

        # if the provenance given is a string, convert to ID
        if prov.dtype.char is np.dtype(np.str).char:
            self.provenance = np.array([_provenance_from_string(s) \
                                         for s in prov])
        else:
            self.provenance = prov.astype(np.int)

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
        self._astype(np.asanyarray(self.type).ravel())

        lt=len(self.time)
        if not lt==len(self.x)==len(self.y)==\
               len(self.value)==len(self.error)==len(self.type):
            # If these lengths are not equal, then there is a serious issue
            raise ValueError("Lengths of observation attributes are not equal.")
        else:
            # For the others, we can pad the information to ensure
            # consistency
            def _resizearr(arr,n):
                if not hasattr(self,arr):
                    return np.zeros(n)
                else:
                    arr=np.asanyarray(getattr(self,arr)).ravel()
                    la=len(arr)
                    if la==n:
                        return arr
                    elif la<n:
                        return np.pad(arr, (0,n-la), 'edge')
                    else:
                        return arr[:n]

            self.z = _resizearr("z", lt)
            self.lat = _resizearr("lat", lt)
            self.lon = _resizearr("lon", lt)
            self.depth = _resizearr("depth", lt)
            self._asprovenance(_resizearr("provenance", lt))
            self.meta = _resizearr("meta", lt)

        # Eliminate bad values
        good_vals=np.logical_and(np.isfinite(self.value), np.isfinite(self.x))
        good_vals=np.logical_and(good_vals,np.isfinite(self.y))
        good_vals=np.logical_and(good_vals,np.isfinite(self.error))
        good_vals=np.logical_and(good_vals,np.isfinite(self.time))
        if np.any(~good_vals):
            self.delete(np.where(good_vals==False))

    def __len__(self):
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
  "{:.2f}, [{:s}<{:s}] ({:.2f},{:.2f},{:.2f}) = {:.4f} +/- {:.4f}".format( \
                    t, obs_types[self.type[n]],
                    obs_provenance[self.provenance[n]],
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
        nc = seapy.roms.ncgen.create_da_obs(filename,
                survey=len(self.survey_time), provenance=obs_provenance,
                title=self.title)
        nc.variables["spherical"][:] = 1
        nc.variables["Nobs"][:] = self.nobs
        nc.variables["survey_time"][:] = self.survey_time
        nc.variables["obs_variance"][:] = np.ones(max(7,np.max(self.type)))*0.1
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


def gridder(grid, raw, dx=1, dy=1, dz=np.arange(0,5000,100), dt=1, nx=1, ny=1,
            threads=2):
    """
    Construct an observations set from raw observations by placing them
    onto a grid.

    Parameters
    ----------
    grid : seapy.model.grid or filename string,
        Grid to place the raw observations onto
    raw : dict,
        dictionary of the raw observations. Dict must have keys:
            "time" : list of floats [days since epoch]
            "lat"  : list of floats [degrees]
            "lon"  : list of floats [degrees]
            "depth": list of floats [m]
        Any other keys are considered data and compared against the
        obs_types dictionary.
    dx : float, optional,
        half-width of the x-dimension box to consider from the
        center of each grid cell in km (e.g., dx=4 would result in an
        8km width)
    dy : float, optional,
        half-width of the y-dimension box to consider from the
        center of each grid cell in km (e.g., dy=4 would result in an
        8km width)
    dz : list of floats, optional,
        The bins of vertical depths to consider for in situ observations.
        Positive values are down and in units of m.
    dt : float, optional,
        The bin size of time for observations to be considered at the
        same time in units of days.
    nx : int, optional,
        stride value of grid cells to consider in the x-direction (e.g.,
        1 is every point, 2 is every other point)
    ny : int, optional,
        stride value of grid cells to consider in the y-direction (e.g.,
        1 is every point, 2 is every other point)
    threads : int, optional,
        number of threads to run to grid the data

    Returns
    -------
    y : obs class
        Resulting observations from the raw data as placed onto grid.
    """
    # Make sure everything is in proper form.
    grid=seapy.model.asgrid(grid)
    raw_vars=[]
    for key in raw:
        raw[key] = np.asanyarray(raw[key]).ravel()
        if key.upper() in list(obs_types.values()):
            raw_vars.append(key)
    if not raw_vars:
        raise AttributeError("The raw data must contain a valid obs data type")

    l=np.where(raw["lon"]>180)
    raw["lon"][l]=raw["lon"][l]-360
    if "depth" in raw:
        raw["depth"]=np.sort(-np.abs(raw["depth"]))[::-1]

    # Construct the grid boxes to search for data
    rng=np.s_[::ny,::nx]
    ocean=np.where(grid.mask_rho[rng]==1)
    lon=grid.lon_rho[rng][ocean]
    lat=grid.lat_rho[rng][ocean]
    dx = dx*1000 / seapy.earth_distance(lon, lat, lon+1, lat)
    dy = dy*1000 / seapy.earth_distance(lon, lat, lon, lat+1)
    dxy = dx * np.sin(grid.angle[rng][ocean]);
    dx  = dx * np.cos(grid.angle[rng][ocean]);
    dyx = dy * np.sin(grid.angle[rng][ocean]);
    dy  = dy * np.cos(grid.angle[rng][ocean]);
    xv = np.vstack(((lon+dx-dyx).ravel(),
                    (lon+dx+dyx).ravel(),
                    (lon-dx+dyx).ravel(),
                    (lon-dx-dyx).ravel(),
                    (lon+dx-dyx).ravel())).T
    yv = np.vstack(((lat+dy+dxy).ravel(),
                    (lat-dy+dxy).ravel(),
                    (lat-dy-dxy).ravel(),
                    (lat+dy-dxy).ravel(),
                    (lat+dy+dxy).ravel())).T
    l=np.where(xv>180)
    xv[l]=xv[l]-360

    # Limit the data and grid boxes to those that are overlapping
    # First, eliminate raw data that are outside of our grid boxes
    region=np.where(np.logical_and( \
                    np.logical_and(raw["lon"] <= np.max(xv),
                                   raw["lon"] >= np.min(xv)),
                    np.logical_and(raw["lat"] <= np.max(yv),
                                   raw["lat"] >= np.min(yv))))
    for key in raw:
        if key is not "depth":
            raw[key]=raw[key][region]

    # Next, eliminate grid boxes that are outside of the raw data
    search=np.where(np.logical_and( \
                    np.logical_and(xv[:,0] <= np.max(raw["lon"]),
                                   xv[:,2] >= np.min(raw["lon"])),
                    np.logical_and(yv[:,0] <= np.max(raw["lat"]),
                                   yv[:,2] >= np.min(raw["lat"]))))
    xv=xv[search[0],:]
    yv=yv[search[0],:]

    # Build arrays that will be used many times
    pts=np.vstack((raw["lon"],raw["lat"])).T
    times=np.unique(raw["time"])
    d=np.diff(times)<dt
    while d.any():
        i=np.min(np.where(d))
        times[i+1]=times[i]
        times=np.unique(times)
        d=np.diff(times)<dt

    # Create an internal function for gridding the data
    def __gridder_thread(xv, yv, raw, pts, dt):
        # Find the points inside the current grid point
        poly=matplotlib.path.Path(list(zip(xv,yv)))
        inside=poly.contains_points(pts)
        if inside.any():
            pu.db
            # Loop over the available times
            time=np.mean(raw["time"][inside])
            lat=np.mean(raw["lat"][inside])
            lon=np.mean(raw["lon"][inside])
            num=np.count_nonzero(inside)
            z=1
            for v in raw_vars:
                data=np.nanmean(raw[v][inside])
                var=np.nanvar(raw[v][inside])
            pu.db
        return

    # Loop over the data, calling the internal gridder
    for i in range(xv.shape[0]):
        # pu.db
        __gridder_thread(xv[i,:],yv[i,:],raw,pts,dt)

    # Put the results together
    pass

 # ndata = np.ma.array(Parallel(n_jobs=threads,verbose=2)\
#                   (delayed(__interp2_thread) (
#    getattr(src_grid,"lon_"+grd), getattr(src_grid,"lat_"+grd),
#    ncsrc.variables[k][i,:,:],
#    getattr(child_grid,"lon_"+grd), getattr(child_grid,"lat_"+grd),
#    pmap["pmap"+grd], weight,
#    nx, ny, getattr(child_grid,"mask_"+grd))
#  for i in records), copy=False)
#
