#!/usr/bin/env python
"""
  genobs.py

  State Estimation and Analysis for PYthon

  Module to process observations:
    obsgen : class to convert from raw to ROMS observations using
             specific subclasses

  Written by Brian Powell on 08/15/15
  Copyright (c)2015 University of Hawaii under the BSD-License.
"""

from __future__ import print_function

import numpy as np
import netCDF4
import seapy
import datetime
from warnings import warn

class obsgen(object):
    def __init__(self, grid, dt, epoch=None):
        """
        class for abstracting the processing of raw observation files
        (satellite, in situ, etc.) into ROMS observations files. All
        processing has commonalities which this class encapsulates, while
        leaving the loading and translation of individual data formats
        to subclasses.

        Parameters
        ----------
        grid: seapy.model.grid or string,
            grid to use for generating observations
        dt: float,
            Model time-step or greater in units of days
        epoch: datetime, optional,
            Time to reference all observations from

        Returns
        -------
        None

        """
        self.grid = seapy.model.asgrid(grid)
        self.dt = dt
        self.epoch = datetime.datetime(2000,1,1) if epoch is None else epoch

    def convert_file(self, file, title=None):
        """
        convert a raw observation file into a ROMS observations structure.
        The subclasses are responsible for the conversion, and this method
        is obsgen is only a stub.

        Parameters
        ----------
        file : string,
            filename of the file to process
        title : string,
            Title to give the new observation structure global attribute

        Returns
        -------
        seapy.roms.obs.obs,
            observation structure from raw obs
        """
        pass

    def batch_files(self, in_files, out_files):
        """
        Given a list of input files, process each one and save each result
        into the given output file.

        Parameters
        ----------
        in_files : list of strings,
            filenames of the files to process
        out_files : list of strings,
            filenames of the files to create for each of the input filenames.
            If a single string is given, the character '#' will be replaced
            by the starting time of the observation (e.g. out_files="out_#.nc"
            will become out_03234.nc)

        Returns
        -------
        None
        """
        import re

        outtime = False
        if isinstance(out_files, str):
            outtime = True
            time = re.compile('\#')

        for n,file in enumerate(in_files):
            try:
                print(file)
                obs = self.convert_file(file)
                if outtime:
                    obs.to_netcdf(time.sub("{:05d}".format(int(obs.time[0])),
                                           out_files))
                else:
                    obs.to_netcdf(out_files[n])
            except:
                warn("WARNING: "+file+" cannot be processed.")
        pass

class aquarius_sss(obsgen):
    """
    class to process Aquarius SSS HDF5 files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """
    def __init__(self, grid, dt, epoch=None, salt_limits=(10,36)):
        self.salt_limits = salt_limits
        super().__init__(grid, dt, epoch)

    def convert_file(self, file, title="AQUARIUS Obs"):
        """
        Load an Aquarius file and convert into an obs structure
        """
        import h5py

        f = h5py.File(file,'r')
        salt = np.ma.masked_equal(np.flipud(f['l3m_data'][:]),
                                  f['l3m_data'].attrs['_FillValue'])
        year = f.attrs['Period End Year']
        day = f.attrs['Period End Day']
        nlat = f.attrs['Northernmost Latitude']-0.5
        slat = f.attrs['Southernmost Latitude']+0.5
        wlon = f.attrs['Westernmost Longitude']+0.5
        elon = f.attrs['Easternmost Longitude']-0.5
        dlat = f.attrs['Latitude Step']
        dlon = f.attrs['Longitude Step']
        f.close()

        [lon, lat] = np.meshgrid(np.arange(wlon,elon+dlon,dlon),
                                 np.arange(slat,nlat+dlat,dlat))
        time = (datetime.datetime(year,1,1) + datetime.timedelta(int(day)) -
               self.epoch).days
        lat = lat.flatten()
        lon = lon.flatten()
        if np.min(self.grid.lon_rho > 0):
            lon[lon<0] += 360

        salt = np.ma.masked_outside(salt.flatten(), self.salt_limits[0],
                                    self.salt_limits[1])
        data = [seapy.roms.obs.raw_data("SALT", "SSS_AQUARIUS",
                                        salt, None, 0.1)]
        # Grid it
        return seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                      data, self.dt, title)
        pass

class argo_ctd(obsgen):
    """
    class to process ARGO CTD netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """
    def __init__(self, grid, dt, epoch=None):
        super().__init__(grid, dt, epoch)

    def convert_file(self, file, title="AVISO Obs"):
        """
        Load an Argo file and convert into an obs structure
        """
        nc = netCDF4.Dataset(file)

        # Load the position of all profiles in the file
        lon = nc.variables["LONGITUDE"][:]
        lat = nc.variables["LATITUDE"][:]
        pro_q = nc.variables["POSITION_QC"][:].astype(int)
        # Find the profiles that are in our area with known locations quality
        if np.min(self.grid.lon_rho < 0):
            lon[lon>180] -= 360
        profile_list = np.where(np.logical_and.reduce((
                    lat>=np.min(self.grid.lat_rho),
                    lat<=np.max(self.grid.lat_rho),
                    lon>=np.min(self.grid.lon_rho),
                    lon<=np.max(self.grid.lon_rho),
                    pro_q==1)))[0]

        # Check which are good profiles
        profile_qc = nc.variables["PROFILE_PRES_QC"][profile_list].astype('<U1')
        profile_list = profile_list[profile_qc == 'A']
        if not len(profile_list):
            return None

        # Load only the data from those in our area
        lon = lon[profile_list]
        lat = lat[profile_list]
        julian_day = nc.variables["JULD_LOCATION"][profile_list]
        argo_epoch = datetime.datetime.strptime(''.join( \
            nc.variables["REFERENCE_DATE_TIME"][:].astype('<U1')),'%Y%m%d%H%M%S')
        temp = nc.variables["TEMP"][profile_list,:]
        temp_qc = nc.variables["TEMP_QC"][profile_list,:]
        salt = nc.variables["PSAL"][profile_list,:]
        salt_qc = nc.variables["PSAL_QC"][profile_list,:]
        pres = nc.variables["PRES"][profile_list,:]
        pres_qc = nc.variables["PRES_QC"][profile_list,:]
        nc.close()

        # Combine the QC codes
        qc = np.mean(np.vstack((temp_qc.compressed(), salt_qc.compressed(),
                                pres_qc.compressed())).astype(int),axis=0)
        good_data = np.where(qc==1)

        # Put everything together into individual observations
        time_delta = (self.epoch - argo_epoch).days
        time = np.resize(julian_day-time_delta,
                         pres.shape[::-1]).T[~temp.mask][good_data]
        lat = np.resize(lat, pres.shape[::-1]).T[~temp.mask][good_data]
        lon = np.resize(lon, pres.shape[::-1]).T[~temp.mask][good_data]
        depth = -seapy.seawater.depth(pres.compressed()[good_data], lat)
        temp = np.ma.masked_outside(temp.compressed()[good_data], 2, 30)
        salt = np.ma.masked_outside(salt.compressed()[good_data], 10, 36)

        data = [ seapy.roms.obs.raw_data("TEMP", "CTD_ARGO", temp, None, 0.25),
                 seapy.roms.obs.raw_data("SALT", "CTD_ARGO", salt, None, 0.1)]

        return seapy.roms.obs.gridder(self.grid, time, lon, lat, depth,
                                      data, self.dt, title)

class aviso_sla_map(obsgen):
    """
    class to process AVISO SLA map netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """
    def __init__(self, grid, dt, epoch=None, sshmean=None):
        if sshmean is not None:
            self.sshmean = seapy.convolve_mask(sshmean, ksize=5, copy=True)
        else:
            self.sshmean = None
        super().__init__(grid, dt, epoch)

    def convert_file(self, file, title="AVISO Obs"):
        """
        Load an AVISO file and convert into an obs structure
        """
        # Load AVISO Data
        nc = netCDF4.Dataset(file)
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]
        dat = np.squeeze(nc.variables["sla"][:])
        err = np.squeeze(nc.variables["err"][:])
        time = netCDF4.num2date(nc.variables["time"][0],
                                nc.variables["time"].units) - self.epoch
        time = time.total_seconds()/86400
        nc.close()
        lon, lat = np.meshgrid(lon, lat)
        lat = lat.flatten()
        lon = lon.flatten()
        if np.min(self.grid.lon_rho < 0):
            lon[lon>180] -= 360
        data = [seapy.roms.obs.raw_data("ZETA", "SSH_AVISO_MAP",
                                        dat.flatten(), err.flatten(), 0.05)]
        # Grid it
        obs = seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                     data, self.dt, title)

        # Apply the model mean ssh to the sla data
        if self.sshmean is not None:
            m,p = seapy.oasurf(self.grid.I, self.grid.J, self.sshmean,
                               obs.x, obs.y, nx=1, ny=1, weight=7)
            obs.value += m
        return obs

class ostia_sst_map(obsgen):
    """
    class to process OSTIA SST map netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """
    def __init__(self, grid, dt, epoch=None):
        super().__init__(grid, dt, epoch)

    def convert_file(self, file, title="OSTIA SST Obs"):
        """
        Load an OSTIA file and convert into an obs structure
        """
        # Load OSTIA Data
        nc = netCDF4.Dataset(file)
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]
        dat = np.squeeze(nc.variables["analysed_sst"][:]) - 273.15
        err = np.squeeze(nc.variables["analysis_error"][:])
        time = netCDF4.num2date(nc.variables["time"][0],
                                nc.variables["time"].units) - self.epoch
        time = time.total_seconds()/86400
        nc.close()
        lon, lat = np.meshgrid(lon, lat)
        lat = lat.flatten()
        lon = lon.flatten()
        if np.min(self.grid.lon_rho < 0):
            lon[lon>180] -= 360
        data = [seapy.roms.obs.raw_data("TEMP", "SST_OSTIA", dat.flatten(),
                                        err.flatten(), 0.4)]
        # Grid it
        return seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                     data, self.dt, title)

class seaglider_profile(obsgen):
    """
    class to process SeaGlider .pro files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """
    def __init__(self, grid, dt, epoch=None, temp_limits=(5,30),
                 salt_limits=(33,35.5), depth_limit=-15):
        self.temp_limits = temp_limits
        self.salt_limits = salt_limits
        self.depth_limit = depth_limit
        super().__init__(grid, dt, epoch)

    def convert_file(self, file, title="SeaGlider Obs"):
        """
        Load a SeaGlider .pro file and convert into an obs structure
        """
        import re

        dtype = { 'names': ('time','pres','depth','temp','cond',
                            'salt','sigma','lat','lon'),
                  'formats': ['f4']*9 }

        # Load the text file. All data goes into the pro dictionary
        # as defined by dtype. The header information needs to be parsed
        with open(file) as myfile:
            header = [ myfile.readline() for i in range(19) ]
            pro = np.loadtxt(myfile, dtype, delimiter=',', comments='%')

        # Parse the header information
        parser = re.compile('^%(\w+): (.*)$')
        params = {}
        for line in header:
            try:
                opt = parser.findall(line)
                params[opt[0][0]] = opt[0][1]
            except:
                pass

        # Determine the needed information from the headers
        glider_name = "GLIDER" if params.get("glider", None) is None else \
                      "GLIDER_SG"+params["glider"]
        provenance = seapy.roms.obs.asprovenance(glider_name)
        try:
            date = [ int(s) for s in re.findall('([\d]{2})\s', params["start"]) ]
            start_time = datetime.datetime.strptime(params["start"].strip(),
                                                    "%m %d 1%y %H %M %S")
            dtime = (start_time - self.epoch).total_seconds()/86400
        except:
            raise ValueError("date format incorrect in file: "+file)

        # Make sure that the GPS fix isn't screwy
        if np.min(self.grid.lon_rho > 0):
            pro["lon"][pro["lon"]<0] += 360
        dist = seapy.earth_distance(pro["lon"][0], pro["lat"][0],
                                    pro["lon"][-1], pro["lat"][-1])
        velocity = dist / pro["time"][-1]
        if velocity > 2:
            warn("WARNING: GPS fix is incorrect for "+file)
            return None

        # Build the data with masked entries
        temp = np.ma.masked_outside(pro["temp"], self.temp_limits[0],
                                    self.temp_limits[1])
        salt = np.ma.masked_outside(pro["salt"], self.salt_limits[0],
                                    salt.temp_limits[1])
        depth = np.ma.masked_greater(-pro["depth"], self.depth_limit)

        # Grid it
        data = [ seapy.roms.obs.raw_data("TEMP", provenance, temp,
                                         None, 0.2),
                 seapy.roms.obs.raw_data("SALT", provenance, salt,
                                         None, 0.05)]
        return seapy.roms.obs.gridder(self.grid, pro["time"]/86400+dtime,
                                      pro["lon"], pro["lat"], depth,
                                      data, self.dt, title)

