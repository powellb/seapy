#!/usr/bin/env python
"""
  obsgen.py

  State Estimation and Analysis for PYthon

  Module to process observations:
    obsgen : class to convert from raw to ROMS observations using
             specific subclasses

  Written by Brian Powell on 08/15/15
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import numpy as np
import netCDF4
import h5py
import seapy
import datetime
from warnings import warn


def error_profile(obs, depth, error, provenance=None):
    """
    Apply a vertical error profile to a given observation structure.
    This allows for error minimums to vary by depth and observation
    type.

    Parameters
    ----------
    obs : seapy.roms.obs.obs or string,
      The observations to enforce the error profile upon.
    depth : ndarray,
      Array of depths for the errors provided
    error : dict,
      Dictionary of the errors, where the key is the type of observation
      (as defined by seapy.roms.obs.obs_types) and the value is
      an ndarray of same length as depth with the error [in squared units]
      of the observation profile.
    provenance : list of int or string, optional,
      The provenance to apply the errors to (ignore other observations
      of the same type, but different instrument)

    Returns
    -------
    None:
      The obs structure is mutable is changed in place

    Examples
    --------
    >>> obs = obs('observation_file.nc')
    >>> depth = [10, 30, 50, 1000, 2000]
    >>> error['temp'] = [0.5, 0.2, 0.4, 0.1, 0.01]
    >>> error_profile(obs, depth, error)

    The resulting 'obs' class will have had its error profiles
    modified.
    """
    from scipy.interpolate import interp1d
    obs = seapy.roms.obs.asobs(obs)
    depth = np.atleast_1d(depth).flatten()
    depth = np.abs(depth)
    pro = seapy.roms.obs.asprovenance(provenance) if provenance else None

    # Loop over all of the profiles in the error dictionary and
    # apply them to the observations
    for var in error:
        typ = seapy.roms.obs.astype(var)
        try:
            fint = interp1d(depth, error[var].flatten(), copy=False)
            if pro.any():
                l = np.where(np.logical_and(obs.type == typ,
                                            np.in1d(obs.provenance, pro)))
            else:
                l = np.where(np.logical_and(obs.type == typ, obs.depth < 0))
            nerr = fint(np.abs(obs.depth[l]))
            obs.error[l] = np.maximum(obs.error[l], nerr)
        except ValueError:
            warn("Error for {:s} is the wrong size".format(var))
            continue
    pass


def add_ssh_tides(obs, tide_file, tide_error, tide_start=None, provenance=None,
                  reftime=seapy.default_epoch):
    """
    Apply predicted barotropic tides to the SSH values of given observations
    using the tide_file given.

    Parameters
    ----------
    obs : seapy.roms.obs.obs or string,
      The observations to enforce the error profile upon.
    tide_file : string,
      The name of the ROMS tidal forcing file to use for predicting the
      barotropic tides.
    tide_error : np.masked_array
      A two dimensional array of the tidal fit errors to apply to
      the ssh errors when adding the tides. This should be the same size
      as the rho-grid. The units of the error must be in meters. If it is
      masked, the mask will be honored and obs that are in the mask will
      be removed. This allows you to filter on regions of high error.
    tide_start : bool, optional,
      If given, the tide_start of the tide file. If not specified,
      will read the attribute of the tidal forcing file
    provenance : list of int or string, optional,
      The provenance to apply the tides to (ignore other observations
      of the same type, but different instrument)
    reftime: datetime,
      Reference time for the observation times

    Returns
    -------
    None:
      The obs structure is mutable is changed in place

    Examples
    --------
    >>> obs = obs('observation_file.nc')
    >>> add_ssh_tides(obs, 'tide_frc.nc', errmap)

    The resulting 'obs' variable will have modified data. To save it:

    >>> obs.to_netcdf()
    """

    # Load tidal file data
    frc = seapy.roms.tide.load_forcing(tide_file)
    if not tide_start:
        tide_start = frc['tide_start']

    # Make sure that the sizes are the same
    if frc['Eamp'].shape[1:] != tide_error.shape:
        raise ValueError(
            "The error array is not the same size as the tidal grid")

    # Gather the observations that need tidal information
    obs = seapy.roms.obs.asobs(obs)
    pro = seapy.roms.obs.asprovenance(provenance) if provenance else None
    if pro:
        l = np.where(np.logical_and(obs.type == 1,
                                    np.in1d(obs.provenance, pro)))
    else:
        l = np.where(obs.type == 1)

    # If we have any, then do tidal predictions and add the signal
    # and error to the observations
    bad = []
    if l[0].any():
        ox = np.rint(obs.x[l]).astype(int)
        oy = np.rint(obs.y[l]).astype(int)
        idx = seapy.unique_rows((ox, oy))
        for cur in seapy.progressbar.progress(idx):
            pts = np.where(np.logical_and(ox == ox[cur], oy == oy[cur]))
            # If this point is masked, remove from the observations
            if not tide_error[oy[cur], ox[cur]]:
                bad.append(l[0][pts].tolist())
            else:
                time = [reftime + datetime.timedelta(t) for t in
                        obs.time[l][pts]]
                amppha = seapy.tide.pack_amp_phase(
                    frc['tides'], frc['Eamp'][:, oy[cur], ox[cur]],
                    frc['Ephase'][:, oy[cur], ox[cur]])
                zpred = seapy.tide.predict(time, amppha,
                                           lat=obs.lat[l][cur],
                                           tide_start=tide_start)
                # Add the information to the observations
                obs.value[l[0][pts]] += zpred
                obs.error[l[0][pts]] = np.maximum(
                    obs.error[l[0][pts]], tide_error[oy[cur], ox[cur]]**2)

    # If any were bad, then remove them
    if bad:
        obs.delete(seapy.flatten(bad))
    pass


class obsgen(object):

    def __init__(self, grid, dt, reftime=seapy.default_epoch):
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
        self.epoch = reftime

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

    def datespan_file(self, file):
        """
        check the given file and return the date span of data that are
        covered by the file.

        Parameters
        ----------
        file : string,
            filename of the file to process

        Returns
        -------
        start : datetime
            starting date and time of the data
        end : datetime
            ending date and time of the data
        """

        return None, None

    def batch_files(self, in_files, out_files, start_time=None,
                    end_time=None, clobber=True):
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
        start_time : datetime, optional
            starting date and time for data to process (ignore files that are
            outside of the time period)
        end_time : datetime, optional
            ending date and time for data to process (ignore files that are
            outside of the time period). If start_time is provided, and
            end_time is not, then a period of one day is assumed.
        clobber : bool, optional
            If TRUE, overwrite any existing output files. If False, the
            file is given a letter suffix.

        Returns
        -------
        None
        """
        import re
        import os

        datecheck = False
        if start_time is not None:
            datecheck = True
            if end_time is None:
                end_time = start_time + datetime.timedelta(1)

        outtime = False
        if isinstance(out_files, str):
            outtime = True
            time = re.compile('\#')

        for n, file in enumerate(in_files):
            try:
                # Check the times if user requested
                print(file, end="")
                if datecheck:
                    st, en = self.datespan_file(file)
                    if (en is not None and en < start_time) or \
                            (st is not None and st > end_time):
                        print(": SKIPPED")
                        continue

                # Convert the file
                obs = self.convert_file(file)
                if obs is None:
                    print(": NO OBS")
                    continue

                # Output the obs to the correct file
                if outtime:
                    ofile = time.sub("{:05d}".format(int(obs.time[0])),
                                     out_files)
                else:
                    ofile = out_files[n]

                if clobber:
                    obs.to_netcdf(ofile, True)
                else:
                    for i in "abcdefgh":
                        if os.path.isfile(ofile):
                            ofile = re.sub("[a-h]{0,1}\.nc", i + ".nc", ofile)
                        else:
                            break
                    obs.to_netcdf(ofile, False)
                print(": SAVED")

            except (BaseException, UserWarning) as e:
                warn("WARNING: {:s} cannot be processed.\nError: {:}".format(
                    file, e.args))
        pass


##############################################################################
#
# REMOTE-SENSING DATA
#
##############################################################################

class aquarius_sss(obsgen):
    """
    class to process Aquarius SSS HDF5 files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, salt_limits=None,
                 salt_error=0.2):
        if salt_limits is None:
            self.salt_limits = (10, 36)
        else:
            self.salt_limits = salt_limits
        self.salt_error = salt_error
        super().__init__(grid, dt, reftime)

    def datespan_file(self, file):
        f = h5py.File(file, 'r')
        try:
            year = f.attrs['Period End Year']
            day = f.attrs['Period End Day']
            st = datetime.datetime(year, 1, 1) + datetime.timedelta(int(day))
            en = st + datetime.timedelta(1)
        except:
            st = en = None
            pass
        finally:
            f.close()
            return st, en

    def convert_file(self, file, title="AQUARIUS Obs"):
        """
        Load an Aquarius file and convert into an obs structure
        """
        f = h5py.File(file, 'r')
        salt = np.ma.masked_equal(np.flipud(f['l3m_data'][:]),
                                  f['l3m_data'].attrs['_FillValue'])
        year = f.attrs['Period End Year']
        day = f.attrs['Period End Day']
        nlat = f.attrs['Northernmost Latitude'] - 0.5
        slat = f.attrs['Southernmost Latitude'] + 0.5
        wlon = f.attrs['Westernmost Longitude'] + 0.5
        elon = f.attrs['Easternmost Longitude'] - 0.5
        dlat = f.attrs['Latitude Step']
        dlon = f.attrs['Longitude Step']
        f.close()

        [lon, lat] = np.meshgrid(np.arange(wlon, elon + dlon, dlon),
                                 np.arange(slat, nlat + dlat, dlat))
        time = (datetime.datetime(year, 1, 1) + datetime.timedelta(int(day)) -
                self.epoch).days
        lat = lat.flatten()
        lon = lon.flatten()
        if self.grid.east():
            lon[lon < 0] += 360

        salt = np.ma.masked_outside(salt.flatten(), self.salt_limits[0],
                                    self.salt_limits[1])
        data = [seapy.roms.obs.raw_data("SALT", "SSS_AQUARIUS",
                                        salt, None, self.salt_error)]
        # Grid it
        return seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                      data, self.dt, title)
        pass


class aviso_sla_map(obsgen):
    """
    class to process AVISO SLA map netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, ssh_mean=None,
                 ssh_error=0.05):
        if ssh_mean is not None:
            self.ssh_mean = seapy.convolve_mask(ssh_mean, ksize=5, copy=True)
        else:
            self.ssh_mean = None
        self.ssh_error = ssh_error
        super().__init__(grid, dt, reftime)

    def datespan_file(self, file):
        nc = seapy.netcdf(file)
        try:
            st = datetime.datetime.strptime(nc.getncattr("time_coverage_start"),
                                            "%Y-%m-%dT%H:%M:%SZ")
            en = datetime.datetime.strptime(nc.getncattr("time_coverage_end"),
                                            "%Y-%m-%dT%H:%M:%SZ")
        except:
            st = en = None
            pass
        finally:
            nc.close()
            return st, en

    def convert_file(self, file, title="AVISO Obs"):
        """
        Load an AVISO file and convert into an obs structure
        """
        # Load AVISO Data
        nc = seapy.netcdf(file)
        lonname = 'lon' if 'lon' in nc.variables.keys() else 'longitude'
        lon = nc.variables[lonname][:]
        latname = 'lat' if 'lat' in nc.variables.keys() else 'latitude'
        lat = nc.variables[latname][:]
        dat = np.squeeze(nc.variables["sla"][:])
        err = np.squeeze(nc.variables["err"][:])
        time = seapy.roms.get_time(
            nc, "time", records=[0], epoch=self.epoch)[0]
        nc.close()
        lon, lat = np.meshgrid(lon, lat)
        lat = lat.flatten()
        lon = lon.flatten()
        if not self.grid.east():
            lon[lon > 180] -= 360
        data = [seapy.roms.obs.raw_data("ZETA", "SSH_AVISO_MAP",
                                        dat.flatten(), err.flatten(), self.ssh_error)]
        # Grid it
        obs = seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                     data, self.dt, title)

        # Apply the model mean ssh to the sla data
        if self.ssh_mean is not None:
            m, p = seapy.oasurf(self.grid.I, self.grid.J, self.ssh_mean,
                                obs.x, obs.y, nx=1, ny=1, weight=7)
            obs.value += m
        return obs


_aviso_sla_errors = {
    "SSH_AVISO_ENVISAT": 0.06,
    "SSH_AVISO_JASON1": 0.05,
    "SSH_AVISO_JASON2": 0.05,
    "SSH_AVISO_JASON3": 0.05,
    "SSH_AVISO_GFO": 0.05,
    "SSH_AVISO_ALTIKA": 0.07,
    "SSH_AVISO_CRYOSAT2": 0.07,
    "SSH_AVISO_HAIYANG": 0.07,
    "SSH_AVISO_ERS1": 0.06,
    "SSH_AVISO_ERS2": 0.06,
    "SSH_AVISO_TOPEX_POSEIDON": 0.05,
    "SSH_AVISO_SENTINEL3A": 0.05
}


class aviso_sla_track(obsgen):
    """
    class to process AVISO SLA track netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data. THIS COVERS ALL SATELLITES/INSTRUMENTS FROM AVISO TRACK:
    al, c2, e1, e2, en, enn, g2, h2, j1, j1g, j1n, j2, tp and tpn.

    Parameters
    ----------
    ssh_mean : ndarray,
      Spatial map of rho-grid shape that contains the model mean SSH
    ssh_error: dict, optional
      Dictionary of the minimum errors for each satellite. The default
      uses the errors defined in _aviso_sla_errors
    repeat: int
      Number of hours to repeat the track before and after its initial
      pass
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, ssh_mean=None,
                 ssh_error=None, repeat=3, provenance="SSH"):
        self.provenance = provenance.upper()
        self.repeat = repeat
        self.ssh_error = ssh_error if ssh_error else _aviso_sla_errors
        if ssh_mean is not None:
            self.ssh_mean = seapy.convolve_mask(ssh_mean, ksize=5, copy=True)
        else:
            self.ssh_mean = None
        super().__init__(grid, dt, reftime)

    def convert_file(self, file, title="AVISO SLA Track Obs"):
        """
        Load an AVISO file and convert into an obs structure
        """
        # Load AVISO Data
        nc = seapy.netcdf(file)
        lon = nc.variables["longitude"][:]
        lat = nc.variables["latitude"][:]
        slaname = 'SLA' if 'SLA' in nc.variables.keys() else 'sla_filtered'
        dat = nc.variables[slaname][:]
        time = seapy.roms.num2date(nc, "time", epoch=self.epoch)
        nc.close()

        # make them into vectors
        lat = lat.ravel()
        lon = lon.ravel()
        dat = dat.ravel()
        err = np.ones(dat.shape) * _aviso_sla_errors.get(self.provenance, 0.1)

        if not self.grid.east():
            lon[lon > 180] -= 360

        good = dat.nonzero()
        data = [seapy.roms.obs.raw_data("ZETA", self.provenance,
                                        dat[good], err[good], err[0])]
        # Grid it
        obs = seapy.roms.obs.gridder(self.grid, time, lon[good], lat[good], None,
                                     data, self.dt, title)

        # Apply the model mean ssh to the sla data
        if self.ssh_mean is not None and obs is not None:
            m, p = seapy.oasurf(self.grid.I, self.grid.J, self.ssh_mean,
                                obs.x, obs.y, nx=1, ny=1, weight=7)
            obs.value += m

        # Duplicate the observations before and after as per the repeat
        # time unless it is zero
        if self.repeat and obs:
            prior = obs.copy()
            after = obs.copy()
            prior.time -= self.repeat / 24
            after.time += self.repeat / 24
            obs.add(prior)
            obs.add(after)

        return obs


class ostia_sst_map(obsgen):
    """
    class to process OSTIA SST map netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, temp_error=0.4,
                 temp_limits=None):
        self.temp_error = temp_error
        if temp_limits is None:
            self.temp_limits = (2, 35)
        else:
            self.temp_limits = temp_limits
        super().__init__(grid, dt, reftime)

    def convert_file(self, file, title="OSTIA SST Obs"):
        """
        Load an OSTIA file and convert into an obs structure
        """
        # Load OSTIA Data
        nc = seapy.netcdf(file)
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]
        dat = np.ma.masked_outside(np.squeeze(
            nc.variables["analysed_sst"][:]) - 273.15,
            self.temp_limits[0], self.temp_limits[1])
        err = np.ma.masked_outside(np.squeeze(
            nc.variables["analysis_error"][:]), 0.01, 2.0)
        dat[err.mask] = np.ma.masked
        time = seapy.roms.num2date(
            nc, "time", records=[0], epoch=self.epoch)[0]
        nc.close()
        if self.grid.east():
            lon[lon < 0] += 360
        lon, lat = np.meshgrid(lon, lat)
        good = dat.nonzero()
        lat = lat[good]
        lon = lon[good]
        data = [seapy.roms.obs.raw_data("TEMP", "SST_OSTIA", dat.compressed(),
                                        err[good], self.temp_error)]
        # Grid it
        return seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                      data, self.dt, title)


class navo_sst_map(obsgen):
    """
    class to process NAVO SST map netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """

    def __init__(self, grid, dt, depth=None, reftime=seapy.default_epoch,
                 temp_error=0.25, temp_limits=None, provenance="SST_NAVO_MAP"):

        self.temp_error = temp_error
        self.provenance = provenance.upper()
        self.temp_limits = (2, 35) if temp_limits is None else temp_limits
        self.depth = 4 if depth is None else np.abs(depth)
        super().__init__(grid, dt, reftime)

    def datespan_file(self, file):
        nc = seapy.netcdf(file)
        try:
            st = datetime.datetime.strptime(nc.getncattr("start_date"),
                                            "%Y-%m-%d UTC")
            en = datetime.datetime.strptime(nc.getncattr("stop_date"),
                                            "%Y-%m-%d UTC")
        except:
            st = en = None
            pass
        finally:
            nc.close()
            return st, en

    def convert_file(self, file, title="NAVO SST Obs"):
        """
        Load a NAVO map file and convert into an obs structure
        """
        import re
        import sys

        nc = seapy.netcdf(file)
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]
        dat = np.ma.masked_outside(np.squeeze(nc.variables["analysed_sst"][:]) - 273.15,
                                   self.temp_limits[0], self.temp_limits[1])
        err = np.ma.array(np.squeeze(
            nc.variables["analysis_error"][:]), mask=dat.mask)

        # this is an analyzed product and provides errors as a function
        # of space and time directly the temperature is the bulk
        # temperature (ie at around 4m depth, below the e-folding depths of
        # sunlight in the ocean so the product does not have a diuranl cycle
        # (ie you don;t have to worry about hourly variations)
        time = seapy.roms.num2date(
            nc, "time", records=[0], epoch=self.epoch)[0]
        nc.close()

        # here we set the depth to be 4 m below the surface
        if self.grid.east():
            lon[lon < 0] += 360
        lon, lat = np.meshgrid(lon, lat)
        good = dat.nonzero()
        lat = lat[good]
        lon = lon[good]
        data = [seapy.roms.obs.raw_data("TEMP", self.provenance, dat.compressed(),
                                        err[good], self.temp_error)]
        # Grid it
        obs = seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                     data, self.dt, depth_adjust=True, title=title)
        obs.z *= 0
        obs.depth = -self.depth * np.ones(len(obs.depth))
        return obs


class modis_sst_map(obsgen):
    """
    class to process MODIS SST map netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, temp_error=0.5,
                 temp_limits=None, provenance="SST_MODIS_AQUA"):

        self.temp_error = temp_error
        self.provenance = provenance.upper()
        if temp_limits is None:
            self.temp_limits = (2, 35)
        else:
            self.temp_limits = temp_limits
        super().__init__(grid, dt, reftime)

    def convert_file(self, file, title="MODIS SST Obs"):
        """
        Load an MODIS file and convert into an obs structure
        """
        # Load MODIS Data
        import re

        nc = seapy.netcdf(file)
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]
        dat = np.ma.masked_outside(nc.variables["sst"][:],
                                   self.temp_limits[0], self.temp_limits[1])
        err = np.ones(dat.shape) * self.temp_error

        time = seapy.date2day(datetime.datetime.strptime(
            re.sub('\.[0-9]+Z$', '', nc.time_coverage_end),
            "%Y-%m-%dT%H:%M:%S"), self.epoch)

        # Check the data flags
        flags = np.ma.masked_not_equal(nc.variables["qual_sst"][:], 0)
        dat[flags.mask] = np.ma.masked

        nc.close()

        if self.grid.east():
            lon[lon < 0] += 360
        lon, lat = np.meshgrid(lon, lat)
        good = dat.nonzero()
        lat = lat[good]
        lon = lon[good]
        data = [seapy.roms.obs.raw_data("TEMP", self.provenance, dat.compressed(),
                                        err[good], self.temp_error)]
        # Grid it
        return seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                      data, self.dt, title)


class remss_swath(obsgen):
    """
    class to process REMSS SST swath netcdf files into ROMS observation
    files. The files may be AMSRE, TMI, etc. This is a subclass of
    seapy.roms.genobs.genobs, and handles the loading of the data.
    """

    def __init__(self, grid, dt, check_qc_flags=True, reftime=seapy.default_epoch, temp_error=0.4,
                 temp_limits=None, provenance="SST_REMSS"):
        self.temp_error = temp_error
        self.provenance = provenance.upper()
        self.check_qc_flags = check_qc_flags
        if temp_limits is None:
            self.temp_limits = (2, 35)
        else:
            self.temp_limits = temp_limits
        super().__init__(grid, dt, reftime)

    def convert_file(self, file, title="REMSS SST Obs"):
        """
        Load an REMSS file and convert into an obs structure
        """
        # Load REMSS Data
        nc = seapy.netcdf(file)
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]
        dat = np.ma.masked_outside(np.squeeze(
            nc.variables["sea_surface_temperature"][:]) - 273.15,
            self.temp_limits[0], self.temp_limits[1])
        err = np.ma.masked_outside(np.squeeze(
            nc.variables["sses_standard_deviation"][:]), 0.01, 2.0)
        dat[err.mask] = np.ma.masked

        # Check the data flags
        if self.check_qc_flags:
            flags = np.ma.masked_not_equal(
                np.squeeze(nc.variables["quality_level"][:]), 5)
            dat[flags.mask] = np.ma.masked
        else:
            dat = np.ma.masked_where(
                np.squeeze(nc.variables["quality_level"][:]).data == 1, dat)

        # Grab the observation time
        time = seapy.roms.num2date(nc, "time", records=[0])[0] - self.epoch
        dtime = nc.variables["sst_dtime"][:]
        time = np.squeeze((time.total_seconds() + dtime) * seapy.secs2day)
        nc.close()
        if self.grid.east():
            lon[lon < 0] += 360
        good = dat.nonzero()
        data = [seapy.roms.obs.raw_data("TEMP", self.provenance,
                                        dat.compressed(),
                                        err[good], self.temp_error)]
        # Grid it
        return seapy.roms.obs.gridder(self.grid, time[good], lon[good], lat[good],
                                      None, data, self.dt, title)


class remss_map(obsgen):
    """
    class to process REMSS SST map netcdf files into ROMS observation
    files. The files may be AMSRE, TMI, etc. This is a subclass of
    seapy.roms.genobs.genobs, and handles the loading of the data.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, temp_error=0.4,
                 temp_limits=None, provenance="SST_REMSS"):
        self.temp_error = temp_error
        self.provenance = provenance.upper()
        if temp_limits is None:
            self.temp_limits = (2, 35)
        else:
            self.temp_limits = temp_limits
        super().__init__(grid, dt, reftime)

    def convert_file(self, file, title="REMSS SST Obs"):
        """
        Load an REMSS file and convert into an obs structure
        """
        # Load REMSS Data
        nc = seapy.netcdf(file)
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]
        dat = np.ma.masked_outside(np.squeeze(
            nc.variables["sea_surface_temperature"][:]) - 273.15,
            self.temp_limits[0], self.temp_limits[1])
        err = np.ma.masked_outside(np.squeeze(
            nc.variables["SSES_standard_deviation_error"][:]), 0.01, 2.0)
        dat[err.mask] = np.ma.masked

        # Check the data flags
        flags = np.ma.masked_not_equal(
            np.squeeze(nc.variables["rejection_flag"][:]), 0)
        dat[flags.mask] = np.ma.masked
        err[flags.mask] = np.ma.masked

        # Grab the observation time
        time = seapy.roms.num2date(nc, "time", epoch=self.epoch)
        sst_time = nc.variables["sst_dtime"][:] * seapy.secs2day
        for n, i in enumerate(time):
            sst_time[n, :, :] += i
        sst_time[dat.mask] = np.ma.masked

        # Set up the coordinate
        lon, lat = np.meshgrid(lon, lat)
        lon = np.ma.masked_where(dat.mask, seapy.adddim(lon, len(time)))
        lat = np.ma.masked_where(dat.mask, seapy.adddim(lat, len(time)))

        nc.close()

        if self.grid.east():
            lon[lon < 0] += 360
        data = [seapy.roms.obs.raw_data("TEMP", self.provenance,
                                        dat.compressed(),
                                        err.compressed(), self.temp_error)]
        # Grid it
        return seapy.roms.obs.gridder(self.grid, sst_time.compressed(),
                                      lon.compressed(), lat.compressed, None,
                                      data, self.dt, title)


class viirs_swath(obsgen):
    """
    class to process VIIRS SST swath netcdf files into ROMS observation
    files.  This is a subclass of
    seapy.roms.obsgen.obsgen, and handles the loading of the data.
    """

    def __init__(self, grid, dt, check_qc_flags=True, reftime=seapy.default_epoch,
                 temp_error=0.4, temp_limits=None, provenance="SST_VIIRS"):
        self.temp_error = temp_error
        self.provenance = provenance.upper()
        self.check_qc_flags = check_qc_flags
        if temp_limits is None:
            self.temp_limits = (2, 35)
        else:
            self.temp_limits = temp_limits
        super().__init__(grid, dt, reftime)

    def convert_file(self, file, title="VIIRS SST Obs"):
        """
        Load a VIIRS file and convert into an obs structure
        """
        # Load VIIRS Data
        nc = seapy.netcdf(file, aggdim="time")
        lon = nc.variables["lon"][:]
        lat = nc.variables["lat"][:]
        dat = np.ma.masked_outside(
            nc.variables["sea_surface_temperature"][:] - 273.15,
            self.temp_limits[0], self.temp_limits[1])
        err = np.ma.masked_outside(
            nc.variables["sses_standard_deviation"][:], 0.01, 2.0)
        dat[err.mask] = np.ma.masked

        # Check the data flags
        if self.check_qc_flags:
            flags = np.ma.masked_not_equal(
                nc.variables["quality_level"][:], 5)
            dat[flags.mask] = np.ma.masked
        else:
            dat = np.ma.masked_where(
                nc.variables["quality_level"][:].data == 1, dat)

        # Grab the observation time
        time = netCDF4.num2date(nc.variables["time"][:],
                                nc.variables["time"].units) - self.epoch
        time = np.asarray([x.total_seconds() for x in time])[
            :, np.newaxis, np.newaxis]
        dtime = nc.variables["sst_dtime"][:]
        time = (time + dtime) * seapy.secs2day
        nc.close()

        # Set up the coordinate
        lon = np.ma.masked_where(dat.mask, seapy.adddim(lon, len(time)))
        lat = np.ma.masked_where(dat.mask, seapy.adddim(lat, len(time)))
        if self.grid.east():
            lon[lon < 0] += 360
        good = dat.nonzero()
        data = [seapy.roms.obs.raw_data("TEMP", self.provenance,
                                        dat.compressed(),
                                        err[good], self.temp_error)]
        # Grid it
        return seapy.roms.obs.gridder(self.grid, time[good], lon[good], lat[good],
                                      None, data, self.dt, title)


##############################################################################
#
# IN SITU DATA
#
##############################################################################


class seaglider_profile(obsgen):
    """
    class to process SeaGlider .pro files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, dtype=None, temp_limits=None,
                 salt_limits=None, depth_limit=-15, temp_error=0.2,
                 salt_error=0.05):
        if temp_limits is None:
            self.temp_limits = (5, 30)
        else:
            self.temp_limits = temp_limits
        if salt_limits is None:
            self.salt_limits = (31, 35.5)
        else:
            self.salt_limits = salt_limits
        if dtype is None:
            self.dtype = {'names': ('time', 'pres', 'depth', 'temp', 'cond',
                                    'salt', 'sigma', 'lat', 'lon'),
                          'formats': ['f4'] * 9}
        else:
            self.dtype = dtype
        self.depth_limit = depth_limit
        self.temp_error = temp_error
        self.salt_error = salt_error
        super().__init__(grid, dt, reftime)

    def convert_file(self, file, title="SeaGlider Obs"):
        """
        Load a SeaGlider .pro file and convert into an obs structure
        """
        import re

        # Load the text file. All data goes into the pro dictionary
        # as defined by dtype. The header information needs to be parsed
        with open(file) as myfile:
            header = [myfile.readline() for i in range(19)]
            pro = np.loadtxt(myfile, self.dtype, delimiter=',', comments='%')

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
                      "GLIDER_SG" + params["glider"]
        provenance = seapy.roms.obs.asprovenance(glider_name)
        try:
            date = [int(s) for s in re.findall('([\d]{2})\s', params["start"])]
            start_time = datetime.datetime.strptime(params["start"].strip(),
                                                    "%m %d 1%y %H %M %S")
            dtime = (start_time - self.epoch).total_seconds() / 86400
        except:
            raise ValueError("date format incorrect in file: " + file)

        # Make sure that the GPS fix isn't screwy
        if self.grid.east():
            pro["lon"][pro["lon"] < 0] += 360
        dist = seapy.earth_distance(pro["lon"][0], pro["lat"][0],
                                    pro["lon"][-1], pro["lat"][-1])
        velocity = dist / pro["time"][-1]
        if velocity > 2:
            warn("WARNING: GPS fix is incorrect for " + file)
            return None

        # Build the data with masked entries
        temp = np.ma.masked_outside(pro["temp"], self.temp_limits[0],
                                    self.temp_limits[1])
        salt = np.ma.masked_outside(pro["salt"], self.salt_limits[0],
                                    self.salt_limits[1])
        depth = np.ma.masked_greater(-pro["depth"], self.depth_limit)
        good = ~np.ma.getmaskarray(depth)

        # Grid it
        data = [seapy.roms.obs.raw_data("TEMP", provenance, temp[good],
                                        None, self.temp_error),
                seapy.roms.obs.raw_data("SALT", provenance, salt[good],
                                        None, self.salt_error)]
        return seapy.roms.obs.gridder(self.grid, pro["time"][good] / 86400 + dtime,
                                      pro["lon"][good],
                                      pro["lat"][good],
                                      depth.compressed(),
                                      data, self.dt, title)


class mooring(obsgen):
    """
    Class to process generic moorings into ROMS observation files. This
    handles temp, salt, u, and v.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, temp_limits=None,
                 salt_limits=None, u_limits=None, v_limits=None,
                 depth_limit=0, temp_error=0.25, salt_error=0.08,
                 u_error=0.08, v_error=0.08, lat=None, lon=None,
                 provenance=None):
        if temp_limits is None:
            self.temp_limits = (5, 35)
        else:
            self.temp_limits = temp_limits
        if salt_limits is None:
            self.salt_limits = (31, 35.5)
        else:
            self.salt_limits = salt_limits
        if u_limits is None:
            self.u_limits = (-3, 3)
        else:
            self.u_limits = u_limits
        if v_limits is None:
            self.v_limits = (-3, 3)
        else:
            self.v_limits = v_limits
        if provenance is None:
            self.provenance = seapy.roms.obs.asprovenance("MOORING")
        else:
            self.provenance = provenance.upper()
        self.depth_limit = depth_limit
        self.temp_error = temp_error
        self.salt_error = salt_error
        self.u_error = u_error
        self.v_error = v_error
        self.lat = np.atleast_1d(lat)
        self.lon = np.atleast_1d(lon)
        super().__init__(grid, dt, reftime)

    def convert_data(self, time, depth, data, error=None, title="Mooring Obs"):
        """
        Given a set of data, process into an observation structure

        Parameters
        ----------
        time : ndarray
          time of observations
        depth : ndarray
          depth of observations. depth is in rows, time in columns.
          If depth does not change with time, it will be replicated in time.
        data : dict
          data to put into observations. A dictionary using seapy.roms.fields
          as keys.
        error : dict, optional
          error of the observations (same keys and sizes as data)
        title : string, optional
          title for obs

        Returns
        -------
        obs: seapy.roms.obs.obs
        """
        # Check that the lat/lon is in the grid
        if self.grid.east():
            self.lon[self.lon <= 0] += 360
        else:
            self.lon[self.lon >= 180] -= 360
        if not np.logical_and.reduce((
                self.lon >= np.min(self.grid.lon_rho),
                self.lon <= np.max(self.grid.lon_rho),
                self.lat >= np.min(self.grid.lat_rho),
                self.lat <= np.max(self.grid.lat_rho))):
            warn("Mooring location is not in grid")
            return
        depth = np.atleast_1d(depth)

        if not error:
            error = {}

        if not data:
            warn("No data is provided")
            return

        # Process the data
        obsdata = []
        for field in data:
            limit = getattr(self, field + '_limits')
            vals = np.ma.masked_outside(data[field], limit[0], limit[1],
                                        copy=False)
            obsdata.append(seapy.roms.obs.raw_data(field, self.provenance,
                                                   vals, getattr(
                                                       error, field, None),
                                                   getattr(self, field + '_error')))

        ndep = depth.size
        nt = len(time)
        lat = np.resize(self.lat, (nt, ndep))
        lon = np.resize(self.lon, (nt, ndep))
        depth = np.resize(depth, (nt, ndep))
        time = np.resize(time, (nt, ndep))
        return seapy.roms.obs.gridder(self.grid, time, lon, lat, depth,
                                      obsdata, self.dt, title)


class tao_mooring(mooring):
    """
    class to process TAO files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, temp_limits=None,
                 salt_limits=None, u_limits=None, v_limits=None,
                 depth_limit=0, temp_error=0.25, salt_error=0.08,
                 u_error=0.08, v_error=0.08):
        super().__init__(grid, dt, reftime)

    def convert_file(self, file, title="TAO Obs"):
        """
        Load a TAO netcdf file and convert into an obs structure
        """
        vals = {"temp": ["T_20", "QT_5020"],
                "salt": ["S_41", "QS_5041"],
                "u": ["U_320", "QS_5300"],
                "v": ["V_321", "QS_5300"]}
        nc = seapy.netcdf(file)
        lat = nc.variables["lat"][:]
        lon = nc.variables["lon"][:]
        if not self.grid.east():
            lon[lon > 180] -= 360
        lat, lon = np.meshgrid(lat, lon)
        time = seapy.roms.num2date(nc, "time", epoch=self.epoch)
        depth = -nc.variables["depth"][:]
        profile_list = np.where(np.logical_and.reduce((
            lon >= np.min(self.grid.lon_rho),
            lon <= np.max(self.grid.lon_rho),
            lat >= np.min(self.grid.lat_rho),
            lat <= np.max(self.grid.lat_rho))))

        # If nothing is in the area, return nothing
        if not profile_list[0].size:
            return None

        # Process each of the variables that are present
        obsdata = []
        for field in vals:
            limit = getattr(self, field + '_limits')
            if vals[field][0] in nc.variables:
                data = nc.variables[vals[field][0]][:]
                data = np.ma.masked_outside(
                    data[profile_list[0], profile_list[1], :, :],
                    limit[0], limit[1], copy=False)
                qc = nc.variables[vals[field][1]][:]
                qc = qc[profile_list[0], profile_list[1], :, :]
                bad = np.where(np.logical_and(qc != 1, qc != 2))
                data[bad] = np.ma.masked
                obsdata.append(seapy.roms.obs.raw_data(field, "TAO_ARRAY",
                                                       data.compressed(), None,
                                                       getattr(self, field + '_error')))
        nc.close()

        # Build the time, lon, lat, and depth arrays of appropriate size
        npts = profile_list[0].size
        ndep = depth.size
        nt = len(time)
        lat = np.resize(lat[profile_list], (nt, ndep, npts))
        lat = np.squeeze(np.transpose(lat, (2, 1, 0)))[~data.mask]
        lon = np.resize(lon[profile_list], (nt, ndep, npts))
        lon = np.squeeze(np.transpose(lon, (2, 1, 0)))[~data.mask]
        depth = np.resize(depth, (npts, nt, ndep))
        depth = np.squeeze(np.transpose(depth, (0, 2, 1)))[~data.mask]
        time = np.squeeze(np.resize(time, (npts, ndep, nt)))[~data.mask]
        return seapy.roms.obs.gridder(self.grid, time, lon, lat, depth,
                                      obsdata, self.dt, title)


class argo_ctd(obsgen):
    """
    class to process ARGO CTD netcdf files into ROMS observation
    files. This is a subclass of seapy.roms.genobs.genobs, and handles
    the loading of the data.
    """

    def __init__(self, grid, dt, reftime=seapy.default_epoch, temp_limits=None,
                 salt_limits=None, temp_error=0.25,
                 salt_error=0.1):
        if temp_limits is None:
            self.temp_limits = (2, 35)
        else:
            self.temp_limits = temp_limits
        if salt_limits is None:

            self.salt_limits = (10, 35.5)
        else:
            self.salt_limits = salt_limits
        self.temp_error = temp_error
        self.salt_error = salt_error
        super().__init__(grid, dt, reftime)

    def datespan_file(self, file):
        """
        return the just the day that this argo file covers
        """
        nc = seapy.netcdf(file)
        try:
            d = netCDF4.num2date(nc.variables['JULD'][0],
                                 nc.variables['JULD'].units)
            st = datetime.datetime(*d.timetuple()[:3])
            en = datetime.datetime(*d.timetuple()[:3] + (23, 59, 59))
        except:
            st = en = None
            pass
        finally:
            nc.close()
            return st, en

    def convert_file(self, file, title="Argo Obs"):
        """
        Load an Argo file and convert into an obs structure
        """
        nc = seapy.netcdf(file, aggdim="N_PROF")

        # Load the position of all profiles in the file
        lon = nc.variables["LONGITUDE"][:]
        lat = nc.variables["LATITUDE"][:]
        pro_q = nc.variables["POSITION_QC"][:].astype(int)
        # Find the profiles that are in our area with known locations quality
        if self.grid.east():
            lon[lon < 0] += 360
        profile_list = np.where(np.logical_and.reduce((
            lat >= np.min(self.grid.lat_rho),
            lat <= np.max(self.grid.lat_rho),
            lon >= np.min(self.grid.lon_rho),
            lon <= np.max(self.grid.lon_rho),
            pro_q == 1)))[0]

        # Check which are good profiles
        profile_qc = nc.variables["PROFILE_PRES_QC"][
            profile_list].astype('<U1')
        profile_list = profile_list[profile_qc == 'A']
        if not profile_list.size:
            return None

        # Load only the data from those in our area
        julian_day = nc.variables["JULD_LOCATION"][profile_list]
        argo_epoch = datetime.datetime.strptime(''.join(
            nc.variables["REFERENCE_DATE_TIME"][:].astype('<U1')), '%Y%m%d%H%M%S')
        time_delta = (self.epoch - argo_epoch).days
        file_stamp = datetime.datetime.strptime(''.join(
            nc.variables["DATE_CREATION"][:].astype('<U1')), '%Y%m%d%H%M%S')

        # Grab data over the previous day
        file_time = np.minimum((file_stamp - argo_epoch).days,
                               int(np.max(julian_day)))
        time_list = np.where(julian_day >= file_time - 1)[0]
        julian_day = julian_day[time_list]
        lon = lon[profile_list[time_list]]
        lat = lat[profile_list[time_list]]
        profile_list = profile_list[time_list]

        # Load the data in our region and time
        temp = nc.variables["TEMP"][profile_list, :]
        temp_qc = nc.variables["TEMP_QC"][profile_list, :]
        salt = nc.variables["PSAL"][profile_list, :]
        salt_qc = nc.variables["PSAL_QC"][profile_list, :]
        pres = nc.variables["PRES"][profile_list, :]
        pres_qc = nc.variables["PRES_QC"][profile_list, :]
        nc.close()

        # Ensure consistency
        full_mask = np.logical_or.reduce((temp.mask, salt.mask, pres.mask))
        temp[full_mask] = np.ma.masked
        temp_qc[full_mask] = np.ma.masked
        salt[full_mask] = np.ma.masked
        salt_qc[full_mask] = np.ma.masked
        pres[full_mask] = np.ma.masked
        pres_qc[full_mask] = np.ma.masked

        # Combine the QC codes
        qc = np.mean(np.vstack((temp_qc.compressed(), salt_qc.compressed(),
                                pres_qc.compressed())).astype(int), axis=0)
        good_data = np.where(qc == 1)

        # Put everything together into individual observations
        time = np.resize(julian_day - time_delta,
                         pres.shape[::-1]).T[~temp.mask][good_data]
        lat = np.resize(lat, pres.shape[::-1]).T[~temp.mask][good_data]
        lon = np.resize(lon, pres.shape[::-1]).T[~temp.mask][good_data]
        depth = -seapy.seawater.depth(pres.compressed()[good_data], lat)

        # Apply the limits
        temp = np.ma.masked_outside(temp.compressed()[good_data],
                                    self.temp_limits[0], self.temp_limits[1])
        salt = np.ma.masked_outside(salt.compressed()[good_data],
                                    self.salt_limits[0], self.salt_limits[1])

        data = [seapy.roms.obs.raw_data("TEMP", "CTD_ARGO", temp,
                                        None, self.temp_error),
                seapy.roms.obs.raw_data("SALT", "CTD_ARGO", salt,
                                        None, self.salt_error)]

        return seapy.roms.obs.gridder(self.grid, time, lon, lat, depth,
                                      data, self.dt, title)
