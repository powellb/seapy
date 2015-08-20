#!/usr/bin/env python
"""
  analysis.py

  Methods to assist in the analysis of ROMS fields

  Written by Brian Powell on 05/24/15
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function
import numpy as np
from joblib import Parallel, delayed
import seapy
import netCDF4

# Create a function to do the multiprocessing vertical interpolation.
# This function is used because we don't want to construct an array of
# tuples from the oa routine.
def __dinterp(x,y,z,dat,fz,pmap):
    ndat, pm = seapy.oavol(x, y, z, dat, x, y, fz, pmap, 5, 1, 1)
    ndat = np.squeeze(ndat)
    return np.ma.masked_where(np.abs(ndat)>9e10, ndat, copy=False)

def constant_depth(field, grid, depth, zeta=None, threads=-2):
    """
    Find the values of a 3-D field at a constant depth for all times given.

    Parameters
    ----------
    field : ndarray,
        ROMS 3-D field to interpolate onto a constant depth level. Can be
        two- or three-dimensional array (first dimension assumed to be time).
    grid : seapy.model.grid or string,
        Grid that defines the depths and stretching for the field given
    depth : float,
        Depth (in meters) to find all values
    zeta : ndarray, optional,
        ROMS zeta field corresponding to field if you wish to apply the SSH
        correction to the depth calculations.
    threads : int, optional,
        Number of threads to use for processing

    Returns
    -------
    nfield : ndarray,
        Values from ROMS field on the given constant depth
    """

    # Make sure our inputs are all valid
    grid = seapy.model.asgrid(grid)
    if np.ndim(field)==3:
        field = seapy.adddim(field)
    if zeta is not None and np.ndim(zeta==2):
        zeta = seapy.adddim(zeta)
    depth = depth if depth < 0 else -depth

    # Set up some arrays
    x, y = np.meshgrid(np.arange(field.shape[-1]),np.arange(field.shape[-2]))
    fz, pmap = seapy.oasurf(x,y,x,x,y,None,5,1,1)
    fz = seapy.adddim(np.ones(x.shape))*depth
    # Loop over all times, generate new field at depth
    nfield = np.ma.array(Parallel(n_jobs=threads,verbose=2)\
                  (delayed(__dinterp)(x,y,grid.depth_rho,
                                      np.squeeze(field[i,:,:,:]),fz,pmap)
               for i in range(field.shape[0])), copy=False)

    return nfield

def gen_std_i(roms_file, std_file, std_window=5, pad=1, skip=30, fields=None):
    """
    Create a std file for the given ocean fields. This std file can be used
    for initial conditions constraint in 4D-Var. This requires a long-term
    model spinup file from which to compute the standard deviation.

    Parameters
    ----------
    roms_file: string,
        The ROMS (history or average) file from which to compute the std
    std_file: string,
        The name of the file to store the standard deviations fields
    std_window: int,
        The size of the window (in number of records) to compute the std over
    pad: int,
        How much to pad each side of the window for overlap. For example,
        std_window=10 and pad=2 would give a total window of 14 with 2 records
        used in the prior window and 2 in the post window as well.
    skip: int,
        How many records to skip at the beginning of the file
    fields: list of str,
        The fields to compute std for. Default is to use the ROMS prognostic
        variables.

    Returns
    -------
        None
    """
    if fields is None:
        fields = list(seapy.roms.fields.keys())

    grid = seapy.model.asgrid(roms_file)
    nc = netCDF4.Dataset(roms_file)
    time_var = seapy.roms.gettimevar(nc)
    epoch = netCDF4.num2date(0, nc.variables[time_var].units)
    time = nc.variables[time_var][:]
    ncout = seapy.roms.ncgen.create_da_ini_std(std_file,
                      eta_rho=grid.ln, xi_rho=grid.lm, s_rho=grid.n,
                      reftime=epoch, title="std from " + roms_file)
    grid.to_netcdf(nc)

    # If there are any fields that are not part of the standard, add them
    # to the output file
    for f in fields:
        if f not in seapy.roms.fields:
            ncout.createVariable(f, np.float32,
                                 ('ocean_time', "s_rho", "eta_rho", "xi_rho"))

    # Loop over the time with the variance window:
    for n, t in enumerate(seapy.progressbar.progress(np.arange(skip + pad,
                                                len(time) - pad, std_window))):
        idx = np.arange(t - pad, t + std_window + pad)
        ncout.variables[time_var][n] = np.mean(time[idx])
        for v in fields:
            dat = np.std(nc.variables[v][idx, :], axis=0)
            dat[dat > 10] = 0.0
            ncout.variables[v][n, :] = dat
        ncout.sync()
    ncout.close()
    nc.close()

def gen_std_f(roms_file, std_file, std_window=5, pad=1, skip=30, fields=None):
    """
    Create a std file for the given atmospheric forcing fields. This std
    file can be used for the forcing constraint in 4D-Var. This requires a
    long-term model spinup file from which to compute the standard deviation.

    Parameters
    ----------
    roms_file: string,
        The ROMS (history or average) file from which to compute the std
    std_file: string,
        The name of the file to store the standard deviations fields
    std_window: int,
        The size of the window (in number of records) to compute the std over
    pad: int,
        How much to pad each side of the window for overlap. For example,
        std_window=10 and pad=2 would give a total window of 14 with 2 records
        used in the prior window and 2 in the post window as well.
    skip: int,
        How many records to skip at the beginning of the file
    fields: list of str,
        The fields to compute std for. Default is to use the ROMS atmospheric
        variables (sustr, svstr, shflux, ssflux).

    Returns
    -------
        None
    """
    if fields is None:
        fields = ["sustr", "svstr", "shflux", "ssflux"]

    grid = seapy.model.asgrid(roms_file)
    nc = netCDF4.Dataset(roms_file)
    time_var = seapy.roms.gettimevar(nc)
    epoch = netCDF4.num2date(0, nc.variables[time_var].units)
    time = nc.variables[time_var][:]
    ncout = seapy.roms.ncgen.create_da_frc_std(std_file,
                      eta_rho=grid.ln, xi_rho=grid.lm, s_rho=grid.n,
                      reftime=epoch, title="std from " + roms_file)
    grid.to_netcdf(nc)

    # If there are any fields that are not part of the standard, add them
    # to the output file
    for f in fields:
        if f not in seapy.roms.fields:
            ncout.createVariable(f, np.float32,
                                 ('ocean_time', "eta_rho", "xi_rho"))

    # Loop over the time with the variance window:
    for n, t in enumerate(seapy.progressbar.progress(np.arange(skip + pad,
                                                len(time) - pad, std_window))):
        idx = np.arange(t - pad, t + std_window + pad)
        ncout.variables[time_var][n] = np.mean(time[idx])
        for v in fields:
            dat = np.std(nc.variables[v][idx, :], axis=0)
            dat[dat > 10] = 0.0
            ncout.variables[v][n, :] = dat
        ncout.sync()
    ncout.close()
    nc.close()
