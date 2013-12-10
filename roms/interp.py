#!/usr/bin/env python
"""
  roms.interp
  
  Methods to interpolate ROMS fields onto other grids

  Written by Brian Powell on 11/02/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np
import netCDF4
import os.path
import re
import seapy
from progressbar import ProgressBar
from joblib import Parallel, delayed
import pudb

def _interp2_thread(rx, ry, data, zx, zy, pmap, weight, nx, ny, mask):
    # Convolve the water over the land
    seapy.convolve_mask(data, 7)
    
    # Interpolate the field and return the result
    res, pm = seapy.oasurf(rx, ry, data, zx, zy, pmap, weight, nx, ny)
    return np.ma.array(res, mask=logical_or(mask==0,res>9e10), copy=False)

def _interp3_thread(rx, ry, rz, data, zx, zy, zz, pmap, weight, nx, ny, mask):
    # Make the mask 3D
    mask = seapy.adddim(mask, zz.shape[0])
    
    # Convolve the water over the land
    seapy.convolve_mask(data, 7)
    
    # Interpolate the field and return the result
    res, pm = seapy.oavol(rx, ry, rz, data, zx, zy, zz, pmap, \
                            weight, nx, ny)
    return np.ma.array(res, mask=np.logical_or(mask==0,res>9e10), copy=False)
    
def to_grid(src_file, dest_file, records=None, threads=1, nx=0, ny=0):
    """
    Given a roms file (average, history, etc.), interpolate the fields
    onto another gridded file. 
    """
    # Load the Source file
    ncsrc = seapy.model.grid(src_file, minimal=False)

    # Load the z-grid
    ncdst = seapy.model.grid(dest_file, minimal=False, nc="a")

    # Generate a file to store the pmap information
    rname = re.search("[^\.]*",os.path.basename(src_file));
    zname = re.search("[^\.]*",os.path.basename(dest_file));
    pmap_file = rname.group() + "_" + zname.group() + "_pmap.npz"

    # Create or load the pmaps depending on if they exist
    weight=1
    nx = np.ceil( np.round( np.mean( (ncsrc.dm / ncdst.dm).flatten() ),1)) \
        if nx==0 else nx
    ny = np.ceil( np.round( np.mean( (ncsrc.dn / ncdst.dn).flatten() ),1)) \
        if ny==0 else ny
    if os.path.isfile(pmap_file):
        pmap = np.load(pmap_file)
    else:
        tmp = np.ones(ncsrc.lat_rho.shape)
        tmp, pmaprho = seapy.oasurf(ncsrc.lon_rho, ncsrc.lat_rho, \
                            tmp, ncdst.lon_rho, ncdst.lat_rho, \
                            weight=weight, nx=nx, ny=ny)
        tmp = np.ones(ncsrc.lat_u.shape)
        tmp, pmapu = seapy.oasurf(ncsrc.lon_u, ncsrc.lat_u, \
                            tmp, ncdst.lon_rho, ncdst.lat_rho, \
                            weight=weight, nx=nx, ny=ny)
        tmp = np.ones(ncsrc.lat_v.shape)
        tmp, pmapv = seapy.oasurf(ncsrc.lon_v, ncsrc.lat_v, \
                            tmp, ncdst.lon_rho, ncdst.lat_rho, \
                            weight=weight, nx=nx, ny=ny)
        np.savez(pmap_file, pmaprho=pmaprho, pmapu=pmapu, pmapv=pmapv)
        pmap = {"pmaprho":pmaprho, "pmapu":pmapu, "pmapv":pmapv}
        
    # Interpolate the fields
    records = np.arange(0, len(ncsrc._nc.variables["ocean_time"][:])) \
                 if records == None else records

    for k in seapy.roms.fields.keys():
        if k not in ncdst._nc.variables:
            continue
        print(k)
        grd = seapy.roms.fields[k]["grid"]
        if seapy.roms.fields[k]["dims"]==2:
            ndata = np.ma.array(Parallel(n_jobs=threads,verbose=2)\
                             (delayed(_interp2_thread) (
              getattr(ncsrc,"lon_"+grd), getattr(ncsrc,"lat_"+grd),
              ncsrc._nc.variables[k][i,:,:],
              getattr(ncdst,"lon_"+grd), getattr(ncdst,"lat_"+grd),
              pmap["pmap"+grd], weight,
              nx, ny, getattr(ncdst,"mask_"+grd)) for i in records))
        else:
            ndata = np.ma.array( Parallel(n_jobs=threads,verbose=2)
                             (delayed(_interp3_thread)( 
              getattr(ncsrc,"lon_"+grd), getattr(ncsrc,"lat_"+grd),
              getattr(ncsrc,"depth_"+grd),
              ncsrc._nc.variables[k][i,:,:,:], 
              getattr(ncdst,"lon_"+grd), getattr(ncdst,"lat_"+grd),
              getattr(ncdst,"depth_"+grd),
              pmap["pmap"+grd], weight,
              nx, ny, getattr(ncdst,"mask_"+grd)) for i in records))
        pass
        # ndata.fill_value = None
        ncdst._nc.variables[k][:] = ndata

def to_zgrid(roms_file, z_file, depth=None, records=None, threads=1):
    """
    Given an existing ROMS history or average file, create (if does not exit)
    a new z-grid file with the same horizontal extent and the specified depths
    and interpolate the ROMS fields onto the z-grid. If an existing z-grid
    file is given, it is interpolated onto the specified.
    """
    if not os.path.isfile(z_file):
        ncroms = seapy.model.grid(roms_file)
        records = np.arange(0, len(ncroms._nc.variables["ocean_time"][:])) \
                 if records == None else records
        lat=ncroms.lat_rho.shape[0]
        lon=ncroms.lat_rho.shape[1]
        origin=seapy.roms.get_timebase(ncroms._nc.variables["ocean_time"])
        if depth==None:
            raise ValueError("depth must be specified")
        ncz=seapy.roms.ncgen.create_zlevel(z_file,lat,lon,len(depth),origin,"ROMS z-level")
        ncz.variables["lat"][:]=ncroms.lat_rho
        ncz.variables["lon"][:]=ncroms.lon_rho
        ncz.variables["depth"][:]=depth
        ncz.variables["mask"][:]=ncroms.mask_rho
        ncz.variables["time"][:]=ncroms._nc.variables["ocean_time"][:]/86400
        ncz.close()
    # Call the interpolation
    to_grid(roms_file, z_file, records=records, threads=threads)
    
def to_child(src_file, dest_file, dest_grid=None, records=None, threads=1):
    """
    Given an existing ROMS history or average file, create (if does not exit) a
    new ROMS file using the given ROMS destination grid and interpolate the
    ROMS fields onto the new grid. If an existing destination file is given, it
    is interpolated onto the specified.
    """
    if not os.path.isfile(dest_file) and dest_grid != None:
        if isinstance(dest_grid,basestring):
            destg = seapy.model.grid(dest_grid)
        else:
            destg = dest_grid
        ncroms = seapy.model.grid(src_file)
        records = np.arange(0, len(ncroms._nc.variables["ocean_time"][:])) \
                 if records == None else records
        origin=seapy.roms.get_timebase(ncroms._nc.variables["ocean_time"])
        # ncdest=seapy.roms.ncgen.create_zlevel(dest_file, 
        #          eta_rho=destg.,xi_rho=destg.,N=destg.n,timebase=origin,
        #          title="interpolated from "+src_file)
        ncdest.variables["lat_rho"][:]=destg.lat_rho
        ncdest.variables["lon_rho"][:]=destg.lon_rho
        ncdest.variables["lat_u"][:]=destg.lat_u
        ncdest.variables["lon_u"][:]=destg.lon_u
        ncdest.variables["lat_v"][:]=destg.lat_v
        ncdest.variables["lon_v"][:]=destg.lon_v
        ncdest.variables["Vtransform"][:]=destg.vtransform
        ncdest.variables["Vstretching"][:]=destg.vstretching
        ncdest.variables["theta_s"][:]=destg.theta_s
        ncdest.variables["theta_b"][:]=destg.theta_b
        ncdest.variables["hc"][:]=destg.hc
        ncdest.variables["Tcline"][:]=destg.tcline
        ncdest.variables["s_rho"][:]=destg.s_rho
        ncdest.variables["Cs_r"][:]=destg.cs_r
        ncdest.variables["h"][:]=destg.h
        ncdest.variables["time"][:]=ncroms._nc.variables["ocean_time"][:]/86400
        ncdest.close()
    else:
        raise AttributeError("you must supply a destination file or a grid to make the file")

    # Call the interpolation
    to_grid(roms_file, dest_file, records=records, threads=threads)
    
pass