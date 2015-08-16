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
import os
import seapy
from seapy.timeout import timeout,TimeoutError
from joblib import Parallel, delayed
from warnings import warn

_up_scaling = {"zeta":1.0, "u":1.0, "v":1.0, "temp":1.0, "salt":1.0}
_down_scaling = {"zeta":1.0, "u":0.95, "v":0.95, "temp":0.98, "salt":1.02}
_ksize_range = (7,15)
_max_memory = 512*1024*1024    # Limit arrays to 512MB in size [bytes]

def __mask_z_grid(z_data, src_depth, z_depth):
    """
    When interpolating to z-grid, we need to apply depth dependent masking
    based on the original ROMS depths
    """
    for k in np.arange(0,z_depth.shape[0]):
        idx=np.nonzero(z_depth[k,:,:]<src_depth)
        if z_data.ndim==4:
            z_data.mask[:,k,idx[0],idx[1]]=True
        elif z_data.ndim==3:
            z_data.mask[k,idx[0],idx[1]]=True

def __interp2_thread(rx, ry, data, zx, zy, pmap, weight, nx, ny, mask):
    """
    internal routine: 2D interpolation thread for parallel interpolation
    """
    data = np.ma.fix_invalid(data, copy=False, fill_value=-999999.0)
    # Convolve the water over the land
    ksize=2*np.round(np.sqrt((nx/np.median(np.diff(rx)))**2 +
                    (ny/np.median(np.diff(ry.T)))**2))+1
    if ksize < _ksize_range[0]:
        warn("nx or ny values are too small for stable OA, {:f}".format(ksize))
        ksize=_ksize_range[0]
    elif ksize > _ksize_range[1]:
        warn("nx or ny values are too large for stable OA, {:f}".format(ksize))
        ksize=_ksize_range[1]
    data = seapy.convolve_mask(data, ksize=ksize, copy=False)

    # Interpolate the field and return the result
    with timeout(minutes=30):
        res, pm = seapy.oasurf(rx, ry, data, zx, zy, pmap, weight, nx, ny)

    return np.ma.masked_where(np.logical_or(mask==0,np.abs(res)>9e10), res,
                       copy=False)

def __interp3_thread(rx, ry, rz, data, zx, zy, zz, pmap,
                    weight, nx, ny, mask, up_factor=1.0, down_factor=1.0):
    """
    internal routine: 3D interpolation thread for parallel interpolation
    """
    # Make the mask 3D
    mask = seapy.adddim(mask, zz.shape[0])
    data = np.ma.fix_invalid(data, copy=False, fill_value=-999999.0)

    # To avoid extrapolation, add a new top and bottom layer that replicates
    # the data of the existing current and top. 1) Determine which way the
    # depth goes and add/subtract new layers, and 2) fill in masked values
    # from the layer above/below.
    gradsrc = (rz[0,1,1]-rz[-1,1,1]) > 0
    graddest = (zz[0,1,1]-zz[-1,1,1]) > 0
    nrz = np.zeros((data.shape[0]+2,data.shape[1],data.shape[2]))
    nrz[1:-1,:,:]=rz
    if not gradsrc:
        # The first level is the bottom
        nrz[0,:,:]=rz[0,:,:]-500
        nrz[-1,:,:]=np.minimum(rz[-1,:,:]+50,0)
        factor=down_factor
        # Fill in missing values where we have them from above (level above)
        for k in np.arange(data.shape[0]-2,-1,-1):
            idx=np.nonzero(np.logical_xor(data.mask[k,:,:],data.mask[k+1,:,:]))
            data.mask[k,idx[0],idx[1]]=data.mask[k+1,idx[0],idx[1]]
            data[k,idx[0],idx[1]]=data[k+1,idx[0],idx[1]]*factor
    else:
        # The first level is the top
        nrz[0,:,:]=np.minimum(rz[0,:,:]+50,0)
        nrz[-1,:,:]=rz[-1,:,:]-500
        factor=up_factor
        # Fill in missing values where we have them from below (level below)
        for k in np.arange(1,data.shape[0]):
            idx=np.nonzero(np.logical_xor(data.mask[k,:,:],data.mask[k-1,:,:]))
            data.mask[k,idx[0],idx[1]]=data.mask[k-1,idx[0],idx[1]]
            data[k,idx[0],idx[1]]=data[k-1,idx[0],idx[1]]*factor

    # Convolve the water over the land
    ksize=2*np.round(np.sqrt((nx/np.median(np.diff(rx)))**2 +
                    (ny/np.median(np.diff(ry.T)))**2))+1
    if ksize < _ksize_range[0]:
        warn("nx or ny values are too small for stable OA, {:f}".format(ksize))
        ksize=_ksize_range[0]
    elif ksize > _ksize_range[1]:
        warn("nx or ny values are too large for stable OA, {:f}".format(ksize))
        ksize=_ksize_range[1]
    data = seapy.convolve_mask(data, ksize=ksize, copy=False)

    # Add upper and lower boundaries
    ndat = np.zeros((data.shape[0]+2,data.shape[1],data.shape[2]))
    ndat[0,:,:]=data[0,:,:].filled(np.nan)*factor
    ndat[1:-1,:,:]=data.filled(np.nan)
    ndat[-1,:,:]=data[-1,:,:].filled(np.nan)*factor

    # Interpolate the field and return the result
    with timeout(minutes=30):
        if gradsrc:
            res, pm = seapy.oavol(rx, ry, \
                        nrz[np.arange(nrz.shape[0]-1,-1,-1),:,:], \
                        ndat[np.arange(nrz.shape[0]-1,-1,-1),:,:], \
                        zx, zy, zz, pmap, \
                        weight, nx, ny)
        else:
            res, pm = seapy.oavol(rx, ry, nrz, ndat, zx, zy, zz, \
                                pmap, weight, nx, ny)

    return np.ma.masked_where(np.logical_or(mask==0,np.abs(res)>9e5), res,
                       copy=False)

def __interp3_vel_thread(rx, ry, rz, ra, u, v, zx, zy, zz, za, pmap,
                    weight, nx, ny, mask):
    """
    internal routine: 3D velocity interpolation thread for parallel interpolation
    """
    # Put on the same grid
    if u.shape != v.shape:
        u = seapy.model.u2rho(u, fill=True)
        v = seapy.model.v2rho(v, fill=True)

    # Rotate the fields (NOTE: ROMS angle is negative relative to "true")
    if ra is not None:
        u, v = seapy.rotate(u, v, ra)

    # Interpolate
    u = __interp3_thread(rx, ry, rz, u, zx, zy, zz, pmap,
                        weight, nx, ny, mask, _up_scaling["u"],
                        _down_scaling["u"])
    v = __interp3_thread(rx, ry, rz, v, zx, zy, zz, pmap,
                        weight, nx, ny, mask, _up_scaling["v"],
                        _down_scaling["v"])

    # Rotate to destination (NOTE: ROMS angle is negative relative to "true")
    if za is not None:
        u, v = seapy.rotate(u, v, -za)

    # Return the masked data
    return u, v

def __interp_grids(src_grid, child_grid, ncout, records=None,
            threads=2, nx=0, ny=0, weight=10, vmap=None, z_mask=False,
            pmap=None):
    """
    internal method:  Given a model file (average, history, etc.),
    interpolate the fields onto another gridded file.

    Parameters
    ----------
    src_grid : seapy.model.grid data source (History, Average, etc. file)
    child_grid : seapy.model.grid output data grid
    ncout : netcdf output file
    [records] : array of the record indices to interpolate
    [threads] : number of processing threads
    [nx] : decorrelation length in grid-cells for x
    [ny] : decorrelation length in grid-cells for y
    [vmap] : variable name mapping
    [z_mask] : mask out depths in z-grids
    [pmap] : use the specified pmap rather than compute it

    Returns
    -------
    None

    """
    # If we don't have a variable map, then do a one-to-one mapping
    if vmap is None:
        vmap=dict()
        for k in seapy.roms.fields:
            vmap[k]=k

    # Generate a file to store the pmap information
    sname = getattr(src_grid, 'name', None)
    cname = getattr(child_grid, 'name', None)
    pmap_file = None if any(v is None for v in (sname,cname)) else \
        sname + "_" + cname + "_pmap.npz"

    # Create or load the pmaps depending on if they exist
    if nx==0:
        if hasattr(src_grid,"dm") and hasattr(child_grid,"dm"):
            nx = np.ceil( np.mean( src_grid.dm ) / np.mean( child_grid.dm) )
        else:
            nx = 5
    if ny==0:
        if hasattr(src_grid,"dn") and hasattr(child_grid,"dn"):
            ny = np.ceil( np.mean( src_grid.dn ) / np.mean( child_grid.dn) )
        else:
            ny = 5

    if pmap is None:
        if pmap_file is not None and os.path.isfile(pmap_file):
            pmap = np.load(pmap_file)
        else:
            tmp = np.ones(src_grid.lat_rho.shape,order="F")
            tmp, pmaprho = seapy.oasurf(src_grid.lon_rho, src_grid.lat_rho, \
                                tmp, child_grid.lon_rho, child_grid.lat_rho, \
                                weight=weight, nx=nx, ny=ny)
            tmp = np.ones(src_grid.lat_u.shape,order="F")
            tmp, pmapu = seapy.oasurf(src_grid.lon_u, src_grid.lat_u, \
                                tmp, child_grid.lon_rho, child_grid.lat_rho, \
                                weight=weight, nx=nx, ny=ny)
            tmp = np.ones(src_grid.lat_v.shape,order="F")
            tmp, pmapv = seapy.oasurf(src_grid.lon_v, src_grid.lat_v, \
                                tmp, child_grid.lon_rho, child_grid.lat_rho, \
                                weight=weight, nx=nx, ny=ny)
            if pmap_file is not None:
                np.savez(pmap_file, pmaprho=pmaprho, pmapu=pmapu, pmapv=pmapv)
            pmap = {"pmaprho":pmaprho, "pmapu":pmapu, "pmapv":pmapv}

    # Get the time field
    ncsrc = netCDF4.Dataset(src_grid.filename)
    time = seapy.roms.get_timevar(ncsrc)

    # Interpolate the depths from the source to final grid
    src_depth=np.min(src_grid.depth_rho,0)
    dst_depth=__interp2_thread(src_grid.lon_rho, src_grid.lat_rho, src_depth,
                    child_grid.lon_rho, child_grid.lat_rho, pmap["pmaprho"],
                    weight, nx, ny, child_grid.mask_rho)
    # Interpolate the scalar fields
    records = np.arange(0, ncsrc.variables[time].shape[0]) \
                 if records is None else np.atleast_1d(records)
    for src in vmap:
        dest = vmap[src]
        # Only interpolate the fields we want in the destination
        if (dest not in ncout.variables ) or \
           ("rotate" in seapy.roms.fields[dest]):
                continue
        if seapy.roms.fields[dest]["dims"]==2:
            # Compute the max number of hold in memory
            maxrecs = np.minimum(len(records),
                np.int(_max_memory/(child_grid.lon_rho.nbytes +
                                    src_grid.lon_rho.nbytes)))
            for rn,recs in enumerate(seapy.chunker(records, maxrecs)):
                outr = np.s_[rn*maxrecs:np.minimum((rn+1)*maxrecs,len(records))]
                ndata = np.ma.array(Parallel(n_jobs=threads, verbose=2) \
                                 (delayed(__interp2_thread) (
                                        src_grid.lon_rho, src_grid.lat_rho,
                                        ncsrc.variables[src][i, :, :],
                                        child_grid.lon_rho, child_grid.lat_rho,
                                        pmap["pmaprho"], weight,
                                        nx, ny, child_grid.mask_rho)
                            for i in recs), copy=False)
                ncout.variables[dest][outr, :, :] = ndata
                ncout.sync()
        else:
            maxrecs = np.minimum(len(records), np.int(_max_memory /
                (child_grid.lon_rho.nbytes * child_grid.n +
                 src_grid.lon_rho.nbytes * src_grid.n)))
            for rn,recs in enumerate(seapy.chunker(records, maxrecs)):
                outr = np.s_[rn*maxrecs:np.minimum((rn+1)*maxrecs,len(records))]
                ndata = np.ma.array(Parallel(n_jobs=threads, verbose=2)
                                 (delayed(__interp3_thread)(
                                    src_grid.lon_rho, src_grid.lat_rho,
                                    src_grid.depth_rho,
                                    ncsrc.variables[src][i, :, :, :],
                                    child_grid.lon_rho, child_grid.lat_rho,
                                    child_grid.depth_rho,
                                    pmap["pmaprho"], weight,
                                    nx, ny, child_grid.mask_rho,
                                    up_factor=_up_scaling.get(dest, 1.0),
                                    down_factor=_down_scaling.get(dest, 1.0))
                            for i in recs), copy=False)
                if z_mask:
                    __mask_z_grid(ndata, dst_depth, child_grid.depth_rho)
                ncout.variables[dest][outr, :, :, :] = ndata
                ncout.sync()

    # Rotate and Interpolate the vector fields. First, determine which
    # are the "u" and the "v" vmap fields
    try:
        velmap = {
            "u": list(vmap.keys())[list(vmap.values()).index("u")],
            "v": list(vmap.keys())[list(vmap.values()).index("v")]}
    except:
        warn("velocity not present in source file")
        return

    srcangle = src_grid.angle if src_grid.cgrid else None
    dstangle = child_grid.angle if child_grid.cgrid else None
    maxrecs = np.minimum(len(records), np.int(_max_memory/
        (2 * (child_grid.lon_rho.nbytes * child_grid.n +
            src_grid.lon_rho.nbytes*src_grid.n))))
    for nr,recs in enumerate(seapy.chunker(records, maxrecs)):
        vel = Parallel(n_jobs=threads, verbose=2) \
                 (delayed(__interp3_vel_thread)(
                    src_grid.lon_rho, src_grid.lat_rho,
                    src_grid.depth_rho, srcangle,
                    ncsrc.variables[velmap["u"]][i, :, :, :],
                    ncsrc.variables[velmap["v"]][i, :, :, :],
                    child_grid.lon_rho, child_grid.lat_rho,
                    child_grid.depth_rho, dstangle,
                    pmap["pmaprho"], weight, nx, ny,
                    child_grid.mask_rho) for i in recs)

        for j in range(len(vel)):
            vel_u = np.ma.array(vel[j][0], copy=False)
            vel_v = np.ma.array(vel[j][1], copy=False)
            if z_mask:
                __mask_z_grid(vel_u, dst_depth, child_grid.depth_rho)
                __mask_z_grid(vel_v, dst_depth, child_grid.depth_rho)

            if child_grid.cgrid:
                vel_u = seapy.model.rho2u(vel_u)
                vel_v = seapy.model.rho2v(vel_v)

            ncout.variables["u"][nr*maxrecs+j, :] = vel_u
            ncout.variables["v"][nr*maxrecs+j, :] = vel_v

            if "ubar" in ncout.variables:
                # Create ubar and vbar
                # depth = seapy.adddim(child_grid.depth_u, vel_u.shape[0])
                ncout.variables["ubar"][nr*maxrecs+j, :] = \
                    np.sum(vel_u * child_grid.depth_u, axis=0) /  \
                    np.sum(child_grid.depth_u, axis=0)

            if "vbar" in ncout.variables:
                # depth = seapy.adddim(child_grid.depth_v, vel_v.shape[0])
                ncout.variables["vbar"][nr*maxrecs+j, :] = \
                    np.sum(vel_v * child_grid.depth_v, axis=0) /  \
                    np.sum(child_grid.depth_v, axis=0)

            ncout.sync()

    # Return the pmap that was used
    return pmap

def to_zgrid(roms_file, z_file, z_grid=None, depth=None, records=None,
             threads=2, nx=0, ny=0, weight=10, vmap=None, cdlfile=None, dims=2,
             pmap=None):
    """
    Given an existing ROMS history or average file, create (if does not exit)
    a new z-grid file. Use the given z_grid or otherwise build one with the
    same horizontal extent and the specified depths and interpolate the
    ROMS fields onto the z-grid.

    Parameters
    ----------
    roms_file  : string,
        File name of src file to interpolate from
    z_file : string,
        Name of desination file to write to
    z_grid: (string or seapy.model.grid), optional:
        Name or instance of output definition
    depth: numpy.ndarray, optional:
        array of depths to use for z-level
    records : numpy.ndarray, optional:
        Record indices to interpolate
    threads : int, optional:
        number of processing threads
    nx : float, optional:
        decorrelation length-scale for OA (same units as source data)
    ny : float, optional:
        decorrelation length-scale for OA (same units as source data)
    weight : int, optional:
        number of points to use in weighting matrix
    vmap : dictionary, optional
        mapping source and destination variables
    cdlfile : string, optional
        cdlfile to use for generating the z-file
    dims : int, optional
        number of dimensions to use for lat/lon arrays (default 2)
    pmap : numpy.ndarray, optional:
        use the specified pmap rather than compute it

    Returns
    -------
    pmap : ndarray
        the weighting matrix computed during the interpolation

    """
    roms_grid = seapy.model.asgrid(roms_file)
    ncroms = netCDF4.Dataset(roms_file)
    src_ref, time = seapy.roms.get_reftime(ncroms)
    records = np.arange(0, ncroms.variables[time].shape[0]) \
        if records is None else np.atleast_1d(records)

    # Load the grid
    if z_grid is not None:
        z_grid = seapy.model.asgrid(z_grid)
    elif os.path.isfile(z_file):
        z_grid = seapy.model.asgrid(z_file)

    if not os.path.isfile(z_file):
        if z_grid is None:
            lat = roms_grid.lat_rho.shape[0]
            lon = roms_grid.lat_rho.shape[1]
            if depth is None:
                raise ValueError("depth must be specified")
            ncout = seapy.roms.ncgen.create_zlevel(z_file, lat, lon,
                                    len(depth), src_ref, "ROMS z-level",
                                    cdlfile=cdlfile, dims=dims)
            if dims==1:
                ncout.variables["lat"][:] = roms_grid.lat_rho[:, 0]
                ncout.variables["lon"][:] = roms_grid.lon_rho[0, :]
            else:
                ncout.variables["lat"][:] = roms_grid.lat_rho
                ncout.variables["lon"][:] = roms_grid.lon_rho
            ncout.variables["depth"][:] = depth
            ncout.variables["mask"][:] = roms_grid.mask_rho
            ncout.sync()
            z_grid = seapy.model.grid(z_file)
        else:
            lat = z_grid.lat_rho.shape[0]
            lon = z_grid.lat_rho.shape[1]
            dims = z_grid.spatial_dims
            ncout = seapy.roms.ncgen.create_zlevel(z_file, lat, lon,
                                len(z_grid.z), src_ref, "ROMS z-level",
                                cdlfile=cdlfile, dims=dims)
            if dims==1:
                ncout.variables["lat"][:] = z_grid.lat_rho[:, 0]
                ncout.variables["lon"][:] = z_grid.lon_rho[0, :]
            else:
                ncout.variables["lat"][:] = z_grid.lat_rho
                ncout.variables["lon"][:] = z_grid.lon_rho
            ncout.variables["depth"][:] = z_grid.z
            ncout.variables["mask"][:] = z_grid.mask_rho
    else:
        ncout = netCDF4.Dataset(z_file, "a")

    ncout.variables["time"][:] = netCDF4.date2num(
            netCDF4.num2date(ncroms.variables[time][records],
                             ncroms.variables[time].units),
            ncout.variables["time"].units)
    ncroms.close()

    # Call the interpolation
    try:
        roms_grid.set_east(z_grid.east())
        pmap = __interp_grids(roms_grid, z_grid, ncout, records=records,
                  threads=threads, nx=nx, ny=ny, vmap=vmap, weight=weight,
                  z_mask=True, pmap=pmap)
    except TimeoutError:
        print("Timeout: process is hung, deleting output.")
        # Delete the output file
        os.remove(z_file)
    finally:
        # Clean up
        ncout.close()

    return pmap

def to_grid(src_file, dest_file, dest_grid=None, records=None, threads=2,
            nx=0, ny=0, weight=10, vmap=None, pmap=None):
    """
    Given an existing model file, create (if does not exit) a
    new ROMS history file using the given ROMS destination grid and
    interpolate the ROMS fields onto the new grid. If an existing destination
    file is given, it is interpolated onto the specified.

    Parameters
    ----------
    src_file  : string,
        Filename of src file to interpolate from
    dest_file : string,
        Name of desination file to write to
    dest_grid: (string or seapy.model.grid), optional:
        Name or instance of output definition
    records : numpy.ndarray, optional:
        Record indices to interpolate
    threads : int, optional:
        number of processing threads
    nx : float, optional:
        decorrelation length-scale for OA (same units as source data)
    ny : float, optional:
        decorrelation length-scale for OA (same units as source data)
    weight : int, optional:
        number of points to use in weighting matrix
    vmap : dictionary, optional
        mapping source and destination variables
    pmap : numpy.ndarray, optional:
        use the specified pmap rather than compute it

    Returns
    -------
    pmap : ndarray
        the weighting matrix computed during the interpolation
    """
    src_grid = seapy.model.asgrid(src_file)
    if dest_grid is not None:
        destg = seapy.model.asgrid(dest_grid)

        if not os.path.isfile(dest_file):
            ncsrc = netCDF4.Dataset(src_file)
            src_ref, time = seapy.roms.get_reftime(ncsrc)
            records = np.arange(0, ncsrc.variables[time].shape[0]) \
                 if records is None else np.atleast_1d(records)
            ncout = seapy.roms.ncgen.create_ini(dest_file,
                     eta_rho=destg.eta_rho, xi_rho=destg.xi_rho, s_rho=destg.n,
                     reftime=src_ref, title="interpolated from "+src_file)
            destg.to_netcdf(ncout)
            ncout.variables["ocean_time"][:] = netCDF4.date2num(
                    netCDF4.num2date(ncsrc.variables[time][records],
                                     ncsrc.variables[time].units),
                    ncout.variables["ocean_time"].units)
            ncsrc.close()

    if os.path.isfile(dest_file):
        ncout = netCDF4.Dataset(dest_file,"a")
        if dest_grid is None:
            destg = seapy.model.asgrid(dest_file)

    # Call the interpolation
    try:
        src_grid.set_east(destg.east())
        pmap = __interp_grids(src_grid, destg, ncout, records=records,
                              threads=threads, nx=nx, ny=ny, weight=weight,
                              vmap=vmap, pmap=pmap)
    except TimeoutError:
        print("Timeout: process is hung, deleting output.")
        # Delete the output file
        os.remove(dest_file)
    finally:
        # Clean up
        ncout.close()

    return pmap

def to_clim(src_file, dest_file, dest_grid=None, records=None, threads=2,
            nx=0, ny=0, weight=10, vmap=None, pmap=None):
    """
    Given an model output file, create (if does not exit) a
    new ROMS climatology file using the given ROMS destination grid and
    interpolate the ROMS fields onto the new grid. If an existing destination
    file is given, it is interpolated onto the specified.

    Parameters
    ----------
    src_file  : string,
        Filename of src file to interpolate from
    dest_file : string,
        Name of desination file to write to
    dest_grid: (string or seapy.model.grid), optional:
        Name or instance of output definition
    records : numpy.ndarray, optional:
        Record indices to interpolate
    threads : int, optional:
        number of processing threads
    nx : float, optional:
        decorrelation length-scale for OA (same units as source data)
    ny : float, optional:
        decorrelation length-scale for OA (same units as source data)
    weight : int, optional:
        number of points to use in weighting matrix
    vmap : dictionary, optional
        mapping source and destination variables
    pmap : numpy.ndarray, optional:
        use the specified pmap rather than compute it

    Returns
    -------
    pmap : ndarray
        the weighting matrix computed during the interpolation
    """
    if dest_grid is not None:
        destg = seapy.model.asgrid(dest_grid)
        src_grid = seapy.model.asgrid(src_file)
        ncsrc = netCDF4.Dataset(src_file)
        src_ref, time = seapy.roms.get_reftime(ncsrc)
        records = np.arange(0, ncsrc.variables[time].shape[0]) \
                 if records is None else np.atleast_1d(records)
        ncout=seapy.roms.ncgen.create_clim(dest_file,
                 eta_rho=destg.ln, xi_rho=destg.lm, s_rho=destg.n,
                 ntimes=records.size, reftime=src_ref,
                 title="interpolated from "+src_file)
        src_time = netCDF4.num2date(ncsrc.variables[time][records],
                                    ncsrc.variables[time].units)
        for dtime in ("zeta", "v2d", "v3d", "temp", "salt"):
            ncout.variables[dtime + "_time"][:] = netCDF4.date2num(src_time,
                                       ncout.variables[dtime + "_time"].units)
        ncsrc.close()
    else:
        raise AttributeError("you must supply a destination file or a grid to make the file")

    # Call the interpolation
    try:
        src_grid.set_east(destg.east())
        pmap = __interp_grids(src_grid, destg, ncout, records=records, threads=threads,
                  nx=nx, ny=ny, vmap=vmap, weight=weight, pmap=pmap)
    except TimeoutError:
        print("Timeout: process is hung, deleting output.")
        # Delete the output file
        os.remove(dest_file)
    finally:
        # Clean up
        ncout.close()
    return pmap
pass
