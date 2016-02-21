#!/usr/bin/env python
"""
  clim.py

  ROMS climatology utilities

  Written by Brian Powell on 08/15/15
  Copyright (c)2016 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import seapy
import numpy as np
import netCDF4

clim_times = ('zeta_time', 'v2d_time', 'v3d_time', 'temp_time', 'salt_time')


def gen_bry_clim(clim_file, grid, bry):
    """
    Taking the results of gen_ncks and interpolation, stitch together
    climatology files that were interpolated using only the boundary regions
    into a single climatology (with no data where interpolation wasn't
    performed).

    Parameters
    ----------
    clim_file: str,
        The name of the output climate file
    grid: seapy.model.grid or str,
        The output ROMS grid
    bry: dict,
        A dictionary prescribing the climatology file interpolated for each
        boundary side.
        {"west":filename, "south":filename}, ...}

    Returns
    -------
        None
    """
    grid = seapy.model.asgrid(grid)

    # Grab the first dictionary record and use it to determine the number
    # of times in the new climatology file
    nc = netCDF4.Dataset(bry[list(bry.keys())[0]])
    reftime, time = seapy.roms.get_reftime(nc)
    times = nc.variables[time][:]
    nc.close()

    # Create the new climatology file
    ncout = seapy.roms.ncgen.create_clim(clim_file,
                                         eta_rho=grid.ln, xi_rho=grid.lm, s_rho=grid.n,
                                         ntimes=len(times), reftime=reftime,
                                         title="stitched from boundary interpolation")
    for dtime in ("zeta", "v2d", "v3d", "temp", "salt"):
        ncout.variables[dtime + "_time"][:] = times

    for side in bry:
        if bry[side] is None:
            continue
        ncin = netCDF4.Dataset(bry[side])
        for fld in seapy.roms.fields:
            idx = [np.s_[:] for i in range(seapy.roms.fields[fld]["dims"] + 1)]
            dat = ncin.variables[fld][:]
            shp = dat.shape
            if side == "west":
                idx[-1] = np.s_[:shp[-1]]
                pass
            elif side == "east":
                idx[-1] = np.s_[-shp[-1]:]
                pass
            elif side == "north":
                idx[-2] = np.s_[-shp[-2]:]
                pass
            elif side == "south":
                idx[-2] = np.s_[:shp[-2]]
                pass
            ncout.variables[fld][idx] = dat
            ncout.sync()
        ncin.close()
    ncout.close()


def nc_concat(clim_files, out_file, maxrecs=20):
    """
    Since climatology files have multiple non-record dimensions, they are
    difficult to work with using standard netcdf operations. This routine
    takes a list of climatology files and creates a new one with the list
    concatenated together.

    Parameters
    ----------
    clim_files: list,
        List of climatology files to concatenate
    out_file: string,
        Name of output file
    maxrecs: int,
        For large records, we can only deal with so many at each chunk.
        Define how many to grab in one pass.

    Returns
    -------
    None
    """
    # Assume that zeta_time is consistent with all times
    clim_time = []
    for f in clim_files:
        nc = netCDF4.Dataset(f)
        t = nc.variables["zeta_time"][:]
        units = nc.variables["zeta_time"].units
        eta_rho = len(nc.dimensions["eta_rho"])
        xi_rho = len(nc.dimensions["xi_rho"])
        s_rho = len(nc.dimensions["s_rho"])
        nc.close()
        clim_time.append(t)
    reftime = netCDF4.num2date(0, units)

    # Figure out the unique records to use from each file so we
    # don't have overlap
    recs = [list(range(len(i))) for i in clim_time]
    filenum = [np.ones(len(i)) * n for n, i in enumerate(clim_time)]
    final_t, idx = np.unique(np.hstack(clim_time), return_index=True)
    recs = np.hstack(recs)[idx]
    filenum = np.hstack(filenum)[idx]

    # Create the output file
    ncout = seapy.roms.ncgen.create_clim(out_file, eta_rho=eta_rho,
                                         xi_rho=xi_rho, s_rho=s_rho, ntimes=len(
                                             final_t),
                                         reftime=reftime, clobber=True,
                                         title='concatenated: ' + ', '.join(clim_files))

    try:
        # Add all of the times (if it isn't in the file, don't
        # let it be fatal)
        out_t = netCDF4.date2num(netCDF4.num2date(final_t, units),
                                 ncout.variables["zeta_time"].units)
        for t in clim_times:
            try:
                ncout.variables[t][:] = out_t
            except KeyError:
                pass

        # Go over all files and concatenate the appropriate records
        for n, f in enumerate(seapy.progressbar.progress(clim_files)):
            nc = netCDF4.Dataset(f)
            for fld in seapy.roms.fields:
                out = np.where(filenum == n)[0]
                try:
                    # These can be very large, so we need to chunk through
                    # the 3D variables
                    if seapy.roms.fields[fld]['dims'] == 3:
                        for j in seapy.chunker(out, maxrecs):
                            ncout.variables[fld][j, :] = \
                                nc.variables[fld][recs[j], :]
                    else:
                        ncout.variables[fld][out, :] = \
                            nc.variables[fld][recs[out], :]
                    ncout.sync()
                except KeyError:
                    pass
            nc.close()
    finally:
        ncout.close()
    pass
