#!/usr/bin/env python
"""
  boundary.py

  ROMS boundary utilities

  Written by Brian Powell on 01/15/14
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import seapy
import numpy as np
import netCDF4
import textwrap
from collections import namedtuple

# Define the sides of ROMS boundaries along with the DA ordering
__side_info = namedtuple("__side_info", "indices order xi")
sides = {"west": __side_info((np.s_[:], 0), 1, False),
         "south": __side_info((0, np.s_[:]), 2, True),
         "east": __side_info((np.s_[:], -1), 3, False),
         "north": __side_info((-1, np.s_[:]), 4, True)}


def from_roms(roms_file, bry_file, grid=None, records=None,
              clobber=False, cdl=None):
    """
    Given a ROMS history, average, or climatology file, generate
    boundary conditions on the same grid.

    Parameters
    ----------
    roms_file : string or list,
        ROMS source (history, average, climatology file)
    bry_file : string,
        output boundary file,
    grid : seapy.model.grid or string, optional,
        ROMS grid for boundaries
    records : array, optional,
        record indices to put into the boundary
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.

    Returns
    -------
    None

    """
    if grid is None:
        grid = seapy.model.asgrid(roms_file)
    else:
        grid = seapy.model.asgrid(grid)
    ncroms = seapy.netcdf(roms_file)
    src_ref, time = seapy.roms.get_reftime(ncroms)
    records = np.arange(0, len(ncroms.variables[time][:])) \
        if records is None else records

    # Create the boundary file and fill up the descriptive data
    ncbry = seapy.roms.ncgen.create_bry(bry_file,
                                        eta_rho=grid.eta_rho, xi_rho=grid.xi_rho,
                                        s_rho=grid.n, reftime=src_ref,
                                        cdl=cdl, clobber=clobber,
                                        title="generated from " + roms_file)
    brytime = seapy.roms.get_timevar(ncbry)
    grid.to_netcdf(ncbry)
    ncbry.variables[brytime][:] = seapy.roms.date2num(
        seapy.roms.num2date(ncroms, time, records), ncbry, brytime)

    for var in seapy.roms.fields:
        if var in ncroms.variables:
            for bry in sides:
                ndim = seapy.roms.fields[var]["dims"]
                if ndim == 3:
                    ncbry.variables["_".join((var, bry))][:] = \
                        ncroms.variables[var][records, :,
                                              sides[bry].indices[0],
                                              sides[bry].indices[1]]
                elif ndim == 2:
                    ncbry.variables["_".join((var, bry))][:] = \
                        ncroms.variables[var][records,
                                              sides[bry].indices[0],
                                              sides[bry].indices[1]]
    ncbry.close()
    pass


def gen_ncks(parent_file, grid, sponge=0, pad=3):
    """
    Create the ncks commands for extracting a section of a global model
    that encompasses each of the boundaries of the given grid. The reason
    for this is that often we end up with massive global files, and we do
    not want to interpolate the entirety onto our regional grids (if our
    regional grid is large also), but only the boundary and nudge/sponge
    region.

    This script simply does the computation and outputs the necessary `ncks`
    commands that the user needs to produce global boundary region and
    boundary "grid" for the regional grid to then use in the interpolation
    (e.g., seapy.interp.to_clim). This is purely to save disk and cpu
    expense, and it is non-trivial to use.

    Parameters
    ----------
    parent_file : seapy.model.grid or string,
        Parent file (HYCOM, etc.)
    grid : seapy.model.grid or string, optional,
        ROMS grid for boundaries
    sponge : int,
        Width to extract along each boundary. If 0, only the boundary itself
        will be extracted.
    pad : int,
        Additional rows/columns to extract from parent around the region

    Returns
    -------
    None
    """
    import re

    parent_grid = seapy.model.asgrid(parent_file)
    child_grid = seapy.model.asgrid(grid)
    fre = re.compile('.nc')
    # Make sure we are on the same coordinates
    if parent_grid.east() != child_grid.east():
        if child_grid.east():
            parent_grid.lon_rho[parent_grid.lon_rho < 0] += 360
        else:
            parent_grid.lon_rho[parent_grid.lon_rho > 180] -= 360

    # Loop over each side of the grid and determine the indices from the
    # parent and child
    for side in sides:
        # Figure out which dimension this side is on and determine all
        # of the indices needed
        idx = sides[side].indices
        if isinstance(idx[0], int):
            pdim = parent_grid.key["lat_rho"]
            cdim = "eta"
            if idx[0] == -1:
                idx = np.s_[-(sponge + 2):, :]
                # pdidx = "{:d},".format(parent_grid.lat_rho.shape[0]-sponge-2)
                cdidx = "{:d},".format(child_grid.eta_rho - sponge - 2)
                pass
            else:
                idx = np.s_[:sponge + 1, :]
                cdidx = ",{:d}".format(sponge + 1)
            l = np.where(np.logical_and(
                parent_grid.lat_rho >= np.min(child_grid.lat_rho[idx]),
                parent_grid.lat_rho <= np.max(child_grid.lat_rho[idx])))
            i = np.maximum(0, np.min(l[0]) - pad)
            j = np.minimum(parent_grid.lat_rho.shape[0],
                           np.max(l[0]) + pad + 1)
        else:
            pdim = parent_grid.key["lon_rho"]
            cdim = "xi"
            if idx[1] == -1:
                idx = np.s_[:, -(sponge + 2):]
                cdidx = "{:d},".format(child_grid.xi_rho - sponge - 2)
            else:
                idx = np.s_[:, :sponge + 1]
                cdidx = ",{:d}".format(sponge + 1)
            l = np.where(np.logical_and(
                parent_grid.lon_rho >= np.min(child_grid.lon_rho[idx]),
                parent_grid.lon_rho <= np.max(child_grid.lon_rho[idx])))
            i = np.maximum(0, np.min(l[1]) - pad)
            j = np.minimum(parent_grid.lon_rho.shape[1],
                           np.max(l[1]) + pad + 1)

        # Put the indices together into strings for extracting out new files
        pdidx = "{:d},{:d}".format(i, j)

        # Display the commands:
        cmd = "ncks"
        pfiles = "{:s} {:s}".format(parent_grid.filename,
                                    fre.sub("_{:s}.nc".format(side),
                                            parent_grid.filename))
        cfiles = "{:s} {:s}".format(child_grid.filename,
                                    fre.sub("_{:s}.nc".format(side),
                                            child_grid.filename))

        grids = ('rho', 'u', 'v', 'psi')
        if parent_grid.cgrid:
            pdim = ' -d'.join(["{:s}_{:s},{:s}".format(pdim, k, pdidx)
                               for k in grids])
        else:
            pdim = "{:s},{:s}".format(pdim, pdidx)

        if child_grid.cgrid:
            cdim = ' -d'.join(["{:s}_{:s},{:s}".format(cdim, k, cdidx)
                               for k in grids])
        else:
            cdim = "{:s},{:s}".format(cdim, cdidx)

        print("-" * 40 + "\n" + side + "\n" + "-" * 40)
        print("{:s} -O -d{:s} {:s}".format(cmd, pdim, pfiles))
        print("{:s} -O -d{:s} {:s}".format(cmd, cdim, cfiles))

    pass


def from_std(std_filename, bry_std_file, fields=None, clobber=False, cdl=None):
    """
    Generate the boundary standard deviations file for 4DVAR from the
    standard deviation of a boundary file. Best to use nco:

    $ ncwa -a bry_time roms_bry_file tmp.nc

    $ ncbo -O -y sub roms_bry_file tmp.nc tmp.nc

    $ ncra -y rmssdn tmp.nc roms_bry_std.nc

    to generate the standard deviations. This method simply takes the
    standard deviations of the boundaries and puts them into the
    proper format for ROMS 4DVAR.

    Parameters
    ----------
    std_filename : string or list,
        Filename of the boundary standard deviation file
    bry_std_file : string,
        Filename of the boundary standard deviations file to create
    fields : list, optional,
        ROMS fields to generate boundaries for. The default are the
        standard fields as defined in seapy.roms.fields
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.

    Returns
    -------
    None
    """
    ncstd = seapy.netcdf(std_filename)
    eta_rho = len(ncstd.dimensions["eta_rho"])
    xi_rho = len(ncstd.dimensions["xi_rho"])
    s_rho = len(ncstd.dimensions["s_rho"])
    ncout = seapy.roms.ncgen.create_da_bry_std(bry_std_file,
                                               eta_rho=eta_rho, xi_rho=xi_rho,
                                               s_rho=s_rho, clobber=clobber, cdl=cdl,
                                               title='STD from ' + std_filename)
    ncout.variables["ocean_time"][:] = ncstd.variables["bry_time"][0]

    if fields is None:
        fields = seapy.roms.fields

    # Go over every side for every field and put it together
    for var in fields:
        vname = var + "_obc"
        if vname not in ncout.variables:
            ncout.createVariable(vname, np.float32,
                                 ('ocean_time', "boundary", "s_rho", "IorJ"))
        ndat = np.zeros(ncout.variables[vname].shape)
        for bry in sides:
            order = sides[bry].order - 1
            dat = ncstd.variables[var + "_" + bry][0, :]
            if dat.ndim == 1:
                ndat[0, order, :len(dat)] = dat
            else:
                ndat[0, order, :, :dat.shape[1]] = dat
        ncout.variables[vname][:] = ndat
        ncout.sync()
    pass


def gen_stations(filename, grid):
    """
    Generate a station file with stations at every boundary location for use in
    nesting one grid within another.

    Parameters
    ----------
    filename: string
        Input name of station file to create
    grid: string or seapy.model.grid
        Input grid to generate station file from. If string, it will open
        the grid file. If grid, it will use the grid information

    Returns
    -------
    None

    """
    grid = seapy.model.asgrid(grid)

    # Put together the boundaries
    lon = np.concatenate([grid.lon_rho[0, :], grid.lon_rho[-1, :],
                          grid.lon_rho[:, 0], grid.lon_rho[:, -1]])
    lat = np.concatenate([grid.lat_rho[0, :], grid.lat_rho[-1, :],
                          grid.lat_rho[:, 0], grid.lat_rho[:, -1]])
    Npts = len(lon)

    header = """\

    ! Switch to control the writing of stations data within nested and/or multiple
    ! connected grids, [1:Ngrids].

       Lstations == T

    ! Logical switches (TRUE/FALSE) to activate writing of fields in STATION
    ! output file, [Sout(:,ng), ng=1, Ngrids].

    Sout(idUvel) == T       ! u                  3D U-velocity
    Sout(idVvel) == T       ! v                  3D V-velocity
    Sout(idWvel) == F       ! w                  3D W-velocity
    Sout(idOvel) == F       ! omega              3D omega vertical velocity
    Sout(idUbar) == T       ! ubar               2D U-velocity
    Sout(idVbar) == T       ! vbar               2D V-velocity
    Sout(idFsur) == T       ! zeta               free-surface
    Sout(idBath) == F       ! bath               time-dependent bathymetry

    Sout(idTvar) == T T     ! temp, salt, ...    all (NT) tracers

    Sout(idUsms) == F       ! sustr              surface U-stress
    Sout(idVsms) == F       ! svstr              surface V-stress
    Sout(idUbms) == F       ! bustr              bottom U-stress
    Sout(idVbms) == F       ! bvstr              bottom V-stress

    Sout(idUbrs) == F       ! bustrc             bottom U-current stress
    Sout(idVbrs) == F       ! bvstrc             bottom V-current stress
    Sout(idUbws) == F       ! bustrw             bottom U-wave stress
    Sout(idVbws) == F       ! bvstrw             bottom V-wave stress
    Sout(idUbcs) == F       ! bustrcwmax         bottom max wave-current U-stress
    Sout(idVbcs) == F       ! bvstrcwmax         bottom max wave-current V-stress

    Sout(idUbot) == F       ! Ubot               bed wave orbital U-velocity
    Sout(idVbot) == F       ! Vbot               bed wave orbital V-velocity
    Sout(idUbur) == F       ! Ur                 bottom U-velocity above bed
    Sout(idVbvr) == F       ! Vr                 bottom V-velocity above bed

    Sout(idW2xx) == F       ! Sxx_bar            2D radiation stress, Sxx component
    Sout(idW2xy) == F       ! Sxy_bar            2D radiation stress, Sxy component
    Sout(idW2yy) == F       ! Syy_bar            2D radiation stress, Syy component
    Sout(idU2rs) == F       ! Ubar_Rstress       2D radiation U-stress
    Sout(idV2rs) == F       ! Vbar_Rstress       2D radiation V-stress
    Sout(idU2Sd) == F       ! ubar_stokes        2D U-Stokes velocity
    Sout(idV2Sd) == F       ! vbar_stokes        2D V-Stokes velocity

    Sout(idW3xx) == F       ! Sxx                3D radiation stress, Sxx component
    Sout(idW3xy) == F       ! Sxy                3D radiation stress, Sxy component
    Sout(idW3yy) == F       ! Syy                3D radiation stress, Syy component
    Sout(idW3zx) == F       ! Szx                3D radiation stress, Szx component
    Sout(idW3zy) == F       ! Szy                3D radiation stress, Szy component
    Sout(idU3rs) == F       ! u_Rstress          3D U-radiation stress
    Sout(idV3rs) == F       ! v_Rstress          3D V-radiation stress
    Sout(idU3Sd) == F       ! u_stokes           3D U-Stokes velocity
    Sout(idV3Sd) == F       ! v_stokes           3D V-Stokes velocity

    Sout(idWamp) == F       ! Hwave              wave height
    Sout(idWlen) == F       ! Lwave              wave length
    Sout(idWdir) == F       ! Dwave              wave direction
    Sout(idWptp) == F       ! Pwave_top          wave surface period
    Sout(idWpbt) == F       ! Pwave_bot          wave bottom period
    Sout(idWorb) == F       ! Ub_swan            wave bottom orbital velocity
    Sout(idWdis) == F       ! Wave_dissip        wave dissipation

    Sout(idPair) == F       ! Pair               surface air pressure
    Sout(idUair) == F       ! Uair               surface U-wind component
    Sout(idVair) == F       ! Vair               surface V-wind component

    Sout(idTsur) == F F     ! shflux, ssflux     surface net heat and salt flux
    Sout(idLhea) == F       ! latent             latent heat flux
    Sout(idShea) == F       ! sensible           sensible heat flux
    Sout(idLrad) == F       ! lwrad              longwave radiation flux
    Sout(idSrad) == F       ! swrad              shortwave radiation flux
    Sout(idEmPf) == F       ! EminusP            E-P flux
    Sout(idevap) == F       ! evaporation        evaporation rate
    Sout(idrain) == F       ! rain               precipitation rate

    Sout(idDano) == F       ! rho                density anomaly
    Sout(idVvis) == F       ! AKv                vertical viscosity
    Sout(idTdif) == F       ! AKt                vertical T-diffusion
    Sout(idSdif) == F       ! AKs                vertical Salinity diffusion
    Sout(idHsbl) == F       ! Hsbl               depth of surface boundary layer
    Sout(idHbbl) == F       ! Hbbl               depth of bottom boundary layer
    Sout(idMtke) == F       ! tke                turbulent kinetic energy
    Sout(idMtls) == F       ! gls                turbulent length scale

    ! Logical switches (TRUE/FALSE) to activate writing of exposed sediment
    ! layer properties into STATIONS output file.  Currently, MBOTP properties
    ! are expected for the bottom boundary layer and/or sediment models:
    !
    ! idBott( 1=isd50)   grain_diameter          mean grain diameter
    ! idBott( 2=idens)   grain_density           mean grain density
    ! idBott( 3=iwsed)   settling_vel            mean settling velocity
    ! idBott( 4=itauc)   erosion_stres           critical erosion stress
    ! idBott( 5=irlen)   ripple_length           ripple length
    ! idBott( 6=irhgt)   ripple_height           ripple height
    ! idBott( 7=ibwav)   bed_wave_amp            wave excursion amplitude
    ! idBott( 8=izdef)   Zo_def                  default bottom roughness
    ! idBott( 9=izapp)   Zo_app                  apparent bottom roughness
    ! idBott(10=izNik)   Zo_Nik                  Nikuradse bottom roughness
    ! idBott(11=izbio)   Zo_bio                  biological bottom roughness
    ! idBott(12=izbfm)   Zo_bedform              bed form bottom roughness
    ! idBott(13=izbld)   Zo_bedload              bed load bottom roughness
    ! idBott(14=izwbl)   Zo_wbl                  wave bottom roughness
    ! idBott(15=iactv)   active_layer_thickness  active layer thickness
    ! idBott(16=ishgt)   saltation               saltation height
    !
    !                                 1 1 1 1 1 1 1
    !               1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6

    Sout(idBott) == F F F F F F F F F F F F F F F F

    ! Number of stations to process in each nested grid.  These values are
    ! essential because the station arrays are dynamically allocated using
    ! these values, [1:Ngrids].

    """
    stations = """
    ! Station locations for all grids in any desired order.  The horizontal
    ! location for a particular station may be specified in terms of fractional
    ! (I,J) grid pairs (FLAG=0) or (longitude,latitude) grid pairs (FLAG=1).
    ! Here, FLAG is a special switch and may be used for multiple purposes.
    ! The GRID column indicates nested grid number to process. This value must
    ! be one in non-nested applications.  The COMMENT section is ignored during
    ! reading and may be used to help documentation.

    POS =  GRID  FLAG      X-POS       Y-POS     COMMENT
    """
    with open(filename, "w") as text_file:
        print("! BOUNDARY STATIONS FOR GRID: {}".format(grid.filename),
              file=text_file)
        print(textwrap.dedent(header), file=text_file)
        print("        NSTATION ==  {}".format(Npts), file=text_file)
        print(textwrap.dedent(stations), file=text_file)
        for i in range(Npts):
            print("        1     1    {0:10.6f}   {1:10.6f}   BRY".format(
                lon[i], lat[i]), file=text_file)

    pass


def from_stations(station_file, bry_file, grid=None):
    """
    Construct a boundary forcing file from a stations file generated by a parent-grid.
    The stations.in file must have been generated by the seapy.roms.gen_stations method;
    otherwise, the order will be incorrect.

    Parameters
    ==========
    station_file : string
        Filename of the stations file that is the source for the boundary data
    bry_file : string
        Filename of the boundary conditions file to generate
    grid : string or seapy.model.grid
        Grid that the boundary conditions are created for

    Returns
    -------
    None

    """
    grid = seapy.model.asgrid(grid)
    ncstation = netCDF4.Dataset(station_file)
    src_ref, time = seapy.roms.get_reftime(ncstation)

    # Create the boundary file and fill up the descriptive data
    ncbry = seapy.roms.ncgen.create_bry(bry_file,
                                        eta_rho=grid.eta_rho, xi_rho=grid.xi_rho,
                                        s_rho=grid.n, reftime=src_ref, clobber=False,
                                        title="generated from " + station_file)
    grid.to_netcdf(ncbry)

    # Load the times: we need to see if the times are duplicated
    # because if using assimilation, they may be duplicated for every
    # outer-loop. Currently, it overwrites the first one each time (this
    # will need to be fixed if ROMS is fixed).
    statime = ncstation.variables[time][:]
    dup = np.where(statime[1:] == statime[0])[0]
    rng = np.s_[:]
    if dup.size > 0:
        rng = np.s_[0:np.min(dup)]
        statime = statime[rng]
    brytime = seapy.roms.get_timevar(ncbry)
    ncbry.variables[brytime][:] = seapy.roms.date2num(
        seapy.roms.num2date(ncstation, time, rng), ncbry, brytime)

    # Set up the indices
    bry = {
        "south": range(0, grid.lm),
        "north": range(grid.lm, 2 * grid.lm),
        "west": range(2 * grid.lm, 2 * grid.lm + grid.ln),
        "east": range(2 * grid.lm + grid.ln, 2 * (grid.lm + grid.ln))
    }

    # Get the information to construct the depths of the station data
    sta_vt = ncstation.variables["Vtransform"][:]
    sta_hc = ncstation.variables["hc"][:]
    sta_s_rho = ncstation.variables["s_rho"][:]
    sta_cs_r = ncstation.variables["Cs_r"][:]
    sta_h = ncstation.variables["h"][:]
    sta_angle = ncstation.variables["angle"][:]
    sta_lon = ncstation.variables["lon_rho"][:]
    sta_lat = ncstation.variables["lat_rho"][:]
    sta_mask = np.ones(sta_lat.shape)
    sta_mask[sta_lon * sta_lat > 1e10] = 0

    # Load the station data as we need to manipulate it
    sta_zeta = np.ma.masked_greater(ncstation.variables["zeta"][rng], 100)
    sta_ubar = np.ma.masked_greater(ncstation.variables["ubar"][rng], 100)
    sta_vbar = np.ma.masked_greater(ncstation.variables["vbar"][rng], 100)
    sta_temp = np.ma.masked_greater(ncstation.variables["temp"][rng], 100)
    sta_salt = np.ma.masked_greater(ncstation.variables["salt"][rng], 100)
    sta_u = np.ma.masked_greater(ncstation.variables["u"][rng], 100)
    sta_v = np.ma.masked_greater(ncstation.variables["v"][rng], 100)
    ncstation.close()

    # Create the true positions and mask
    grid_h = np.concatenate([grid.h[0, :], grid.h[-1, :],
                             grid.h[:, 0], grid.h[:, -1]])
    grid_lon = np.concatenate([grid.lon_rho[0, :], grid.lon_rho[-1, :],
                               grid.lon_rho[:, 0], grid.lon_rho[:, -1]])
    grid_lat = np.concatenate([grid.lat_rho[0, :], grid.lat_rho[-1, :],
                               grid.lat_rho[:, 0], grid.lat_rho[:, -1]])
    grid_mask = np.concatenate([grid.mask_rho[0, :], grid.mask_rho[-1, :],
                                grid.mask_rho[:, 0], grid.mask_rho[:, -1]])
    grid_angle = np.concatenate([grid.angle[0, :], grid.angle[-1, :],
                                 grid.angle[:, 0], grid.angle[:, -1]])

    # Search for bad stations due to child grid overlaying parent mask.
    # Unfortunately, ROMS will give points that are not at the locations
    # you specify if those points conflict with the mask. So, these points
    # are simply replaced with the nearest.
    dist = np.sqrt((sta_lon - grid_lon)**2 + (sta_lat - grid_lat)**2)
    bad_pts = np.where(np.logical_and(dist > 0.001, grid_mask == 1))[0]
    good_pts = np.where(np.logical_and(dist < 0.001, grid_mask == 1))[0]
    for i in bad_pts:
        didx = np.sqrt((sta_lon[i] - sta_lon[good_pts])**2 +
                       (sta_lat[i] - sta_lat[good_pts])**2).argmin()
        index = good_pts[didx]
        sta_h[i] = sta_h[index]
        sta_angle[i] = sta_angle[index]
        sta_lon[i] = sta_lon[index]
        sta_lat[i] = sta_lat[index]
        sta_zeta[:, i] = sta_zeta[:, index]
        sta_ubar[:, i] = sta_ubar[:, index]
        sta_vbar[:, i] = sta_vbar[:, index]
        sta_temp[:, i, :] = sta_temp[:, index, :]
        sta_salt[:, i, :] = sta_salt[:, index, :]
        sta_u[:, i, :] = sta_u[:, index, :]
        sta_v[:, i, :] = sta_v[:, index, :]

    # Construct the boundaries: a dictionary of boundary side and two element
    # array whether the u[0] or v[1] dimensions need to be averaged
    sides = {"north": [True, False], "south": [True, False],
             "east": [False, True], "west": [False, True]}
    delta_angle = sta_angle - grid_angle
    sta_ubar, sta_vbar = seapy.rotate(sta_ubar, sta_vbar, delta_angle)
    sta_u, sta_v = seapy.rotate(sta_u, sta_v, np.tile(delta_angle,
                                                      (sta_u.shape[-1], 1)).T)

    # Set up the parameters for depth-interpolated
    wght = 5
    nx = 3
    ny = 9

    # Build a non-extrapolating field to interpolate. Generate the
    # position and depth
    def __expand_field(x):
        shp = x.shape
        y = np.zeros((shp[0] + 2, shp[1] + 2))
        y[1:-1, 1:-1] = x
        y[1:-1, 0] = x[:, 0]
        y[1:-1, -1] = x[:, -1]
        y[0, :] = y[1, :]
        y[-1, :] = y[-2, :]
        return y

    for side in sides:
        print(side)

        # Masks
        sta_ocean = np.where(sta_mask[bry[side]] == 1)[0]
        ocean = np.where(grid_mask[bry[side]] == 1)[0]

        # If we have a masked boundary, skip it
        if not np.any(ocean):
            continue

        # 1) Zeta
        ncbry.variables["zeta_" + side][:,
                                        ocean] = sta_zeta[:, bry[side]][:, ocean]

        # 2) Ubar
        if sides[side][0]:
            ncbry.variables["ubar_" + side][:] = 0.5 * (
                sta_ubar[:, bry[side][0:-1]] + sta_ubar[:, bry[side][1:]])
        else:
            ncbry.variables["ubar_" + side][:] = sta_ubar[:, bry[side]]

        # 3) Vbar
        if sides[side][1]:
            ncbry.variables["vbar_" + side][:] = 0.5 * (
                sta_vbar[:, bry[side][0:-1]] + sta_vbar[:, bry[side][1:]])
        else:
            ncbry.variables["vbar_" + side][:] = sta_vbar[:, bry[side]]

        # For 3D variables, we need to loop through time and interpolate
        # onto the child grid. Construct the distances
        x = np.zeros(len(bry[side]))
        x[1:] = np.cumsum(seapy.earth_distance(grid_lon[bry[side][0:-1]],
                                               grid_lat[bry[side][0:-1]],
                                               grid_lon[bry[side][1:]],
                                               grid_lat[bry[side][1:]]))
        sta_x = seapy.adddim(x, len(sta_s_rho))
        x = seapy.adddim(x, len(grid.s_rho))

        for n, t in seapy.progressbar.progress(enumerate(statime), statime.size):
            sta_depth = seapy.roms.depth(sta_vt, sta_h[bry[side]], sta_hc,
                                         sta_s_rho, sta_cs_r, sta_zeta[n, bry[side]])
            depth = seapy.roms.depth(grid.vtransform, grid_h[bry[side]],
                                     grid.hc, grid.s_rho, grid.cs_r, sta_zeta[n, bry[side]])

            in_x = __expand_field(sta_x[:, sta_ocean])
            in_x[:, 0] = in_x[:, 0] - 3600
            in_x[:, -1] = in_x[:, -1] + 3600
            in_depth = __expand_field(sta_depth[:, sta_ocean])
            in_depth[0, :] = in_depth[0, :] - 1000
            in_depth[-1, :] = in_depth[-1, :] + 10

            # 4) Temp
            in_data = __expand_field(np.transpose(
                sta_temp[n, bry[side], :][sta_ocean, :]))
            ncbry.variables["temp_" + side][n, :] = 0.0
            ncbry.variables["temp_" + side][n, :, ocean], pmap = seapy.oa.oasurf(
                in_x, in_depth, in_data,
                x[:, ocean], depth[:, ocean], nx=nx, ny=ny, weight=wght)

            # 5) Salt
            in_data = __expand_field(np.transpose(
                sta_salt[n, bry[side], :][sta_ocean, :]))
            ncbry.variables["salt_" + side][n, :] = 0.0
            ncbry.variables["salt_" + side][n, :, ocean], pmap = seapy.oa.oasurf(
                in_x, in_depth, in_data,
                x[:, ocean], depth[:, ocean], pmap=pmap, nx=nx, ny=ny, weight=wght)

            # 6) U
            in_data = __expand_field(np.transpose(
                sta_u[n, bry[side], :][sta_ocean, :]))
            data = np.zeros(x.shape)
            data[:, ocean], pmap = seapy.oa.oasurf(in_x, in_depth, in_data,
                                                   x[:, ocean],
                                                   depth[:, ocean],
                                                   pmap=pmap, nx=nx, ny=ny, weight=wght)
            if sides[side][0]:
                ncbry.variables["u_" + side][n, :] = 0.5 * (
                    data[:, 0:-1] + data[:, 1:])
            else:
                ncbry.variables["u_" + side][n, :] = data

            # 7) V
            in_data = __expand_field(np.transpose(
                sta_v[n, bry[side], :][sta_ocean, :]))
            data = data * 0
            data[:, ocean], pmap = seapy.oa.oasurf(in_x, in_depth, in_data,
                                                   x[:, ocean],
                                                   depth[:, ocean],
                                                   pmap=pmap, nx=nx, ny=ny,
                                                   weight=wght)
            if sides[side][1]:
                ncbry.variables["v_" + side][n, :] = 0.5 * (
                    data[:, 0:-1] + data[:, 1:])
            else:
                ncbry.variables["v_" + side][n, :] = data
            ncbry.sync()
    ncbry.close()
    pass


def detide(grid, bryfile, tidefile, tides=None, tide_start=None):
    """
    Given a boundary file, detide the barotropic components and create tidal
    forcing file for the grid. This method will update the given boundary file.

    Parameters
    ----------
    grid : seapy.model.grid or string,
       The grid that defines the boundaries shape and mask
    bryfile : string,
       The boundary file to detide
    tidefile : string,
       The output tidal forcing file with the tide spectral forcing
    tides : string array, optional
       Array of strings defining which tides to extract. Defaults to the
       standard 11 constituents.
    tide_start : datetime, optional
       The reference date to use for the tide forcing. If None, the
       center of the time period is used.

    Returns
    -------
    None

    Examples
    --------
    Make a long time-series boundary conditions from a group of boundary files,
    skipping the last record of each file to prevent overlap (if there are 100 records
    in each file). Afterwards, detide the resulting file.

    >>> !ncrcat -dbry_time,0,,100,99 bry_out_*nc bry_detide.nc
    >>> seapy.roms.boundary.detide("mygrid.nc", "bry_detide.nc", "tide_out.nc")

    """
    import datetime

    if not tides:
        tides = seapy.tide.default_tides
    else:
        tides = np.atleast_1d(tides)

    # Load Files
    grid = seapy.model.asgrid(grid)
    bry = netCDF4.Dataset(bryfile, "a")

    # Get the definitions of the boundary file
    epoch, timevar = seapy.roms.get_reftime(bry)
    time = seapy.roms.num2date(bry, timevar)

    # Pick the time for the tide file reference
    if not tide_start:
        tide_start = time[0] + (time[-1] - time[0]) / 2
        tide_start = datetime.datetime(
            tide_start.year, tide_start.month, tide_start.day)

    try:
        s_rho = len(bry.dimensions['s_rho'])
    except:
        s_rho = grid.n

    # Set variables to detide
    detide_vars = ['zeta', 'ubar', 'vbar']

    # Create the tide forcing file
    bry.detide = "Detided to generate tide forcing: {:s}".format(tidefile)

    # Detide the free-surface
    eamp = np.zeros((len(tides), grid.eta_rho, grid.xi_rho))
    epha = np.zeros((len(tides), grid.eta_rho, grid.xi_rho))
    cmin = np.zeros((len(tides), grid.eta_rho, grid.xi_rho))
    cmax = np.zeros((len(tides), grid.eta_rho, grid.xi_rho))
    cang = np.zeros((len(tides), grid.eta_rho, grid.xi_rho))
    cpha = np.zeros((len(tides), grid.eta_rho, grid.xi_rho))

    for side in sides:
        lvar = "zeta_" + side
        idx = sides[side].indices
        lat = grid.lat_rho[idx[0], idx[1]]
        size = grid.xi_rho if sides[side].xi else grid.eta_rho
        if lvar in bry.variables:
            print(lvar)
            zeta = np.ma.array(bry.variables[lvar][:])
            mask = np.ma.getmaskarray(zeta)
            # Detide
            for i in seapy.progressbar.progress(range(size)):
                if np.any(mask[:, i]):
                    continue
                out = seapy.tide.fit(time, zeta[:, i], tides=tides, lat=lat[i],
                                     tide_start=tide_start)
                zeta[:, i] -= out['fit'].data

                # Save the amp/phase in the tide file
                for n, t in enumerate(tides):
                    if sides[side].xi:
                        eamp[n, idx[0], i] = out['major'][t].amp
                        epha[n, idx[0], i] = np.mod(
                            out['major'][t].phase, 2 * np.pi)
                    else:
                        eamp[n, i, idx[1]] = out['major'][t].amp
                        epha[n, i, idx[1]] = np.mod(
                            out['major'][t].phase, 2 * np.pi)

            # Save out the detided information
            bry.variables[lvar][:] = zeta
            zeta = [0]
            bry.sync()

        # Detide the barotropic velocity
        uvar = "ubar_" + side
        vvar = "vbar_" + side
        if uvar in bry.variables and vvar in bry.variables:
            print(uvar, vvar)
            ubar = np.zeros((len(time), size))
            vbar = np.zeros((len(time), size))

            # Load data, put onto rho-grid, and rotate
            bubar = np.ma.array(bry.variables[uvar][:]).filled(0)
            bvbar = np.ma.array(bry.variables[vvar][:]).filled(0)
            if sides[side].xi:
                ubar[:, 1:-1] = 0.5 * (bubar[:, 1:] + bubar[:, :-1])
                ubar[:, 0] = bubar[:, 1]
                ubar[:, -1] = bubar[:, -2]
                vbar = bvbar.copy()
            else:
                vbar[:, 1:-1] = 0.5 * (bvbar[:, 1:] + bvbar[:, :-1])
                vbar[:, 0] = bvbar[:, 1]
                vbar[:, -1] = bvbar[:, -2]
                ubar = bubar.copy()
            ubar, vbar = seapy.rotate(ubar, vbar, grid.angle[idx[0], idx[1]])
            bubar = bvbar = []

            # Detide
            for i in seapy.progressbar.progress(range(size)):
                if np.any(mask[:, i]):
                    continue
                out = seapy.tide.fit(
                    time, ubar[:, i] + 1j * vbar[:, i], tides=tides, lat=lat[i],
                    tide_start=tide_start)
                ubar[:, i] -= np.real(out['fit'])
                vbar[:, i] -= np.imag(out['fit'])

                # Save the amp/phase in the tide file
                for n, t in enumerate(tides):
                    if sides[side].xi:
                        cmax[n, idx[0], i] = out['major'][t].amp
                        cmin[n, idx[0], i] = out['minor'][t].amp
                        cpha[n, idx[0], i] = out['major'][t].phase
                        cang[n, idx[0], i] = out['minor'][t].phase
                    else:
                        cmax[n, i, idx[1]] = out['major'][t].amp
                        cmin[n, i, idx[1]] = out['minor'][t].amp
                        cpha[n, i, idx[1]] = out['major'][t].phase
                        cang[n, i, idx[1]] = out['minor'][t].phase

            ubar, vbar = seapy.rotate(ubar, vbar, -grid.angle[idx[0], idx[1]])
            if sides[side].xi:
                bry.variables[uvar][:] = 0.5 * (ubar[:, 1:] + ubar[:, :-1])
                bry.variables[vvar][:] = vbar
            else:
                bry.variables[vvar][:] = 0.5 * (vbar[:, 1:] + vbar[:, :-1])
                bry.variables[uvar][:] = ubar
            bry.sync()
            ubar = vbar = []

    # Have to duplicate the boundary tide info into the inner row/column
    eamp[:, 1:-1, 1] = eamp[:, 1:-1, 0]
    eamp[:, 1:-1, -2] = eamp[:, 1:-1, -1]
    eamp[:, 1, 1:-1] = eamp[:, 0, 1:-1]
    eamp[:, -2, 1:-1] = eamp[:, -1, 1:-1]
    epha[:, 1:-1, 1] = epha[:, 1:-1, 0]
    epha[:, 1:-1, -2] = epha[:, 1:-1, -1]
    epha[:, 1, 1:-1] = epha[:, 0, 1:-1]
    epha[:, -2, 1:-1] = epha[:, -1, 1:-1]
    cmax[:, 1:-1, 1] = cmax[:, 1:-1, 0]
    cmax[:, 1:-1, -2] = cmax[:, 1:-1, -1]
    cmax[:, 1, 1:-1] = cmax[:, 0, 1:-1]
    cmax[:, -2, 1:-1] = cmax[:, -1, 1:-1]
    cmin[:, 1:-1, 1] = cmin[:, 1:-1, 0]
    cmin[:, 1:-1, -2] = cmin[:, 1:-1, -1]
    cmin[:, 1, 1:-1] = cmin[:, 0, 1:-1]
    cmin[:, -2, 1:-1] = cmin[:, -1, 1:-1]
    cpha[:, 1:-1, 1] = cpha[:, 1:-1, 0]
    cpha[:, 1:-1, -2] = cpha[:, 1:-1, -1]
    cpha[:, 1, 1:-1] = cpha[:, 0, 1:-1]
    cpha[:, -2, 1:-1] = cpha[:, -1, 1:-1]
    cang[:, 1:-1, 1] = cang[:, 1:-1, 0]
    cang[:, 1:-1, -2] = cang[:, 1:-1, -1]
    cang[:, 1, 1:-1] = cang[:, 0, 1:-1]
    cang[:, -2, 1:-1] = cang[:, -1, 1:-1]

    # Set the tide reference
    tideout = {}
    tideout['tides'] = tides
    tideout['tide_start'] = tide_start
    tideout['Eamp'] = eamp
    tideout['Ephase'] = epha
    tideout['Cmajor'] = cmax
    tideout['Cminor'] = cmin
    tideout['Cphase'] = cpha
    tideout['Cangle'] = cang

    seapy.roms.tide.create_forcing(tidefile, tideout,
                                   title="Tides from " + bryfile, epoch=epoch)
    bry.close()

    pass
