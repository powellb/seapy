#!/usr/bin/env python
"""
  boundary.py
  
  ROMS boundary utilities

  Written by Brian Powell on 01/15/14
  Copyright (c)2014 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import seapy
import os
import numpy as np
import netCDF4
import netcdftime
import textwrap
from scipy import interpolate

import pudb

def from_roms(roms_file, bry_file, grid=None, records=None):
    """
    Given a ROMS history, average, or climatology file, generate 
    boundary conditions on the same grid.
    
    Parameters
    ----------
    roms_file : string,
        ROMS source (history, average, climatology file)
    bry_file : string,
        output boundary file
    grid : seapy.model.grid or string, optional,
        ROMS grid for boundaries
    [records] : array, optional 
        record indices to put into the boundary
    
    Returns
    -------
    None
    
    """
    grid = seapy.model.asgrid(grid)
    ncroms = netCDF4.Dataset(roms_file)
    time = seapy.roms.get_timevar(ncroms)
    records = np.arange(0, len(ncroms.variables[time][:])) \
             if records is None else records
    try:
        src_time=netcdftime.utime(ncroms.variables[time].units)
    except AttributeError:
        src_time=netcdftime.utime(seapy.roms.default_epoch)

    # Create the boundary file and fill up the descriptive data
    if os.path.isfile(bry_file):
        ncbry=netCDF4.Dataset(bry_file,"a")
    else:
        ncbry=seapy.roms.ncgen.create_bry(bry_file, 
             eta_rho=grid.eta_rho,xi_rho=grid.xi_rho,s_rho=grid.n,
             timebase=src_time.origin,title="generated from "+roms_file)
    brytime = seapy.roms.get_timevar(ncbry)
    bry_time=netcdftime.utime(ncbry.variables[brytime].units)
    grid.to_netcdf(ncbry)
    ncbry.variables["bry_time"][:]=bry_time.date2num(
        src_time.num2date(ncroms.variables[time][records]))

    # Go over the variables for each side
    sides={"north":[-2,"-1"], "south":[-2,"0"], "east":[-1,"-1"], "west":[-1,"0"]}
    index=["records",":",":",":"]
    for var in seapy.roms.fields.keys():
        if var in ncroms.variables:
            for side in sides: 
                outvar=var+"_"+side
                ndim=seapy.roms.fields[var]["dims"]
                if ndim==3:
                    myindex=list(index)
                elif ndim==2:
                    myindex=list(index[:-1])
                myindex[sides[side][0]]=sides[side][1]
                indexstr=",".join(myindex)
                ncbry.variables[outvar][:]=eval("ncroms.variables[var]["+indexstr+"]")
    ncbry.close()
    
    pass

def gen_stations(filename, grid=None):
    """
    Generate a station file with stations at every boundary location for use in 
    nesting one grid within another.
    
    Parameters
    ----------
    filename: string
        Input name of station file to create
    grid: string or seapy.model.grid
        Input grid to generate station file from. If string, it will open the grid file.
        If grid, it will use the grid information

    Returns
    -------
    None
    
    """
    grid = seapy.model.asgrid(grid)
    
    # Put together the boundaries
    lon=np.concatenate([ grid.lon_rho[0,:], grid.lon_rho[-1,:], 
                         grid.lon_rho[:,0], grid.lon_rho[:,-1]])
    lat=np.concatenate([ grid.lat_rho[0,:], grid.lat_rho[-1,:], 
                         grid.lat_rho[:,0], grid.lat_rho[:,-1]])
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
        print(" BOUNDARY STATIONS FOR GRID: {}".format(grid.filename), 
              file=text_file)
        print(textwrap.dedent(header), file=text_file)
        print("        NSTATION ==  {}".format(Npts), file=text_file)
        print(textwrap.dedent(stations), file=text_file)
        for i in range(Npts):
            print("        1     1    {0:10.6f}   {1:10.6f}   BRY".format( \
                lon[i],lat[i]), file=text_file)
            
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
    time = seapy.roms.get_timevar(ncstation)
    try:
        src_time=netcdftime.utime(ncstation.variables[time].units)
    except AttributeError:
        src_time=netcdftime.utime(seapy.roms.default_epoch)

    # Create the boundary file and fill up the descriptive data
    if os.path.isfile(bry_file):
        ncbry=netCDF4.Dataset(bry_file,"a")
    else:
        ncbry=seapy.roms.ncgen.create_bry(bry_file, 
             eta_rho=grid.eta_rho,xi_rho=grid.xi_rho,s_rho=grid.n,
             timebase=src_time.origin,title="generated from "+station_file)
    brytime = seapy.roms.get_timevar(ncbry)
    bry_time=netcdftime.utime(ncbry.variables[brytime].units)
    grid.to_netcdf(ncbry)
    ncbry.variables["bry_time"][:]=bry_time.date2num(
        src_time.num2date(ncstation.variables[time][:]))

    # Set up the indices
    bry = {
        "south":range(0,grid.lm),
        "north":range(grid.lm,2*grid.lm),
        "west":range(2*grid.lm,2*grid.lm+grid.ln),
        "east":range(2*grid.lm+grid.ln,2*(grid.lm+grid.ln))
    }

    # Get the information to construct the depths of the station data
    sta_vt=ncstation.variables["Vtransform"][:]
    sta_hc=ncstation.variables["hc"][:]
    sta_s_rho=ncstation.variables["s_rho"][:]
    sta_cs_r=ncstation.variables["Cs_r"][:]
    sta_h=ncstation.variables["h"][:]
    sta_angle=ncstation.variables["angle"][:]
    sta_lon=ncstation.variables["lon_rho"][:]
    sta_lat=ncstation.variables["lat_rho"][:]
    sta_mask=np.ones(sta_lat.shape)
    sta_mask[np.where(sta_lon*sta_lat>1e10)[0]]=0

    # Load the station data as we need to manipulate it
    sta_zeta=np.ma.masked_greater(ncstation.variables["zeta"][:],100)
    sta_ubar=np.ma.masked_greater(ncstation.variables["ubar"][:],100)
    sta_vbar=np.ma.masked_greater(ncstation.variables["vbar"][:],100)
    sta_temp=np.ma.masked_greater(ncstation.variables["temp"][:],100)
    sta_salt=np.ma.masked_greater(ncstation.variables["salt"][:],100)
    sta_u=np.ma.masked_greater(ncstation.variables["u"][:],100)
    sta_v=np.ma.masked_greater(ncstation.variables["v"][:],100)

    
    # Create the true positions and mask
    grid_h=np.concatenate([ grid.h[0,:], grid.h[-1,:], 
                         grid.h[:,0], grid.h[:,-1]])
    grid_lon=np.concatenate([ grid.lon_rho[0,:], grid.lon_rho[-1,:], 
                         grid.lon_rho[:,0], grid.lon_rho[:,-1]])
    grid_lat=np.concatenate([ grid.lat_rho[0,:], grid.lat_rho[-1,:], 
                         grid.lat_rho[:,0], grid.lat_rho[:,-1]])
    grid_mask=np.concatenate([ grid.mask_rho[0,:], grid.mask_rho[-1,:], 
                              grid.mask_rho[:,0], grid.mask_rho[:,-1]])
    grid_angle=np.concatenate([ grid.angle[0,:], grid.angle[-1,:], 
                              grid.angle[:,0], grid.angle[:,-1]])

    # Search for bad stations due to child grid overlaying parent mask.
    # Unfortunately, ROMS will give points that are not at the locations
    # you specify if those points conflict with the mask. So, these points
    # are simply replaced with the nearest.
    dist = np.sqrt( (sta_lon-grid_lon)**2 + (sta_lat-grid_lat)**2 )
    bad_pts = np.where( np.logical_and(dist > 0.001, grid_mask == 1) )[0]
    good_pts = np.where( np.logical_and(dist < 0.001, grid_mask == 1) )[0]
    for i in bad_pts:
        dist = np.sqrt( (sta_lon[i]-sta_lon[good_pts])**2 + 
                        (sta_lat[i]-sta_lat[good_pts])**2 )
        index = good_pts[ np.where(dist==np.min(dist))[0] ]
        sta_h[i] = sta_h[index]
        sta_angle[i] = sta_angle[index]
        sta_lon[i] = lon[index]
        sta_lat[i] = lat[index]
        sta_zeta[:,i] = sta_zeta[:,index]
        sta_ubar[:,i] = sta_ubar[:,index]
        sta_vbar[:,i] = sta_vbar[:,index]
        sta_temp[:,i,:] = sta_temp[:,index,:]
        sta_salt[:,i,:] = sta_salt[:,index,:]
        sta_u[:,i,:] = sta_u[:,index,:]
        sta_v[:,i,:] = sta_v[:,index,:]


    # Construct the boundaries: a dictionary of boundary side and two element
    # array whether the u[0] or v[1] dimensions need to be averaged
    sides={"north":[True,False], "south":[True,False], 
           "east":[False,True], "west":[False,True]}
    delta_angle = sta_angle-grid_angle
    sta_ubar, sta_vbar = seapy.rotate(sta_ubar,sta_vbar,delta_angle)
    sta_u, sta_v = seapy.rotate(sta_u,sta_v,delta_angle.T)
    for side in sides.keys():
        print(side)
        # 1) Zeta
        ncbry.variables["zeta_"+side][:] = sta_zeta[:,bry[side]]

        # 2) Ubar
        if sides[side][0]:
            ncbry.variables["ubar_"+side][:] = 0.5 * ( \
                sta_ubar[:,bry[side][0:-1]]+sta_ubar[:,bry[side][1:]])
        else:
            ncbry.variables["ubar_"+side][:] = sta_ubar[:,bry[side]]
            
        # 3) Vbar
        if sides[side][1]:
            ncbry.variables["vbar_"+side][:] = 0.5 * ( \
                sta_vbar[:,bry[side][0:-1]]+sta_vbar[:,bry[side][1:]])
        else:
            ncbry.variables["vbar_"+side][:] = sta_vbar[:,bry[side]]
        
        # For 3D variables, we need to loop through time and interpolate
        # onto the child grid. Construct the distances
        x = bry[side]*0
        x[1:] = np.cumsum(seapy.earth_distance(grid_lon[bry[side][0:-1]],
                                     grid_lat[bry[side][0:-1]],
                                     grid_lon[bry[side][1:]],
                                     grid_lat[bry[side][1:]]))
        sta_x = seapy.adddim(x,len(sta_s_rho))
        x = seapy.adddim(x,len(grid.s_rho))
        for n,t in seapy.progress(enumerate(ncstation.variables[time][:])):
            sta_depth = seapy.roms.depth(sta_vt, sta_h[bry[side]], sta_hc, 
                            sta_s_rho, sta_cs_r, sta_zeta[n,bry[side]])
            depth = seapy.roms.depth(grid.vtransform, grid_h[bry[side]], 
                        grid.hc, grid.s_rho, grid.cs_r, sta_zeta[n,bry[side]])
            sta_ocean = np.where(sta_mask[bry[side]]==1)[0]
            ocean = np.where(grid_mask[bry[side]]==1)[0]
            # 4) Temp
            ncbry.variables["temp_"+side][n,:]=0.0
            ncbry.variables["temp_"+side][n,:,ocean],pmap = seapy.oa.oasurf( \
                sta_x[:,ocean], sta_depth[:,ocean], 
                np.transpose(sta_temp[n,bry[side],:][ocean,:]), \
                x[:,ocean], depth[:,ocean],nx=0, ny=2, weight=2)
            # pu.db
        
            # 5) Salt
            ncbry.variables["salt_"+side][n,:]=0.0
            ncbry.variables["salt_"+side][n,:,ocean],pmap = seapy.oa.oasurf( \
                sta_x[:,ocean], sta_depth[:,ocean], 
                np.transpose(sta_salt[n,bry[side],:][ocean,:]), \
                x[:,ocean], depth[:,ocean], pmap=pmap, nx=0, ny=2, weight=2)

            # 6) U
            data = np.zeros(x.shape)
            data[:,ocean],pmap = seapy.oa.oasurf( \
                sta_x[:,ocean], sta_depth[:,ocean], 
                np.transpose(sta_u[n,bry[side],:][ocean,:]), \
                x[:,ocean], depth[:,ocean], pmap=pmap, nx=0, ny=2, weight=2)
            if sides[side][0]:
                ncbry.variables["u_"+side][n,:] = 0.5 * ( \
                    data[:,0:-1]+data[:,1:])
            else:
                ncbry.variables["u_"+side][n,:] = data

            # 7) V
            data = data * 0
            data[:,ocean],pmap = seapy.oa.oasurf( \
                sta_x[:,ocean], sta_depth[:,ocean], 
                np.transpose(sta_v[n,bry[side],:][ocean,:]), \
                x[:,ocean], depth[:,ocean], pmap=pmap, nx=0, ny=2, weight=2)
            if sides[side][1]:
                ncbry.variables["v_"+side][n,:] = 0.5 * ( \
                    data[:,0:-1]+data[:,1:])
            else:
                ncbry.variables["v_"+side][n,:] = data
            
