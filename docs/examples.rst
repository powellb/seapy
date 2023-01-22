Examples
========

Many of the time-saving features are in generating fields for running the ROMS model.

1. To load the meta information about a model (ROMS, HYCOM, MITgcm, POM, SODA), load an output file (history, average, climatology, grid, etc.) via:

::

    >> mygrid = seapy.model.asgrid(filename)

    >> mygrid
    C-Grid: 32x194x294

    >> print(mygrid)
    filename
    32x194x294: C-Grid with S-level
    Available: I,J,_isroms,_nc,angle,cgrid,cs_r,depth_rho,depth_u,depth_v,dm,dn,eta_rho,eta_u,eta_v,f,filename,h,hc,lat_rho,lat_u,lat_v,lm,ln,lon_rho,lon_u,lon_v,mask_rho,mask_u,mask_v,n,name,pm,pn,s_rho,shape,spatial_dims,tcline,theta_b,theta_s,thick_rho,thick_u,thick_v,vstretching,vtransform,xi_rho,xi_u,xi_v

2. Most methods available in SEAPY require a grid, which can be specified as a "filename" or as a grid object.

3. Find out how to download global HYCOM data that will span my grid from 1/1/2015 through 5/1/2015:

::


    >> seapy.model.hycom.load_history("hycom_file.nc", start_time=datetime(2015,1,1),
                                     end_time=datetime(2015,5,1),
                                     grid=mygrid, load_data=False)
    ncks -v water_temp,salinity,surf_el,water_v,water_u -d time,352,352 -d lat,1204,1309 -d lon,2438,2603 http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1 hycom_file.nc

This will display the 'ncks' command necessary to download the data. If you want to have SEAPY download it (not recommended due to server-speed), use 'load_data=True'.

4. Once you have HYCOM data, interpolate it to your grid

::

    >> seapy.roms.interp.to_clim("hycom_file.nc", "my_clim.nc",
                          dest_grid=mygrid, nx=1/6, ny=1/6,
                          vmap={"surf_el":"zeta", "water_temp":"temp",
                          "water_u":"u", "water_v":"v", "salinity":"salt"})

5. Generate boundary conditions from the climatology

::

    >> seapy.roms.boundary.from_roms("my_clim.nc", "my_bry.nc")

6. Generate initial conditions from the climatology

::

    >> seapy.roms.initial.from_roms("my_clim.nc", "my_ini.nc")

7. You now have what you need to run your model

8. To set up data assimilation, download the raw observations (e.g., aviso_map_day1.nc, aviso_map_day2.nc, argo_day1.nc ). You can then process the data:

::

    >> dt = 400/86400       # time-step of the model in days
    >> aviso_gen = seapy.roms.obsgen.aviso_sla_map(mygrid, dt)
    >> aviso_gen.batch_files(seapy.list_files('.','aviso.*nc'), 'aviso_roms_#.nc')
    >> argo_gen = seapy.roms.obsgen.argo_ctd(mygrid, dt)
    >> obs = argo_gen.convert_file("argo_day1.nc")
    >> obs.to_netcdf("argo_roms_1.nc")

9. Put all of the processed observations files together into a file for a given assimilation window

::

    >> seapy.roms.obs.merge_files(seapy.list_files('.*roms_[0-9]+.nc'), 'roms_obs_#.nc', np.arange([0, 10.1, 5]))

There are many more things that can be done, but these show some of the power available via simple commands.



