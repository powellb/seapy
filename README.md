# State Estimation and Analysis in PYthon (SEAPY)

[![Test](https://github.com/rsoutelino/seapy/actions/workflows/test.yml/badge.svg)](https://github.com/rsoutelino/seapy/actions/workflows/test.yml)
[![Build](https://github.com/rsoutelino/seapy/actions/workflows/build.yml/badge.svg)](https://github.com/rsoutelino/seapy/actions/workflows/build.yml)



Tools for working with ocean models and data.

SEAPY requires: basemap, h5py, joblib, netcdf4, numpy, numpy_groupies, rich and scipy.


## Installation

### Install with pip

Pip installation is based on [meson-python](https://mesonbuild.com/meson-python/), following example from [pyfort](https://github.com/larsbuntemeyer/pyfort).

- create a python virtual environment using your preferred tool, activate the environment and run `pip install .`

Simple example using virtualenvwrapper assuming you are in the root directory of cloned
version of this repository:

```bash
mkvirtualenv seapy
workon seapy # if not activated automatically
pip install .
```

### Install from Conda-Forge

Install from [conda-forge](https://conda-forge.org/) with the Conda package manager:
```
$ conda install -c conda-forge seapy
```

You should also consider making conda-forge your default channel. See the [conda-forge tips and tricks page](https://conda-forge.org/docs/user/tipsandtricks.html).

The Conda-Forge [SEAPY feedstock](https://github.com/conda-forge/seapy-feedstock) is maintained by Filipe Fernandes, [ocefpaf](https://github.com/ocefpaf/). As of February 2021 there are binary packages on all the platforms that Conda-Forge supports: Python 3.6 through 3.9 on Linux, Windows and Mac OSX (all 64-bit).

To remove seapy:
```
$ conda remove seapy
```

## Install from PyPI

Install from [PyPI](https://pypi.org/) with PIP:
```
$ pip install seapy-ocean
```

Note that on PyPI (but nowhere else) the package name is seapy-ocean to avoid a name clash with another package. The module name is still seapy.

SEAPY packages on PyPI have been built and uploaded by Mark Hadfield, [hadfieldnz](https://pypi.org/user/hadfieldnz/). There is a source distribution that should build with no problems on Linux (and Mac OSX, but we haven't tested it). In the pst there have been binary distributions for Windows (64-bit), but these have now been deleted as binary builds with PIP are no longer supported.

In a Conda environment, it is quite possible to install with PIP, but dependency handling and updating will be cleaner if you use the Conda package.

To remove seapy-ocean
```
$ pip uninstall seapy-ocean
```

## Install from source code on GitHub.com

The SEAPY source code is maintained by Brian Powell, (powellb)[https://github.com/powellb]. Releases are made on the [master branch](https://github.com/powellb/seapy/tree/master)

Install from [GitHub.com](https://github.com/) with PIP:
```
$ pip install git+git://github.com/powellb/seapy@master
```

OR clone a copy of the source and install in editable mode, eg:
```
$ git clone https://github.com/powellb/seapy.git
$ pip install -e seapy
```

With an [editable-mode](https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs) installation, changes you make to your copy of the source code will take effect when you import the module.

In principle it is possible to build from source on Windows--and success with this has been achieved in the past--but the process tends to break with changes in the environment or Python version, so we don't recommend it or support it.

## Contributing


If you've installed from source in editable mode, then you should definitely consider forking your own copy of the repository. This allows you to keep your changes under revision control on GitHub.com and potentially contribute them to the main project. You should follow the procedures described in this [Git Workflow](https://www.asmeurer.com/git-workflow/) document.

[Forking on GitHub.com](https://help.github.com/en/articles/fork-a-repo) is a lightweight process that won't complicate your workflow and keeps the relationship between your work and the original project clear, so it is strongly advised to do it early. However the immutable and unique nature of Git commits means that you can create and populate a fork later if you want to, as long as you have saved your work somewhere in Git format. To create a fork you will need a [GitHub.com user account](https://help.github.com/en/articles/signing-up-for-a-new-github-account).


All your changes should be committed to a branch other than "master", which is reserved for the master branch in Brian Powell's repository (or copies thereof). A common practice in the [existing SEAPY forks](https://github.com/powellb/seapy/network/members) is to use a branch name matching your user name for your own work. However if you are developing a specific feature or bug fix to be pulled into master, it may be sensible to name the branch after that feature or bug fix.


## Examples

Many of the time-saving features are in generating fields for running the ROMS model.

1. To load the meta information about a model (ROMS, HYCOM, MITgcm, POM, SODA), load an output file (history, average, climatology, grid, etc.) via:

        >> mygrid = seapy.model.asgrid(filename)

        >> mygrid
        C-Grid: 32x194x294

        >> print(mygrid)
        filename
        32x194x294: C-Grid with S-level
        Available: I,J,_isroms,_nc,angle,cgrid,cs_r,depth_rho,depth_u,depth_v,dm,dn,eta_rho,eta_u,eta_v,f,filename,h,hc,lat_rho,lat_u,lat_v,lm,ln,lon_rho,lon_u,lon_v,mask_rho,mask_u,mask_v,n,name,pm,pn,s_rho,shape,spatial_dims,tcline,theta_b,theta_s,thick_rho,thick_u,thick_v,vstretching,vtransform,xi_rho,xi_u,xi_v


2. Most methods available in SEAPY require a grid, which can be specified as a "filename" or as a grid object.

3. Find out how to download global HYCOM data that will span my grid from 1/1/2015 through 5/1/2015:


        >> seapy.model.hycom.load_history("hycom_file.nc", start_time=datetime(2015,1,1),
                                         end_time=datetime(2015,5,1),
                                         grid=mygrid, load_data=False)
        ncks -v water_temp,salinity,surf_el,water_v,water_u -d time,352,352 -d lat,1204,1309 -d lon,2438,2603 http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1 hycom_file.nc

This will display the 'ncks' command necessary to download the data. If you want to have SEAPY download it (not recommended due to server-speed), use `'load_data=True'`.

4. Once you have HYCOM data, interpolate it to your grid

        >> seapy.roms.interp.to_clim("hycom_file.nc", "my_clim.nc",
                          dest_grid=mygrid, nx=1/6, ny=1/6,
                          vmap={"surf_el":"zeta", "water_temp":"temp",
                          "water_u":"u", "water_v":"v", "salinity":"salt"})

5. Generate boundary conditions from the climatology

        >> seapy.roms.boundary.from_roms("my_clim.nc", "my_bry.nc")

6. Generate initial conditions from the climatology

        >> seapy.roms.initial.from_roms("my_clim.nc", "my_ini.nc")

7. You now have what you need to run your model

8. To set up data assimilation, download the raw observations (e.g., `aviso_map_day1.nc`, `aviso_map_day2.nc`, `argo_day1.nc` ). You can then process the data:

        >> dt = 400/86400       # time-step of the model in days
        >> aviso_gen = seapy.roms.obsgen.aviso_sla_map(mygrid, dt)
        >> aviso_gen.batch_files(seapy.list_files('./aviso.*nc'), 'aviso_roms_#.nc')
        >> argo_gen = seapy.roms.obsgen.argo_ctd(mygrid, dt)
        >> obs = argo_gen.convert_file("argo_day1.nc")
        >> obs.to_netcdf("argo_roms_1.nc")

9. Put all of the processed observations files together into a file for a given assimilation window

        >> seapy.roms.obs.merge_files(seapy.list_files('.*roms_[0-9]+.nc'), 'roms_obs_#.nc', np.arange([0, 10.1, 5]))

There are many more things that can be done, but these show some of the power available via simple commands.
