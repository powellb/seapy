project( 'seapy', 'c',
  # Note that the git commit hash cannot be added dynamically here
  version: '0.1.0',
  license: 'MIT',
  meson_version: '>= 0.63.0',
  default_options: [
    'buildtype=debugoptimized',
    # TODO: the below -Wno flags are all needed to silence warnings in
    # f2py-generated code. This should be fixed in f2py itself.
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types',
    'fortran_args=-Wno-conversion',  # silence "conversion from REAL(8) to INTEGER(4)"
    'fortran_std=legacy',
  ],
)

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)


fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()
message(py3.full_path())
message(py3.get_install_dir())

# Determine whether this is being compiled with numpy 1.X or 2.X
numpy_ver = run_command(py3,
  ['-c', 'import numpy; print(numpy.__version__.split(".")[0])'],
  check : true
).stdout().strip()
is_numpytwo = numpy_ver.to_int() >= 2

message('The Numpy major version is: ', numpy_ver)


# Don't use the deprecated NumPy C API. Define this to a fixed version instead
# of NPY_API_VERSION in order not to break compilation for released versions
# when NumPy introduces a new deprecation. Use in a meson.build file::
#
#   py3.extension_module('_name',
#     'source_fname',
#     numpy_nodepr_api)
#
if is_numpytwo
  numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_2_0_API_VERSION'
else
  numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'
endif

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
  ['-c', 'import os; os.chdir(".."); from numpy import f2py; print(f2py.get_include())'],
  check : true
).stdout().strip()

oalib_source = custom_target('oalibmodule.c',
  input : ['src/oalib.f', 'src/hindices.f'],  # .f so no F90 wrappers
  output : ['oalibmodule.c', 'oalib-f2pywrappers.f'],
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'oalib', '--lower']
)

inc_dirs = include_directories(incdir_numpy, incdir_f2py)

# Platform detection to set more flags for Windows systems
is_windows = host_machine.system() == 'windows'
is_mac = host_machine.system() == 'darwin'
is_mingw = is_windows and cc.get_id() == 'gcc'

# Declare the sources
py3.install_sources([
  'seapy/__init__.py',
  'seapy/cdl_parser.py',
  'seapy/debug.py',
  'seapy/environ.py',
  'seapy/filt.py',
  'seapy/lib.py',
  'seapy/mapping.py',
  'seapy/oa.py',
  'seapy/plot.py',
  'seapy/qserver.py',
  'seapy/river.py',
  'seapy/tidal_energy.py',
  'seapy/tide.py',
  'seapy/timeout.py',
  'seapy/constituents.npz',
  'seapy/model/__init__.py',
  'seapy/model/grid.py',
  'seapy/model/hycom.py',
  'seapy/model/lib.py',
  'seapy/model/soda.py',
  'seapy/roms/__init__.py',
  'seapy/roms/analysis.py',
  'seapy/roms/boundary.py',
  'seapy/roms/clim.py',
  'seapy/roms/ezgrid.py',
  'seapy/roms/forcing.py',
  'seapy/roms/initial.py',
  'seapy/roms/interp.py',
  'seapy/roms/lib.py',
  'seapy/roms/ncgen.py',
  'seapy/roms/obs.py',
  'seapy/roms/obsgen.py',
  'seapy/roms/psource.py',
  'seapy/roms/tide.py',
  'seapy/roms/cdl/adsen.cdl',
  'seapy/roms/cdl/bry_unlimit.cdl',
  'seapy/roms/cdl/clm_ts.cdl',
  'seapy/roms/cdl/frc_bulk.cdl',
  'seapy/roms/cdl/frc_direct.cdl',
  'seapy/roms/cdl/frc_fluxclm.cdl',
  'seapy/roms/cdl/frc_qcorr.cdl',
  'seapy/roms/cdl/frc_rivers.cdl',
  'seapy/roms/cdl/frc_srelax.cdl',
  'seapy/roms/cdl/frc_tides.cdl',
  'seapy/roms/cdl/frc_wave.cdl',
  'seapy/roms/cdl/frc_windstress.cdl',
  'seapy/roms/cdl/ini_hydro.cdl',
  'seapy/roms/cdl/nudge_coef.cdl',
  'seapy/roms/cdl/roms_grid.cdl',
  'seapy/roms/cdl/s4dvar_obs.cdl',
  'seapy/roms/cdl/s4dvar_obs_rays.cdl',
  'seapy/roms/cdl/s4dvar_std_b.cdl',
  'seapy/roms/cdl/s4dvar_std_f.cdl',
  'seapy/roms/cdl/s4dvar_std_i.cdl',
  'seapy/roms/cdl/s4dvar_std_m.cdl',
  'seapy/roms/cdl/zlevel_1d.cdl',
  'seapy/roms/cdl/zlevel_1d_grid.cdl',
  'seapy/roms/cdl/zlevel_2d.cdl',
  'seapy/roms/cdl/zlevel_2d_grid.cdl',
  'seapy/roms/cdl/zlevel_microbe_1d.cdl',
  'seapy/roms/cdl/zlevel_microbe_2d.cdl',
  'seapy/roms/cobalt/__init__.py',
  'seapy/roms/cobalt/bry.cdl',
  'seapy/roms/cobalt/clim.cdl',
  'seapy/roms/cobalt/cobalt.py',
  'seapy/roms/cobalt/frc_bulk.cdl',
  'seapy/roms/cobalt/frc_cobalt.cdl',
  'seapy/roms/cobalt/his-all.cdl',
  'seapy/roms/cobalt/his.cdl',
  'seapy/roms/cobalt/ini.cdl',
  'seapy/roms/cobalt/nudge.cdl',
  'seapy/roms/cobalt/zevel_dia.cdl',
  'seapy/roms/cobalt/zlevel-cobalt-2d.cdl',
  'seapy/roms/cobalt/zlevel.cdl',
  'seapy/seawater/__init__.py', 
  'seapy/seawater/density.py', 
  'seapy/seawater/heat.py', 
  'seapy/seawater/misc.py', 
  'seapy/seawater/salinity.py',
  'seapy/external/__init__.py',
  'seapy/hawaii_coast/hawaii.dbf',
  'seapy/hawaii_coast/hawaii.shp',
  'seapy/hawaii_coast/hawaii.shx',
  ],
  pure: false,
  subdir: 'seapy')

py3.extension_module('oalib',
  ['src/oalib.f', 'src/hindices.f', oalib_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_dirs,
  dependencies: py3_dep,
  install: true,
  install_dir: join_paths(py3.get_install_dir(), 'seapy/external')
)
