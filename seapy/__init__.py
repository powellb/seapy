"""
  __init__.py

  State Estimation and Analysis for PYthon

    Module for working with oceanographic data and models

  Copyright (c)2020 University of Hawaii under the MIT-License.

  Requires the following packages: joblib

  Import classes include:

  - :class:`~seapy.environ.opt`
  - :class:`~seapy.progressbar.ProgressBar`
  - :class:`~seapy.tidal_energy.energetics`

  Imported functions include:

  - :func:`~seapy.lib.adddim`
  - :func:`~seapy.lib.chunker`
  - :func:`~seapy.lib.convolve_mask`
  - :func:`~seapy.lib.day2date`
  - :func:`~seapy.lib.date2day`
  - :func:`~seapy.lib.earth_angle`
  - :func:`~seapy.lib.earth_distance`
  - :func:`~seapy.lib.flatten`
  - :func:`~seapy.lib.list_files`
  - :func:`~seapy.lib.netcdf`
  - :func:`~seapy.lib.rotate`
  - :func:`~seapy.lib.today2day`
  - :func:`~seapy.lib.unique_rows`
  - :func:`~seapy.lib.vecfind`
  - :func:`~seapy.oa.oasurf`
  - :func:`~seapy.oa.oavol`
  - :func:`~seapy.tidal_energy.tidal_energy`
  - :func:`~seapy.progressbar.progress`

"""

from .lib import *
from . import roms
from . import model
from . import qserver
from . import mapping
from . import filt
from . import plot
from . import progressbar
from . import seawater
from . import tide
from .tidal_energy import tidal_energy
from .environ import opt
from .hawaii import hawaii
from .oa import *
