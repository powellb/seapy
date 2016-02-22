"""
  __init__.py

  State Estimation and Analysis for PYthon

    Module for working with oceanographic data and models

  Copyright (c)2016 University of Hawaii under the BSD-License.

  Requires the following packages: joblib

  Import classes include:

  - :class:`~seapy.environ.opt`
  - :class:`~seapy.progressbar.ProgressBar`

  Imported functions include:

  - :func:`~seapy.lib.adddim`
  - :func:`~seapy.lib.convolve_mask`
  - :func:`~seapy.lib.earth_distance`
  - :func:`~seapy.lib.rotate`
  - :func:`~seapy.lib.vecfind`
  - :func:`~seapy.lib.list_files`
  - :func:`~seapy.lib.day2date`
  - :func:`~seapy.lib.date2day`
  - :func:`~seapy.lib.today2day`
  - :func:`~seapy.oa.oasurf`
  - :func:`~seapy.oa.oavol`
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
from .environ import opt
from .hawaii import hawaii
from .oa import *
