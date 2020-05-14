#!/usr/bin/env python
"""
   cobalt.py

   Define fields and utility functions for working with GFDL COBALT
   and ROMS

   Author: Brian Powell <powellb@hawaii.edu>
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""
import seapy
import os

# Define the COBALT fields that are used
fields = {"alk": {"grid": "rho", "dims": 3},
          "cadet_arag": {"grid": "rho", "dims": 3},
          "cadet_calc": {"grid": "rho", "dims": 3},
          "dic": {"grid": "rho", "dims": 3},
          "fed": {"grid": "rho", "dims": 3},
          "fedet": {"grid": "rho", "dims": 3},
          "fedi": {"grid": "rho", "dims": 3},
          "felg": {"grid": "rho", "dims": 3},
          "fesm": {"grid": "rho", "dims": 3},
          "ldon": {"grid": "rho", "dims": 3},
          "ldop": {"grid": "rho", "dims": 3},
          "lith": {"grid": "rho", "dims": 3},
          "lithdet": {"grid": "rho", "dims": 3},
          "nbact": {"grid": "rho", "dims": 3},
          "ndet": {"grid": "rho", "dims": 3},
          "ndi": {"grid": "rho", "dims": 3},
          "nlg": {"grid": "rho", "dims": 3},
          "nsm": {"grid": "rho", "dims": 3},
          "nh4": {"grid": "rho", "dims": 3},
          "no3": {"grid": "rho", "dims": 3},
          "o2": {"grid": "rho", "dims": 3},
          "pdet": {"grid": "rho", "dims": 3},
          "po4": {"grid": "rho", "dims": 3},
          "srdon": {"grid": "rho", "dims": 3},
          "srdop": {"grid": "rho", "dims": 3},
          "sldon": {"grid": "rho", "dims": 3},
          "sldop": {"grid": "rho", "dims": 3},
          "sidet": {"grid": "rho", "dims": 3},
          "silg": {"grid": "rho", "dims": 3},
          "sio4": {"grid": "rho", "dims": 3},
          "nsmz": {"grid": "rho", "dims": 3},
          "nmdz": {"grid": "rho", "dims": 3},
          "nlgz": {"grid": "rho", "dims": 3}}

# Extra aggregate fields that are required to initialize, but are not
# provided by the COBALT output
ini_fields = {"chl": {"grid": "rho", "dims": 3},
              "irr_mem": {"grid": "rho", "dims": 3},
              "htotal": {"grid": "rho", "dims": 3},
              "co3_ion": {"grid": "rho", "dims": 3},
              "mu_mem_di": {"grid": "rho", "dims": 3},
              "mu_mem_sm": {"grid": "rho", "dims": 3},
              "mu_mem_lg": {"grid": "rho", "dims": 3}}

# Diagnostic fields that are saved in diagnostic files that are useful
dia_fields = {"chl": {"grid": "rho", "dims": 3},
              "irr_mem": {"grid": "rho", "dims": 3},
              "htotal": {"grid": "rho", "dims": 3},
              "co3_ion": {"grid": "rho", "dims": 3},
              "fe_bulk_flx": {"grid": "rho", "dims": 3},
              "omega_cadet_calc": {"grid": "rho", "dims": 3},
              "omega_cadet_arag": {"grid": "rho", "dims": 3}}

# Extra fields that are required in the atmospheric forcing
frc_fields = {"atmCO2": {"grid": "rho", "dims": 2},
              "ironsed": {"grid": "rho", "dims": 2},
              "fecoast": {"grid": "rho", "dims": 2},
              "solublefe": {"grid": "rho", "dims": 2},
              "mineralfe": {"grid": "rho", "dims": 2}}

# Define the vmap for nesting. This simply forms a one-to-one correspondence
# between the fields.
vmap = {k: k for k in fields}

# Create a dictionary of CDL files
_cdl_dir = os.path.dirname(__file__)
cdl = {"ini": _cdl_dir + "/ini.cdl",
       "his": _cdl_dir + "/his.cdl",
       "nudge": _cdl_dir + "/nudge.cdl",
       "clim": _cdl_dir + "/clim.cdl",
       "bry": _cdl_dir + "/bry.cdl",
       "frc": _cdl_dir + "/frc_bulk.cdl",
       "zlevel": _cdl_dir + "/zlevel.cdl"}


# Keep track of original ROMS fields
roms_fields = dict(seapy.roms.fields)

# Helper functions to enable/disable COBALT fields


def add_psource_var(nc, name):
    """
    Add a list of COBALT variables to an existing netCDF file that
    has been opened for writing.

    Input
    -----
      nc: netCDF4,
          The netCDF4 river file object to write to.
          NOTE: it must be opened as writeable.
      name: string or list of strings,
          The COBALT variable names to add as point sources to the
          river file.

    Output
    ------
       None
    """
    from numpy import asarray
    name = asarray(name)
    for n in name:
        var = {"name": f"river_{n}",
               "type": "float",
               "dims": "river_time, s_rho, river",
               "attr": {"units": "mol/kg",
                        "long_name": f"River runoff for {n}",
                        "fields": f"river {n}, scalar, series"}}
        seapy.roms.ncgen.add_variable(nc, var)


def enable():
    """
    Switch seapy to use all fields from ROMS hydrodynamics and COBALT
    """
    seapy.roms.fields.update(fields)


def disable():
    """
    Switch seapy to use only fields from ROMS hydrodynamics
    """
    seapy.roms.fields = dict(roms_fields)


def enable_ini():
    """
    Switch seapy to use all fields from ROMS hydrodynamics and COBALT ini fields
    """
    enable()
    seapy.roms.fields.update(ini_fields)


def disable_ini():
    """
    Switch seapy to use ROMS and COBALT fields without the ini fields
    """
    disable()
    enable()


enable()
