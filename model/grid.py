#!/usr/bin/env python
"""
  grid
  
  This module handles general model grid information, whether from ROMS or
  other models; however, it is mostly geared towards ROMS

  Written by Brian Powell on 10/09/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import netCDF4
import datetime
import os
import re
import seapy
import numpy as np

class grid:
    def __init__(self, file=None, lat=None, lon=None, depth=None,
                 minimal=True, cgrid=False):
        """
            Class to wrap around a numerical model grid for oceanography. 
            It attempts to track latitude, longitude, depth, and other
            parameters.
        """
        self.file = file
        self.cgrid = cgrid

        if self.file != None:
            self._initfile()
            self._isroms = True if \
              (len(list(set(("s_rho","pm","pn","theta_s","theta_b", \
                            "vtransform", "vstretching")).intersection( \
                            set(self.__dict__)))) > 0) \
              else False
            self.cgrid = True if self._isroms else self.cgrid
        else:
            self._nc = None
            self.lat_rho = lat
            self.lon_rho = lon
            self.depth = depth
            self.cgrid = False
        self._verify_shape()
        if not minimal:
            self.set_dims()
            self.set_depth()
        if self._nc is not None:
            self._nc.close()
            self._nc = None
            
    def _initfile(self):
        """
            Using an input file, try to load as much information
            as can be found in the given file.
        """
        # Define a dictionary to go through and convert netcdf variables
        # to internal class attributes
        gvars = {"lat_rho": ["lat_rho", "lat", "latitude"],
                 "lon_rho": ["lon_rho", "lon", "longitude"],
                 "lat_u": ["lat_u"],
                 "lon_u": ["lon_u"],
                 "lat_v": ["lat_v"],
                 "lon_v": ["lon_v"],
                 "mask_rho": ["mask_rho", "mask"],
                 "mask_u": ["mask_u"],
                 "mask_v": ["mask_v"],
                 "angle": ["angle"],
                 "h": ["h"],
                 "n": ["N"],
                 "theta_s": ["theta_s"],
                 "theta_b": ["theta_b"],
                 "tcline": ["Tcline"],
                 "hc": ["hc"],
                 "vtransform": ["Vtransform"],
                 "vstretching": ["Vstretching"],
                 "s_rho": ["s_rho"],
                 "cs_r": ["Cs_r"],
                 "depth": ["depth", "lev"],
                 "f": ["f"],
                 "pm": ["pm"],
                 "pn": ["pn"],
                }

        # Open the file
        self._nc = netCDF4.Dataset(self.file,"r")
        self.name = re.search("[^\.]*",os.path.basename(self.file)).group();
        for var in gvars.keys():
            for inp in gvars[var]:
                if inp in self._nc.variables:
                    self.__dict__[var] = self._nc.variables[inp][:]

    def _verify_shape(self):
        """
            Verify the dimensionality of the system, create variables that
            can be generated from the others if they aren't already loaded
        """
        # Check that we have the minimum required data
        if ("lat_rho" or "lon_rho") not in self.__dict__:
            raise AttributeError("grid does not have attribute lat_rho or lon_rho")
        
        # Check that it is formatted into 2-D
        if self.lat_rho.ndim==1 and self.lon_rho.ndim==1:
            [self.lon_rho, self.lat_rho] = np.meshgrid(self.lon_rho, 
                                                       self.lat_rho)
        
        # Compute the dimensions
        self.ln = self.lat_rho.shape[0]
        self.lm = self.lat_rho.shape[1]
        
        pass
        
    def set_dims(self):
        """
          Using the available information, compute the remaining fields that
          we may be missing for this grid.
        """  
        # Set the number of layers
        if "n" not in self.__dict__:
            if "s_rho" in self.__dict__:
                self.n = self.s_rho.size
            elif "depths" in self.__dict__:
                self.n = self.depths.size
            
        # Generate the u- and v-grids
        if ("lat_u" or "lon_u") not in self.__dict__:
            if self.cgrid is True:
                self.lat_u = 0.5 * ( self.lat_rho[:,1:] - self.lat_rho[:,0:-1])
                self.lon_u = 0.5 * ( self.lon_rho[:,1:] - self.lon_rho[:,0:-1])
            else:
                self.lat_u = self.lat_rho
                self.lon_u = self.lon_rho
        if ("lat_v" or "lon_v") not in self.__dict__:
            if self.cgrid is True:
                self.lat_v = 0.5 * ( self.lat_rho[1:,:] - self.lat_rho[0:-1,:])
                self.lon_v = 0.5 * ( self.lon_rho[1:,:] - self.lon_rho[0:-1,:])
            else:
                self.lat_v = self.lat_rho
                self.lon_v = self.lon_rho
        if "mask_rho" in self.__dict__:
            if "mask_u" not in self.__dict__:
                if self.cgrid is True:
                    self.mask_u = self.mask_rho[:,1:] * self.mask_rho[:,0:-1]
                else:
                    self.mask_u = self.mask_rho
            if "mask_v" not in self.__dict__:
                if self.cgrid is True:
                    self.mask_v = self.mask_rho[1:,:] * self.mask_rho[0:-1,:]
                else:
                    self.mask_v = self.mask_rho
        
        # Compute the resolution
        if "pm" in self.__dict__:
            self.dm = 1.0/self.pm
        else:
            self.dm = np.ones(self.lon_rho.shape)
            self.dm[:,0:-1] = seapy.earth_distance( self.lon_rho[:,1:], \
                                              self.lat_rho[:,1:], \
                                              self.lon_rho[:,0:-1], \
                                              self.lat_rho[:,0:-1] )
            self.dm[:,-1] = self.dm[:,-2]
        if "pn" in self.__dict__:
            self.dn = 1.0/self.pn
        else:
            self.dn = np.ones(self.lat_rho.shape)
            self.dn[0:-1,:] = seapy.earth_distance( self.lon_rho[1:,:], \
                                              self.lat_rho[1:,:], \
                                              self.lon_rho[0:-1,:], \
                                              self.lat_rho[0:-1,:] )
            self.dn[-1,:] = self.dn[-2,:]
        
        # Compute the Coriolis
        if "f" not in self.__dict__:
            omega=2*np.pi/86400;
            self.f=2*omega*np.sin(np.radians(self.lat_rho))

        pass
        
    def set_mask_h(self, fld=None, bad=None):
        """
          Given a 3D field for a z-level model, compute a mask and an h
          based on where the values exist
        """
        if fld == None and self._nc != None:
            # Try to load a field from the file
            for f in ["temp","temperature"]:
                if f in self._nc.variables:
                    fld = self._nc.variables[f][0,:,:,:]
                    fld = np.ma.array( fld, mask=np.isnan(fld) )
                    break
        
        # If we don't have a field to examine, then we cannot compute the
        # mask and bathymetry
        if fld == None:
            raise AttributeError("Missing 3D field to evaluate")

        # Next, we go over the field to examine the depths and mask
        self.h = np.zeros(self.lat_rho.shape)
        self.mask_rho = np.zeros(self.lat_rho.shape)
        for k in range(0,self.depth.size):
            water = np.nonzero( np.logical_not(fld.mask[k,:,:]) )
            self.h[water] = self.depth[k]
            if k==0:
                self.mask_rho[water] = 1.0
        pass
        
        
    def set_depth(self):
        """
          Create a depth array for the model grid.
        """
        if self._isroms:
            if "s_rho" not in self.__dict__:
                self.s_rho, self.cs_r = seapy.roms.stretching(self.vtransform,
                   self.theta_s, self.theta_b, self.hc, self.n)
            self.depth_rho = seapy.roms.depth(self.vtransform, 
                 self.h, self.hc, self.s_rho, self.cs_r)
            self.depth_u=seapy.model.rho2u(self.depth_rho)
            self.depth_v=seapy.model.rho2v(self.depth_rho)
        else:
            d = self.depth.copy()
            l = np.nonzero(d>0)
            d[l]=-d[l]
            self.depth_rho = np.kron( np.kron( d,
                                    np.ones(self.lon_rho.shape[1])),
                                    np.ones(self.lon_rho.shape[0])).reshape(
                                    [self.depth.size,self.lon_rho.shape[0],
                                     self.lon_rho.shape[1]])
            if self.cgrid is True:
                self.depth_u=seapy.model.rho2u(self.depth_rho)
                self.depth_v=seapy.model.rho2v(self.depth_rho)
            else:
                self.depth_u = self.depth_rho
                self.depth_v = self.depth_rho

        pass
