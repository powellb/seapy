#!/usr/bin/env python
"""
  lib.py
  
  General ROMS utils

  Written by Brian Powell on 05/24/13
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function
import numpy as np
import netcdftime

fields={ # "zeta":{"grid":"rho","dims":2},
           # "ubar":{"grid":"u","dims":2},
           # "vbar":{"grid":"v","dims":2},
           "u":{"grid":"u","dims":3}}
           # "v":{"grid":"v","dims":3},
           # "temp":{"grid":"rho","dims":3},
           # "salt":{"grid":"rho","dims":3}}


def stretching(vstretching=2, theta_s=2, theta_b=0.1, hc=100, N=10,
               w_grid=False):
    """
      Compute the stretching function for ROMS
    """
    ds=1.0/N
    if w_grid:
        lev=np.arange(0,N+1)
    else:
        lev=np.arange(1,N+1)-0.5
    s=(lev-N)*ds
    
    if vstretching == 1:
        if theta_s > 0:
            ptheta=np.sinh(theta_s*s)/np.sinh(theta_s)
            rtheta=np.tanh(theta_s*(s+0.5))/(2.0*np.tanh(0.5*theta_s)) - 0.5
            cs=(1.0-theta_b)*ptheta+theta_b*rtheta
        else:
            cs=s
        pass
    elif vstretching == 2:
        if theta_s > 0:
            csur = (1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
            if theta_b > 0:
                cbot=-1.0+np.sinh(theta_b*(s+1.0))/np.sinh(theta_b)
                weight=(s+1.0)*(1.0-s)
                cs=weight*csur+(1.0-weight)*cbot
            else:
                cs=csur
        else:
            cs=s
        pass
    elif vstretching == 3:
        if theta_s > 0:
            exp_s=theta_s
            exp_b=theta_b
            alpha=3
            cbot=np.log(np.cosh(alpha*(s+1.0)**exp_b)) / \
                     np.log(np.cosh(alpha))-1.0
            csur=-np.log(np.cosh(alpha*np.abs(s)**exp_s)) / \
                     np.log(np.cosh(alpha))
            weight=(1.0-np.tanh( alpha*(s+.5)))/2.0
            cs=weight*cbot+(1.0-weight)*csur
        else:
            cs=s
        pass
    elif vstretching == 4:
        if theta_s > 0:
            csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
        else:
            csur=-(s*s)
        pass
        if theta_b > 0:
            cs=(np.exp(theta_b*csur)-1.0)/(1.0-np.exp(-theta_b))
        else:
            cs=s
        pass
    elif vstretching == 5:
        s = -(lev*lev-2*lev*N+lev+N*N-N)/(1.0*N*N-N) - \
                0.01*(lev*lev-lev*N)/(1.0-N)
        if theta_s > 0:
            csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1)
        else:
            csur=-(s*s)
        if theta_b > 0:
            cs=(np.exp(theta_b*(csur+1.0))-1.0)/(np.exp(theta_b)-1.0)-1.0
        else:
            cs=csur;
        pass
    else:
        raise ValueError("stretching value must be between 1 and 5")

    return s, cs

def depth(vtransform=1, h=None, hc=100, scoord=None,
          stretching=None, zeta=0, w_grid=False):
    """
        Given the transform method, the bathyemtry, the sea surface height
        [optional], the critical depth, the s-coordinates, and stretching
        function (returned by the roms.utils.stretching function), determine
        the depth of the given bathymetry.
    """
    if h == None or scoord == None or stretching == None:
        raise AttributeError("you must supply h, scoord, and stretching")
    if scoord.size != stretching.size:
        raise ValueError("the stretching and scoord arrays must be the same size")
    N=scoord.size
    hinv = 1/h
    if w_grid:
        z = np.zeros([N+1,h.shape[0],h.shape[1]])
        z[0,:,:] = -h
        r = np.arange(1,N+1)
    else:
        z = np.zeros([N,h.shape[0],h.shape[1]])
        r = np.arange(0,N)

    if vtransform == 1:
        cff = hc*(scoord-stretching)
        for k in r:
            z0 = cff[k] + stretching[k]*h
            z[k,:,:] = z0+zeta * (1.0 + z0*hinv)
    elif vtransform == 2:
        cff = 1/(hc + h)
        for k in r:
            cff1 = hc*scoord[k] + h*stretching[k]
            z[k,:,:] = zeta + ( zeta + h )*cff*cff1
    else:
        raise ValueError("transform value must be between 1 and 2")

    return z

def get_timebase(time):
    """
    Give a netCDF4 time record from a ROMS file, compute the timebase
    for the file
    """
    return netcdftime.utime(time.units).origin

    
pass

