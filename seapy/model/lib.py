#!/usr/bin/env python
"""
  lib.py

  Library of utilities for ocean models, imported into the namespace
  when importing the model module

  Written by Brian Powell on 10/18/13
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""
import numpy as np
import seapy
import scipy.constants

# Define reference constants
_R0 = 999.83
_R0a = 5.053e-3
_R0b = 0.048e-6


def _cgrid_rho_vel(rho, dim, fill):
    """
    Private Method: Compute the u- or v-grid velocity from a rho-field for a c-grid
    """
    rho = np.ma.array(rho, copy=False)
    if fill:
        rho = seapy.convolve_mask(rho, copy=True)
    shp = np.array(rho.shape)
    fore = np.product([shp[i] for i in np.arange(0, dim)]).astype(int)
    aft = np.product([shp[i]
                      for i in np.arange(dim + 1, rho.ndim)]).astype(int)
    nfld = 0.5 * (rho.reshape([fore, shp[dim], aft])[:, 0:-1, :].filled(np.nan) +
                  rho.reshape([fore, shp[dim], aft])[:, 1:, :].filled(np.nan))
    shp[dim] = shp[dim] - 1
    return np.ma.fix_invalid(nfld.reshape(shp), copy=False, fill_value=1e+37)


def rho2u(rho, fill=False):
    """
    Put the rho field onto the u field for the c-grid

    Parameters
    ----------
    rho : masked array like
        Input rho field
    fill : bool, optional
        Fill the masked data before moving grids

    Returns
    -------
    u : masked array
    """
    return _cgrid_rho_vel(rho, rho.ndim - 1, fill)


def rho2v(rho, fill=False):
    """
    Put the rho field onto the v field for the c-grid

    Parameters
    ----------
    rho : masked array like
        Input rho field
    fill : bool, optional
        Fill the masked data before moving grids

    Returns
    -------
    v : masked array
    """
    return _cgrid_rho_vel(rho, rho.ndim - 2, fill)


def u2rho(u, fill=False):
    """
    Put the u field onto the rho field for the c-grid

    Parameters
    ----------
    u : masked array like
        Input u field
    fill : bool, optional
        Fill the masked data before moving grids

    Returns
    -------
    rho : masked array
    """
    u = np.ma.array(u, copy=False)
    if fill:
        u = seapy.convolve_mask(u, copy=True)
    shp = np.array(u.shape)
    nshp = shp.copy()
    nshp[-1] = nshp[-1] + 1
    fore = np.product([shp[i] for i in np.arange(0, u.ndim - 1)]).astype(int)
    nfld = np.ones([fore, nshp[-1]])
    nfld[:, 1:-1] = 0.5 * \
        (u.reshape([fore, shp[-1]])[:, 0:-1].filled(np.nan) +
         u.reshape([fore, shp[-1]])[:, 1:].filled(np.nan))
    nfld[:, 0] = nfld[:, 1] + (nfld[:, 2] - nfld[:, 3])
    nfld[:, -1] = nfld[:, -2] + (nfld[:, -2] - nfld[:, -3])
    return np.ma.fix_invalid(nfld.reshape(nshp), copy=False, fill_value=1e+37)


def v2rho(v, fill=False):
    """
    Put the v field onto the rho field for the c-grid

    Parameters
    ----------
    v : masked array like
        Input v field
    fill : bool, optional
        Fill the masked data before moving grids

    Returns
    -------
    rho : masked array
    """
    v = np.ma.array(v, copy=False)
    if fill:
        v = seapy.convolve_mask(v, copy=True)
    shp = np.array(v.shape)
    nshp = shp.copy()
    nshp[-2] = nshp[-2] + 1
    fore = np.product([shp[i] for i in np.arange(0, v.ndim - 2)]).astype(int)
    nfld = np.ones([fore, nshp[-2], nshp[-1]])
    nfld[:, 1:-1, :] = 0.5 * \
        (v.reshape([fore, shp[-2], shp[-1]])[:, 0:-1, :].filled(np.nan) +
         v.reshape([fore, shp[-2], shp[-1]])[:, 1:, :].filled(np.nan))
    nfld[:, 0, :] = nfld[:, 1, :] + (nfld[:, 2, :] - nfld[:, 3, :])
    nfld[:, -1, :] = nfld[:, -2, :] + (nfld[:, -2, :] - nfld[:, -3, :])
    return np.ma.fix_invalid(nfld.reshape(nshp), copy=False, fill_value=1e+37)


def density(depth, temp, salt):
    """
    Calcuate density and adjoint-forcing terms of density using the quadratic
    equation of state.

    Parameters
    ----------
    depth: ndarray
      Depth of temperature and salinity values
    temp: ndarray,
      Model temperature
    salt: ndarray,
      Model salt

    Returns
    -------
    rho: ndarray,
      Density of field
    drhodt : ndarray,
      Variation of density with respect to temperature
    drhods : ndarray,
      Variation of density with respect to salt
    """
    Ba = 0.808
    Bb = 0.0085e-3
    Aa = 0.0708
    Ab = 0.351e-3
    Ac = 0.068
    Ad = 0.0683e-3
    Ga = 0.003
    Gb = 0.059e-3
    Gc = 0.012
    Gd = 0.064e-3

    Z = np.abs(np.ma.asarray(depth))
    T = np.ma.asarray(temp)
    S = np.ma.asarray(salt)

    # Calculate the density
    rho0 = _R0 + _R0a * Z - _R0b * Z * Z
    alpha = Aa * (1 + Ab * Z + Ac * (1 - Ad * Z) * T)
    beta = Ba - Bb * Z
    gamma = Ga * (1 - Gb * Z - Gc * (1 - Gd * Z) * T)
    rho = rho0 + beta * S - alpha * T - gamma * (35 - S) * T

    # Calculate drhodt
    rho0 = 0
    alpha = Aa * Ab * Z - 2 * Aa * Ac * Ad * T * Z + 2 * Aa * Ac * T + Aa
    beta = 0
    gamma = Ga * Gb * S * Z - 35 * Ga * Gb * Z - 2 * Ga * Gc * Gd * S * T * Z + \
        2 * 35 * Ga * Gc * Gd * T * Z + 2 * Ga * Gc * S * T - 2 * 35 * Ga * Gc * T - \
        Ga * S + 35 * Ga
    drhodt = rho0 + beta - alpha - gamma

    # Calculate drhods
    rho0 = 0
    alpha = 0
    beta = Ba - Bb * Z
    gamma = Ga * Gb * T * Z - Ga * Gc * Gd * T * T * Z + Ga * Gc * T * T - Ga * T
    drhods = rho0 + beta - alpha - gamma

    return rho, drhodt, drhods


def w(grid, u, v):
    """
    Compute the vertical velocity for the grid from the u and v velocity.

    For a standard, z-level model the formulation is:

    w_ij = u_ij * delta Hz / delta x + v_ij * delta Hz / delta y

    Parameters
    ----------
    grid: seapy.model.grid or string
      The grid to use to compute the vertical velocity
    u   : ndarray
      u-component of velocity
    v   : ndarray
      v-component of velocity

    Returns
    -------
    w: ndarray,
      vertical velocity as [m s**-1]
    """
    grid = seapy.model.asgrid(grid)
    u = np.ma.array(u)
    v = np.ma.array(v)

    if grid.cgrid:
        w = u2rho(0.5 * u * (grid.thick_rho[:, :, 1:] -
                             grid.thick_rho[:, :, :-1]) *
                  (grid.pm[:, 1:] + grid.pm[:, :-1])) + \
            v2rho(0.5 * v * (grid.thick_rho[:, 1:, :] -
                             grid.thick_rho[:, :-1, :]) *
                  (grid.pn[1:, :] + grid.pn[:-1, :]))
    else:
        w = np.ma.array(np.zeros(u.shape))
        w = u * grid.thick_rho * grid.pm + \
            v * grid.thick_rho * grid.pn

    return w


def pressure(hz, rho, axis=None, drhodt=None, drhods=None):
    """
    Calculate the water pressure and the adjoint forcing of pressure
    (if adjoint forcing of density is provided).

    Parameters
    ----------
    hz : ndarray
      Thickness of rho values
    rho : ndarray
      Density of water
    axis : int, optional
      Axis to integrate for pressure. If not specified, try to
      estimate the best dimension to integrate
    drhodt: ndarray, optional
      Variation of rho with temperature. Only required to compute
      the variation of pressure
    drhods: ndarray, optional
      Variation of rho with salt. Only required to compute
      the variation of pressure

    Returns
    -------
    pres: ndarray
      Pressure
    dpresdt : ndarray
      Variation of pressure with respect to temperature
    dpresds : ndarray
      Variation of pressure with respect to salt
    """
    hz = np.ma.array(hz)
    rho = np.ma.array(rho)
    dpresdt = None
    dpresds = None

    if not axis:
        # Figure out the dimensionss. If rho has fewer dimension than depth,
        # then time is the mismatched. If they are the same, we will integrate
        # the first dimension
        if hz.ndim != rho.ndim:
            nd = rho.ndim - hz.ndim
            if hz.shape == rho.shape[nd:]:
                axis = nd
            elif hz.shape == rho.shape[:-nd]:
                axis = -nd
        else:
            axis = 1

    flip = [slice(None, None, None) for i in range(rho.ndim)]
    flip[axis] = slice(None, None, -1)

    pres = np.cumsum(scipy.constants.g *
                     (rho * hz)[flip], axis=axis)[flip] * 1e-4

    # Calculate adjoint forcing
    if drhodt is not None:
        dpresdt = np.cumsum(scipy.constants.g * (drhodt * hz)[flip],
                            axis=axis)[flip] * 1e-4
    if drhods is not None:
        dpresds = np.cumsum(scipy.constants.g * (drhods * hz)[flip],
                            axis=axis)[flip] * 1e-4

    return pres, dpresdt, dpresds


def sound(depth, temp, salt):
    """
    Calcuate speed of sound and adjoint-forcing terms of sound speed using the
    quadratic equation of state.

    Parameters
    ----------
    depth: ndarray
      Depth of temperature and salinity values
    temp: ndarray,
      Model temperature
    salt: ndarray,
      Model salt

    Returns
    -------
    c: ndarray,
      Speed of Sound field
    dcdt : ndarray,
      Variation of sound speed with respect to temperature
    dcds : ndarray,
      Variation of sound speed with respect to salt
    """
    Ca = 1448.96
    Cb = 4.591
    Cc = 0.05304
    Cd = 2.374e-4
    Ce = 1.34
    Cf = 0.01630
    Cg = 1.675e-7
    Ch = 1.025e-2
    Ci = 7.139e-13

    Z = np.abs(np.ma.asarray(depth))
    T = np.ma.asarray(temp)
    S = np.ma.asarray(salt)

    c = Ca + Cb * T - Cc * T * T + Cd * T**3 + Ce * S - Ce * 35 + Cf * Z + \
        Cg * Z * Z - Ch * T * S + Ch * 35 * T - Ci * T * Z**3

    dcdt = Cb - 2 * Cc * T + 3 * Cd * T * T - Ch * S + Ch * 35 - Ci * Z**3
    dcds = Ce - Ch * T

    return c, dcdt, dcds


def bvf(depth, rho, axis=None, drhodt=None, drhods=None):
    """
    Calculate the Brunt-Vaisala Frequency of the density

    Parameters
    ----------
    depth: ndarray
      Depth of temperature and salinity values
    rho : ndarray
      Density of water
    axis : int, optional
      Axis to integrate for pressure. If not specified, try to
      estimate the best dimension to integrate
    drhodt: ndarray, optional
      Variation of rho with temperature. Only required to compute
      the variation of pressure
    drhods: ndarray, optional
      Variation of rho with salt. Only required to compute
      the variation of pressure

    Returns
    -------
    bvf : ndarray
      Brunt-Vaisala Frequency
    dbvfdt : ndarray
      Variation of BVF with respect to temperature
    dbvfds : ndarray
      Variation of BVF with respect to salt
    """
    Z = np.abs(np.ma.asarray(depth))
    rho = np.ma.asarray(rho)
    dbvfdt = None
    dbvfds = None

    if not axis:
        # Figure out the dimensionss. If rho has fewer dimension than depth,
        # then time is the mismatched. If they are the same, we will integrate
        # the first dimension
        if depth.ndim != rho.ndim:
            nd = rho.ndim - depth.ndim
            if depth.shape == rho.shape[nd:]:
                axis = nd
            elif depth.shape == rho.shape[:-nd]:
                axis = -nd
        else:
            axis = 1

    fill = [slice(None, None, None) for i in range(rho.ndim)]
    fill[axis] = slice(None, -1, None)

    rho0 = _R0 + _R0a * Z - _R0b * Z * Z
    bvf = np.zeros(rho.shape)
    bvf[fill] = - scipy.constants.g / rho0[::-1, ...][:-1, ...] * \
        np.diff(rho, axis=axis) / -np.diff(Z, axis=0)
    if drhodt:
        dbvfdt[fill] = - scipy.constants.g / rho0[::-1, ...][:-1, ...] * \
            np.diff(drhodt, axis=axis) / -np.diff(Z, axis=0)
    if drhods:
        dbvfds[fill] = - scipy.constants.g / rho0[::-1, ...][:-1, ...] * \
            np.diff(drhods, axis=axis) / -np.diff(Z, axis=0)

    return bvf, dbvfdt, dbvfds
