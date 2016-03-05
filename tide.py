#!/usr/bin/env python
"""
  tide.py

  Functions for working with tidal time-series

  Written by Glenn Carter and Dale Partridge
  Copyright (c)2016 University of Hawaii under the BSD-License.
"""
import numpy as np
import datetime
from collections import namedtuple
import os

amp_pha = namedtuple('amp_pha', 'amp pha')
__cformat = namedtuple('__cformat', 'freq doodson semi sat')
__satinfo = namedtuple('__satinfo', 'deldood phcorr amprat ilatfac')
__vuf_vals = namedtuple('__vuf_vals', 'v u f')

# Load the constituent data when the module is imported
__reftime = datetime.datetime(1899, 12, 31, 12, 0, 0)

with np.load(os.path.dirname(__file__) + "/constituents.npz") as data:
    __const_file = data['__const'][()]

__const = {}
for f in list(__const_file.keys()):
    __sat = __satinfo(__const_file[f]['sat']['deldood'], __const_file[f]['sat'][
        'phcorr'], __const_file[f]['sat']['amprat'], __const_file[f]['sat']['ilatfac'])
    __const[f] = __cformat(np.float(__const_file[f]['freq']), __const_file[f][
        'doodson'], np.float(__const_file[f]['semi']), __sat)

default_tides = ['M4', 'K2', 'S2', 'M2', 'N2',
                 'K1', 'P1', 'O1', 'Q1', 'MF', 'MM']


def __set_tides(tides=None):
    """
    Private method: make sure that tides passed in are a proper array
    and in uppercase
    """
    return default_tides if tides is None else np.array([t.upper() for t in np.atleast_1d(tides)])


def frequency(tides=None):
    """
    Returns array of the frequency in cycles per hour for the requested 
    tidal constituents.

    Parameters
    ----------
    tides: string or list of strings,
        List of consituents to obtain the frequency for

    Returns
    -------
    ndarray:
        List of frequencies matching the tidal names

    Examples
    --------
    Display the period (in hours) of M2 and S2

    >>> 1 / frequency(['M2', 'S2'])
    array([ 12.4206012,  12.       ])

    """
    return np.array([__const[x].freq for x in __set_tides(tides)])


def __astron(ctime):
    """
    PRIVATE METHOD
    --------------
        This subroutine calculates the following five ephermides of the sun 
        and moon following t_tide's t_astron.
                h = mean longitude of the sum
                pp = mean longitude of the solar perigee
                s = mean longitude of the moon
                p = mean longitude of the lunar perigee
                np = negative of the longitude of the mean ascending node
        and their rates of change. Units for the ephermides are cycles and
        for their derivatives are cycles/365 days
        The formulae for calculating this ephermides were taken from pages 98 
        and 107 of the Explanatory Supplement to the Astronomical  
        Ephermeris and the American Ephermis and Nautical Almanac (1961).
    """
    # Compute number of days from epoch of 12:00 UT Dec 31, 1899.
    dt = ctime - __reftime
    d = dt.days + dt.seconds / 86400

    # Coefficients used in t_tide
    sc = np.array([270.434164, 13.1763965268, -0.0000850, 0.000000039])
    hc = np.array([279.696678, 0.9856473354, 0.00002267, 0.000000000])
    pc = np.array([334.329556, 0.1114040803, -0.0007739, -0.00000026])
    npc = np.array([-259.183275, 0.0529539222, -0.0001557, -0.000000050])
    ppc = np.array([281.220844, 0.0000470684, 0.0000339, 0.000000070])
    coefs = np.vstack((sc, hc, pc, npc, ppc))

    # output variables
    astro = np.empty(6,)
    ader = np.empty(6,)

    # Compute astronomical constants at ctime; we only need the fractional
    # part of the cycle.
    D = d * 1e-4
    args = np.vstack((1, d, D**2, D**3))
    astro[1:] = (np.remainder(np.dot(coefs, args) / 360, 1)).ravel()

    # Compute lunar time tau, based on fractional part of solar day.
    # Add the hour angle to the longitude of the sun and subtract the
    # longitude of the moon.
    astro[0] = (ctime - datetime.datetime(ctime.year, ctime.month, ctime.day)) / \
        datetime.timedelta(1) + astro[2] - astro[1]

    # Compute rates of change.
    dargs = np.vstack((0, 1, 2e-4 * D, 3e-4 * D**2))

    ader[1:] = (np.dot(coefs, dargs) / 360).ravel()
    ader[0] = 1.0 + ader[2] - ader[1]

    return astro, ader


def _vuf(time, tides=None, lat=5.0):
    """
    Returns the astronomical phase V, the nodal phase modulation U, and 
    the nodal amplitude correction F at ctime for the requested 
    constituents and latitude.

    Note, there are differences compared to t_tide:
        - That V and U are returned as radians, not cycles. 
        - Shallow water constituents are not supported.

    Parameters
    ----------
    time : datetime,
        Time for the nodal correction
    tides : string or list of strings, optional,
        Names of tidal constituents
    lat : float
        Latitude for nodal correction

    Returns
    -------
    dict : 
        dictionary of vuf corrections with the tide names as keys
    """
    # If the user tries to pass a list, we only want the mid-point;
    # otherwise, just use what was sent
    try:
        time = time[np.floor(len(time) / 2)]
    except:
        pass

    tides = __set_tides(tides)

    # Calculate astronomical arguments at the requested time (mid-point of
    # timeseries).
    astro, ader = __astron(time)

    # According to t_tide, the second order terms in the tidal potential
    # go to zero at the equator, but the third order terms do not. Hence,
    # when trying to infer the third-order terms from the second-order
    # terms the nodal correction factors blow up. In order to prevent
    # this it is assumed that the equatorial forcing at the equator is
    # due to second-order forcing off the equator from about the 5 degree
    # location. Latitudes are hence (somewhat arbitrarily) forced to be
    # no closer than 5 degrees to the equator.
    if(np.abs(lat) < 5.0):
        lat = 5 * np.sign(lat + np.finfo('float').eps)
    slat = np.sin(lat * np.pi / 180.)

    # Setup output dictionaries
    vufs = {}

    for c in tides:
        v = np.fmod(np.dot(__const[c].doodson, astro) +
                    __const[c].semi, 1) * (2 * np.pi)
        if np.isnan(v):
            v = 0
        N = len(__const[c].sat.ilatfac)
        if N == 0:
            f = 1
            u = 0
        else:
            fsum = 1
            for k in range(N):
                uu = np.fmod(np.dot(__const[c].sat.deldood[k, :], astro[3:])
                             + __const[c].sat.phcorr[k], 1)
                rr = {
                    0: __const[c].sat.amprat[k],
                    1: __const[c].sat.amprat[k] * 0.36309 *
                    (1.0 - 5.0 * slat**2) / slat,
                    2: __const[c].sat.amprat[k] * 2.59808 * slat
                }.get(__const[c].sat.ilatfac[k], 0)
                fsum += rr * np.exp(1j * 2 * np.pi * uu)
            f = np.absolute(fsum)
            u = np.angle(fsum)
        vufs[c] = __vuf_vals(v, u, f)

    return vufs


def predict(times, tide, lat=None, tide_start=None):
    """
    Generate a tidal time-series for the given tides. Nodal correction
    is applied for the time as well as the given latitude (if specified).

    Parameters
    ----------
    times : datetime array,
        The times of the predicted tide(s)
    tide : dict, 
        Dictionary of the tides to predict with the constituent name as 
        the key, and the value is an amp_ph namedtuple.
    lat : float optional,
        latitude of the nodal correction
    tide_start : datetime optional,
        If specified, the nodal correction are applied to adjust from
        the tide_start to the center of the times record

    Returns
    -------
    ndarray : 
        Values of the total tidal prediction of all constituents provided

    Examples
    --------
    Generate a time series every hour for the month of April, 2014 for
    a specific phase and amplitude of M2 and K1.

    >>> times = [ datetime(2014,4,1) + timedelta(t/24) for t in range(31*24) ]
    >>> tide = {'M2': amp_pha(2.3, np.radians(22)), 
                'K1': amp_pha(0.4, np.radians(227))}
    >>> z = predict(times, tide, 23)


    """
    times = np.atleast_1d(times)
    clist = list(tide.keys())
    freq = frequency(clist)

    # Calculate midpoint of time series
    ctime = times[0] + (times[-1] - times[0]) / 2

    # If we have a start time, then we have to adjust for that
    if tide_start:
        # vufs = _vuf(__reftime + (ctime - tide_start), clist, lat)
        vufs = _vuf(tide_start, clist, lat)
        for ap in tide:
            tide[ap] = amp_pha(tide[ap].amp / vufs[ap].f,
                               -tide[ap].pha + (vufs[ap].v + vufs[ap].u))

    # Calculate astronomical and nodal values
    vufs = _vuf(ctime, clist, lat)

    # Time series as hours from ctime
    hours = np.array([(t - ctime).total_seconds() / 3600.0 for t in times])

    # Calulate time series
    ts = np.zeros(len(times))
    for i, ap in enumerate(tide):
        c = tide[ap]
        ap = ap.upper()
        ts += c.amp * vufs[ap].f * np.cos(2.0 * np.pi * np.dot(freq[i], hours)
                                          + (vufs[ap].v + vufs[ap].u) - c.pha)

    return ts


def fit(times, xin, tides=None, lat=None, use_reftime=False):
    """
    Perform a harmonic fit of tidal constituents to a time-series of data. The
    series can be unevenly spaced in time, but every time must be specified.

    Parameters
    ----------
    times : datetime array,
        List of times matching each value of the time-series input
    xin : ndarray,
        Data values to perform the harmonic fit upon
    tides : string or list of strings, optional,
        list of tidal constituents to fit to the data provided. If not
        specified, the dominant 11 consituents are used:
        ['M4', 'K2', 'S2', 'M2', 'N2', 'K1', 'P1', 'O1', 'Q1', 'MF', 'MM']
    lat : float, optional,
        latitude of the nodal correction
    use_reftime : bool, optional,
        If True, apply nodal corrections to the fit to put the phases
        back to the reference time. If False [Default], the phases
        are relative to the returned 'tide_start'.

    Returns
    -------
    dict :
        A dictionary of the results are returned with the following keys:
        'tide_start': Reference date used as phase reference in fit
        'fit': the time-series of the tidal fit
        'percent': percentage of the signal explained by the tides
        'major': dictionary of the major axis fit comprised of:
            tidename: amp_ph namedtuple
            Providing the amplitude and phase of each fit constituent
        'minor': dictionary of the minor axis fit comprised of:
            tidename: amp_ph namedtuple
            Providing the amplitude and phase of each fit constituent

    Examples
    --------
    Fit the 'M2' and 'K1' to the given time-series, x.

    >>> data = np.load('water_level.npz')
    >>> out = fit(data.times, data.level, ['M2','K1'], lat=28)
    >>> plt.plot(data.times, out['fit'])
    """

    xin = np.atleast_1d(xin).flatten()
    if len(xin) != len(times):
        raise ValueError("The times and input data must be of same size.")
    tides = __set_tides(tides)

    # The tide_start is the beginning of the time record
    ctime = times[0]
    tide_start = times[0]

    # Nodal Corrections values
    if use_reftime:
        # Calculate midpoint of time series, but ensure 00:00 hour
        tide_start = __reftime
        ctime = times[0] + (times[-1] - times[0]) / 2
        ctime = datetime.datetime(ctime.year, ctime.month, ctime.day)
        vufs = _vuf(ctime, tides, lat)

    # time series as hours from ctime
    hours = np.array([(t - ctime).total_seconds() / 3600.0 for t in times])

    # get frequencies
    freq = frequency(tides)

    # Generate cosines and sines for all the requested constitutents.
    A = np.hstack([np.ones((len(hours), 1)),
                   np.cos(2 * np.pi * np.outer(hours, freq)),
                   np.sin(2 * np.pi * np.outer(hours, freq))])

    # Calculate coefficients
    coef = np.linalg.lstsq(A, xin)[0]
    xout = np.dot(A[:, 1:], coef[1:])

    # Explained variance
    var_exp = 100 * (np.cov(np.real(xout)) + np.cov(np.imag(xout))) / \
        (np.cov(np.real(xin)) + np.cov(np.imag(xin)))

    # Calculate amplitude & phase
    num = len(tides)
    ap = np.squeeze((coef[1:1 + num] - 1j * coef[1 + num:]) / 2.0)
    am = np.squeeze((coef[1:1 + num] + 1j * coef[1 + num:]) / 2.0)

    # Compute major/minor axis amplitude and phase
    maj_amp = np.empty((len(tides),))
    maj_pha = maj_amp.copy()
    min_amp = maj_amp.copy()
    min_pha = maj_amp.copy()
    for i, c in enumerate(tides):
        maj_amp[i] = (np.abs(ap[i]) + np.abs(am[i]))
        min_amp[i] = (np.abs(ap[i]) - np.abs(am[i]))
        min_pha[i] = np.mod(
            ((np.angle(ap[i]) + np.angle(am[i])) / 2), np.pi)
        maj_pha[i] = np.mod(np.mod(np.angle(ap[i]), 2.0 * np.pi) + min_pha[i],
                            2.0 * np.pi)
        if use_reftime:
            maj_amp[i] /= vufs[c].f
            min_amp[i] /= vufs[c].f
            maj_pha[i] = np.mod(np.mod(vufs[c].v + vufs[c].u -
                                       np.angle(ap[i]), 2.0 * np.pi) + min_pha[i],
                                2.0 * np.pi)

    return {
        'tide_start': tide_start,
        'fit': xout,
        'percent': var_exp,
        'major': make_amp_ph(tides, maj_amp, maj_pha),
        'minor': make_amp_ph(tides, min_amp, min_pha)
    }


def make_amp_ph(tides, amp, pha):
    """
    Makes a dictionary of amp_pha named tuples given arrays of names, amplitudes and phases
    Inputs:-
        tides - List of constituent strings
        amp - array of amplitudes for each constituent 
        pha - array of phases (rad) for each constituent
    Outputs:-
        amp_ph - dictionary of amp_ph named tuples, with constiutents as keys
    """
    amp_ph = {}
    for i, c in enumerate(tides):
        amp_ph[c] = amp_pha(amp[i], pha[i])

    return amp_ph
