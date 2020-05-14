#!/usr/bin/env python
"""
  tide.py

  Functions for working with tidal time-series

  Written by Glenn Carter and Dale Partridge
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""
import numpy as np
import datetime
from collections import namedtuple
import os
from warnings import warn

amp_phase = namedtuple('amp_phase', 'amp phase')
tellipse = namedtuple('tellipse', 'major minor angle phase')
vuf_vals = namedtuple('vuf_vals', 'v u f')
__cformat = namedtuple('__cformat', 'freq doodson semi sat shallow')
__satinfo = namedtuple('__satinfo', 'deldood phcorr amprat ilatfac')
__shallowinfo = namedtuple('__shallowinfo', 'isshallow iname coef')
# Load the constituent data when the module is imported
__reftime = datetime.datetime(1899, 12, 31, 12, 0, 0)

with np.load(os.path.dirname(__file__) + "/constituents.npz", allow_pickle=True) as data:
    __const_file = data['__const'][()]

__const = {}
for f in list(__const_file.keys()):
    __shallow = __shallowinfo(__const_file[f]['shallow']['isshallow'],
                              __const_file[f]['shallow']['iname'], __const_file[f]['shallow']['coef'])
    __sat = __satinfo(__const_file[f]['sat']['deldood'],
                      __const_file[f]['sat']['phcorr'],
                      __const_file[f]['sat']['amprat'],
                      __const_file[f]['sat']['ilatfac'])
    __const[f] = __cformat(np.float(__const_file[f]['freq']),
                           __const_file[f]['doodson'],
                           np.float(__const_file[f]['semi']),
                           __sat, __shallow)

default_tides = ['M4', 'K2', 'S2', 'M2', 'N2',
                 'K1', 'P1', 'O1', 'Q1', 'MF', 'MM']


def _set_tides(tides=None):
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
    return np.array([__const[x].freq for x in _set_tides(tides)])


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


def vuf(time, tides=None, lat=55):
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

    tides = _set_tides(tides)

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
        shtides = __const[c].shallow.iname if __const[
            c].shallow.isshallow else [c]
        v, u, f = 0, 0, 1
        for i, s in enumerate(shtides):
            vtemp = np.fmod(np.dot(__const[s].doodson, astro) +
                            __const[s].semi, 1) * (2 * np.pi)

            fsum = 1
            for k in range(len(__const[s].sat.ilatfac)):
                uu = np.fmod(np.dot(__const[s].sat.deldood[k, :], astro[3:])
                             + __const[s].sat.phcorr[k], 1)
                rr = {
                    0: __const[s].sat.amprat[k],
                    1: __const[s].sat.amprat[k] * 0.36309 *
                    (1.0 - 5.0 * slat**2) / slat,
                    2: __const[s].sat.amprat[k] * 2.59808 * slat
                }.get(__const[s].sat.ilatfac[k], 0)
                fsum += rr * np.exp(1j * 2 * np.pi * uu)
            ftemp = np.absolute(fsum)
            utemp = np.angle(fsum)

            if __const[c].shallow.isshallow:
                v += vtemp * __const[c].shallow.coef[i]
                u += utemp * __const[c].shallow.coef[i]
                f *= ftemp**np.abs(__const[c].shallow.coef[i])
            else:
                v, u, f = vtemp, utemp, ftemp
        vufs[c] = vuf_vals(v, u, f)

    return vufs


def vel_ellipse(u, v):
    """
    Generate ellipse parameters from the U and V amplitude and phases

    Parameters
    ----------
    u : dict,
       Dictionary of the u-component of velocity tides with the
       constituent name as a key, and the value is an amp_phase namedtuple.
    v : dict,
       Dictionary of the v-component of velocity tides with the
       constituent name as a key, and the value is an amp_phase namedtuple.

    Returns
    -------
    ellipse: dict,
       Dictionary of the tidal ellipses with the constituent name as the
       keys and the values are the parameters of the tellipse namedtuple.

    Examples
    --------
    With a set of u and v amplitude and phases, generate the ellipse
    parameters:
    >>> u = {"M2":amp_phase(4.3, np.radians(231.4))}
    >>> v = {"M2":amp_phase(0.7, np.radians(10.1))}
    >>> ell = vel_ellipse(u, v)
    >>> print(ell)
    {'M2': ellipse(major=4.3324053381635519, minor=0.45854551121121889,
     angle=6.1601050372480319, phase=4.0255995338808006)}
    """
    ell = {}
    for c in u:
        # Compute the parameters of the tide
        au = u[c].amp * np.exp(-1j * u[c].phase)
        av = v[c].amp * np.exp(-1j * v[c].phase)
        rccw = (au + 1j * av) / 2.0
        rcw = ((au - 1j * av) / 2.0).conjugate()
        theta_ccw = np.angle(rccw)
        theta_cw = np.angle(rcw)
        rccw = np.abs(rccw)
        rcw = np.abs(rcw)

        # Set the ellipse parameters
        major = rccw + rcw
        minor = rccw - rcw
        phase = (theta_cw - theta_ccw) / 2.0
        angle = (theta_cw + theta_ccw) / 2.0
        phase = np.mod(phase, 2 * np.pi) if phase > 0 else phase + 2 * np.pi
        angle = np.mod(angle, 2 * np.pi) if angle > 0 else angle + 2 * np.pi

        # Store the result
        ell[c.upper()] = tellipse(major, minor, angle, phase)

    return ell


def predict(times, tide, tide_minor=None, lat=55, tide_start=None):
    """
    Generate a tidal time-series for the given tides. Nodal correction
    is applied for the time as well as the given latitude (if specified).

    Parameters
    ----------
    times : datetime array,
        The times of the predicted tide(s)
    tide : dict,
        Dictionary of the tides to predict with the constituent name as
        the key, and the value is an amp_phase namedtuple.
    tide_minor : dict optional,
        Dictionary of the minor axis amplitude and angle to predict with
        the constituent name as the key, and the value is an amp_phase namedtuple.
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
    >>> tide = {'M2': amp_phase(2.3, np.radians(22)),
                'K1': amp_phase(0.4, np.radians(227))}
    >>> z = predict(times, tide, lat=23)


    """
    times = np.atleast_1d(times)
    clist = list(tide.keys())
    freq = frequency(clist)

    # If given a tide_start, then the phase is relative to that datetime,
    # and no corrections need to be applied; furthermore, the times to predict
    # are relative to tide_start.
    #
    # If no tide_start is given, then everything is done as per standard.
    if tide_start:
        vufs = dict((ap.upper(), vuf_vals(0, 0, 1)) for ap in clist)
        hours = np.array(
            [(t - tide_start).total_seconds() / 3600.0 for t in times])
    else:
        ctime = times[0] + (times[-1] - times[0]) / 2
        vufs = vuf(ctime, clist, lat)
        hours = np.array([(t - ctime).total_seconds() / 3600.0 for t in times])

    # Calculate time series
    ts = np.zeros(len(times))
    if tide_minor:
        ts = np.zeros(len(times), dtype=np.complex)
    for i, ap in enumerate(tide):
        c = tide[ap]
        ap = ap.upper()
        if tide_minor:
            m = tide_minor[ap]

            ts += np.exp(1j * m.phase) * (
                c.amp * vufs[ap].f * np.cos(2.0 * np.pi * np.dot(freq[i], hours)
                                            + (vufs[ap].v + vufs[ap].u) - c.phase)
                + m.amp * vufs[ap].f * np.sin(2.0 * np.pi * np.dot(freq[i], hours)
                                            + (vufs[ap].v + vufs[ap].u) - c.phase))
        else:
            ts += c.amp * vufs[ap].f * np.cos(2.0 * np.pi * np.dot(freq[i], hours)
                                              + (vufs[ap].v + vufs[ap].u) - c.phase)

    return ts


def fit(times, xin, tides=None, lat=55, tide_start=None, trend=True):
    """
    Perform a harmonic fit of tidal constituents to a time-series of data. The
    series can be unevenly spaced in time, but every time must be specified.
    Note that returned amplitudes and phases are zero if time series is not long
    enough to resolve specific tides

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
    tide_start : datetime, optional,
        If specified, the phases of the fit will be relative to the
        given tide_start day.
    trend : boolean, optional,
        If True (default), adjust the fit for a linear time trend

    Returns
    -------
    dict :
        A dictionary of the results are returned with the following keys:
        'tide_start': Reference date used as phase reference in fit
        'fit': the time-series of the tidal fit
        'percent': percentage of the signal explained by the tides
        'major': dictionary of the major axis fit comprised of:
            tidename: amp_phase namedtuple
            Providing the amplitude and phase of each fit constituent
        'minor': dictionary of the minor axis fit comprised of:
            tidename: amp_phase namedtuple
            Providing the minimum axis amplitude and angle of each fit constituent

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
    tides = _set_tides(tides)

    # Exclude long period tides if time series not long enough
    freq = frequency(tides)
    total_tides = len(tides)
    total_hours = (times[-1] - times[0]).total_seconds() / 3600
    invalid_tides = [t for t, f in zip(tides, 1 / freq) if 2 * f > total_hours]
    tides = [t for t in tides if t not in invalid_tides]
    freq = frequency(tides)

    # time series as hours from ctime
    ctime = times[0] + (times[-1] - times[0]) / 2
    hours = np.array([(t - ctime).total_seconds() / 3600.0 for t in times])

    # Generate cosines and sines for all the requested constitutents.
    if trend:
        A = np.hstack([np.cos(2 * np.pi * np.outer(hours, freq)),
                       np.sin(2 * np.pi * np.outer(hours, freq)),
                       np.atleast_2d(hours).T,
                       np.ones((len(hours), 1))])
    else:
        A = np.hstack([np.cos(2 * np.pi * np.outer(hours, freq)),
                       np.sin(2 * np.pi * np.outer(hours, freq)),
                       np.ones((len(hours), 1))])

    # Calculate coefficients
    ntides = len(tides)
    coef = np.linalg.lstsq(A, xin, rcond=None)[0]
    xout = np.dot(A[:, :2 * ntides], coef[:2 * ntides])

    # Explained variance
    var_exp = 100 * (np.cov(np.real(xout)) + np.cov(np.imag(xout))) / \
        (np.cov(np.real(xin)) + np.cov(np.imag(xin)))

    # Calculate amplitude & phase
    ap = (coef[:ntides] - 1j * coef[ntides:2 * ntides]) / 2.0
    am = (coef[:ntides] + 1j * coef[ntides:2 * ntides]) / 2.0

    # Nodal Corrections values
    vufs = vuf(ctime, tides, lat)
    if tide_start:
        vufs_ref = vuf(tide_start, tides, lat)
        for v in vufs:
            vufs[v] = vuf_vals(np.mod(vufs[v].v - vufs_ref[v].v, 2.0 * np.pi),
                               np.mod(vufs[v].u - vufs_ref[v].u, 2.0 * np.pi),
                               vufs[v].f / vufs_ref[v].f)

    # Compute major/minor axis amplitude and phase
    maj_amp = np.zeros((total_tides,))
    maj_pha = maj_amp.copy()
    min_amp = maj_amp.copy()
    min_pha = maj_amp.copy()
    for i, c in enumerate(tides):
        maj_amp[i] = (np.abs(ap[i]) + np.abs(am[i])) / vufs[c].f
        min_amp[i] = (np.abs(ap[i]) - np.abs(am[i])) / vufs[c].f
        min_pha[i] = np.mod(
            ((np.angle(ap[i]) + np.angle(am[i])) / 2), np.pi)
        maj_pha[i] = np.mod(vufs[c].v + vufs[c].u - np.angle(ap[i]) + min_pha[i],
                            2.0 * np.pi)

    return {
        'tide_start': tide_start,
        'fit': xout,
        'percent': var_exp,
        'major': pack_amp_phase(tides + invalid_tides, maj_amp, maj_pha),
        'minor': pack_amp_phase(tides + invalid_tides, min_amp, min_pha)
    }


def pack_amp_phase(tides, amp, phase):
    """
    Pack together the tides, amplitudes, and phases into a dictionary of
    tide names as keys and amp_phase namedtuples as values

    Parameters
    ----------
    tides : ndarray of strings,
       List of constituent strings
    amp : ndarray,
       Amplitudes for each constituent
    phase : ndarray,
       phases (rad) for each constituent

    Returns
    -------
       dict:  amplitudes and phases
        constituent name as key and amp_phase namedtuple as value
    """
    tides = _set_tides(tides)
    amp = np.atleast_1d(amp)
    phase = np.atleast_1d(phase)
    if np.any(phase > 2 * np.pi):
        warn("Phases appear to be degrees. Beware of results.")

    amp_ph = {}
    for i, c in enumerate(tides):
        amp_ph[c] = amp_phase(amp[i], phase[i])

    return amp_ph


def unpack_amp_phase(amp_ph, tides=None):
    """
    Given a dictionary of amplitude and phase parameters, unpack the dictionary
    and return an amp_phase namedtuple with arrays for each element that
    are in the order of the tides array.

    Parameters
    ----------
    amp_ph : dict,
      amp_phase dictionary returned by fit or pack_amp_phase
    tides: list of strings,
      List of strings to provide the order of the arrays

    Returns
    -------
      amp_phase : namedtuple

    """
    tides = _set_tides(tides)
    am = np.zeros((len(tides),))
    ph = am.copy()

    for n, t in enumerate(tides):
        try:
            am[n] = amp_ph[t].amp
            ph[n] = amp_ph[t].phase
        except KeyError:
            continue
    return amp_phase(am, ph)


def unpack_vuf(vuf, tides=None):
    """
    Given a dictionary of vuf parameters, unpack the dictionary
    and return a vuf_vals namedtuple with arrays for each element that
    are in the order of the tides array.

    Parameters
    ----------
    vuf : dict,
      vuf_vals dictionary returned by vuf
    tides: list of strings,
      List of strings to provide the order of the arrays

    Returns
    -------
      vuf_vals : namedtuple

    """
    tides = _set_tides(tides)
    v = np.zeros((len(tides),))
    u = v.copy()
    f = v.copy()

    for n, t in enumerate(tides):
        try:
            v[n] = vuf[t].v
            u[n] = vuf[t].u
            f[n] = vuf[t].f
        except KeyError:
            continue
    return vuf_vals(v, u, f)


def unpack_ellipse(ellipse, tides=None):
    """
    Given a dictionary of tidal ellipse parameters, unpack the dictionary
    and return an ellipse namedtuple with arrays for each element that
    are in the order of the tides array.

    Parameters
    ----------
    ellipse : dict,
      Ellipse dictionary returned by tidal_ellipse
    tides: list of strings,
      List of strings to provide the order of the arrays

    Returns
    -------
      tellipse : namedtuple

    """
    tides = _set_tides(tides)
    mj = np.zeros((len(tides),))
    mn = mj.copy()
    ph = mj.copy()
    ag = mj.copy()

    for n, t in enumerate(tides):
        try:
            mj[n] = ellipse[t].major
            mn[n] = ellipse[t].minor
            ph[n] = ellipse[t].phase
            ag[n] = ellipse[t].angle
        except KeyError:
            continue
    return tellipse(mj, mn, ag, ph)
