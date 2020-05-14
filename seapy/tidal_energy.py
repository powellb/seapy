#!/usr/bin/env python
"""
tidal_energy.py

State Estimation and Analysis for PYthon

Module to compute tidal energy from a column of data.

Written by Brian Powell on 03/30/16
Copyright (c)2020 University of Hawaii under the MIT-License.

Notes
-----

Barotropic to Baroclinic conversion is given by:

.. math::

    C=1 / T_t \int_0^T_t P'_t * wbar_t * dt,                    (1)

where, T_t is the tidal period for consituent, t, P' is the pressure perturbation,
wbar is the vertical velocity. Hence, conversion is the time average of the
vertical motion of the bottom pressure perturbation. We can do it spectrally if
we represent P'_t and wbar_t as waves:

.. math::

    P'_t = Amp_P'_t * sin( 2 * pi * t / T_t + Pha_P'_t )                   (2) \\
    wbar_t = Amp_wbar_t * sin( 2 * pi * t / T_t + Pha_wbar_t )             (3)

If we substitute (2) and (3) into (1) using trig. identity and integrate over
the period (recall that integrating a wave over one period is zero):

.. math::

    Conversion = 0.5 * Amp_P'_t * Amp_wbar_t * cos( Pha_P'_t - Pha_wbar_t )(4)

Energy Flux is given by:

.. math::

    Flux_u = 1 / T_t * \int_0^T_t u'_t * P'_t * dt,                        (5) \\
    Flux_v = 1 / T_t * \int_0^T_t v'_t * P'_t * dt,                        (6)

where u' and v' are the velocity anomalies for the constituent, t. As per
above, we can express as waves to yield:

.. math::

    Flux_u = 0.5 * Amp_u'_t * Amp_P'_t * cos( Pha_u'_t - Pha_P'_t )        (7) \\
    Flux_v = 0.5 * Amp_v'_t * Amp_P'_t * cos( Pha_v'_t - Pha_P'_t )        (8)

Displacement is given by:

.. math::

    Displace = \int_0^T_t/2 g * rho'_t / ( rho0 * N_t**2 ) * dt,           (9)

where rho' is the density anomaly and N**2 is the Brunt-Vaisala. NOTE:
this is integrated over one-half period becuase (by definition), it would
integrate to zero. However, if we know the tidal vertical velocity, then
we can integrate it for one-half period for the todal displacement:

.. math::

    Displace = \int_0^T_t/2 w_t * dt \\
        = \int_0^T_t/2 Amp_w_t * sin( 2 * pi * t / T_t )               (10) \\
        = Amp_w_t * T_t / pi

Horizontal Kinetic Energy is given by:

.. math::

    HKE = 0.5 * rho0 * 1 / T_t * \int_0^T_t (u'_t**2 + v'_t**2) * dt        (11)

substitute u' and v' as waveforms and integrate over a period,

.. math::

    HKE = 0.5 * rho0 * 0.5 * ( Amp_u'_t**2 _ Amp_v'_t**2 )                  (12)

Available Potential Energy is given by:

.. math::

    APE = 0.5 * rho0 * 1 / T_t * \int_0^T_t N_t**2 * Displace_t**2 * dt     (13)

For this, we will use the time-average N**2 (not at the specific tidal
frequency) and use (10); hence, it becomes:

.. math::

    APE = 0.5 * rho0 * (Amp_w_t * T_t / pi)**2 * 1/T_t \int_0^T_t N**2 * dt (14)
"""

import numpy as np
import seapy

_rho0 = 1000


class energetics():
    """
      This class is a container for the energetics produced by the tidal_energy
      calculation to simplify access to the resulting data.
    """

    def __init__(self, tides, energy, integrals, ellipse):
        try:
            self.tides = tides.tolist()
        except AttributeError:
            self.tides = tides
        self.energy = energy
        if len(tides) != energy.shape[0]:
            raise ValueError(
                "The number of tides and energy values are inconsistent")
        self.integrals = integrals
        self.ellipse = ellipse
        pass

    def __getitem__(self, key):
        """
        Return the energetics for a tidal constituent
        """
        t = self.tides.index(key.upper())
        return {"conversion": self.integrals[t, 2],
                "flux_u": self.energy[t, :, 0],
                "flux_v": self.energy[t, :, 1],
                "disp": self.energy[t, :, 2],
                "hke": self.energy[t, :, 3],
                "ape": self.energy[t, :, 4],
                "total_flux_u": self.integrals[t, 0],
                "total_flux_v": self.integrals[t, 1],
                "total_hke": self.integrals[t, 3],
                "total_ape": self.integrals[t, 4],
                "ellipse": self.ellipse[key.upper()]}
        pass


def tidal_energy(time, hz, u, v, w, pressure, bvf=None, tides=None,
                 ubar=None, vbar=None, wbar=None):
    """
    Calculate aspects of tidal energy from the given data: baroclinic energy flux,
    HKE, APE, displacement, and conversion.

    This only works for a single depth profile, and the arrays are to be 2D with
    dimensions of [time, depth] with depth index 0 as the bottom and inded -1 as
    the surface. Likewise, the hz field is oriented the same.

    Parameters
    ----------
    time : list of datetime,
      times of data
    hz : ndarray,
      Thickness of the water column represented by 3D quantities [m]
    u : ndarray,
      u-component of 3D velocity [m s**-1]
    v : ndarray,
      v-component of 3D velocity  [m s**-1]
    w : ndarray,
      w-component of 3D velocity  [m s**-1]
    pressure : ndarray,
      pressure of the 3D field [dbar]
    bvf : ndarray, optional
      Brunt-Vaisala Frequency of the 3D field [s**-1]. If not specified
      the APE will not be computed
    tides: list of strings, optional
      The names of the tides to use for analysis. If none
      provided, use the defaults from seapy.tide
    ubar : ndarray, optional
      u-component of barotropic velocity [m s**-1]. If none
      provided, compute from u
    vbar : ndarray, optional
      v-component of barotropic velocity [m s**-1]. If none
      provided, compute from v
    wbar : ndarray, optional
      w-component of barotropic velocity [m s**-1]. If none
      provided, compute from w

    Returns
    -------
    energetics : class,
      The energetics for each tidal consituent as well as the
      vertically integrated properties. The energetics class
      provides various methods for accessing the data
    """

    # Ensure arrays as needed
    u = np.ma.array(u)
    v = np.ma.array(v)
    w = np.ma.array(w)
    pressure = np.ma.array(pressure)

    # Setup the thicknesses in time
    hz = np.ma.array(hz)
    if hz.ndims == 1:
        hz = np.tile(hz, (u.shape[0]))
        total_h = np.sum(hz, axis=1)
        ndep = hz.shape[1]

    # If BVF not given, set to zero
    if bvf:
        bvf = np.ma.array(bvf).mean(axis=0)
    else:
        bvf = np.zeros(ndep)

    # Setup the tides
    tides = seapy.tide._set_tides(tides)
    ntides = len(tides)
    period = 3600 / seapy.tide.frequency(tides)

    # Setup the barotropic velocities
    if ubar and vbar:
        ubar = np.ma.array(ubar)
        vbar = np.ma.array(vbar)
        wbar = np.ma.array(wbar)
    else:
        ubar = np.sum(hz * u, axis=1) / total_h
        vbar = np.sum(hz * v, axis=1) / total_h
        wbar = np.sum(hz * w, axis=1) / total_h

    # Calculate Pressure Anomalies
    p_prime = pressure - pressure.mean(axis=0)
    # Apply baroclinicity
    p_prime -= np.sum(p_prime * hz) / np.sum(hz)

    # Calculate tides
    tidal_vel = seapy.tide.fit(time, ubar + 1j * vbar, tides)
    wbar = seapy.tide.fit(time, wbar, tides)

    # Store the tidal ellipse
    ellipse = {}
    for t in tides:
        ellipse[t] = seapy.tide.tellipse(tidal_vel['major'][t].amp,
                                         tidal_vel['minor'][t].amp,
                                         tidal_vel['minor'][t].phase,
                                         tidal_vel['major'][t].phase)

    # Velocity anomalies
    u_prime = u - u.mean(axis=0) - np.real(tidal_vel['fit'])
    v_prime = v - v.mean(axis=0) - np.imag(tidal_vel['fit'])
    w_prime = w - v.mean(axis=0) - wbar['fit']
    wbar = wbar['major']

    # Set the results structure: for each tide, and for each
    # depth, we will store five values (flux_u, flux_v,
    # displacement, HKE, APE)
    energy = np.zeros((ntides, ndep, 5))

    # For vertically integrated, we will also have five values:
    # (flux_u, flux_v, conversion, HKE, APE)
    integrals = np.zeros((ntides, 5))

    # Compute over all depths
    for d in seapy.progressbar.progress(np.arange(ndep)):
        # Generate the tidal components
        t_pres = seapy.tide.fit(time, p_prime[:, d], tides)['major']

        # velocity
        t_u = seapy.tide.fit(time, u_prime[:, d], tides)['major']
        t_v = seapy.tide.fit(time, v_prime[:, d], tides)['major']
        t_w = seapy.tide.fit(time, w_prime[:, d], tides)['major']

        # Compute each term for each tide
        for n, t in enumerate(tides):
            # If this is the bottom, generate the conversion
            if d == 0:
                integrals[n, 2] = 0.5 * t_pres[t].amp * wbar[t].amp * \
                    np.cos(t_pres[t].phase - wbar[t].phase)

            # Calculate Energy Flux
            energy[n, d, 0] = 0.5 * t_u[t].amp * \
                t_pres[t].amp * np.cos(t_u[t].phase - t_pres[t].phase)
            energy[n, d, 1] = 0.5 * t_v[t].amp * \
                t_pres[t].amp * np.cos(t_v[t].phase - t_pres[t].phase)

            # Calculate Displacement
            energy[n, d, 2] = t_w[t].amp * 3600 * period[t] / np.pi

            # Calculate HKE and APE
            energy[n, d, 3] = 0.5 * _rho0 * bvf[d] * displace
            energy[n, d, 4] = 0.25 * _rho0 * (t_u[t].amp + t_v[t].amp)

    # Vertically Integrate
    for n, t in enumerate(tides):
        for i in [0, 1, 3, 4]:
            integrals[n, i] = np.sum(energy[n, :, i] * hz, axis=1) / total_h

    # Put it all together to return
    return energetics(tides, energy, integrals, ellipse)
