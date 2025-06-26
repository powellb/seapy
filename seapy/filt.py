#!/usr/bin/env python
"""
  filt.py

  Functions for quickly filtering time-series of data or analyzing
  spectra of data.


  Written by Brian Powell on 02/09/16
  Copyright (c)2010--2025 University of Hawaii under the MIT-License.
"""
import numpy as np
import seapy
import scipy.signal
import scipy.stats as ss

def powerspectra(x, dt=1, nsmooth=5):
    """
    Calculate the variance-preserving power spectra density from a
    time-series of data. The

    Based on previous implementation by E. Firing

    Parameters
    ----------
    x : ndarray,
        The time-series of data to calculate the power spectral density.
    dt : float,
        The time-step between the values in x.
    nsmooth : int,
        The number of time-steps to smooth the spectra over. The smoothing
        preserves the variance of the spectra.

    Returns
    -------
    freq : ndarray,
        The frequencies of the spectra.
    ps : ndarray,
        The power spectra at each frequency.
    psd : ndarray,
        The power spectral density at each frequency.
    conf : ndarray,
        The 95% confidence interval that can be multiplied to each
        power spectral density to get the range of confidence.

    Examples
    --------
    With a time-series of 'data' sampled every 3 hours, calculate the
    power-spectral density and smooth it via variance preserving at the
    24 hour period (8 dt steps).

      >>>  freq, ps, psd, conf = powerspectra(data, dt=3, nsmooth=8)
      >>>  fig, ax = plt.subplots()
      >>>  ax.loglog(freq, psd * freq)
      >>>  ax.set_title('Variance Preserving Spectra')
      >>>  ax.set_xlabel('Cycles per hour')
      >>>  ax.set_ylabel('Variance per CPH');
      >>>  plot_period = lambda x: 1/x
      >>>  tax = ax.secondary_xaxis('top', functions=(plot_period, plot_period))
      >>>  tax.set_xlabel('Period [hours]')

    """
    x = np.atleast_1d(x).copy()
    nx = len(x)
    t = np.arange(nx)

    # Remove linear trend and mean
    f = np.polyfit(t, x, 1)
    x -= np.polyval(f, t)

    # Smooth the data with simple boxcar
    xweight = 1 - ((t - 0.5 * nx) / (0.5 * nx)) ** 2
    x *= xweight

    # Generate the spectrum
    npositive = nx//2
    plus = np.s_[1:npositive]
    freqs = np.fft.fftfreq(nx, d=dt)[plus]
    ft = np.fft.fft(x)[plus]
    ps = 2 * np.abs(ft) ** 2 / (nx*nx)
    # Convert PS to Power Spectral Density
    psd = ps * dt * nx

    # Smooth the spectra with boxcar
    if nsmooth:
        pweight = np.ones(nsmooth, dtype=float) / nsmooth
        ps = np.convolve(ps, pweight, mode='valid')
        psd = np.convolve(psd, pweight, mode='valid')
        freqs = np.convolve(freqs, pweight, mode='valid')

    # Compensate for the energy suppressed by the boxcar
    psd *= nx / (xweight**2).sum()
    ps *= nx*nx / xweight.sum()**2

    # Confidence intervals
    df = nsmooth*2
    conf = np.array(df / ss.chi2.ppf([0.025, 0.975], df))

    return freqs, ps, psd, conf

def average_err(x, window=5):
    """
    Generate a moving (boxcar) average and variance of the time-series, x, using
    the specified window size.

    Parameters
    ----------
    x: ndarray,
        The time-series of data to bandpass filter.
    cutoff: float,
        The period at which the bandpass filter will apply the cutoff.
        Units are same as the time-step of the signal provided (e.g., if the
        data are provided every hour, then a cutoff=12 would be a 12 hour
        cutoff.)

    Returns
    -------
    x, variance: ndarray, ndarray
        Returns the moving average of x with the moving variance of the
        average window

    Examples
    --------
    Create data every 30 minutes for 3 days with time in days:

      >>>  t = np.linspace(0, 3.0, 2 * 24 * 3.0, endpoint=False)
      >>>  x = 0.1 * np.sin(2 * np.pi / .008 * t)
      >>>  x += 0.2 * np.cos(2 * np.pi / 0.6 * t + 0.1)
      >>>  x += 0.2 * np.cos(2 * np.pi / 1.6 * t + .11)
      >>>  x += 1 * np.cos(2 * np.pi / 10 * t + 11)

    Average the data over 6 hour period

      >>>  nx, err = average_err(x, window=12)
      >>>  plt.plot(t, x, 'k', t, nx, 'r', label=['Raw', 'Average'])
      >>>  plt.figure()
      >>>  plt.plot(t, err, 'g', label='Variance')

    """
    x = np.atleast_1d(x).flatten()
    nx = np.ma.masked_all(x.shape)
    err = np.ma.masked_all(x.shape)
    filt = np.ones(window) / window
    padlen = window * 4

    # Go over all contiguous regions
    regions = seapy.contiguous(x)
    for r in regions:
        if ((r.stop - r.start) >= padlen):
            nx[r] = scipy.signal.filtfilt(
                filt, [1.0], x[r], padlen=padlen, axis=0)
            err[r] = scipy.signal.filtfilt(
                filt, [1.0], (nx[r] - x[r])**2, padlen=padlen, axis=0)

    return nx, err


def bandpass(x, dt, low_cutoff=None, hi_cutoff=None, order=7):
    """
    Perform a bandpass filter at the cutoff period (same units as the
    time-series step).

    Parameters
    ----------
    x  : ndarray,
        The time-series of data to bandpass filter.
    dt : float,
        The time-step between the values in x. Units must be consistent
        with the cutoff period.
    low_cutoff: float,
        The period at which the bandpass filter will apply the lowpass filter.
        Units are same as the time-step of the signal provided (e.g., if the
        data are provided every hour, then a cutoff=12 would be a 12 hour
        cutoff.) Everything that has a longer period will remain. If you only
        want a hi-pass filter, this value should be None.
    hi_cutoff: float,
        The period at which the bandpass filter will apply the high-pass filter.
        Units are same as the time-step of the signal provided (e.g., if the
        data are provided every hour, then a cutoff=12 would be a 12 hour
        cutoff.) Everything that has a shorter period will remain. If you only
        want a low-pass filter, this value should be None.
    order: int optional,
        The order of the filter to apply

    Returns
    -------
    x : ndarray
        The bandpass filtered time-series

    Examples
    --------
    Create data every 30 minutes for 3 days with time in days:

      >>>  t = np.linspace(0, 3.0, 2 * 24 * 3.0, endpoint=False)
      >>>  x = 0.1 * np.sin(2 * np.pi / .008 * t)
      >>>  x += 0.2 * np.cos(2 * np.pi / 0.6 * t + 0.1)
      >>>  x += 0.2 * np.cos(2 * np.pi / 1.6 * t + .11)
      >>>  x += 1 * np.cos(2 * np.pi / 10 * t + 11)

    Filter the data to low-pass everything longer than the 1 day period

      >>>  nx = bandpass(x, dt=0.5, low_cutoff=24 )
      >>>  plt.plot(t, x, 'k', t, nx, 'r', label=['Raw', 'Filter'])

    Filter the data to low-pass everything longer the 2 day period

      >>>  nx = bandpass(x, dt=0.5, low_cutoff=48 )
      >>>  plt.plot(t, x, 'k', t, nx, 'r', label=['Raw', 'Filter'])

    Filter the data to band-pass everything shorter the 2 day period
    and longer the 1 hour period

      >>>  nx = bandpass(x, dt=0.5, low_cutoff=48, hi_cutoff=1 )
      >>>  plt.plot(t, x, 'k', t, nx, 'r', label=['Raw', 'Filter'])
    """
    x = np.ma.array(np.atleast_1d(x).flatten(), copy=False)
    nx = np.ma.masked_all(x.shape)

    if low_cutoff and hi_cutoff:
        freq = 2.0 * dt / np.array([hi_cutoff, low_cutoff])
        btype = 'bandpass'
    elif low_cutoff:
        freq = 2.0 * dt / low_cutoff
        btype = 'lowpass'
    elif hi_cutoff:
        freq = 2.0 * dt / hi_cutoff
        btype = 'highpass'
    else:
        raise AttributeError("You must specify either low or hi cutoff.")
    b, a = scipy.signal.butter(order, freq, btype=btype)
    padlen = max(len(a), len(b))

    # Go over all contiguous regions
    regions = seapy.contiguous(x)
    for r in regions:
        if ((r.stop - r.start) >= 5 * padlen):
            nx[r] = scipy.signal.filtfilt(
                b, a, x[r], padlen=5 * padlen, axis=0)
    return nx
