#!/usr/bin/env python
"""
  lib.py

  State Estimation and Analysis for PYthon

  Library of utilities for general seapy module, imported into the namespace
  when importing the seapy module

  Written by Brian Powell on 10/18/13
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import numpy as np
from scipy import ndimage
import os
import re
import datetime
import itertools


secs2day = 1.0 / 86400.0
default_epoch = datetime.datetime(2000, 1, 1)
_default_timeref = "days since " + default_epoch.strftime("%Y-%m-%d %H:%M:%S")


def adddim(fld, size=1):
    """
    replicate a field and add a new first dimension with given size

    Parameters
    ----------
    fld : array_like
        Input field.
    size : int, optional
        Size of additional first dimension

    Returns
    -------
    fld : array

    Examples
    --------
    >>> a=np.array([4, 5, 6, 7])
    >>> a.shape
    (4,)
    >>> b = seapy.adddim(a, 2)
    >>> b.shape
    (2, 4)
    >>> b
    array([[4, 5, 6, 7],
           [4, 5, 6, 7]])

    """
    fld = np.atleast_1d(fld)
    s = np.ones(fld.ndim + 1).astype(int)
    s[0] = int(size)
    return np.tile(fld, s)


def fill(x, max_gap=None, kind='linear'):
    """
    Fill missing data from a 1-D vector. When data are missing from a
    vector, this method will interpolate to fill gaps that are less than
    the specified max (or ignored).

    Parameters
    ----------
    x : array
      The array to be filled. It will be cast as a masked array for
      invalid values. If already a masked array, then that mask will
      persist.
    max_gap : int, optional
      The maximum number of continuous values to interpolate (e.g.,
      if this value is 10 and there are 12 continuous missing values,
      they will be left unfilled). Default is to fill everything.
    kind : str, optional
      The kind of interpolant to use (see scipy.interpolate.interp1d).
      Default is 'linear'

    Returns
    -------
    x : array
      The filled array
    """
    from scipy.interpolate import interp1d
    x = np.ma.masked_invalid(np.atleast_1d(x).flatten(), copy=False)
    # If no gaps or empty data, do nothing
    if not np.any(x.mask) or len(x.compressed()) < 3:
        return x
    f = interp1d(x.nonzero()[0], x.compressed())
    nx = x.copy()
    if max_gap is not None:
        regions = contiguous(x)
        for r in regions:
            print(f"{r}")
            if ((r.stop - r.start) <= max_gap) and \
                    (r.stop < f.x.max()) and (r.start > f.x.min()):
                nx[r] = f(np.arange(r.start, r.stop))
                print(f"fill {r}: {nx[r]}")
    else:
        bad = np.nonzero(x.mask)[0]
        bad = np.delete(bad, np.nonzero(
            np.logical_or(bad <= f.x.min(), bad >= f.x.max())))
        nx[bad] = f(bad)
    return nx


def contiguous(x):
    """
    Find the indices that provide contiguous regions of a numpy.masked_array.
    This will find all regions of valid data. NOTE: this casts as 1-D.

    Parameters
    ----------
    x : np.array or np.ma.array
      The data to find the contiguous regions

    Returns
    -------
    idx : array of slices
      Array of slices for each contiguous region

    Examples
    --------
    >>> a = np.array([4, 3, 2, np.nan, 6, 7, 2])
    >>> r = contiguous(a)
    [slice(0, 2, None), slice(4, 6, None)]

    If no contiguous regions are available, an empty array is returned.

    """
    x = np.ma.masked_invalid(np.atleast_1d(x).flatten(), copy=False)
    idx = x.nonzero()[0]
    try:
        d = np.diff(idx) - 1
        dl = idx[np.nonzero(d)[0] + 1]
        d = d[np.nonzero(d)]
        return np.array([np.s_[r[0]:r[1]] for r in
                         zip(np.hstack((idx.min(), dl)),
                             np.hstack((dl - d - 1, idx.max())))])
    except:
        return []


def chunker(seq, size):
    """
    Iterate over an iterable in 'chunks' of a given size

    Parameters
    ----------
    seq : iterable,
        The sequence to iterate over
    size : int,
        The number of items to be returned in each 'chunk'

    Returns
    -------
    chunk : seq,
        The items of the chunk to be iterated

    Examples
    --------
    >>> x = [0,3,4,7,9,10,12,14]
    >>> for i in chunker(x, 3):
    >>>     print(i)
    [0, 3, 4]
    [7, 9, 10]
    [12, 14]

    """
    return (seq[pos:pos + size] for pos in range(0, len(seq), max(1, size)))


def smooth(data, ksize=3, kernel=None, copy=True):
    """
    Smooth the data field using a specified convolution kernel
    or a default averaging kernel.

    Parameters
    ----------
    data : masked array_like
        Input field.
    ksize : int, optional
        Size of square kernel
    kernel : ndarray, optional
        Define a convolution kernel. Default is averaging
    copy : bool, optional
        If true, a copy of input array is made

    Returns
    -------
    fld : masked array
    """
    fld = np.ma.array(data, copy=copy)
    mask = np.ma.getmaskarray(fld).copy()

    # Make sure ksize is odd
    ksize = int(ksize + 1) if int(ksize) % 2 == 0 else int(ksize)
    if fld.ndim > 3 or fld.ndim < 2:
        raise AttributeError("Can only convolve 2- or 3-D fields")
    if ksize < 3:
        raise ValueError("ksize must be greater than or equal to 3")

    if kernel is None:
        kernel = np.ones((ksize, ksize)) / (ksize * ksize)
    else:
        ksize = kernel.shape[0]

    # First, convole over any masked values
    fld = convolve_mask(fld, ksize=ksize, copy=False)

    # Next, perform the convolution
    if fld.ndim == 2:
        fld = ndimage.convolve(fld.data, kernel,
                               mode="reflect", cval=0.0)
    else:
        kernel = kernel[:, :, np.newaxis]
        fld = np.transpose(ndimage.convolve(
            fld.filled(0).transpose(1, 2, 0), kernel,
            mode="reflect", cval=0.0), (2, 0, 1))

    # Apply the initial mask
    return np.ma.array(fld, mask=mask)


def convolve(data, ksize=3, kernel=None, copy=True, only_mask=False):
    """
    Convolve the kernel across the data to smooth or highlight
    the field across the masked region.

    Parameters
    ----------
    data : masked array_like
        Input field.
    ksize : int, optional
        Size of square kernel
    kernel : ndarray, optional
        Define a convolution kernel. Default is averaging
    copy : bool, optional
        If true, a copy of input array is made
    only_mask : bool, optional
        If true, only consider the smoothing over the masked
        region

    Returns
    -------
    fld : masked array
    """
    fld = np.ma.array(data, copy=copy)
    if not copy:
        fld._sharedmask = False

    # Make sure ksize is odd
    ksize = int(ksize + 1) if int(ksize) % 2 == 0 else int(ksize)
    if fld.ndim > 3 or fld.ndim < 2:
        raise AttributeError("Can only convolve 2- or 3-D fields")
    if ksize < 3:
        raise ValueError("ksize must be greater than or equal to 3")

    if kernel is None:
        center = np.round(ksize / 2).astype(int)
        kernel = np.ones([ksize, ksize])
        kernel[center, center] = 0.0
    else:
        ksize = kernel.shape[0]

    # Convolve the mask
    msk = np.ma.getmaskarray(fld)
    if fld.ndim == 2:
        count = ndimage.convolve((~msk).view(np.int8), kernel,
                                 mode="constant", cval=0.0)
        nfld = ndimage.convolve(fld.data * (~msk).view(np.int8), kernel,
                                mode="constant", cval=0.0)
    else:
        kernel = kernel[:, :, np.newaxis]
        count = np.transpose(ndimage.convolve(
            (~msk).view(np.int8).transpose(1, 2, 0), kernel,
            mode="constant", cval=0.0), (2, 0, 1))
        nfld = np.transpose(ndimage.convolve(
            (fld.data * (~msk).view(np.int8)).transpose(1, 2, 0), kernel,
            mode="constant", cval=0.0), (2, 0, 1))

    if only_mask:
        lst = np.nonzero(np.logical_and(msk, count > 0))
        fld[lst] = np.ma.nomask
        fld[lst] = nfld[lst] / count[lst]
    else:
        lst = np.nonzero(~msk)
        fld[lst] = nfld[lst] / count[lst]
    return fld


def convolve_mask(data, ksize=3, kernel=None, copy=True):
    """
    Convolve data over the missing regions of a mask

    Parameters
    ----------
    data : masked array_like
        Input field.
    ksize : int, optional
        Size of square kernel
    kernel : ndarray, optional
        Define a convolution kernel. Default is averaging
    copy : bool, optional
        If true, a copy of input array is made

    Returns
    -------
    fld : masked array
    """
    return convolve(data, ksize, kernel, copy, True)


def matlab2date(daynum):
    """
    Given a day number from matlab, convert into a datetime

    Parameters
    ----------
    daynum: float
      Scalar or array of matlab day numbers

    Returns
    -------
    datetime : list
    """
    daynum = np.atleast_1d(daynum)
    return [datetime.datetime.fromordinal(d.astype(np.int)) +
            datetime.timedelta(days=(d % 1 - 366)) for d in daynum]


def date2day(date=default_epoch, epoch=default_epoch):
    """
    Compute the fractional number of days elapsed since the epoch to the date
    given.

    Parameters
    ----------
    date : datetime
        Input date
    epoch : datetime
        Date of epoch

    Returns
    -------
    numdays : list
    """
    date = np.atleast_1d(date)
    return [(t - epoch).total_seconds() * secs2day for t in date]


def day2date(day=0, epoch=default_epoch):
    """
    Return a datetime object from the number of days since the epoch

    Parameters
    ----------
    day : scalar
        Input day number
    epoch : datetime
        Date of epoch

    Returns
    -------
    date : list of datetime(s)
    """
    day = np.atleast_1d(day)
    return [epoch + datetime.timedelta(days=float(t)) for t in day]


def matlab2date(daynum=0):
    """
    Return a datetime object from a Matlab datenum value

    Parameters
    ----------
    daynum : scalar
        Input Matlab day number

    Returns
    -------
    date : list of datetime(s)
    """
    daynum = np.atleast_1d(daynum)

    return np.array([datetime.datetime.fromordinal(int(d)) +
                     datetime.timedelta(days=d % 1) -
                     datetime.timedelta(days=366) for d in daynum])


def _distq(lon1, lat1, lon2, lat2):
    """
    Compute the geodesic distance between lat/lon points. This code is
    taken from the dist.f routine and the Matlab version distg.m passed
    around WHOI and APL. This was stripped down to use the WGS84 ellipsoid.

    Parameters
    ----------
    lon1 : array_like or scalar
        Input array of source longitude(s)
    lat1 : array_like or scalar
        Input array of source latitude(s)
    lon2 : array_like or scalar
        Input array of destination longitude(s)
    lat2 : array_like or scalar
        Input array of destination latitude(s)

    Returns
    -------
    distance : array or scalar of distance in meters
    angle: array or scalar of angle in radians

    """
    lon1 = np.asanyarray(np.radians(lon1))
    lat1 = np.asanyarray(np.radians(lat1))
    lon2 = np.asanyarray(np.radians(lon2))
    lat2 = np.asanyarray(np.radians(lat2))

    # # If one of the points is a singleton and the other is an
    # array, make them the same size
    if lon1.size == 1 and lon2.size > 1:
        lon1 = lon1.repeat(lon2.size)
        lat1 = lat1.repeat(lat2.size)
    if lon2.size == 1 and lon1.size > 1:
        lon2 = lon2.repeat(lon1.size)
        lat2 = lat2.repeat(lat1.size)

    # Set the WGS84 parameters
    A = 6378137.
    E = 0.081819191
    B = np.sqrt(A * A - (A * E)**2)
    EPS = E * E / (1.0 - E * E)

    # Move any latitudes off of the equator
    lat1[lat1 == 0] = np.finfo(float).eps
    lat2[lat2 == 0] = -np.finfo(float).eps

    # COMPUTE THE RADIUS OF CURVATURE IN THE PRIME VERTICAL FOR EACH POINT
    xnu1 = A / np.sqrt(1.0 - (E * np.sin(lat1))**2)
    xnu2 = A / np.sqrt(1.0 - (E * np.sin(lat2))**2)

    TPSI2 = (1.0 - E * E) * np.tan(lat2) + E * E * xnu1 * np.sin(lat1) / \
        (xnu2 * np.cos(lat2))
    PSI2 = np.arctan(TPSI2)

    DPHI2 = lat2 - PSI2
    DLAM = (lon2 - lon1) + np.finfo(float).eps
    CTA12 = np.sin(DLAM) / (np.cos(lat1) * TPSI2 - np.sin(lat1) * np.cos(DLAM))
    A12 = np.arctan(CTA12)
    CTA21P = np.sin(DLAM) / (np.sin(PSI2) * np.cos(DLAM) -
                             np.cos(PSI2) * np.tan(lat1))
    A21P = np.arctan(CTA21P)

    # C    GET THE QUADRANT RIGHT
    DLAM2 = (np.abs(DLAM) < np.pi).astype(int) * DLAM + \
        (DLAM >= np.pi).astype(int) * (-2 * np.pi + DLAM) + \
        (DLAM <= -np.pi).astype(int) * (2 * np.pi + DLAM)
    A12 = A12 + (A12 < -np.pi).astype(int) * 2 * np.pi - \
        (A12 >= np.pi).astype(int) * 2 * np.pi
    A12 = A12 + np.pi * np.sign(-A12) * \
        (np.sign(A12).astype(int) != np.sign(DLAM2))
    A21P = A21P + (A21P < -np.pi).astype(int) * 2 * np.pi - \
        (A21P >= np.pi).astype(int) * 2 * np.pi
    A21P = A21P + np.pi * np.sign(-A21P) * \
        (np.sign(A21P).astype(int) != np.sign(-DLAM2))

    SSIG = np.sin(DLAM) * np.cos(PSI2) / np.sin(A12)

    dd1 = np.array([np.cos(lon1) * np.cos(lat1),
                    np.sin(lon1) * np.cos(lat1), np.sin(lat1)])
    dd2 = np.array([np.cos(lon2) * np.cos(lat2),
                    np.sin(lon2) * np.cos(lat2), np.sin(lat2)])
    dd2 = np.sum((dd2 - dd1)**2, axis=0)
    bigbrnch = (dd2 > 2).astype(int)

    SIG = np.arcsin(SSIG) * (bigbrnch == 0).astype(int) + \
        (np.pi - np.arcsin(SSIG)) * bigbrnch

    SSIGC = -np.sin(DLAM) * np.cos(lat1) / np.sin(A21P)
    SIGC = np.arcsin(SSIGC)
    A21 = A21P - DPHI2 * np.sin(A21P) * np.tan(SIG / 2.0)

    # C   COMPUTE RANGE
    G2 = EPS * (np.sin(lat1))**2
    G = np.sqrt(G2)
    H2 = EPS * (np.cos(lat1) * np.cos(A12))**2
    H = np.sqrt(H2)
    SIG2 = SIG * SIG
    TERM1 = -H2 * (1.0 - H2) / 6.0
    TERM2 = G * H * (1.0 - 2.0 * H2) / 8.0
    TERM3 = (H2 * (4.0 - 7.0 * H2) - 3.0 * G2 * (1.0 - 7.0 * H2)) / 120.0
    TERM4 = -G * H / 48.0
    rng = xnu1 * SIG * (1.0 + SIG2 * (TERM1 + SIG * TERM2 + SIG2 * TERM3 +
                                      SIG2 * SIG * TERM4))

    return rng, A12


def earth_distance(lon1, lat1, lon2, lat2):
    """
    Compute the geodesic distance between lat/lon points.

    Parameters
    ----------
    lon1 : array_like or scalar
        Input array of source longitude(s)
    lat1 : array_like or scalar
        Input array of source latitude(s)
    lon2 : array_like or scalar
        Input array of destination longitude(s)
    lat2 : array_like or scalar
        Input array of destination latitude(s)

    Returns
    -------
    distance : array or scalar of distance in meters

    """
    rng, _ = _distq(lon1, lat1, lon2, lat2)
    return rng


def earth_angle(lon1, lat1, lon2, lat2):
    """
    Compute the angle between lat/lon points. NOTE: The bearing angle
    is computed, but then converted to geometric (counter-clockwise)
    angle to be returned.

    Parameters
    ----------
    lon1 : array_like or scalar
        Input array of source longitude(s)
    lat1 : array_like or scalar
        Input array of source latitude(s)
    lon2 : array_like or scalar
        Input array of destination longitude(s)
    lat2 : array_like or scalar
        Input array of destination latitude(s)

    Returns
    -------
    angle : array or scalar of bearing in radians

    """
    _, angle = _distq(lon1, lat1, lon2, lat2)
    return (np.pi / 2.0 - angle)


def flatten(l, ltypes=(list, tuple, set)):
    """
    Flatten a list or tuple that contains additional lists or tuples. Like
    the numpy flatten, but for python types.

    Parameters
    ----------
    l: tuple or list,
        The data that is to be flattened
    ltypes: tuple,
        Data types to attempt to flatten

    Returns
    -------
    list

    See Also
    --------
    numpy.flatten()

    Notes
    -----
    This code was taken from:
    <http://rightfootin.blogspot.com.au/2006/09/more-on-python-flatten.html>

    Examples
    --------
    >>> a=[[1,3,4,1], ('test', 'this'), [5,2]]
    >>> flatten(a)
    [1, 3, 4, 1, 'test', 'this', 5, 2]
    """
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)


def list_files(path=".", regex=None, full_path=True):
    """
    list all sorted file names in the given path that conform to the regular
    expression pattern. This is not a generator function because it sorts
    the files in alphabetic/numeric order.

    Parameters
    ----------
    path : string
        Search for the given matches
    regex : string, optional
        Input regular expression string to filter filenames
    full_path : bool, optional
        If True, return the full path for each found object. If false,
        return only the filename

    Returns
    -------
    files : array

    Examples
    --------
    >>> files = seapy.list_files('/path/to/dir/test_.*txt')
    >>> print(files)
    ['/path/to/dir/test_001.txt', '/path/to/dir/test_002.txt']

    NOTE: this is equivalent for separating:
    >>> files = seapy.list_files('/path/to/dir', 'test_.*txt')
    """
    # If only one parameter is given, parse into its components
    if regex is None:
        regex = os.path.basename(path)
        path = os.path.dirname(path)
    if not path:
        path = './'
    elif path[-1] != '/':
        path += '/'
    files = []
    prog = re.compile(regex)
    for file in os.listdir(path):
        if prog.search(file) is not None:
            if full_path:
                files.append(path + file)
            else:
                files.append(file)
    files.sort()
    return files


def netcdf(file, aggdim=None):
    """
    Wrapper around netCDF4 to open a file as either a Dataset or an
    MFDataset.

    Parameters
    ----------
    file : string or list,
        Filename(s) to open. If the string has wildcards or is a list,
        this attempts to open an MFDataset
    aggdim : string,
        Name of dimension to concatenate along if loading a set of files.
        A value of None (default) uses the unlimited dimension.

    Returns
    -------
    netCDF4 Dataset or MFDataset
    """
    import netCDF4
    try:
        nc = netCDF4.Dataset(file)
    except (OSError, RuntimeError):
        try:
            nc = netCDF4.MFDataset(file, aggdim=aggdim)
        except IndexError:
            raise FileNotFoundError("{:s} cannot be found.".format(file))
    return nc


def primes(number):
    """
    Return a list of primes less than or equal to a given value.

    Parameters
    ----------
    number : int
        Find prime values up to this value

    Returns
    -------
    primes : ndarray

    Notes
    -----
    This code was taken from "Cooking with Python, Part 2" by Martelli, et al.

    <http://archive.oreilly.com/pub/a/python/excerpt/pythonckbk_chap1/index1.html?page=last>

    """
    def __erat2():
        D = {}
        yield 2
        for q in itertools.islice(itertools.count(3), 0, None, 2):
            p = D.pop(q, None)
            if p is None:
                D[q * q] = q
                yield q
            else:
                x = p + q
                while x in D or not (x & 1):
                    x += p
                D[x] = p

    return np.array(list(itertools.takewhile(lambda p: p < number, __erat2())))


def rotate(u, v, angle):
    """
    Rotate a vector field by the given angle

    Parameters
    ----------
    u : array like
        Input u component
    v : array like
        Input v component
    angle : array like
        Input angle of rotation in radians

    Returns
    -------
    rotated_u, rotated_v : array
    """
    u = np.asanyarray(u)
    v = np.asanyarray(v)
    angle = np.asanyarray(angle)
    sa = np.sin(angle)
    ca = np.cos(angle)

    return u * ca - v * sa, u * sa + v * ca


def today2day(epoch=default_epoch):
    """
    Return the day number of today (UTC time) since the epoch.

    Parameters
    ----------
    epoch : datetime
        Date of epoch

    Returns
    -------
    numdays : scalar
    """
    return date2day(datetime.datetime.utcnow(), epoch)


def unique_rows(x):
    """
    Convert rows into godelnumbers and find the rows that are unique using
    np.unique

    Parameters
    ----------
    x : ndarray or tuple,
        array of elements to find unique value. If columns are greater
        than 1, then the columns are combined into a single Godel number.
        If a tuple of arrays are passed, they are combined.

    Returns
    -------
    idx : ndarray,
        Indices of the unique values

    Examples
    --------
    >>> a = np.array([3, 3, 5, 5, 6])
    >>> b = np.array([2, 3, 3, 3, 3])
    >>> idx = unique_rows((a, b))
    >>> idx
    array([0, 1, 2, 4])
    """
    if isinstance(x, tuple):
        x = np.vstack(x).T
    else:
        x = np.atleast_1d(x)
    vals, idx = np.unique(godelnumber(x), return_index=True)

    return idx


def vecfind(a, b, tolerance=None):
    """
    Find all occurences of b in a within the given tolerance and return
    the sorted indices of a and b that yield the corresponding values.
    The indices are of equal length, such that

    Written by Eric Firing, University of Hawaii.

    Parameters
    ----------
    a : array
        Input vector
    b : array
        Input vector
    tolerance : same type as stored values of a and b, optional
        Input tolerance for how close a is to b. If not specified,
        then elements of a and b must be equal.

    Returns
    -------
    index_a, index_b : arrays of indices for each vector where values are equal,
            such that a[index_a] == b[index_b]

    Examples
    --------
    >>> a = np.array([3,4,1,8,9])
    >>> b = np.array([4,7,1])
    >>> ia, ib = vecfind(a, b)

    By definition,

    >>> len(ia) == len(ib)
    True
    >>> a[ia] == b[ib]
    True

    """
    a = np.asanyarray(a).flatten()
    b = np.asanyarray(b).flatten()

    # if no tolerance, compute a zero distance  the proper type
    if tolerance is None:
        tolerance = a[0] - a[0]

    _, uniq_a = np.unique(a, return_index=True)
    _, uniq_b = np.unique(b, return_index=True)
    na = len(uniq_a)
    t = np.hstack((a[uniq_a], b[uniq_b]))
    is_a = np.zeros(t.shape, dtype=np.int8)
    is_a[:na] = 1
    isorted = np.argsort(t)
    tsorted = t[isorted]
    is_a_sorted = is_a[isorted]

    dt = np.diff(tsorted)
    mixed = np.abs(np.diff(is_a_sorted)) == 1
    ipair = np.nonzero((np.abs(dt) <= tolerance) & mixed)[0]

    # Now ipair should be the indices of the first elements
    # of consecutive pairs in tsorted for which the two items
    # are from different arrays, and differ by less than tolerance.
    # The problem is that they could be in either order.

    iswap = np.nonzero(is_a_sorted[ipair] == 0)[0]  # b is first, so swap

    temp = isorted[ipair[iswap] + 1]
    isorted[ipair[iswap] + 1] = isorted[ipair[iswap]]
    isorted[ipair[iswap]] = temp

    isorted_a = isorted[ipair]
    isorted_b = isorted[ipair + 1] - na

    return uniq_a[isorted_a], uniq_b[isorted_b]


def godelnumber(x):
    """
    Convert the columns of x into godel numbers. If x is MxN, return an Mx1
    vector. The Godel number is prime**x

    Parameters
    ----------
    x : ndarray,
        Values to convert into Godel number(s)

    Returns
    -------
    godel : ndarray
    """
    x = np.atleast_2d(x.astype(int))
    if x.ndim > 1:
        primevals = primes(x.shape[1] * 10)[:x.shape[1]].astype(float)
        return(np.prod(primevals**x, axis=1))
    else:
        return 2.0**x


pass
