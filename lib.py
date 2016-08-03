#!/usr/bin/env python
"""
  lib.py

  State Estimation and Analysis for PYthon

  Library of utilities for general seapy module, imported into the namespace
  when importing the seapy module

  Written by Brian Powell on 10/18/13
  Copyright (c)2016 University of Hawaii under the BSD-License.
"""


import numpy as np
from scipy import ndimage
import os
import re
import datetime
import itertools


secs2day = 1.0 / 86400.0
default_epoch = datetime.datetime(2000, 1, 1)
_default_timeref = "days since " + default_epoch.strftime("%Y-%m-%m %H:%M:%S")


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
    s = np.ones(fld.ndim + 1)
    s[0] = size
    return np.tile(fld, s)


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
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


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
        kernel = np.expand_dims(kernel, axis=3)
        fld = np.transpose(ndimage.convolve(
            fld.filled(0).transpose(1, 2, 0), kernel,
            mode="reflect", cval=0.0), (2, 0, 1))

    # Apply the initial mask
    return np.ma.array(fld, mask=mask)


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
    fld = np.ma.array(data, copy=copy)
    if not copy: fld._sharedmask = False
    
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
        kernel = np.expand_dims(kernel, axis=3)
        count = np.transpose(ndimage.convolve(
            (~msk).view(np.int8).transpose(1, 2, 0), kernel,
            mode="constant", cval=0.0), (2, 0, 1))
        nfld = np.transpose(ndimage.convolve(
            (fld.data * (~msk).view(np.int8)).transpose(1, 2, 0), kernel,
            mode="constant", cval=0.0), (2, 0, 1))

    lst = np.nonzero(np.logical_and(msk, count > 0))
    fld[lst] = np.ma.nomask
    fld[lst] = nfld[lst] / count[lst]
    return fld


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
    lon2 : array_like or scalar
        Input array of destination longitude(s)

    Returns
    -------
    angle : array or scalar of bearing in radians

    """
    lon1 = np.asanyarray(lon1)
    lon2 = np.asanyarray(lon2)
    lat1 = np.radians(np.asanyarray(lat1))
    lat2 = np.radians(np.asanyarray(lat2))

    dlon = np.radians(lon2 - lon1)
    angle = np.arctan2(np.sin(dlon) * np.cos(lat2),
                       np.cos(lat1) * np.sin(lat2) -
                       np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
    return (np.pi / 2.0 - angle)


def earth_distance(lon1, lat1, lon2, lat2):
    """
    Compute the distance between lat/lon points

    Parameters
    ----------
    lon1 : array_like or scalar
        Input array of source longitude(s)
    lat1 : array_like or scalar
        Input array of source latitude(s)
    lon2 : array_like or scalar
        Input array of destination longitude(s)
    lon2 : array_like or scalar
        Input array of destination longitude(s)

    Returns
    -------
    distance : array or scalar of distance in meters

    """
    epsilon = 0.99664718940443  # This is Sqrt(1-epsilon^2)
    radius = 6378137  # Radius in meters
    d2r = np.pi / 180.0

    lon1 = np.asanyarray(lon1)
    lat1 = np.asanyarray(lat1)
    lon2 = np.asanyarray(lon2)
    lat2 = np.asanyarray(lat2)

    # Using trig identities of tan(atan(b)), cos(atan(b)), sin(atan(b)) for
    # working with geocentric where lat_gc = atan(epsilon * tan(lat))
    tan_lat = epsilon * np.tan(d2r * lat1.astype(np.float64))
    cos_lat = 1.0 / np.sqrt(1.0 + tan_lat**2)
    sin_lat = tan_lat / np.sqrt(1.0 + tan_lat**2)
    tan_lat = epsilon * np.tan(d2r * lat2.astype(np.float64))
    cos_latj = 1.0 / np.sqrt(1.0 + tan_lat**2)
    sin_latj = tan_lat / np.sqrt(1.0 + tan_lat**2)

    return radius * np.sqrt(2.0 * (1.0 - cos_lat * cos_latj *
                                   np.cos(d2r * (lon1 - lon2)) - sin_lat * sin_latj))


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


def netcdf(file):
    """
    Wrapper around netCDF4 to open a file as either a Dataset or an
    MFDataset.

    Parameters
    ----------
    file : string or list,
        Filename(s) to open. If the string has wildcards or is a list,
        this attempts to open an MFDataset

    Returns
    -------
    netCDF4 Dataset or MFDataset
    """
    import netCDF4
    try:
        nc = netCDF4.Dataset(file)
    except RuntimeError:
        try:
            nc = netCDF4.MFDataset(file)
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
        Input angle of rotation

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


def vecfind(a, b, tolerance=0):
    """
    Find all occurences of b in a within the given tolerance and return
    the indices into a and b that correspond.

    Parameters
    ----------
    a : array
        Input vector
    b : array
        Input vector
    tolerance : float, optional
        Input tolerance for how close a==b

    Returns
    -------
    index_a, index_b : arrays of indices for each vector where values are equal,
            such that a[index_a] == b[index_b]

    """
    a = np.asanyarray(a)
    b = np.asanyarray(b)

    # If the tolerance is zero, we can use built-in methods
    if tolerance == 0:
        index_a = np.where(np.in1d(a, b))[0]
        index_b = np.where(np.in1d(b, a))[0]
    # Otherwise, slow method
    else:
        index_a = []
        index_b = []
        for i, c in enumerate(b):
            d = np.abs(a - c)
            x = np.nonzero(d < tolerance)[0]
            if np.any(x):
                index_a.append(np.where(d == np.min(d))[0])
                index_b.append(i)
        index_a = np.array(index_a).flatten()
        index_b = np.array(index_b)
    return index_a, index_b


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
