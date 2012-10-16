import _diptst
import numpy as np
from scipy import stats as st

def multidimensional(func1d):
    """apply 1d function along specified axis
    """
    def _decorated(data, axis = 0):
        if data.ndim <= 1:
            return func1d(data)
        else:
            return np.apply_along_axis(func1d, axis, data)
    _decorated.__doc__ = func1d.__doc__

    return _decorated

def unsqueeze(data, axis):
    """inserts new axis to data at position `axis`.
    
    This is very useful when one wants to do operations which support
    broadcasting without using np.newaxis every time.

    Parameters
    ----------
    data : array
        input array
    axis : int
        axis to be inserted

    Returns
    -------
    out : array
        array which has data.ndim+1 dimensions. Additional dimension
        has length 1
    
    Example
    -------
    >>> data = [[1,2,3], [4,5,6]]
    >>> m = np.mean(data, 1)
    >>> data -= unsqueeze(m, 1)
    >>> data.mean(1)
    array([ 0.,  0.])
    
    """
    shape = data.shape
    shape = np.insert(shape, axis, 1)
    return np.reshape(data, shape)

def std_r(data, axis = 0):
    """Computes robust estimate of standard deviation
    (Quiroga et al, 2004)

    Parameters
    ----------
    data : array
        (n_vecs, n_obs) input data array

    Returns
    -------
    data : array
        (n_vecs) row-wise std values
    """
    median = unsqueeze(np.median(data, axis), axis)
    std_r = np.median(np.abs(data - median), axis)/0.6745
    return std_r

@multidimensional
def dip(data):
    """Computes DIP statistic (Hartigan & Hartigan 1985)

    Parameters
    ----------
    data : array
        input data array
    axis : int
        axis along which to compute dip

    Returns
    -------
    data : float or array
        DIP statistic. If the input data is flat, returns float
    """
    sdata = np.sort(data)
    return _diptst.diptst1(sdata)[0]

@multidimensional
def ks(data):
    """Computes Kolmogorov-Smirnov statistic (Lilliefors modification)

    Parameters
    ----------
    data : array
        (n_vecs) input data array
    axis : int
        axis along which to compute ks

    Returns
    -------
    data : array or float
        KS statistic. If the input data is flat, returns float
    """
    mr = np.median(data)
    stdr = np.median(np.abs(data - mr))/0.6745

    # avoid zero-variance
    if stdr == 0:
        return 0.

    return st.kstest(data, st.norm(loc=mr, scale=stdr).cdf)[0]

std = np.std








