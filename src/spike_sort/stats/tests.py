#from spike_sort.stats._diptst import diptst1
import _diptst
import numpy as np
from scipy import stats as st

def mean_r(data):
    """Computes robust estimate of mean (Quiroga et al, 2004) row-wise

    Parameters
    ----------
    data : array
        (n_vecs, n_obs) input data array

    Returns
    -------
    data : array
        (n_vecs) row-wise mean values
    """
    return np.median(data, 1)

def std_r(data):
    """Computes robust estimate of standard deviation
    (Quiroga et al, 2004) row-wise

    Parameters
    ----------
    data : array
        (n_vecs, n_obs) input data array

    Returns
    -------
    data : array
        (n_vecs) row-wise std values
    """
    return np.median(np.abs(data - mean_r(data)[:, None]), 1)/0.6745

def dip1d(data):
    """Computes DIP statistic (Hartigan & Hartigan 1985)

    Parameters
    ----------
    data : array
        (n_vecs) input data array

    Returns
    -------
    data : float
        DIP statistic
    """
    sdata = np.sort(data)
    return _diptst.diptst1(sdata)[0]

def ks1d(data):
    """Computes Kolmogorov-Smirnov statistic (Lilliefors modification)

    Parameters
    ----------
    data : array
        (n_vecs) input data array

    Returns
    -------
    data : float
        KS statistic
    """
    mr = np.median(data)
    stdr = np.median(np.abs(data - mr))/0.6745

    # avoid zero-variance
    if stdr == 0:
        return 0.

    return st.kstest(data, st.norm(loc=mr, scale=stdr).cdf)[0]

def dip(data):
    """Computes DIP statistic row-wise (see also dip1d)

    Parameters
    ----------
    data : array
        (n_vecs, n_obs) input data array

    Returns
    -------
    data : array
        DIP statistic row-wise
    """
    return np.array([dip1d(x) for x in data])

def ks(data):
    """Computes KS statistic row-wise (see also ks1d)

    Parameters
    ----------
    data : array
        (n_vecs, n_obs) input data array

    Returns
    -------
    data : array
        KS statistic row-wise
    """
    return np.array([ks1d(x) for x in data])

def std(data):
    """Computes standard deviation row-wise

    Parameters
    ----------
    data : array
        (n_vecs, n_obs) input data array

    Returns
    -------
    data : array
        std row-wise
    """
    return np.std(data, 1)
