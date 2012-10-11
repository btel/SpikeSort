#from spike_sort.stats._diptst import diptst1
import _diptst
import numpy as np
from scipy import stats as st

def mean_r(data):
    return np.median(data, 1)

def std_r(data):
    return np.median(np.abs(data - mean_r(data)[:, None]), 1)/0.6745

def dip1d(data):
    sdata = np.sort(data)
    return _diptst.diptst1(sdata)[0]

def ks1d(data):
    mr = np.median(data)
    stdr = np.median(np.abs(data - mr))/0.6745

    # avoid zero-variance
    if stdr == 0:
        return 0.

    return st.kstest(data, st.norm(loc=mr, scale=stdr).cdf)[0]

def dip(data):
    return np.array([dip1d(x) for x in data])

def ks(data):
    return np.array([ks1d(x) for x in data])

def std(data):
    return np.std(data, 1)
