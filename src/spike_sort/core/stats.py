#import _diptst
import numpy as np
from scipy import stats as st

def mean_r(data):
    return np.median(data, 1)[:, None]

def std_r(data):
    np.median(np.abs(data - mean_r(data)), 1)[:, None]/0.6745

#def dip1d(data):
    #sdata = np.sort(data)
    #return _diptst.diptst1(sdata)[0]

def ks1d(data):
    mr = np.median(data)
    stdr = np.median(np.abs(data - mr))/0.6745

    return st.kstest(data, st.norm(loc=mr, scale=stdr).cdf)[0]

#def dip(data):
    #return np.array([dip1d(x) for x in data])[:, None]

def ks(data):
    return np.array([ks1d(x) for x in data])[:, None]

def var(data):
    return std_r(data).T
