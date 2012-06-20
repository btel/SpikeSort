import tables
import tempfile
import numpy as np
from scipy import signal

class ZeroPhaseFilter(object):
    """IIR Filter with zero phase delay"""
    def __init__(self, ftype, fband, tw=200., stop=20):
        self.gstop = stop
        self.gpass = 1
        self.fband = fband
        self.tw = tw
        self.ftype = ftype
        self._coefs_cache = {}

    def _design_filter(self, FS):

        if not FS in self._coefs_cache:
            wp = np.array(self.fband)
            ws = wp + np.array([-self.tw, self.tw])
            wp, ws = wp * 2.0 / FS, ws * 2.0 / FS
            b, a = signal.iirdesign(wp=wp,
                                    ws=ws,
                                    gstop=self.gstop,
                                    gpass=self.gpass,
                                    ftype=self.ftype)
            self._coefs_cache[FS] = (b, a)
        else:
            b, a = self._coefs_cache[FS]
        return b, a

    def __call__(self, x, FS):
        b, a = self._design_filter(FS)
        return signal.filtfilt(b, a, x)


class FilterFir(object):
    """FIR filter with zero phase delay

    Attributes
    ----------
    f_pass : float
             normalised low-cutoff frequency

    f_stop : float
             normalised high-cutoff frequency

    order : int
            filter order

    """
    def __init__(self, f_pass, f_stop, order):
        self._coefs_cache = {}
        self.fp = f_pass
        self.fs = f_stop
        self.order = order

    def _design_filter(self, FS):
        if not FS in self._coefs_cache:
            bands = [0, min(self.fs, self.fp), max(self.fs, self.fp),  FS / 2]
            gains = [int(self.fp < self.fs), int(self.fp > self.fs)]
            b, a = signal.remez(self.order, bands, gains, Hz=FS), [1]
            self._coefs_cache[FS] = (b, a)
        else:
            b, a = self._coefs_cache[FS]
        return b, a

    def __call__(self, x, FS):
        b, a = self._design_filter(FS)
        return signal.filtfilt(b, a, x)


class Filter(object):
    def __init__(self, fpass, fstop, gpass=1, gstop=10, ftype='butter'):
        self.ftype = ftype
        self.fp = np.asarray(fpass)
        self.fs = np.asarray(fstop)
        self._coefs_cache = {}
        self.gstop = gstop
        self.gpass = gpass

    def _design_filter(self, FS):
        if not FS in self._coefs_cache:
            wp, ws = self.fp * 2 / FS, self.fs * 2 / FS
            b, a = signal.iirdesign(wp=wp,
                                    ws=ws,
                                    gstop=self.gstop,
                                    gpass=self.gpass,
                                    ftype=self.ftype)
            self._coefs_cache[FS] = (b, a)
        else:
            b, a = self._coefs_cache[FS]
        return b, a

    def __call__(self, x, FS):
        b, a = self._design_filter(FS)
        return signal.filtfilt(b, a, x)


def filter_proxy(spikes, filter_obj, chunksize=1E6):
    """Proxy object to read filtered data

    Parameters
    ----------
    spikes : dict
        unfiltered raw recording
    filter_object : object
        Filter to filter the data
    chunksize : int
        size of segments in which data is filtered

    Returns
    -------
    sp_dict : dict
        filtered recordings
    """
    data = spikes['data']
    sp_dict = spikes.copy()

    if filter_obj is None:
        return spikes

    tmp_file = tempfile.NamedTemporaryFile(mode='w')
    filename = tmp_file.name
    atom = tables.Atom.from_dtype(np.dtype('float64'))
    shape = data.shape
    h5f = tables.openFile(filename, 'w')
    carray = h5f.createCArray('/', 'test', atom, shape)

    chunksize = int(chunksize)
    n_chunks = int(np.ceil(shape[1] * 1.0 / chunksize))
    for i in range(shape[0]):
        for j in range(n_chunks):
            stop = int(np.min(((j + 1) * chunksize, shape[1])))
            carray[i, j * chunksize:stop] = filter_obj(
                data[i, j * chunksize:stop], sp_dict['FS'])
    sp_dict['data'] = carray
    return sp_dict

def fltLinearIIR(signal, fpass, fstop, gpass=1, gstop=10, ftype='butter'): 
    """An IIR acausal linear filter. Works through
    spike_sort.core.filters.filter_proxy method
    
    Parameters
    ----------
    signal : dict
        input [raw] signal
    fpass, fstop : float
        Passband and stopband edge frequencies [Hz]
        For more details see scipy.signal.iirdesign
    gpass : float
        The maximum loss in the passband (dB).
    gstop : float
        The minimum attenuation in the stopband (dB).
    ftype : str, optional
        The type of IIR filter to design:

            - elliptic    : 'ellip'
            - Butterworth : 'butter',
            - Chebyshev I : 'cheby1',
            - Chebyshev II: 'cheby2',
            - Bessel :      'bessel'
    """
    filt = Filter(fpass, fstop, gpass, gstop, ftype)
    return filter_proxy(signal, filt)
