#!/usr/bin/env python
#coding=utf-8

import tempfile
from warnings import warn
import operator

import numpy as np
from scipy import interpolate, signal
import tables


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


def split_cells(spikes, idx, which='all'):
    """Return the spike features splitted into separate cells
    """

    if which == 'all':
        classes = np.unique(idx)
    else:
        classes = which

    data = spikes['data']
    time = spikes['time']
    spikes_dict = dict([(cl, {'data': data[:, idx == cl, :], 'time': time})
                        for cl in classes])
    return spikes_dict


def remove_spikes(spt_dict, remove_dict, tolerance):
    """Remove spikes with given spike times from the spike time
    structure """
    spt_data = spt_dict['data']
    spt_remove = remove_dict['data']

    mn, mx = tolerance

    for t in spt_remove:
        spt_data = spt_data[(spt_data > (t + mx)) | (spt_data < (t + mn))]

    spt_ret = spt_dict.copy()
    spt_ret['data'] = spt_data
    return spt_ret


def detect_spikes(spike_data, thresh='auto', edge="rising",
                  contact=0, filter=None):
    r"""Detects spikes in extracellular data using amplitude thresholding.

    Parameters
    ----------
    spike_data : dict
        extracellular waveforms
    thresh : float or 'auto'
        threshold for detection. if thresh is 'auto' it will be
        estimated from the data.
    edge : {'rising', 'falling'}
        which edge to trigger on
    contact : int, optional
        index of tetrode contact to use for detection, defaults to
        first contact
    filter : object, optional
        filter used for spike detection; defaults to no filtering

    Returns
    -------
    spt_dict : dict
        dictionary with 'data' key which contains detected threshold
        crossing in miliseconds

    """

    sp_data = spike_data['data'][contact, :]
    n_contacts = spike_data['n_contacts']

    if filter is not None:
        sp_data = filter(sp_data, spike_data['FS'])
    #if n_contacts>1:
    #    sp_data = sp_data[:,contact]

    FS = spike_data['FS']

    if isinstance(thresh, basestring):
        if thresh == 'auto':
            thresh_frac = 8
        else:
            thresh_frac = float(thresh)

        thresh = thresh_frac * np.sqrt(float(np.var(sp_data[:int(10 * FS)])))

        if edge == 'falling' or edge == "min":
            thresh = -thresh

    edges = ('rising', 'max', 'falling', 'min')
    if edge not in edges:
        raise TypeError("'edge' parameter must be 'rising' or 'falling'")
    
    op1, op2 = operator.lt, operator.gt

    if edge in edges[2:]:
        op1, op2 = op2, op1

    i, = np.where(op1(sp_data[:-1], thresh) & op2(sp_data[1:], thresh))      

    spt = i * 1000.0 / FS
    spt_dict = {'data': spt, 'thresh': thresh, 'contact': contact}
    return spt_dict


def filter_spt(spike_data, spt_dict, sp_win):
    spt = spt_dict['data']
    sp_data = spike_data['data']
    FS = spike_data['FS']

    try:
        n_pts = sp_data.shape[1]
    except IndexError:
        n_pts = len(sp_data)
    max_time = (n_pts) * 1000.0 / FS

    t_min = np.max((-sp_win[0], 0))
    t_max = np.min((max_time, max_time - sp_win[1]))
    idx, = np.nonzero((spt >= t_min) & (spt <= t_max))
    return idx


def extract_spikes(spike_data, spt_dict, sp_win, resample=1,
                   contacts='all'):
    """Extract spikes from recording.

    Parameters
    ----------
    spike_data : dict
       extracellular data (see :ref:`raw_recording`)
    spt : dict
       spike times structure (see :ref:`spike_times`)
    sp_win : list of int
       temporal extent of the wave shape

    Returns
    -------
    wavedict : dict
       spike waveforms structure (see :ref:`spike_wave`)


    """
    sp_data = spike_data['data']
    n_contacts = spike_data['n_contacts']

    if contacts == "all":
        contacts = np.arange(n_contacts)
    elif type(contacts) is int:
        contacts = np.array([contacts])
    else:
        contacts = np.asarray(contacts)

    FS = spike_data['FS']
    spt = spt_dict['data']
    idx = np.arange(len(spt))
    inner_idx = filter_spt(spike_data, spt_dict, sp_win)
    outer_idx = idx[np.in1d(idx, inner_idx) == False]

    indices = (spt / 1000.0 * FS).astype(np.int32)
    win = (np.asarray(sp_win) / 1000.0 * FS).astype(np.int32)
    time = np.arange(win[1] - win[0]) * 1000.0 / FS + sp_win[0]
    n_contacts, n_pts = sp_data.shape

    # auxiliary function to find a valid spike window within data range
    minmax = lambda x: np.max([np.min([n_pts, x]), 0])
    spWave = np.zeros((len(time), len(spt), len(contacts)),
                      dtype=np.float32)

    for i in inner_idx:
        sp = indices[i]
        spWave[:, i, :] = np.atleast_2d(sp_data[contacts,
                                                sp + win[0]:sp + win[1]]).T

    for i in outer_idx:
        sp = indices[i]
        l, r = map(minmax, sp + win)
        if l != r:
            spWave[(l - sp) - win[0]:(r - sp) - win[0], i, :] = sp_data[contacts, l:r].T

    wavedict = {"data": spWave, "time": time, "FS": FS}

    if len(idx) != len(inner_idx):
        is_valid = np.zeros(len(spt), dtype=np.bool)
        is_valid[inner_idx] = True
        wavedict['is_valid'] = is_valid

    if resample != 1:
        warn("resample argument is deprecated."
             "Please update your code to use function"
             "resample_spikes", DeprecationWarning)
        wavedict = resample_spikes(wavedict, FS * resample)
    return wavedict


def resample_spikes(spikes_dict, FS_new):
    """Upsample spike waveforms using spline interpolation"""

    sp_waves = spikes_dict['data']
    time = spikes_dict['time']
    FS = spikes_dict['FS']

    resamp_time = np.arange(time[0], time[-1], 1000.0 / FS_new)
    n_pts, n_spikes, n_contacts = sp_waves.shape

    spike_resamp = np.empty((len(resamp_time), n_spikes, n_contacts))

    for i in range(n_spikes):
        for contact in range(n_contacts):
            tck = interpolate.splrep(time, sp_waves[:, i, contact], s=0)
            spike_resamp[:, i, contact] = interpolate.splev(resamp_time, tck, der=0)

    return {"data": spike_resamp, "time": resamp_time, "FS": FS}


def align_spikes(spike_data, spt_dict, sp_win, type="max", resample=1,
                contact=0, remove=True):
    """Aligns spike waves and returns corrected spike times

    Parameters
    ----------
    spike_data : dict
    spt_dict : dict
    sp_win : list of int
    type : {'max', 'min'}, optional
    resample : int, optional
    contact : int, optional
    remove : bool, optiona

    Returns
    -------
    ret_dict : dict
       spike times of aligned spikes

    """

    spt = spt_dict['data'].copy()

    idx_align = np.arange(len(spt))
    #spt_align = {'data': spt}

    #go in a loop until all spikes are correctly aligned
    iter_id = 0
    while len(idx_align) > 0:
        spt_align = {'data': spt[idx_align]}
        spt_inbound = filter_spt(spike_data, spt_align, sp_win)
        idx_align = idx_align[spt_inbound]
        #spt_align = {'data': spt[idx_align]}
        sp_waves_dict = extract_spikes(spike_data, spt_align, sp_win,
                                       resample=resample, contacts=contact)

        sp_waves = sp_waves_dict['data'][:, spt_inbound, 0]
        #if sp_waves_dict.has_key('is_valid'):
        #    sp_waves = sp_waves[:, sp_waves_dict['is_valid']]
        time = sp_waves_dict['time']

        if type == "max":
            i = sp_waves.argmax(0)
        elif type == "min":
            i = sp_waves.argmin(0)

        #move spike markers
        shift = time[i]
        spt[idx_align] += shift

        #if spike maximum/minimum was at the edge we have to extract it at the
        # new marker and repeat the alignment
        tol = 0.1
        idx_align = idx_align[(shift < (sp_win[0] + tol)) | (shift > (sp_win[1] - tol))]
        iter_id += 1
        #print "Align. iteration %d, remaining idx %d" % (iter_id, len(idx_align))
        #print shift

    ret_dict = {'data': spt}

    if remove:
        #remove double spikes
        FS = spike_data['FS']
        ret_dict = remove_doubles(ret_dict, 1000.0 / FS)

    return ret_dict


def remove_doubles(spt_dict, tol):
    new_dict = spt_dict.copy()
    spt = spt_dict['data']

    if len(spt) > 0:
        spt = spt[np.concatenate(([True], np.diff(spt) > tol))]

    new_dict['data'] = spt
    return new_dict


def merge_spikes(spike_waves1, spike_waves2):
    """Merges two sets of spike waves

    Parameters
    ----------

    spike_waves1 : dict
    spike_waves2 : dict
        spike wavefroms to merge; both spike wave sets must be defined
        within the same time window and with the same sampling
        frequency

    Returns
    -------

    spike_waves : dict
        merged spike waveshapes

    clust_idx : array
        labels denoting to which set the given spike originally belonged to
    """

    sp_data1 = spike_waves1['data']
    sp_data2 = spike_waves2['data']

    sp_data = np.hstack((sp_data1, sp_data2))
    spike_waves = spike_waves1.copy()
    spike_waves['data'] = sp_data

    clust_idx = np.concatenate((np.ones(sp_data1.shape[1]),
                               np.zeros(sp_data2.shape[1])))

    return spike_waves, clust_idx


def merge_spiketimes(spt1, spt2, sort=True):
    """Merges two sets of spike times

    Parameters
    ----------
    spt1 : dict
    spt2 : dict
    sort : bool, optional
        False if you don't want to be the spike times sorted.

    Returns
    -------
    spt : dict
        dictionary with merged spike time arrrays under data key
    clust_idx : array
        labels denoting to which set the given spike originally belonged
        to

    """

    spt_data1 = spt1['data']
    spt_data2 = spt2['data']

    spt_data = np.concatenate((spt_data1, spt_data2))
    i = spt_data.argsort()
    spt_data = spt_data[i]

    clust_idx = np.concatenate((np.ones(spt_data1.shape[0]),
                               np.zeros(spt_data2.shape[0])))
    clust_idx = clust_idx[i]
    spt_dict = {"data": spt_data}

    return spt_dict, clust_idx
