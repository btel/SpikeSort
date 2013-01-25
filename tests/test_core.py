import numpy as np
import spike_sort as ss

from nose.tools import ok_, eq_, raises
from numpy.testing import assert_array_almost_equal as almost_equal

import os

class TestFilter(object):
    def __init__(self):
        self.n_spikes = 100

        self.FS = 25E3
        self.period = 1000.0 / self.FS * 100
        self.time = np.arange(0, self.period * self.n_spikes, 1000.0 / self.FS)
        self.spikes = np.sin(2 * np.pi / self.period * self.time)[np.newaxis, :]
        self.spk_data = {"data": self.spikes, "n_contacts": 1, "FS": self.FS}

    def test_filter_proxy(self):
        sp_freq = 1000.0 / self.period
        filter = ss.filters.Filter(sp_freq * 0.5, sp_freq * 0.4, 1, 10, 'ellip')
        spk_filt = ss.filters.filter_proxy(self.spk_data, filter)
        ok_(self.spk_data['data'].shape == spk_filt['data'].shape)

    def test_remove_proxy_files_after_exit(self):
        filter_func = lambda x, f: x
        spk_filt = ss.filters.filter_proxy(self.spk_data, filter_func)
        fname = spk_filt['data']._v_file.filename
        ss.filters.clean_after_exit()
        assert not os.path.isfile(fname)
        

    def test_LiearIIR_detect(self):
        n_spikes = self.n_spikes
        period = self.period
        threshold = 0.5
        sp_freq = 1000.0 / period

        self.spk_data['data'] += 2
        spk_filt = ss.filters.fltLinearIIR(self.spk_data, sp_freq * 0.5, sp_freq * 0.4, 1, 10, 'ellip')

        spt = ss.extract.detect_spikes(spk_filt, thresh=threshold)
        ok_(len(spt['data']) == n_spikes)

class TestExtract(object):
    def __init__(self):
        self.n_spikes = 100

        self.FS = 25E3
        self.period = 1000.0 / self.FS * 100
        self.time = np.arange(0, self.period * self.n_spikes, 1000.0 / self.FS)
        self.spikes = np.sin(2 * np.pi / self.period * self.time)[np.newaxis, :]
        self.spk_data = {"data": self.spikes, "n_contacts": 1, "FS": self.FS}
    
    def test_detect(self):
        n_spikes = self.n_spikes
        period = self.period
        FS = self.FS
        time = self.time
        spikes = self.spikes
        threshold = 0.5
        crossings_real = period / 12.0 + np.arange(n_spikes) * period
        spt = ss.extract.detect_spikes(self.spk_data, thresh=threshold)
        ok_((np.abs(spt['data'] - crossings_real) <= 1000.0 / FS).all())
    
    def test_align(self):
        #check whether spikes are correctly aligned to maxima
        maxima_idx = self.period * (1 / 4.0 + np.arange(self.n_spikes))
        thr_crossings = self.period * (1 / 6.0 + np.arange(self.n_spikes))
        spt_dict = {"data": thr_crossings}
        sp_win = [-self.period / 6.0, self.period / 3.0]
        spt = ss.extract.align_spikes(self.spk_data, spt_dict, sp_win)
        ok_((np.abs(spt['data'] - maxima_idx) <= 1000.0 / self.FS).all())

    def test_align_short_win(self):
        #test spike alignment with windows shorter than total spike duration
        maxima_idx = self.period * (1 / 4.0 + np.arange(self.n_spikes))
        thr_crossings = self.period * (1 / 6.0 + np.arange(self.n_spikes))
        spt_dict = {"data": thr_crossings}
        sp_win = [-self.period / 24.0, self.period / 12.0]
        spt = ss.extract.align_spikes(self.spk_data, spt_dict, sp_win)
        ok_((np.abs(spt['data'] - maxima_idx) <= 1000.0 / self.FS).all())

    def test_align_edge(self):
        #???
        spikes = np.sin(2 * np.pi / self.period * self.time + np.pi / 2.0)[np.newaxis, :]
        maxima_idx = self.period * (np.arange(1, self.n_spikes + 1))
        thr_crossings = self.period * (-1 / 6.0 + np.arange(1, self.n_spikes + 1))
        spt_dict = {"data": thr_crossings}
        sp_win = [-self.period / 24.0, self.period / 12.0]
        spk_data = {"data": spikes, "n_contacts": 1, "FS": self.FS}
        spt = ss.extract.align_spikes(spk_data, spt_dict, sp_win)
        last = spt['data'][-1]
        ok_((last >= (self.time[-1] - sp_win[1])) & (last <= self.time[-1]))

    def test_align_double_spikes(self):
        #double detections of the same spike should be removed
        maxima_idx = self.period * (1 / 4.0 + np.arange(self.n_spikes))
        thr_crossings = self.period * (1 / 6.0 + np.arange(0, self.n_spikes, 0.5))
        spt_dict = {"data": thr_crossings}
        sp_win = [-self.period / 24.0, self.period / 12.0]
        spt = ss.extract.align_spikes(self.spk_data, spt_dict, sp_win)
        ok_((np.abs(spt['data'] - maxima_idx) <= 1000.0 / self.FS).all())

    def test_remove_doubles_roundoff(self):
        # remove_doubles should account for slight variations around
        # the 'tolerance' which may occur due to some round-off errors
        # during spike detection/alignment. This won't affect any useful
        # information, because the tolerance is one sample large in this
        # test.
        
        tol = 1000.0 / self.FS  # [ms], one sample
        data = [1.0, 1.0 + tol + tol * 0.01] 
        spt_dict = {"data": np.array(data)}
        clean_spt_dict = ss.extract.remove_doubles(spt_dict, tol)
        ok_(len(clean_spt_dict['data']) == 1) # duplicate removed

    def test_extract(self):
        zero_crossing = self.period * np.arange(self.n_spikes)
        zero_crossing += 1000.0 / self.FS / 2.0  # move by half a sample to avoid round-off errors
        spt_dict = {"data": zero_crossing}
        sp_win = [0, self.period]
        sp_waves = ss.extract.extract_spikes(self.spk_data, spt_dict, sp_win)
        ref_sp = np.sin(2 * np.pi / self.period * sp_waves['time'])
        ok_((np.abs(ref_sp[:, np.newaxis] - sp_waves['data'][:, :, 0]) < 1E-6).all())
        #ok_((np.abs(sp_waves['data'][:,:,0].mean(1)-ref_sp)<2*1000*np.pi/(self.FS*self.period)).all())
        #ok_(np.abs(np.sum(sp_waves['data'][:,:,0].mean(1)-ref_sp))<1E-6)

    #def test_extract_resample_deprecation(self):
    #    zero_crossing = self.period*np.arange(self.n_spikes)
    #    spt_dict = {"data":zero_crossing}
    #    sp_win = [0, self.period]
    #    with warnings.catch_warnings(True) as w:
    #        sp_waves = ss.extract.extract_spikes(self.spk_data, spt_dict, sp_win,
    #                                             resample=2.)
    #        ok_(len(w)>=1)

    def test_extract_and_resample(self):
        zero_crossing = self.period * np.arange(self.n_spikes)
        zero_crossing += 1000.0 / self.FS / 2.0  # move by half a sample to avoid round-off errors
        spt_dict = {"data": zero_crossing}
        sp_win = [0, self.period]
        sp_waves = ss.extract.extract_spikes(self.spk_data, spt_dict, sp_win)
        sp_resamp = ss.extract.resample_spikes(sp_waves, self.FS * 2)
        ref_sp = np.sin(2 * np.pi / self.period * sp_resamp['time'])
        ok_((np.abs(ref_sp[:, np.newaxis] - sp_resamp['data'][:, :, 0]) < 1E-6).all())

    def test_mask_of_truncated_spikes(self):
        zero_crossing = self.period * np.arange(self.n_spikes + 1)
        spt_dict = {"data": zero_crossing}
        sp_win = [0, self.period]
        sp_waves = ss.extract.extract_spikes(self.spk_data, spt_dict, sp_win)
        correct_mask = np.ones(self.n_spikes + 1, np.bool)
        correct_mask[-1] = False
        ok_((sp_waves['is_valid'] == correct_mask).all())
        #ok_(np.abs(np.sum(sp_waves['data'][:,:,0].mean(1)-ref_sp))<1E-6)

    def test_extract_truncated_spike_end(self):
        zero_crossing = np.array([self.period * (self.n_spikes - 0.5)])
        spt_dict = {"data": zero_crossing}
        sp_win = [0, self.period]
        sp_waves = ss.extract.extract_spikes(self.spk_data, spt_dict, sp_win)
        ref_sp = -np.sin(2 * np.pi / self.period * sp_waves['time'])
        ref_sp[len(ref_sp) / 2:] = 0
        almost_equal(sp_waves['data'][:, 0, 0], ref_sp)

    def test_extract_truncated_spike_end(self):
        zero_crossing = np.array([-self.period * 0.5])
        spt_dict = {"data": zero_crossing}
        sp_win = [0, self.period]
        sp_waves = ss.extract.extract_spikes(self.spk_data, spt_dict, sp_win)
        ref_sp = -np.sin(2 * np.pi / self.period * sp_waves['time'])
        ref_sp[:len(ref_sp) / 2] = 0
        almost_equal(sp_waves['data'][:, 0, 0], ref_sp)

    def test_filter_spt(self):
        #out of band spikes  should be removed
        zero_crossing = self.period * (np.arange(self.n_spikes))
        spt_dict = {"data": zero_crossing}
        sp_win = [0, self.period]
        spt_filt = ss.extract.filter_spt(self.spk_data, spt_dict, sp_win)
        ok_(len(spt_filt) == self.n_spikes)

    def test_filter_spt_shorten_left(self):
        #remove out-of-band spikes from the beginning of the train
        zero_crossing = self.period * (np.arange(self.n_spikes))
        spt_dict = {"data": zero_crossing}
        sp_win = [-self.period / 8, self.period / 8.0]
        spt_filt = ss.extract.filter_spt(self.spk_data, spt_dict, sp_win)
        ok_(len(spt_filt) == (self.n_spikes - 1))

    def test_filter_spt_shorten_right(self):
        #remove out-of-band spikes from the end of the train
        zero_crossing = self.period * (np.arange(self.n_spikes))
        spt_dict = {"data": zero_crossing}
        sp_win = [0, self.period + 1000.0 / self.FS]
        spt_filt = ss.extract.filter_spt(self.spk_data, spt_dict, sp_win)
        ok_(len(spt_filt) == (self.n_spikes - 1))

class TestCluster(object):
    """test clustering algorithms"""

    def _cmp_bin_partitions(self, cl1, cl2):
        return (~np.logical_xor(cl1, cl2)).all() or (np.logical_xor(cl1, cl2)).all()

    def setup(self):

        self.K = 2

        n_dim = 2
        pts_in_clust = 100
        np.random.seed(1234)
        data = np.vstack((np.random.rand(pts_in_clust, n_dim),
                  np.random.rand(pts_in_clust, n_dim) + 2 * np.ones(n_dim)))
        self.labels = np.concatenate((np.zeros(pts_in_clust, dtype=int),
                                            np.ones(pts_in_clust, dtype=int)))
        feature_labels = ["feat%d" % i for i in range(n_dim)]
        self.features = {"data": data, "names": feature_labels}

    def test_k_means(self):
        """test own k-means algorithm"""

        cl = ss.cluster.cluster('k_means', self.features, self.K)
        ok_(self._cmp_bin_partitions(cl, self.labels))

    def test_k_means_plus(self):
        """test scikits k-means plus algorithm"""

        cl = ss.cluster.cluster('k_means_plus', self.features, self.K)
        ok_(self._cmp_bin_partitions(cl, self.labels))

    def test_gmm(self):
        """test gmm clustering algorithm"""

        cl = ss.cluster.cluster('gmm', self.features, self.K)
        ok_(self._cmp_bin_partitions(cl, self.labels))

    def test_random(self):
        cl = np.random.rand(len(self.labels)) > 0.5
        ok_(~self._cmp_bin_partitions(cl, self.labels))

    @raises(NotImplementedError)
    def test_method_notimplemented(self):
        cl = ss.cluster.cluster("notimplemented", self.features)
