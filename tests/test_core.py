import numpy as np
import matplotlib.pyplot as plt
import spike_sort as ss

from nose.tools import ok_, eq_, raises
from numpy.testing import assert_array_almost_equal as almost_equal

class TestExtract:
    
    def __init__(self):
        self.n_spikes = 100
        
        self.FS = 25E3
        self.period = 1000./self.FS*100
        self.time = np.arange(0, self.period*self.n_spikes, 1000./self.FS)
        self.spikes = np.sin(2*np.pi/self.period*self.time)[np.newaxis,:]
        self.spk_data ={"data":self.spikes,"n_contacts":1, "FS":self.FS}
        
    def test_filter_proxy(self):
        sp_freq = 1000./self.period
        filter = ss.extract.Filter(sp_freq*0.5, sp_freq*0.4, 1, 10, 'ellip')
        spk_filt = ss.extract.filter_proxy(self.spk_data, filter)
        ok_(self.spk_data['data'].shape == spk_filt['data'].shape)
        
    
    def test_detect(self):
        n_spikes = self.n_spikes
        period = self.period
        FS = self.FS
        time = self.time
        spikes = self.spikes
        threshold = 0.5
        crossings_real = period/12.+np.arange(n_spikes)*period 
        spt = ss.extract.detect_spikes(self.spk_data, thresh=threshold)
        ok_((np.abs(spt['data']-crossings_real)<=1000./FS).all())
        
    def test_filter_detect(self):
        n_spikes = self.n_spikes
        period = self.period
        FS = self.FS
        time = self.time
        threshold = 0.5
        sp_freq = 1000./self.period
        self.spk_data['data']+=2
        filter = ss.extract.Filter(sp_freq*0.5, sp_freq*0.4, 1, 10, 'ellip')
        spt = ss.extract.detect_spikes(self.spk_data, thresh=threshold, filter=filter)
        ok_(len(spt['data'])==n_spikes)
        
    def test_align(self):
        #check whether spikes are correctly aligned to maxima
        maxima_idx = self.period*(1/4.+np.arange(self.n_spikes))
        thr_crossings = self.period*(1/6. + np.arange(self.n_spikes))
        spt_dict = {"data":thr_crossings}
        sp_win = [-self.period/6., self.period/3.]
        spt = ss.extract.align_spikes(self.spk_data, spt_dict, sp_win)
        ok_((np.abs(spt['data']-maxima_idx)<=1000./self.FS).all())
    
    def test_align_short_win(self):
        #test spike alignment with windows shorter than total spike duration
        maxima_idx = self.period*(1/4.+np.arange(self.n_spikes))
        thr_crossings = self.period*(1/6. + np.arange(self.n_spikes))
        spt_dict = {"data":thr_crossings}
        sp_win = [-self.period/24., self.period/12.]
        spt = ss.extract.align_spikes(self.spk_data, spt_dict, sp_win)
        ok_((np.abs(spt['data']-maxima_idx)<=1000./self.FS).all())
        
    def test_align_edge(self):
        #???
        spikes = np.sin(2*np.pi/self.period*self.time+np.pi/2.)[np.newaxis,:]
        maxima_idx = self.period*(np.arange(1, self.n_spikes+1))
        thr_crossings = self.period*(-1/6.+np.arange(1, self.n_spikes+1))
        spt_dict = {"data":thr_crossings}
        sp_win = [-self.period/24., self.period/12.]
        spk_data ={"data":spikes,"n_contacts":1, "FS":self.FS}
        spt = ss.extract.align_spikes(spk_data, spt_dict, sp_win)
        last = spt['data'][-1]     
        ok_((last>=(self.time[-1]-sp_win[1])) & (last<=self.time[-1]))
    
    def test_align_double_spikes(self):
        #double detections of the same spike should be removed
        maxima_idx = self.period*(1/4.+np.arange(self.n_spikes))
        thr_crossings = self.period*(1/6. + np.arange(0,self.n_spikes, 0.5))
        spt_dict = {"data":thr_crossings}
        sp_win = [-self.period/24., self.period/12.]
        spt = ss.extract.align_spikes(self.spk_data, spt_dict, sp_win)
        ok_((np.abs(spt['data']-maxima_idx)<=1000./self.FS).all())
        
    
    def test_extract(self):
        zero_crossing = self.period*np.arange(self.n_spikes)
        spt_dict = {"data":zero_crossing}
        sp_win = [0, self.period]
        sp_waves = ss.extract.extract_spikes(self.spk_data, spt_dict, sp_win)
        ref_sp = np.sin(2*np.pi/self.period*sp_waves['time'])
        ok_((np.abs(sp_waves['data'][:,:,0].mean(1)-ref_sp)<2*1000*np.pi/(self.FS*self.period)).all())
        #ok_(np.abs(np.sum(sp_waves['data'][:,:,0].mean(1)-ref_sp))<1E-6)
        
    def test_filter_spt(self):
        #out of band spikes  should be removed
        zero_crossing = self.period*(np.arange(self.n_spikes))
        spt_dict = {"data":zero_crossing}
        sp_win = [0, self.period]
        spt_filt = ss.extract.filter_spt(self.spk_data, spt_dict, sp_win)
        ok_(len(spt_filt)==self.n_spikes)
        
    def test_filter_spt_shorten_left(self):
        #remove out-of-band spikes from the beginning of the train
        zero_crossing = self.period*(np.arange(self.n_spikes))
        spt_dict = {"data":zero_crossing}
        sp_win = [-self.period/8, self.period/8.]
        spt_filt = ss.extract.filter_spt(self.spk_data, spt_dict, sp_win)
        ok_(len(spt_filt)==(self.n_spikes-1))
    
    def test_filter_spt_shorten_right(self):
        #remove out-of-band spikes from the end of the train
        zero_crossing = self.period*(np.arange(self.n_spikes))
        spt_dict = {"data":zero_crossing}
        sp_win = [0, self.period+1000./self.FS]
        spt_filt = ss.extract.filter_spt(self.spk_data, spt_dict, sp_win)
        ok_(len(spt_filt)==(self.n_spikes-1))
        
    
        
    
   

class TestFeatures:

    def setup(self):
        self.gain = 1000
        
        n_pts = 100
        FS = 5E3
        time = np.arange(n_pts, dtype=np.float32)/n_pts - 0.5
        cells = self.gain*np.vstack(
            (np.sin(time*2*np.pi)/2.,
             np.abs(time)*2 - 1
            )
        )
        self.cells = cells.astype(np.int)
        
        #define data interface (dictionary) 
        self.spikes_dict = {"data":self.cells.T[:, :, np.newaxis], "time":time, "FS":FS}

    def test_fetP2P(self):
        spikes_dict = self.spikes_dict.copy()
        
        n_spikes = 200
        amps = np.random.randint(1, 100, n_spikes)
        amps = amps[:, np.newaxis]
        spikes = amps*self.cells[0,:]
        spikes_dict['data'] = spikes.T[:,:,np.newaxis]
        p2p = ss.features.fetP2P(spikes_dict)
        
        ok_((p2p['data']==amps*self.gain).all())

    def test_PCA(self):
        n_dim = 2
        n_obs = 100
        raw_data = np.random.randn(n_dim, n_obs)
        mixing = np.array([[-1,1],[1,1]])

        mix_data = np.dot(mixing, raw_data)
        mix_data = mix_data - mix_data.mean(1)[:,np.newaxis]

        evals, evecs, score = ss.features.PCA(mix_data, n_dim)
        
        proj_cov = np.cov(score)
        error = np.mean((proj_cov-np.eye(n_dim))**2)
        ok_(error<0.01)


    def test_fetPC(self):
        spikes_dict = self.spikes_dict.copy()
        n_spikes = 200
        n_cells = 2
        amp_var_fact = 0.4
        amp_var = 1+amp_var_fact*np.random.rand(n_spikes,1)
        _amps = np.random.rand(n_spikes)>0.5
        amps = amp_var*np.vstack((_amps>=0.5,_amps<0.5)).T
        
        spikes = np.dot(amps, self.cells).T
        spikes = spikes.astype(np.float32)
        spikes_dict['data'] = spikes[:,:,np.newaxis]

        pcs = ss.features.fetPCs(spikes_dict, ncomps=1)
        pcs = pcs['data']
        compare = ~np.logical_xor(pcs[:,0].astype(int)+1, _amps)
        correct = np.sum(compare)

        eq_(n_spikes, correct)
        
    def test_getSpProjection(self):
        spikes_dict = self.spikes_dict.copy()
        cells = spikes_dict['data']
        spikes_dict['data'] = np.repeat(cells, 10, 1)
        labels = np.repeat([0,1], 10,0)
        
        feat = ss.features.fetSpProjection(spikes_dict, labels)
        ok_(((feat['data'][:,0]>0.5) == labels).all())
        
class TestCluster:
    """test clustering algorithms"""
    
    def _cmp_bin_partitions(self, cl1, cl2):
        return (~np.logical_xor(cl1, cl2)).all() or (np.logical_xor(cl1,cl2)).all()
        
    def setup(self):

        self.K = 2
        
        n_dim = 2
        pts_in_clust = 100
        np.random.seed(1234)
        data = np.vstack((np.random.rand(pts_in_clust,n_dim), 
                  np.random.rand(pts_in_clust,n_dim)+2*np.ones(n_dim)))
        self.labels = np.concatenate((np.zeros(pts_in_clust, dtype=int),
                                            np.ones(pts_in_clust, dtype=int)))
        feature_labels = ["feat%d" % i for i in range(n_dim)]
        self.features = {"data":data, "names":feature_labels}
        
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
        cl = np.random.rand(len(self.labels))>0.5
        ok_(~self._cmp_bin_partitions(cl, self.labels))
        
    @raises(NotImplementedError)
    def test_method_notimplemented(self):
        cl = ss.cluster.cluster("notimplemented", self.features)
        
        
        
