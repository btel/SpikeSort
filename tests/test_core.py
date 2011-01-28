import numpy as np
import matplotlib.pyplot as plt
import spike_sort as ss

from nose.tools import ok_, eq_, raises

class TestFeatures:

    def setup(self):
        self.gain = 1000
        
        n_pts = 100
        FS = 5E3
        n_contacts = 4
        time = np.arange(n_pts, dtype=np.float32)/n_pts - 0.5
        cells = self.gain*np.vstack(
            (np.sin(time*2*np.pi)/2.,
             np.abs(time)*2 - 1
            )
        )
        self.cells = cells.astype(np.int)
        
        #define data interface (dictionary) 
        self.spikes_dict = {"data":self.cells, "time":time, "FS":FS}

    def test_fetP2P(self):
        spikes_dict = self.spikes_dict.copy()
        
        n_spikes = 200
        amps = np.random.randint(1, 100, n_spikes)
        amps = amps[:, np.newaxis]
        spikes = amps*self.cells[0,:]
        spikes_dict['data'] = spikes.T[:,:,np.newaxis]
        p2p, names = ss.features.fetP2P(spikes_dict)
        
        ok_((p2p==amps*self.gain).all())

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

        pcs, names = ss.features.fetPCs(spikes_dict, ncomps=1)
        compare = ~np.logical_xor(pcs[:,0].astype(int)+1, _amps)
        correct = np.sum(compare)

        eq_(n_spikes, correct)
        
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
        """test own k-means algorithm"""
        
        cl = ss.cluster.cluster('k_means_plus', self.features, self.K)
        ok_(self._cmp_bin_partitions(cl, self.labels))
        
    def test_gmm(self):
        """test own k-means algorithm"""
        
        cl = ss.cluster.cluster('gmm', self.features, self.K)
        ok_(self._cmp_bin_partitions(cl, self.labels))
    
    def test_random(self):
        cl = np.random.rand(len(self.labels))>0.5
        ok_(~self._cmp_bin_partitions(cl, self.labels))
        
    @raises(NotImplementedError)
    def test_method_notimplemented(self):
        cl = ss.cluster.cluster("notimplemented", self.features)
        
        
        