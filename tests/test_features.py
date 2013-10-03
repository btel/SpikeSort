import numpy as np
import spike_sort as ss

from nose.tools import ok_, eq_, raises
from numpy.testing import assert_array_almost_equal as almost_equal


class TestFeatures(object):
    def setup(self):
        self.gain = 1000

        n_pts = 100
        FS = 5E3
        time = np.arange(n_pts, dtype=np.float32) / n_pts - 0.5
        cells = self.gain * np.vstack(
            (np.sin(time * 2 * np.pi) / 2.0,
             np.abs(time) * 2 - 1
            )
        )
        self.cells = cells.astype(np.int)

        #define data interface (dictionary)
        self.spikes_dict = {"data": self.cells.T[:, :, np.newaxis],
                            "time": time, "FS": FS}

    def test_fetP2P(self):
        spikes_dict = self.spikes_dict.copy()

        n_spikes = 200
        amps = np.random.randint(1, 100, n_spikes)
        amps = amps[:, np.newaxis]
        spikes = amps * self.cells[0, :]
        spikes_dict['data'] = spikes.T[:, :, np.newaxis]
        p2p = ss.features.fetP2P(spikes_dict)

        ok_((p2p['data'] == amps * self.gain).all())

    def test_PCA(self):
        n_dim = 2
        n_obs = 100
        raw_data = np.random.randn(n_dim, n_obs)
        mixing = np.array([[-1, 1], [1, 1]])

        mix_data = np.dot(mixing, raw_data)
        mix_data = mix_data - mix_data.mean(1)[:, np.newaxis]

        evals, evecs, score = ss.features.PCA(mix_data, n_dim)

        proj_cov = np.cov(score)
        error = np.mean((proj_cov - np.eye(n_dim)) ** 2)
        ok_(error < 0.01)

    def test_fetPCA(self):
        spikes_dict = self.spikes_dict.copy()
        n_spikes = 200
        n_cells = 2
        amp_var_fact = 0.4
        amp_var = 1 + amp_var_fact * np.random.rand(n_spikes, 1)
        _amps = np.random.rand(n_spikes) > 0.5
        amps = amp_var * np.vstack((_amps >= 0.5, _amps < 0.5)).T

        spikes = np.dot(amps, self.cells).T
        spikes = spikes.astype(np.float32)
        spikes_dict['data'] = spikes[:, :, np.newaxis]

        pcs = ss.features.fetPCA(spikes_dict, ncomps=1)
        pcs = pcs['data']
        compare = ~np.logical_xor(pcs[:, 0].astype(int) + 1, _amps)
        correct = np.sum(compare)
        eq_(n_spikes, correct)

    def test_getSpProjection(self):
        spikes_dict = self.spikes_dict.copy()
        cells = spikes_dict['data']
        spikes_dict['data'] = np.repeat(cells, 10, 1)
        labels = np.repeat([0, 1], 10, 0)

        feat = ss.features.fetSpProjection(spikes_dict, labels)
        ok_(((feat['data'][:, 0] > 0.5) == labels).all())

    def test_fetMarkers(self):
        spikes_dict = self.spikes_dict.copy()

        n_spikes = 200
        n_pts = self.cells.shape[1]
        amps = np.random.randint(1, 100, n_spikes)
        amps = amps[:, np.newaxis]
        spikes = amps * self.cells[0, :]
        spikes_dict['data'] = spikes.T[:, :, np.newaxis]
        time = spikes_dict['time']
        indices = np.array([0, n_pts/2, n_pts-1])
        values = ss.features.fetMarkers(spikes_dict, time[indices])

        ok_((values['data'] == spikes[:, indices]).all())

    def test_WT(self):
        "simple test for linearity of wavelet transform"
        spike1, spike2 = self.cells
        spike3 = 0.1 * spike1 + 0.7 * spike2
        spikes = np.vstack((spike1, spike2, spike3)).T
        spikes = spikes[:, :, np.newaxis] # WT only accepts 3D arrays

        wavelet = 'db3'
        wt = ss.features.WT(spikes, wavelet)
        wt1, wt2, wt3 = wt.squeeze().T

        almost_equal(wt3, 0.1 * wt1 + 0.7 * wt2)

    def test_fetWT_math(self):
        n_samples = 256

        # upsampled haar wavelet
        spike1 = np.hstack((np.ones(n_samples / 2), -1 * np.ones(n_samples / 2)))
        # upsampled haar scaling function
        spike2 = np.ones(n_samples)
        
        spikes = np.vstack((spike1, spike2)).T
        spikes = spikes[:, :, np.newaxis]
        spikes_dict = {'data' : spikes}

        features = ss.features.fetWT(spikes_dict, n_samples, wavelet='haar', select_method=None)
        idx = np.nonzero(features['data']) # nonzero indices
        
        # if nonzero elements are ONLY at (0,1) and (1,0),
        # this should be eye(2)
        eye = np.fliplr(np.vstack(idx))

        ok_((eye == np.eye(2)).all())

    def test_fetWT_selection(self):
        n_samples = 30
        n_channels = 2
        n_spikes = 50
        n_features = 10
        methods = [None, 'std', 'std_r', 'ks', 'dip', 'ksPCA', 'dipPCA']

        spikes = np.random.randn(n_samples, n_spikes, n_channels)
        spikes_dict = {'data' : spikes}

        shapes = [(n_spikes, n_features * n_channels)]

        for met in methods:
            wt = ss.features.fetWT(spikes_dict, n_features, wavelet='haar', select_method=met)
            shapes.append(wt['data'].shape)

        equal = lambda x, y: x == y and y or False
        success = bool(reduce(equal, shapes)) # returned shapes for all methods are correct

        ok_(success)

    def test_add_mask_decorator(self):
        spikes_dict = {'data': np.zeros((10, 2)),
                       'is_valid': np.zeros(2,)}

        fetIdentity = lambda x: {'data': x['data'], 'names': 'Identity'}
        deco_fet = ss.features.add_mask(fetIdentity)
        features = deco_fet(spikes_dict)

        ok_((features['is_valid'] == spikes_dict['is_valid']).all())

    def test_combine_features_without_mask(self):
        feature1 = {'data': np.random.uniform(size=(5, 1)), 'names': ['feature1']}
        feature2 = {'data': np.random.uniform(size=(5, 1)), 'names': ['feature2']}
        combined = ss.features.combine((feature1, feature2))
        ok_('is_valid' not in combined)

    def test_combine_features_with_one_mask(self):
        feature1 = {'data': np.random.uniform(size=(5, 1)), 'names': ['feature1']}
        feature2 = {'data': np.random.uniform(size=(5, 1)), 'names': ['feature2'],
                    'is_valid': np.ones(5, dtype=np.bool)}
        combined = ss.features.combine((feature1, feature2))
        ok_((combined['is_valid'] == feature2['is_valid']).all())

    def test_combine_features_with_different_masks(self):
        mask1 = np.ones(5, dtype=np.bool)
        mask1[-1] = False
        mask2 = mask1.copy()
        mask2[:2] = False
        feature1 = {'data': np.random.uniform(size=(5, 1)), 'names': ['feature1'],
                    'is_valid': mask1}
        feature2 = {'data': np.random.uniform(size=(5, 1)), 'names': ['feature2'],
                    'is_valid': mask2}
        combined = ss.features.combine((feature1, feature2))
        ok_((combined['is_valid'] == (mask1 & mask2)).all())

    def test_add_method_suffix(self):
        names = ['methA', 'methB', 'methA_1', 'methC', 'methA_2', 'methD']

        test1 = ss.features._add_method_suffix('methE', names) == 'methE'
        test2 = ss.features._add_method_suffix('methB', names) == 'methB_1'
        test3 = ss.features._add_method_suffix('methA', names) == 'methA_3'

        ok_(test1 and test2 and test3)
