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

    def test_fetPC(self):
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

        pcs = ss.features.fetPCs(spikes_dict, ncomps=1)
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
