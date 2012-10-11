import numpy as np
from spike_sort import stats

from nose.tools import ok_

class TestStats(object):
    def setup(self):
        n_samples = 1000
        
        # unimodal
        normal = np.random.randn(n_samples)
        self.normal = (normal - normal.mean()) / normal.std()
        uniform = np.random.uniform(size = n_samples)
        self.uniform = (uniform - uniform.mean()) / uniform.std()
        laplace = np.random.laplace(size = n_samples)
        self.laplace = (laplace - laplace.mean()) / laplace.std()

        # multimodal
        multimodal = np.hstack((self.normal, self.uniform + 4, self.laplace + 8))
        np.random.shuffle(multimodal)
        multimodal = multimodal[:n_samples]
        self.multimodal = (multimodal - multimodal.mean()) / multimodal.std()

    def test_multimodality_detection(self):
        data = [self.normal, self.uniform, self.laplace]
        tests = [stats.dip1d, stats.ks1d]

        detected = [test(self.multimodal) > test(dist) for dist in data for test in tests]

        ok_(np.array(detected).all())
