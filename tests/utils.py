from spike_beans import base
import numpy as np

spike_dur = 5.
spike_amp = 100.
FS = 25E3
period = 100
n_spikes = 100

class DummySignalSource(base.Component):

    def __init__(self):

        self.period = period
        self.n_spikes = n_spikes
        self.f_filter = None
        super(DummySignalSource, self).__init__()

        self._generate_data()

    def _generate_data(self):
        n_pts = int(self.n_spikes*self.period/1000.*FS)
        sp_idx = (np.arange(1,self.n_spikes-1)*self.period*FS/1000).astype(int)
        spikes = np.zeros(n_pts)[np.newaxis,:]
        spikes[0,sp_idx]=spike_amp

        n = int(spike_dur/1000.*FS) #spike length
        spikes[0,:] = np.convolve(spikes[0,:], np.ones(n), 'full')[:n_pts]
        self.spt = (sp_idx+0.5)*1000./FS
        self.FS = FS
        self._spikes = spikes

    def read_signal(self):

        #in milisecs

        spk_data ={"data":self._spikes,"n_contacts":1, "FS":self.FS}
        return spk_data

    def _update(self):
        self.period = period*2
        self.n_spikes = n_spikes/2
        self._generate_data()

    signal = property(read_signal)

class DummySpikeDetector(base.Component):
    def __init__(self):

        self.threshold = 500
        self.type = 'min'
        self.contact = 0
        self.sp_win = [-0.6, 0.8]
        super(DummySpikeDetector, self).__init__()
        self._generate_data()

    def _generate_data(self):
        n_pts = int(n_spikes*period/1000.*FS)
        sp_idx = (np.arange(1,n_spikes-1)*period*FS/1000).astype(int)
        spt = (sp_idx+0.5)*1000./FS
        self._spt_data = {'data':spt}

    def read_events(self):
        return self._spt_data

    events = property(read_events)

class DummyLabelSource(base.Component):
    def __init__(self):
        self.labels = np.random.randint(0,5, n_spikes-2)

class DummyFeatureExtractor(base.Component):

    def __init__(self):
        n_feats=2
        features = np.vstack((
            np.zeros((n_spikes, n_feats)),
            np.ones((n_spikes, n_feats))
        ))
        names = ["Fet{0}".format(i) for i in range(n_feats)]

        self._features = {"data": features, "names":names}

        super(DummyFeatureExtractor, self).__init__()

    def read_features(self):

        return self._features

    def add_feature(self, name):
        ''' adds random feature with specifid name '''

        self._features['data'] = np.hstack((self._features['data'], np.random.randn(n_spikes * 2, 1)))
        self._features['names'].append(name)

    def add_spikes(self, num = 10):
        '''appends `num` random values to each feature'''

        n_features = self._features['data'].shape[1]
        self._features['data'] = np.vstack((self._features['data'], np.random.randn(num, n_features)))

    features = property(read_features)

class DummySpikeSource(base.Component):

    def __init__(self):
        n_pts = 100
        spike_shape = np.zeros(n_pts)
        spike_shape[n_pts/2] = 1.
        data = spike_shape[:,np.newaxis, np.newaxis]*np.ones(n_spikes-2)[np.newaxis,:,np.newaxis]
        self._sp_waves = {'data':data, 'time':np.ones(n_pts)*1000./FS}
    def read_spikes(self):

        return self._sp_waves

    spikes = property(read_spikes)

class RandomFeatures(base.Component):
    def read_features(self):
        n_feats=2
        features = np.random.randn(n_spikes, n_feats)
        names = ["Fet{0}".format(i) for i in range(n_feats)]
        return {"data": features, "names":names}
    features = property(read_features)

