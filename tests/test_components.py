from SpikeBeans import base, components
from SpikeBeans.base import features
from nose.tools import ok_,raises
from nose import with_setup
import numpy as np

def setup():
    "set up test fixtures"

def teardown():
    "tear down test fixtures"
    base.features = base.FeatureBroker()
    

class DummySignalSource(base.Component):
    def read_signal(self):
        n_spikes = 100
        
        FS = 25E3
        
        period = 100 #in milisecs
        n_pts = int(n_spikes*period/1000.*FS)
        sp_idx = (np.arange(1,n_spikes-1)*period*FS/1000).astype(int)
        spikes = np.zeros(n_pts)[np.newaxis,:]
        spikes[0,sp_idx]=100.
        self.spt = sp_idx*1000./FS
        self.FS = FS
        spk_data ={"data":spikes,"n_contacts":1, "FS":FS}
        return spk_data


def test_spike_detection():
    features.Provide("SignalSource", DummySignalSource())
    detector = components.SpikeDetector(thresh=50.)
    spt = detector.read_events()
    
    source = features['SignalSource']
    ok_((np.abs(spt['data']-source.spt)<=1000./source.FS).all())
    
    