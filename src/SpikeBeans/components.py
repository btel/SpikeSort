'''
Created on Feb 9, 2011

@author: bartosz
'''

import spike_sort as sort

import base
from spike_sort.io.filters import BakerlabFilter, PyTablesFilter

class BakerlabSource(base.Component, BakerlabFilter):
    
    def __init__(self, conf_file, dataset):
        BakerlabFilter.__init__(self, conf_file)
        self.dataset = dataset
        self.signal = None
        self.events = None
        
    def read_signal(self):
        if self.signal is None:
            self.signal = self.read_sp(self.dataset)
        return self.signal
    
    def read_events(self):
        if self.events is None:
            self.events = self.read_spt(self.dataset)
        return self.events
    

class SpikeDetector(base.Component):
    """Detect Spikes with alignment"""
    waveform_src = base.RequiredFeature("SignalSource", 
                                        base.HasAttributes("read_signal"))

    def __init__(self, thresh='auto', 
                 contact=0, 
                 type='max', 
                 resample=1, 
                 sp_win=[-0.2, 0.8]):
        self.thresh = thresh
        self.contact = contact
        self.type = "max"
        self.resample = resample
        self.sp_win = sp_win
        self.sp_times = None
        
    def _detect(self):
        sp = self.waveform_src.read_signal()
        spt = sort.extract.detect_spikes(sp,  contact=self.contact,
                                               thresh=self.thresh)
        self.sp_times = sort.extract.align_spikes(sp, spt, 
                                                  self.sp_win, 
                                                  type=self.type, 
                                                  resample=self.resample)
        
    def read_events(self):
        if self.sp_times is None:
            self._detect()
        return self.sp_times
    
class SpikeExtractor(base.Component):
    waveform_src = base.RequiredFeature("SignalSource", 
                                    base.HasAttributes("read_signal"))
    spike_times = base.RequiredFeature("SpikeMarkerSource", 
                                    base.HasAttributes("read_events"))
    
    def __init__(self, sp_win=[-0.2,0.8]):
        self.sp_shapes = None
        self.sp_win = sp_win
    
    def _extract_spikes(self):
        sp = self.waveform_src.read_signal()
        spt = self.spike_times.read_events()
        self.sp_shapes = sort.extract.extract_spikes(sp, spt, self.sp_win)
    
    def read_spikes(self):
        if self.sp_shapes is None:
            self._extract_spikes()
        return self.sp_shapes
            