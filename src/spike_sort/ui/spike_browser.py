
import numpy as np
from numpy import arange, sin, pi, float, size

import matplotlib

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from matplotlib.widgets import Button

import wx
import time

class SpikeBrowserFrame(wx.Frame):
    def __init__(self, parent, id):
        
        self.sp_win = [-0.8, 1]
        self.spike_collection = None
        
        wx.Frame.__init__(self,parent, id, 'scrollable plot',
                style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER,
                size=(1200, 400))
        self.panel = wx.Panel(self, -1)

        self.fig = Figure((5, 4), 75)
        self.canvas = FigureCanvasWxAgg(self.panel, -1, self.fig)
        self._mpl_init()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, -1, wx.EXPAND)

        self.panel.SetSizer(sizer)
        self.panel.Fit()

        self.canvas.Bind(wx.EVT_SCROLLWIN, self.OnScrollEvt)
        
        self.canvas.mpl_connect('key_press_event', self._on_key)
        

        
    def _mpl_init(self):
        self.fig.clf()
        self.axes = self.fig.add_axes([0.05, 0.1, 0.95,0.9])
        self.ax_prev = self.fig.add_axes([0.8, 0.0, 0.1,0.05])
        self.ax_next = self.fig.add_axes([0.9, 0.0, 0.1,0.05])
        
        self.b_next = Button(self.ax_next, 'Next')
        self.b_prev = Button(self.ax_prev, "Prev")

        self.b_next.on_clicked(self._next_spike)
        self.b_prev.on_clicked(self._prev_spike)
        self.i_spike = 0
        
    def _next_spike(self, event):
        try:
            t = (self.i_start+self.i_window/2.)*1000./self.FS
            t_spk = self.spt[self.spt>t][0]
            self.i_start = int(np.ceil(t_spk/1000.*self.FS-self.i_window/2.))
            self.i_end = self.i_start + self.i_window
            self.canvas.SetScrollPos(wx.HORIZONTAL, self.i_start)
            self.draw_plot()
        except IndexError:
            pass
        
    def _prev_spike(self, event):
        try:
            t = (self.i_start+self.i_window/2.)*1000./self.FS
            t_spk = self.spt[self.spt<t][-1]
            self.i_start = int(t_spk/1000.*self.FS-self.i_window/2.)
            self.i_end = self.i_start + self.i_window
            self.canvas.SetScrollPos(wx.HORIZONTAL, self.i_start)
            self.draw_plot()
        except IndexError:
            pass

    def _on_key(self, event):
        if event.key=='+' or event.key=='=':
            self.ylims/=2.
        elif event.key == '-':
            self.ylims*=2.
        else:
            return
        offset = self.ylims[1]-self.ylims[0]
        self.offsets = np.arange(self.n_chans)*offset
        self.draw_plot()

    def init_data(self, data, spk_idx=None):

        self.x = data['data']
        self.FS = data['FS']
        n_chans, n_pts = self.x.shape
        if spk_idx:
            self.spt = spk_idx['data']
        else:
            self.spt = None
            self.ax_next.set_visible(False)
            self.ax_prev.set_visible(False)


        self.i_window = int(self.winsz/1000.*self.FS)
        # Extents of data sequence: 
        self.i_min = 0
        self.i_max = n_pts - self.i_window
        self.n_chans = n_chans

        self.canvas.SetScrollbar(wx.HORIZONTAL, 0, 5, self.i_max)      

        # Indices of data interval to be plotted:
        self.i_start = 0
        self.i_end = self.i_start + self.i_window
        
        
        self.time = np.arange(self.i_start,self.i_end)*1./self.FS
        
        self.segs = np.empty((n_chans, self.i_window, 2))
        self.segs[:,:,0] = self.time[np.newaxis,:]
        self.segs[:,:,1] = self.x[:,self.i_start:self.i_end]
         
        ylims = (self.segs[:,:,1].min(), self.segs[:,:,1].max())
        offset = ylims[1]-ylims[0]
        self.offsets = np.arange(n_chans)*offset
        self.segs[:,:,1] += self.offsets[:,np.newaxis]
        
        self.ylims = np.array(ylims)

        self.line_collection = LineCollection(self.segs,
                                              offsets=None,
                                              transform=self.axes.transData)

        self.axes.add_collection(self.line_collection)
        self.axes.set_xlim((self.time[0], self.time[-1]))
        self.axes.set_ylim((self.ylims[0]+self.offsets.min(), 
                            self.ylims[1]+self.offsets.max()))
        start = time.time()
        self.canvas.draw()

    def draw_plot(self):

        self.time = np.arange(self.i_start,self.i_end)*1./self.FS
        self.segs[:,:,0] = self.time[np.newaxis,:]
        self.segs[:,:,1] = self.x[:,self.i_start:self.i_end]+self.offsets[:,np.newaxis]
        self.line_collection.set_segments(self.segs)

        # Adjust plot limits:
        self.axes.set_xlim((self.time[0], self.time[-1]))
        self.axes.set_ylim((self.ylims[0]+self.offsets.min(), 
                            self.ylims[1]+self.offsets.max()))
        
        if self.spt is not None:
            self.draw_spikes()
        # Redraw:                  
        self.canvas.draw()
        
    def draw_spikes(self):
        if self.spike_collection is not None:
            self.spike_collection.remove()
            self.spike_collection = None
        sp_win = self.sp_win 
        time = self.segs[0,:,0]*1000.
        t_min, t_max = time[0]-sp_win[0], time[-1]-sp_win[1]
        spt = self.spt[(self.spt>t_min) & (self.spt<t_max)]
        if len(spt)>0:
            n_pts = int((sp_win[1]-sp_win[0])/1000.*self.FS)
            sp_segs = np.empty((len(spt), self.n_chans, n_pts, 2))
            for i in range(len(spt)):
                start, = np.nonzero(time>=(spt[i]+sp_win[0]))
                start = start[0]
                stop  = start+n_pts
                sp_segs[i,:,:,0] = (time[np.newaxis,start:stop]/1000.)
                sp_segs[i,:,:,1] = self.segs[:, start:stop, 1]
            sp_segs = sp_segs.reshape(-1, n_pts, 2)
            self.spike_collection = LineCollection(sp_segs,
                                                  offsets=None,
                                                  color='r',
                                                  transform=self.axes.transData)
            self.axes.add_collection(self.spike_collection)
            
        
    def OnScrollEvt(self, event):

        # Update the indices of the plot:
        self.i_start = self.i_min + event.GetPosition()
        self.i_end = self.i_min + self.i_window + event.GetPosition()
        self.draw_plot()

class SpikeBrowserApp(wx.App):

    def OnInit(self):
        self.frame = SpikeBrowserFrame(parent=None,id=-1)
        self.frame.Show(True)
        self.SetTopWindow(self.frame)
        return True

    def init_data(self, data, spk_idx, winsz):
        self.frame.winsz = winsz
        self.frame.init_data(data, spk_idx)

def browse_data(data, spk_idx=None, win=100):
    # Generate some data to plot:
    
    app = SpikeBrowserApp()
    app.init_data(data, spk_idx, win)
    app.MainLoop()


