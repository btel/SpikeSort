
import numpy as np
from numpy import arange, sin, pi, float, size

import matplotlib

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from matplotlib.widgets import Button

import wx
import Tkinter as Tk

from spike_sort.ui import label_color

class PlotWithScrollBarTk(object):
    def __init__(self):
        
        self.max = 0
        self.cur_pos = 0
        self.page_sz = 0
        self.root = Tk.Tk()
    def get_canvas(self, fig):
        
        self.canvas = FigureCanvasTkAgg(fig, master=self.root)
        self.canvas.show()
        self.scrollbar = Tk.Scrollbar(self.root, orient=Tk.HORIZONTAL)
        self.scrollbar.pack(side=Tk.BOTTOM, fill=Tk.BOTH)
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.scrollbar.config(command=self._callback)
        return self.canvas
    
    def _callback(self, mode, *args):
        if mode == 'moveto':
            pos, = args
            self.cur_pos = int(float(pos) * self.max)
        elif mode=="scroll":
            delta, units = args
            delta = float(delta)
            if units == 'units':
                self.cur_pos += delta*self.page_sz/5.
            elif units == 'pages':
                self.cur_pos += delta*self.page_sz
        self.set_scroll_pos(self.cur_pos)
        self.handler(self.cur_pos)
            
    def set_scroll_handler(self, handler):
        self.handler = handler
    
    def set_scroll_pos(self, pos):
        min, max = str(pos*1./self.max), str((pos+self.page_sz)*1./self.max)
        self.scrollbar.set(min, max)
    
    def set_scroll_max(self, max, page_size):
        self.page_sz = page_size 
        self.max = max
    
class PlotWithScrollBarWx(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self,parent, id, 'scrollable plot',
                          style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER,
                           size=(1200, 400))
    
    def get_canvas(self, fig):
        self.panel = wx.Panel(self, -1)
        self.canvas = FigureCanvasWxAgg(self.panel, -1, fig)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, -1, wx.EXPAND)
        self.panel.SetSizer(sizer)
        self.panel.Fit()
        self.canvas.Bind(wx.EVT_SCROLLWIN, self._callback)
        
        return self.canvas

    def _callback(self, event):
        pos = event.GetPosition()
        self.handler(pos)
        
    def set_scroll_handler(self, handler):
        self.handler = handler
    
    def set_scroll_pos(self, pos):
        self.canvas.SetScrollPos(wx.HORIZONTAL, pos)
        
    def set_scroll_max(self, max, page_size):
        self.canvas.SetScrollbar(wx.HORIZONTAL, 0, 5, max)  
        
        

class SpikeBrowserUI(object):
    def __init__(self, window):
        self.window = window
        self.sp_win = [-0.8, 1]
        self.spike_collection = None
        
        self.fig = Figure((5, 4), 75)
        
        self.canvas = window.get_canvas(self.fig)
        
        self._mpl_init()
        self.canvas.mpl_connect('key_press_event', self._on_key)
        self.window.set_scroll_handler(self.OnScrollEvt)
        
        

        
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
        self.i_start = 0
        self.line_collection = None
        
    def _next_spike(self, event):
        try:
            t = (self.i_start+self.i_window/2.)*1000./self.FS
            t_spk = self.spt[self.spt>t][0]
            self.i_start = int(np.ceil(t_spk/1000.*self.FS-self.i_window/2.))
            self.i_end = self.i_start + self.i_window
            self.window.set_scroll_pos(self.i_start)
            
            self.draw_plot()
        except IndexError:
            pass
        
    def _prev_spike(self, event):
        try:
            t = (self.i_start+self.i_window/2.)*1000./self.FS
            t_spk = self.spt[self.spt<t][-1]
            self.i_start = int(t_spk/1000.*self.FS-self.i_window/2.)
            self.i_end = self.i_start + self.i_window
            self.window.set_scroll_pos(self.i_start)
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

    def set_data(self, data, spk_idx=None, labels=None):

        self.x = data['data']
        self.FS = data['FS']
        n_chans, n_pts = self.x.shape
        if spk_idx:
            self.spt = spk_idx['data']
            if labels is not None:
                self.labels = labels
                self.color_func = label_color(np.unique(labels))
            else:
                self.labels = None
        else:
            self.spt = None
            self.ax_next.set_visible(False)
            self.ax_prev.set_visible(False)


        self.i_window = int(self.winsz/1000.*self.FS)
        # Extents of data sequence: 
        self.i_min = 0
        self.i_max = n_pts - self.i_window
        self.n_chans = n_chans


        self.window.set_scroll_max(self.i_max, self.i_window)
    

        # Indices of data interval to be plotted:
        
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
        
        if self.line_collection:
            self.line_collection.remove()

        self.line_collection = LineCollection(self.segs,
                                              offsets=None,
                                              transform=self.axes.transData)

        self.axes.add_collection(self.line_collection)
        self.axes.set_xlim((self.time[0], self.time[-1]))
        self.axes.set_ylim((self.ylims[0]+self.offsets.min(), 
                            self.ylims[1]+self.offsets.max()))

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
            if self.labels is not None:
                labs = self.labels[(self.spt>t_min) & (self.spt<t_max)]
                colors = np.repeat(self.color_func(labs), self.n_chans, 0)
                #import pdb; pdb.set_trace()
            else:
                colors = 'r'
            self.spike_collection = LineCollection(sp_segs,
                                                  offsets=None,
                                                  color=colors,
                                                  transform=self.axes.transData)
            self.axes.add_collection(self.spike_collection)
            
        
    def OnScrollEvt(self, pos):

        # Update the indices of the plot:
        self.i_start = self.i_min + pos
        self.i_end = self.i_min + self.i_window + pos
        self.draw_plot()


def browse_data(data, spk_idx=None, win=100):
    # Generate some data to plot:
    class SpikeBrowserApp(wx.App):

        def OnInit(self):
            self.frame = PlotWithScrollBarWx(parent=None,id=-1) 
            self.browser = SpikeBrowserUI(self.frame)
            self.frame.Show(True)
            self.SetTopWindow(self.frame)
            return True
    
        def init_data(self, data, spk_idx, winsz):
            self.browser.winsz = winsz
            self.browser.set_data(data, spk_idx)

    app = SpikeBrowserApp()
    app.init_data(data, spk_idx, win)
    app.MainLoop()

def browse_data_tk(data, spk_idx=None, win=100):
    
    frame = PlotWithScrollBarTk() 
    browser = SpikeBrowserUI(frame)
    browser.winsz = win
    browser.set_data(data, spk_idx)
    Tk.mainloop()
