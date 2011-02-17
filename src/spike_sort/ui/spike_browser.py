#!/usr/bin/env python
#coding=utf-8
import numpy as np
from numpy import arange, sin, pi, float, size

import matplotlib
matplotlib.use("WxAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection

import wx
import time

class MyFrame(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self,parent, id, 'scrollable plot',
                style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER,
                size=(800, 400))
        self.panel = wx.Panel(self, -1)

        self.fig = Figure((5, 4), 75)
        self.canvas = FigureCanvasWxAgg(self.panel, -1, self.fig)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, -1, wx.EXPAND)

        self.panel.SetSizer(sizer)
        self.panel.Fit()

        self.canvas.Bind(wx.EVT_SCROLLWIN, self.OnScrollEvt)

    def init_data(self, data):

        self.x = data['data']
        self.FS = data['FS']
        n_chans, n_pts = self.x.shape


        self.i_window = int(self.winsz/1000.*self.FS)
        # Extents of data sequence: 
        self.i_min = 0
        self.i_max = n_pts - self.i_window

        self.canvas.SetScrollbar(wx.HORIZONTAL, 0, 5, self.i_max)
        # Size of plot window:       

        # Indices of data interval to be plotted:
        self.i_start = 0
        self.i_end = self.i_start + self.i_window
        
        self.fig.clf()
        self.axes = self.fig.add_subplot(111)
        
        self.time = np.arange(self.i_start,self.i_end)/self.FS
        
        self.segs = np.empty((n_chans, self.i_window, 2))
        self.segs[:,:,0] = self.time[np.newaxis,:]
        self.segs[:,:,1] = self.x[:,self.i_start:self.i_end]
        ylims = (self.segs[:,:,1].min(), self.segs[:,:,1].max())
        offset = ylims[1]-ylims[0]
        self.offsets = np.arange(n_chans)*offset
        self.segs[:,:,1] += self.offsets[:,np.newaxis]
        
        self.ylims = (ylims[0], ylims[1] + offset*(n_chans-1))

        self.line_collection = LineCollection(self.segs,
                                              offsets=None,
                                              transform=self.axes.transData)

        self.axes.add_collection(self.line_collection)
        self.axes.set_xlim((self.time[0], self.time[-1]))
        self.axes.set_ylim(self.ylims)
        start = time.time()
        self.canvas.draw()

    def draw_plot(self):

        self.time = np.arange(self.i_start,self.i_end)/self.FS
        self.segs[:,:,0] = self.time[np.newaxis,:]
        self.segs[:,:,1] = self.x[:,self.i_start:self.i_end]+self.offsets[:,np.newaxis]
        self.line_collection.set_segments(self.segs)

        # Adjust plot limits:
        self.axes.set_xlim((self.time[0], self.time[-1]))
        #self.axes.set_ylim(self.ylims)

        # Redraw:                  
        self.canvas.draw()

    def OnScrollEvt(self, event):

        # Update the indices of the plot:
        self.i_start = self.i_min + event.GetPosition()
        self.i_end = self.i_min + self.i_window + event.GetPosition()
        self.draw_plot()

class MyApp(wx.App):

    def OnInit(self):
        self.frame = MyFrame(parent=None,id=-1)
        self.frame.Show(True)
        self.SetTopWindow(self.frame)
        return True

    def init_data(self, data, winsz):
        self.frame.winsz = winsz
        self.frame.init_data(data)

def plot_data(data):
    # Generate some data to plot:
    
    app = MyApp()
    app.init_data(data, 50000)
    app.MainLoop()

def main():
    dt = 0.001
    t = np.arange(1E6)*dt
    x = np.sin(2*np.pi*t)*np.sin(2*np.pi*t*0.01)
    x_col = np.repeat(x[np.newaxis,:], 4,0)
    data = {'data':x_col, "FS":1./dt}
    plot_data(data)


if __name__ == '__main__':
    import cProfile
    cProfile.run('main()',"wxprof")
    


