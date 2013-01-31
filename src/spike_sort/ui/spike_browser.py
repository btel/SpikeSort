
import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from matplotlib.widgets import Button
from matplotlib.spines import Spine

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
        elif mode == "scroll":
            delta, units = args
            delta = float(delta)
            if units == 'units':
                self.cur_pos += delta * self.page_sz / 5.
            elif units == 'pages':
                self.cur_pos += delta * self.page_sz
        self.set_scroll_pos(self.cur_pos)
        self.handler(self.cur_pos)

    def set_scroll_handler(self, handler):
        self.handler = handler

    def set_scroll_pos(self, pos):
        min, max = str(
            pos * 1. / self.max), str((pos + self.page_sz) * 1. / self.max)
        self.scrollbar.set(min, max)

    def set_scroll_max(self, max, page_size):
        self.page_sz = page_size
        self.max = max


class SpikeBrowserUI(object):
    def __init__(self, window):
        self.window = window
        self.sp_win = [-0.8, 1]
        self.spike_collection = None

        self.fig = Figure((9, 5), 75)

        self.canvas = window.get_canvas(self.fig)

        self._mpl_init()
        self.canvas.mpl_connect('key_press_event', self._zoom_key_handler)
        self.canvas.mpl_connect(
            'key_press_event', self._browse_spikes_key_handler)
        self.window.set_scroll_handler(self.OnScrollEvt)

    def _mpl_init(self):
        self.fig.clf()
        self.axes = self.fig.add_axes([0.02, 0.1, 0.96, 0.85])
        self.fancyyaxis = FancyYAxis(self, 0.05)
        self.ax_prev = self.fig.add_axes([0.8, 0.0, 0.1, 0.05])
        self.ax_next = self.fig.add_axes([0.9, 0.0, 0.1, 0.05])

        self.b_next = Button(self.ax_next, 'Next')
        self.b_prev = Button(self.ax_prev, "Prev")

        self.b_next.on_clicked(self._next_spike)
        self.b_prev.on_clicked(self._prev_spike)
        self.i_spike = 0
        self.i_start = 0
        self.line_collection = None

    def _next_spike(self, event):
        try:
            if self.i_spike < len(self.spt) - 1:
                self.i_spike += 1
            t_spk = self.spt[self.i_spike]
            i_start = int(
                np.ceil(t_spk / 1000. * self.FS - self.i_window / 2.))
            i_start = np.maximum(self.i_min, i_start)
            i_start = np.minimum(self.i_max, i_start)
            self.i_start = i_start
            self.i_end = self.i_start + self.i_window
            self.window.set_scroll_pos(self.i_start)
            self.draw_plot()
        except IndexError:
            pass

    def _prev_spike(self, event):
        try:
            if self.i_spike > 0:
                self.i_spike -= 1
            t_spk = self.spt[self.i_spike]
            i_start = int(
                np.ceil(t_spk / 1000. * self.FS - self.i_window / 2.))
            i_start = np.maximum(self.i_min, i_start)
            i_start = np.minimum(self.i_max, i_start)
            self.i_start = i_start
            self.i_end = self.i_start + self.i_window
            self.window.set_scroll_pos(self.i_start)
            self.draw_plot()
        except IndexError:
            pass

    def _zoom_key_handler(self, event):
        if event.key == '+' or event.key == '=':
            self.scale_y(0.5)
        elif event.key == '-':
            self.scale_y(2)
        elif event.key == 'ctrl++' or event.key == 'ctrl+=':
            self.scale_x(0.5)
        elif event.key == 'ctrl+-':
            self.scale_x(2)

        self.draw_plot()

    def _browse_spikes_key_handler(self, event):
        if event.key == 'right':
            self._next_spike(event)
        elif event.key == 'left':
            self._prev_spike(event)
        else:
            return

    def set_spiketimes(self, spk_idx, labels=None, all_labels=None):
        if spk_idx:
            self.spt = spk_idx['data']
            if labels is not None:
                self.labels = labels
                if all_labels is None:
                    self.color_func = label_color(np.unique(labels))
                else:
                    self.color_func = label_color(all_labels)
            else:
                self.labels = None

            self.ax_next.set_visible(True)
            self.ax_prev.set_visible(True)

            self.update_i_sipke()

        else:
            self.spt = None
            self.ax_next.set_visible(False)
            self.ax_prev.set_visible(False)

    def set_data(self, data):
        self.x = data['data']
        self.FS = data['FS']
        n_chans, n_pts = self.x.shape

        # reset spike times data/hide buttons
        self.set_spiketimes(None)

        self.i_window = int(self.winsz / 1000. * self.FS)
        # Extents of data sequence:
        self.i_min = 0
        self.i_max = n_pts - self.i_window
        self.n_chans = n_chans

        self.window.set_scroll_max(self.i_max, self.i_window)

        # Indices of data interval to be plotted:
        self.i_end = self.i_start + self.i_window
        curr_slice = self.x[:, self.i_start:self.i_end]
        ylims = (curr_slice.min(), curr_slice.max())
        offset = ylims[1] - ylims[0]

        self.ylims = np.array(ylims)
        self.offsets = np.arange(n_chans) * offset
        
        # will be filled in draw_plot
        self.segs = np.empty((n_chans, self.i_window, 2))

        if self.line_collection:
            self.line_collection.remove()

        self.line_collection = LineCollection(self.segs,
                                              offsets=None,
                                              transform=self.axes.transData,
                                              color='k')

        self.axes.add_collection(self.line_collection)
        self.fancyyaxis.reset()

        self.draw_plot()

    def draw_plot(self):
        self.time = np.arange(self.i_start, self.i_end) * 1. / self.FS
        self.segs[:, :, 0] = self.time[np.newaxis, :]
        y_signal = self.x[:, self.i_start:self.i_end]
        y_signal = y_signal - np.mean(y_signal, 1)[:, None]
        self.segs[:, :, 1] = y_signal + self.offsets[:, np.newaxis]
        self.line_collection.set_segments(self.segs)

        # Adjust plot limits:
        self.axes.set_xlim((self.time[0], self.time[-1]))

        ygap = np.max(np.abs(self.ylims))
        self.axes.set_ylim((- ygap + self.offsets.min(),
                            ygap + self.offsets.max()))
        self.fancyyaxis.update()

        if self.spt is not None:
            self.draw_spikes()
        # Redraw:
        self.canvas.draw()


    def draw_spikes(self):
        if self.spike_collection is not None:
            self.spike_collection.remove()
            self.spike_collection = None
        sp_win = self.sp_win
        time = self.segs[0, :, 0] * 1000.
        t_min, t_max = time[0] - sp_win[0], time[-1] - sp_win[1]
        spt = self.spt[(self.spt > t_min) & (self.spt < t_max)]
        if len(spt) > 0:
            n_pts = int((sp_win[1] - sp_win[0]) / 1000. * self.FS)
            sp_segs = np.empty((len(spt), self.n_chans, n_pts, 2))
            for i in range(len(spt)):
                start, = np.nonzero(time >= (spt[i] + sp_win[0]))
                start = start[0]
                stop = start + n_pts
                sp_segs[i, :, :, 0] = (time[np.newaxis, start:stop] / 1000.)
                sp_segs[i, :, :, 1] = self.segs[:, start:stop, 1]
            sp_segs = sp_segs.reshape(-1, n_pts, 2)
            if self.labels is not None:
                labs = self.labels[(self.spt > t_min) & (self.spt < t_max)]
                colors = np.repeat(self.color_func(labs), self.n_chans, 0)
            else:
                colors = 'r'
            self.spike_collection = LineCollection(sp_segs,
                                                   offsets=None,
                                                   color=colors,
                                                   transform=self.axes.transData)
            self.axes.add_collection(self.spike_collection)

    def scale_y(self, factor):
        self.ylims *= factor
        offset = self.ylims[1] - self.ylims[0]
        self.offsets = np.arange(self.n_chans) * offset

    def scale_x(self, factor):
        i_center = self.i_start + self.i_window / 2
        self.i_window = int(self.i_window * factor)
        self.i_start = i_center - self.i_window / 2
        self.i_start = self.i_start >= 0 and self.i_start or 0
        self.i_end = self.i_start + self.i_window

        self.segs = np.empty((self.n_chans, self.i_window, 2))


    def OnScrollEvt(self, pos):

        # Update the indices of the plot:
        self.i_start = self.i_min + pos
        self.i_end = self.i_min + self.i_window + pos

        self.update_i_sipke()
        self.draw_plot()

    def update_i_sipke(self):
        '''
        Finds a spike index which is near or inside the current data
        window (between i_start and i_end). The i_spike variable is then
        updated with this index.
        '''

        t_center = (self.i_start + self.i_window / 2.) * 1000. / self.FS
        idx, = np.where(self.spt < t_center)
        if len(idx) > 0:
            self.i_spike = idx[-1]
        else:
            self.i_spike = 0


class FancyYAxis:
    def __init__(self, browser,
            data_xoffset = 0.1,
            spine_xoffset = 0.7,
            text_xoffset = 0.4):
        self.b = browser
        self.data_xoffset = data_xoffset
        self.spine_xoffset = spine_xoffset
        self.text_xoffset = text_xoffset
        self.annotations = {}

    def reset(self):
        """Initialize or Reset"""
        # delete references to previous annotations
        for obj in self.annotations.values():
            obj.remove()

        spine_pos = self.data_xoffset * self.spine_xoffset

        for chan in range(self.b.n_chans):
            # add spines for channel unit bars
            sp = Spine(self.b.axes, 'left', self.b.axes.spines['left']._path)
            sp.set_position(('axes', spine_pos))
            self.b.axes.spines['ch%s' % chan] = sp

            # add channel names
            self.annotations['ch%s' % chan] = self.b.axes.text(
                    0, 0,  # dummy location
                    'Ch%s' % chan,
                    ha='center',
                    va='center',
                    rotation='vertical')

        # hide unneeded spines and move ticks
        self.b.axes.spines['left'].set_color('none')
        self.b.axes.spines['right'].set_color('none')
        self.b.axes.spines['left'].set_position(('axes', spine_pos))

    def update(self):
        """Update"""

        # initialize if needed
        if not self.annotations:
            self.reset()

        # create space for axis drawing
        xlim = self.b.axes.get_xlim()
        xrang = xlim[1] - xlim[0]
        self.b.axes.set_xlim((xlim[0] - xrang * self.data_xoffset, xlim[1]))

        # compute ticks and labels for each bar
        yrange = self.b.ylims[1]-self.b.ylims[0]
        half_range = self.closest_nice_float(yrange * 0.4)
        rel_ticks = np.array([-half_range, 0, half_range])
        ticks = rel_ticks[np.newaxis, :] + self.b.offsets[:, np.newaxis]

        labels = []
        text_offset = xlim[0] \
                - xrang * self.data_xoffset * (1 - self.text_xoffset)

        for chan in range(self.b.n_chans):
            # update vertical ranges and labels for bars
            labels.extend([str(rel_ticks[0]), "", str(rel_ticks[2])])
            self.b.axes.spines['ch%s' % chan].set_bounds(ticks[chan, 0],
                    ticks[chan, -1])

            # update horizontal text position
            self.annotations['ch%s' % chan].set_position((text_offset,
                    self.b.offsets[chan]))

        self.b.axes.set_yticks(ticks.flatten())
        self.b.axes.set_yticklabels(labels)

        # put Y tick labels inside the plot
        for tick in self.b.axes.yaxis.get_major_ticks():
            tick.set_pad(-6)
            tick.label.set_horizontalalignment('left')

    @staticmethod
    def closest_nice_float(value, proximity = 2):
        """
        Find closest value to `value`, who's mantissa has `proximity` digits
        after comma
        """
        exp = int(np.log10(value))
        mnt = value / 10 ** exp
        return round(mnt, proximity) * 10 ** exp


def browse_data_tk(data, spk_idx=None, labels=None, win=100):

    frame = PlotWithScrollBarTk()
    browser = SpikeBrowserUI(frame)
    browser.winsz = win
    browser.set_data(data)
    browser.set_spiketimes(spk_idx, labels)
    Tk.mainloop()


browse_data = browse_data_tk
