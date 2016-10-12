#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: plotGC.py
#          Desc: plot gc stripes and interactively adjust the lm
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-11 15:55:22
#       History:
# =============================================================================
'''


import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, LassoSelector
from matplotlib.path import Path

import numpy as np


class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool highlights
    selected points by fading them out (i.e., reducing their alpha values).
    If your collection has alpha < 1, this tool will permanently alter them.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    """

    def __init__(self, ax, collection, alpha_other=0.001):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, self.Npts).reshape(self.Npts, -1)

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero([path.contains_point(xy) for xy in self.xys])[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 0.2
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 0.5
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


class GCStripePlot():

    """plot the gc stripe and interactively adjust the lm line"""

    def __init__(self, segments, n):
        """init
        segments:  all the segments of data, need to be sampled
        n:  the sample number

        return none
        """

        self.segments = segments
        self.n = n

        # parameters for plot
        self.m0  = 0
        self.c0  = 0
        self.alpha0 = 0.02
        self.area0 = 10

        self.m = 0
        self.c = 0
        self.alpha = 1
        self.area = 10


    def _selectSamples(self, x, y):
        """
        :returns: TODO

        """

        subplot_kw = dict(
            xlim=(
                min(x), max(x)), ylim=(
                min(y), max(y)), autoscale_on=False)
        fig, ax = plt.subplots(subplot_kw=subplot_kw)

        pts = ax.scatter(x, y, s=10)
        selector = SelectFromCollection(ax, pts, alpha_other=0.1)

        plt.draw()
        selector.disconnect()

    def plot(self):
        """plot the graph
        :returns:


        """

        sampledSegs = np.random.choice(self.segments, self.n)
        x = np.array(map(lambda seg: seg.gc, sampledSegs))
        y = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                         np.log(seg.normal_reads_num + 1), sampledSegs))

        fig, ax = plt.subplots()

        pts = ax.scatter(x, y, s=self.area0, alpha=self.alpha0)
        plt.subplots_adjust(left=0.25, bottom=0.35)

        A = np.vstack([x, np.ones(len(x))]).T
        self.m0, self.c0 = np.linalg.lstsq(A, y)[0]
        fl, = ax.plot(x, self.m0*x + self.c0, 'r', label='Fitted line')
        hl = ax.axhline(y = np.median(y), linewidth=1, color='black')




        axcolor = 'lightgoldenrodyellow'
        axm = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        axc = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
        axalpha = plt.axes([0.25, 0.2, 0.65, 0.03], axisbg=axcolor)
        axarea = plt.axes([0.25, 0.25, 0.65, 0.03], axisbg=axcolor)

        sm = Slider(axm, 'slope', -2, 2, valinit=self.m0)
        sc = Slider(axc, 'interception', min(y), max(y), valinit=self.c0)
        salpha = Slider(axalpha, 'alpha', 0, 0.8, valinit=self.alpha0)
        sarea = Slider(axarea, 'area', 1, 50, valinit=self.area0)

        def update_m(val):
            """
            :returns: TODO

            """
            self.m = sm.val
            fl.set_ydata(self.m*x + self.c)
            fig.canvas.draw_idle()

        def update_c(val):
            """
            :returns: TODO

            """
            self.c = sc.val
            fl.set_ydata(self.m*x + self.c)
            fig.canvas.draw_idle()

        def update_alpha(val):
            """
            :returns: TODO

            """
            self.alpha = salpha.val
            pts.set_alpha(self.alpha)
            fig.canvas.draw_idle()

        def update_area(val):
            """
            :returns: TODO

            """
            self.area = sarea.val
            pts.set_sizes(np.ones(len(x)) * self.area)
            fig.canvas.draw_idle()


        sm.on_changed(update_m)
        sc.on_changed(update_c)
        salpha.on_changed(update_alpha)
        sarea.on_changed(update_area)


        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            sm.reset()
            sc.reset()
            salpha.reset()
            sarea.reset()
            plt.ion()

        button.on_clicked(reset)


        def onselect(verts):
            print "lasso selected!"
            fig.canvas.draw_idle()

        lasso = LassoSelector(ax, onselect)

        plt.show()


def main():
    """
    :returns: TODO

    """
    from data import Segment

    segments = []

    num = 200000

    for i in range(num):
        seg = Segment()
        seg.tumor_reads_num = np.random.randint(200, 300)
        seg.normal_reads_num = np.random.randint(300, 400)
        seg.gc = np.random.rand()
        segments.append(seg)

    gsp = GCStripePlot(segments, 20000)
    gsp.plot()


if __name__ == "__main__":
    main()
