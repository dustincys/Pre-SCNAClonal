#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: plotBaseline.py
#          Desc: plot baseline and interactively select the baseline points
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-20 21:51:40
#       History:
# =============================================================================
'''

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, LassoSelector
from matplotlib.path import Path
from matplotlib.colors import colorConverter

import numpy  as np

class BaselinePlot():

    """plot baseline interactively and select points as baseline"""

    def __init__(self, data, max_copynumber, subclone_num):
        """load data

        :data: TODO

        """

        self.data = data
        self.segments = self.data.segments

        self.max_copynumber = max_copynumber
        self.subclone_num = subclone_num

        self.color_back = colorConverter.to_rgba('black')
        self.color_NLOH_NAPM = colorConverter.to_rgba('blue')
        self.color_LOH_APM = colorConverter.to_rgba('green')
        self.color_NLOH_APM = colorConverter.to_rgba('green')
        self.color_baseline = colorConverter.to_rgba('red')

        self.area_back = 10
        self.area_NLOH_NAPM = 20
        self.area_LOH_APM = 20
        self.area_NLOH_APM = 20
        self.area_baseline = 20

        self.alpha_back = 10
        self.alpha_NLOH_NAPM = 20
        self.alpha_LOH_APM = 20
        self.alpha_NLOH_APM = 20
        self.alpha_baseline = 20

        self.LOH_THRED = 0.025
        self.APM_THRED = 0.18


    def plot(self):
        """
        :returns: TODO

        """

        segs_back = self._load_Background()
        segs_NLOH_NAPM = self._load_NLOH_NAPM()
        segs_LOH_APM = self._load_LOH_APM()

        fig, ax = plt.subplots()
        plt.subplots_adjust(left = 0.4, bottom=0.35)

        # LOH and NAPM
        x_back, y_back = self._getXYSegs(segs_back)
        pts_back = ax.scatter(x_back, y_back, s=self.area_back,
                              alpha=self.alpha_back, color=self.color_back)

        # NLOH and NAPM
        x_NLOH_NAPM, y_NLOH_NAPM = self._getXYSegs(segs_NLOH_NAPM)
        pts_NLOH_NAPM = ax.scatter(x_NLOH_NAPM, y_NLOH_NAPM,
                                   s=self.area_NLOH_NAPM,
                                   alpha=self.alpha_NLOH_NAPM,
                                   color=self.color_NLOH_NAPM)

        # LOH and APM
        x_LOH_APM, y_LOH_APM = self._getXYSegs(segs_LOH_APM)
        pts_LOH_APM = ax.scatter(x_LOH_APM, y_LOH_APM,
                                  s=self.area_LOH_APM,
                                  alpha=self.alpha_LOH_APM,
                                  color=self.color_LOH_APM)

        # NLOH and APM
        x_NLOH_APM, y_NLOH_APM = self._getXYSegs(segs_NLOH_APM)


        color_NLOH_APM = np.array([self.color_NLOH_APM if
                                   item.baseline_label == 'FALSE' else
                                   self.color_baseline for item in
                                   segs_NLOH_APM])
        alpha_NLOH_APM = np.array([self.alpha_NLOH_APM if
                                   item.baseline_label == 'FALSE' else
                                   self.alpha_baseline for item in
                                   segs_NLOH_APM])
        area_NLOH_APM = np.array([self.area_NLOH_APM if
                                   item.baseline_label == 'FALSE' else
                                   self.area_baseline for item in
                                   segs_NLOH_APM])

        pts_NLOH_APM = ax.scatter(x_NLOH_APM, y_NLOH_APM,
                                  s=area_NLOH_APM,
                                  alpha=alpha_NLOH_APM,
                                  color=color_NLOH_APM)

        axalpha_back = plt.axes([0.1, 0.8, 0.4, 0.03], axisbg=axcolor)
        axarea_back = plt.axes([0.1, 0.75, 0.4, 0.03], axisbg=axcolor)

        axalpha_NN = plt.axes([0.1, 0.7, 0.4, 0.03], axisbg=axcolor)
        axarea_NN = plt.axes([0.1, 0.4, 0.4, 0.03], axisbg=axcolor)

        axalpha_NT = plt.axes([0.1, 0.6, 0.4, 0.03], axisbg=axcolor)
        axarea_NT = plt.axes([0.1, 0.55, 0.4, 0.03], axisbg=axcolor)

        axalpha_TT = plt.axes([0.1, 0.5, 0.4, 0.03], axisbg=axcolor)
        axarea_TT = plt.axes([0.1, 0.45, 0.4, 0.03], axisbg=axcolor)

        axalpha_baseline = plt.axes([0.1, 0.4, 0.4, 0.03], axisbg=axcolor)
        axarea_baseline = plt.axes([0.1, 0.35, 0.4, 0.03], axisbg=axcolor)


        #salpha_TN = Slider(axalpha, 'alpha', 0, 0.8, valinit=self.alpha0)
        salpha_back = Slider(axalpha_back, 'alpha', 0, 0.8, valinit=self.alpha_back)
        salpha_NN = Slider(axalpha_NN, 'alpha', 0, 0.8, valinit=self.alpha_NN)
        salpha_NT = Slider(axalpha_NT, 'alpha', 0, 0.8, valinit=self.alpha_NT)
        salpha_TT = Slider(axalpha_TT, 'alpha', 0, 0.8, valinit=self.alpha_TT)
        salpha_baseline = Slider(axalpha_baseline, 'alpha', 0, 0.8,
                                 valinit=self.alpha_baseline)
        sarea_back = Slider(axarea_back, 'area', 1, 50, valinit=self.area_back)
        sarea_NN = Slider(axarea_NN, 'area', 1, 50, valinit=self.area_NN)
        sarea_NT = Slider(axarea_NT, 'area', 1, 50, valinit=self.area_NT)
        sarea_TT = Slider(axarea_TT, 'area', 1, 50, valinit=self.area_TT)
        sarea_baseline = Slider(axarea_baseline, 'area', 1, 50,
                                valinit=self.area_baseline)

        def update_alpha_back(val):
            self.alpha_back = salpha_back.val
            pts_back.set_alpha(self.alpha_back)
            fig.canvas.draw_idle()

        def update_alpha_NN(val):
            self.alpha_NLOH_NAPM = salpha_NN.val
            pts_NLOH_NAPM.set_alpha(self.alpha_NLOH_NAPM)
            fig.canvas.draw_idle()

        def update_alpha_NT(val):
            self.alpha_NLOH_APM = salpha_NT.val

            alpha_NLOH_APM = np.array([self.alpha_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.alpha_baseline for item in
                                    segs_NLOH_APM  ])

            pts_NLOH_APM.set_alpha(alpha_LOH_APM)
            fig.canvas.draw_idle()

        def update_alpha_TT(val):
            self.alpha_LOH_APM = salpha_TT.val
            pts_LOH_APM.set_alpha(self.alpha_LOH_APM)
            fig.canvas.draw_idle()

        def update_alpha_baseline(val):
# Here it might be wrong, the index of alpha_LOH_APM may be not right
            self.alpha_baseline = salpha_baseline.val

            alpha_NLOH_APM = np.array([self.alpha_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.alpha_baseline for item in
                                    segs_NLOH_APM  ])

            pts_NLOH_APM.set_alpha(alpha_LOH_APM)
            fig.canvas.draw_idle()

        def update_area_back(val):
            self.area_back = sarea_back.val
            pts_back.set_area(self.area_back)
            fig.canvas.draw_idle()

        def update_area_NN(val):
            self.area_NLOH_NAPM = sarea_NN.val
            pts_NLOH_NAPM.set_area(self.area_NLOH_NAPM)
            fig.canvas.draw_idle()

        def update_area_NT(val):
            self.area_NLOH_APM = sarea_NT.val

            area_NLOH_APM = np.array([self.area_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.area_baseline for item in
                                    segs_NLOH_APM  ])

            pts_baseline.set_area(area_NLOH_APM)
            fig.canvas.draw_idle()

        def update_area_TT(val):
            self.area_LOH_APM = sarea_TT.val
            pts_LOH_APM.set_area(self.area_LOH_APM)
            fig.canvas.draw_idle()

        def update_area_baseline(val):
            self.area_baseline = sarea_baseline.val

            area_NLOH_APM = np.array([self.area_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.area_baseline for item in
                                    segs_NLOH_APM  ])
            pts_baseline.set_area(area_NLOH_APM)
            fig.canvas.draw_idle()

        salpha_back.on_changed(update_alpha_back)
        salpha_NN.on_changed(update_alpha_NN)
        salpha_NT.on_changed(update_alpha_NT)
        salpha_TT.on_changed(update_alpha_TT)
        salpha_baseline.on_changed(update_alpha_baseline)

        sarea_back.on_changed(update_area_back)
        sarea_NN.on_changed(update_area_NN)
        sarea_NT.on_changed(update_area_NT)
        sarea_TT.on_changed(update_area_TT)
        sarea_baseline.on_changed(update_area_baseline)


        axLOH_t = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        axAPM_t = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
        sLOH_t = Slider(axLOH_t, 'LOH', 0, 1, valinit=self.LOH_THRED)
        sAPM_t = Slider(axAPM_t, 'APM', 0, 1, valinit=self.APM_THRED)

        def update_LOH_t(val):
            self.LOH_THRED = sLOH_t.val

            self.data.update_LOHAPM_t(self.LOH_THRED, self.APM_THRED,
                                      self.max_copynumber,
                                      self.subclone_num)

            segs_back = self._load_Background()
            segs_NLOH_NAPM = self._load_NLOH_NAPM()
            segs_NLOH_APM = self._load_NLOH_APM()

            x_back, y_back = self._getXYSegs(segs_back)
            pts_back.set_offsets(zip(x_back, y_back))

            # NLOH and NAPM
            x_NLOH_NAPM, y_NLOH_NAPM = self._getXYSegs(segs_NLOH_NAPM)
            pts_NLOH_NAPM.set_offsets(zip(x_NLOH_NAPM, y_NLOH_NAPM))

            # LOH and APM
            x_LOH_APM, y_LOH_APM = self._getXYSegs(segs_LOH_APM)
            pts_LOH_APM.set_offsets(zip(x_LOH_APM, y_LOH_APM))

            # NLOH and APM
            x_NLOH_APM, y_NLOH_APM = self._getXYSegs(segs_NLOH_APM)
            color_NLOH_APM = np.array([self.color_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.color_baseline for item in
                                    segs_NLOH_APM  ])
            alpha_NLOH_APM = np.array([self.alpha_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.alpha_baseline for item in
                                    segs_NLOH_APM  ])
            area_NLOH_APM = np.array([self.area_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.area_baseline for item in
                                    segs_NLOH_APM  ])

            pts_NLOH_APM.set_offsets(zip(x_NLOH_APM, y_NLOH_APM))
            pts_NLOH_APM.set_alpha(alpha_NLOH_APM)
            pts_NLOH_APM.set_area(area_NLOH_APM)
            fig.canvas.draw_idle()
            pass

        def update_APM_t(val):
            self.APM_THRED = sAPM_t.val

            self.data.update_LOHAPM_t(self.LOH_THRED, self.APM_THRED,
                                      self.max_copynumber,
                                      self.subclone_num)

            segs_back = self._load_Background()
            segs_NLOH_NAPM = self._load_NLOH_NAPM()
            segs_NLOH_APM = self._load_NLOH_APM()

            x_back, y_back = self._getXYSegs(segs_back)
            pts_back.set_offsets(zip(x_back, y_back))

            # NLOH and NAPM
            x_NLOH_NAPM, y_NLOH_NAPM = self._getXYSegs(segs_NLOH_NAPM)
            pts_NLOH_NAPM.set_offsets(zip(x_NLOH_NAPM, y_NLOH_NAPM))

            # LOH and APM
            x_LOH_APM, y_LOH_APM = self._getXYSegs(segs_LOH_APM)
            pts_LOH_APM.set_offsets(zip(x_LOH_APM, y_LOH_APM))

            # NLOH and APM
            x_NLOH_APM, y_NLOH_APM = self._getXYSegs(segs_NLOH_APM)
            color_NLOH_APM = np.array([self.color_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.color_baseline for item in
                                    segs_NLOH_APM  ])
            alpha_NLOH_APM = np.array([self.alpha_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.alpha_baseline for item in
                                    segs_NLOH_APM  ])
            area_NLOH_APM = np.array([self.area_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.area_baseline for item in
                                    segs_NLOH_APM  ])

            pts_NLOH_APM.set_offsets(zip(x_NLOH_APM, y_NLOH_APM))
            pts_NLOH_APM.set_alpha(alpha_NLOH_APM)
            pts_NLOH_APM.set_area(area_NLOH_APM)

            fig.canvas.draw_idle()
            pass

        sLOH_t.on_changed(update_LOH_t)
        sAPM_t.on_changed(update_APM_t)

        def onselect(verts):
            path = Path(verts)
            xys = pts_baseline.get_offsets()

            ind = np.nonzero([path.contains_point(xy) for xy in xys])[0]
            self.x = np.array(select_col(xys[ind], 0))
            self.y = np.array(select_col(xys[ind], 1))

            self.colorin = colorConverter.to_rgba('red', self.alpha0)
            self.colorout = colorConverter.to_rgba('blue', self.alpha0)
            facecolors = np.array([self.colorout for xy in xys])
            facecolors[ind] = self.colorin
            pts.set_facecolors(facecolors)

            fig.canvas.draw_idle()
            pass

        lasso = LassoSelector(ax, onselect)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            salpha_NN.reset()
            salpha_NT.reset()
            salpha_TT.reset()
            salpha_back.reset()
            salpha_baseline.reset()
            sLOH_t.reset()
            sAPM_t.reset()

            segs_back = self._load_Background()
            segs_NLOH_NAPM = self._load_NLOH_NAPM()
            segs_NLOH_APM = self._load_NLOH_APM()

            x_back, y_back = self._getXYSegs(segs_back)
            pts_back.set_offsets(zip(x_back, y_back))

            # NLOH and NAPM
            x_NLOH_NAPM, y_NLOH_NAPM = self._getXYSegs(segs_NLOH_NAPM)
            pts_NLOH_NAPM.set_offsets(zip(x_NLOH_NAPM, y_NLOH_NAPM))

            # LOH and APM
            x_LOH_APM, y_LOH_APM = self._getXYSegs(segs_LOH_APM)
            pts_LOH_APM.set_offsets(zip(x_LOH_APM, y_LOH_APM))

            # NLOH and APM
            x_NLOH_APM, y_NLOH_APM = self._getXYSegs(segs_NLOH_APM)
            color_NLOH_APM = np.array([self.color_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.color_baseline for item in
                                    segs_NLOH_APM  ])
            alpha_NLOH_APM = np.array([self.alpha_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.alpha_baseline for item in
                                    segs_NLOH_APM  ])
            area_NLOH_APM = np.array([self.area_NLOH_APM if
                                    item.baseline_label == 'FALSE' else
                                    self.area_baseline for item in
                                    segs_NLOH_APM  ])

            pts_NLOH_APM.set_offsets(zip(x_NLOH_APM, y_NLOH_APM))
            pts_NLOH_APM.set_alpha(alpha_NLOH_APM)
            pts_NLOH_APM.set_area(area_NLOH_APM)

        button.on_clicked(reset)

        resetax = plt.axes([0.5, 0.025, 0.1, 0.04])
        button_exit = Button(resetax, 'Ok', color=axcolor, hovercolor='0.975')

        def ok_exit(event):
            plt.close()

        button_exit.on_clicked(ok_exit)
        plt.show()

    def _getXYSegs(self, segs):

        x = np.array(map(lambda seg: seg.gc, segs))
        y = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                         np.log(seg.normal_reads_num + 1), segs))

        return x, y

    def _load_Background(self):
        """
        :returns: return the sampled background segments, in order to speed up

        """
        back_segs = filter(lambda seg: seg.LOH_status == 'TRUE' and
                           seg.APM_status == 'FALSE', self.segments)

        sampled_back_segs = np.random.choice(back_segs, self.n)

        return sampled_back_segs

    def _load_NLOH_NAPM(self):
        """TODO: Docstring for _loadNLOHAPM.

        :arg1: TODO
        :returns: TODO

        """
        return filter(lambda seg: seg.LOH_status == 'FALSE' and
                      seg.APM_status == 'FALSE', self.segments)

    def _load_LOH_APM(self):
        """TODO: Docstring for _loadNLOHAPM.

        :arg1: TODO
        :returns: TODO

        """
        return filter(lambda seg: seg.LOH_status == 'TRUE' and
                      seg.APM_status == 'TRUE', self.segments)

    def _load_NLOH_APM(self):
        """TODO: Docstring for _loadNLOHAPM.

        :arg1: TODO
        :returns: TODO

        """
        return filter(lambda seg: seg.LOH_status == 'FALSE' and
                      seg.APM_status == 'TRUE', self.segments)

    def _load_Baseline(self):
        """
        :returns: TODO

        """
        return filter(lambda seg: seg.baseline_label == 'TRUE', self.segments)




def main():
    """
    :returns:

    """
    import data



if __name__ == "__main__":
    main()
