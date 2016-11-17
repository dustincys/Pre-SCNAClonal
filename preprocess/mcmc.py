#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: mcmc.py
#          Desc: The MCMC model for linear regression
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-31 13:15:44
#       History:
# =============================================================================
'''
import sys

import pymc
import numpy as np
from scipy.stats import gaussian_kde

import constants

from matplotlib import pyplot as plt


class MCMCLM(object):

    """The MCMC model for linear regression, return the slope and inlier"""

    def __init__(self, data, n, percentile, prob_threshold):
        """Initialize the MCMCLM model

        :data: the segment data object
        :n: the sampling number
        :percentile: the percentile for calculating standard
            divation for inlier and outlier
        :prob_threshold: the prob_threshold for determine inliner and outliner

        """
        self._data = data
        self._n = n
        self._percentile = percentile
        self._prob_threshold = prob_threshold

    def _getSlopeMC(self, y, x):
        x_30 = np.percentile(x, 30)
        x_70 = np.percentile(x, 70)
        y_up = y[x > x_70]
        y_down = y[x < x_30]

        y_up_y = np.percentile(y_up, 50)
        x_up_x = np.percentile(x, 85)

        y_down_y = np.percentile(y_down, 50)
        x_down_x = np.percentile(x, 15)

        m = (y_up_y - y_down_y) * 1.0 / (x_up_x - x_down_x)
        c = y_down_y - m * x_down_x

        return m, c

    def _correctY(self, y, x, slope, intercept):
        K = np.percentile(y, 50)
        A = slope * x + intercept
        return y - A + K

    def run(self):
        x, y_with_outlier = self._getSampledData()
        m, c = self._getSlopeMC(y_with_outlier, x)
        print "____>>> getSlopeMC: m, c____"
        print m, c
        print "_________end getSlopeMC:m, c______________"
        slope = pymc.Uniform('slope', m-5.0, m+5.0)

        def log_posterior_likelihood_of_slope(
                y_with_outlier,
                slope
                ):
            y_corrected = self._correctY(y_with_outlier, x, slope, 0)
            y_density = gaussian_kde(y_corrected)

            y_30 = np.percentile(y_with_outlier, 30)
            y_70 = np.percentile(y_with_outlier, 70)
            y_xs = np.linspace(y_30, y_70, 20)

            y_ys = y_density(y_xs)
            index = np.argmax(y_ys)

            prob = y_ys[index]

            return prob

        slope_distribution = pymc.stochastic_from_dist(
            'slope_distribution',
            logp=log_posterior_likelihood_of_slope,
            dtype=np.float,
            mv=True)

        slope_dist = slope_distribution('slope_dist',
                                            slope=slope,
                                            observed=True,
                                            value=y_with_outlier)

        model = dict(slope_dist=slope_dist,
                     slope=slope
                     )

        M = pymc.MAP(model)
        M.fit()
        print "M.slope.value = {}".format(M.slope.value)


    def _getSampledData(self):
        if 0 == self._n:
            sampledSegs = self._data.segments
        else:
            sampledSegs = np.random.choice(self._data.segments, self._n)

        print "all sample: {}".format(len(sampledSegs))

        x0 = np.array(map(lambda seg: seg.gc, sampledSegs))
        y0 = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                          np.log(seg.normal_reads_num + 1), sampledSegs))
        l = sorted(zip(y0, x0), reverse=True)
        y0, x0 = [list(t) for t in zip(*l)]

        return np.array(x0), np.array(y0)

    def getInterceptParameters(self, y, x):
        inlier_proportion_prior = constants.INLIER_PROPORTION_PRIOR
        inlier_left_margin = constants.INLIER_LEFT_MARGIN * 100

        y_temp = y[x < inlier_left_margin]
        mu = np.percentile(y_temp, inlier_proportion_prior)
        rho = np.std(y_temp)

        return mu, rho
