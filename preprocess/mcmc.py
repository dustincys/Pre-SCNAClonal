#!/usr/bin/env python
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

import pymc
import numpy as np

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

    def run(self):
        x, y_with_outlier = self._getSampledData()

############## test data ###################################
        #a = -2
        #x1 = np.arange(10, 50, .5)
# Add so#me noise
        #x1 += np.random.normal(size=x1.size)

        #y1 = a * x1 + 10
        #y2 = a * x1 + 22
        #y3 = a * x1 + 24
        #y4 = a * x1 + 34
        #y5 = a * x1 + 54

        #x1 = np.array(x1)
        #y1 = np.array(y1)
        #y2 = np.array(y2)
        #y3 = np.array(y3)

        #x = np.concatenate((x1, x1, x1, x1, x1), axis=0)
        #y_with_outlier = np.concatenate((y1, y2, y3, y4, y5), axis=0)

# Add so#me outliers
        # for ii in np.arange(len(x)/7, len(x), len(x)/7.):
        #    y_with_outlier[ii] = 20*(np.random.random() + 1)

        #y_with_outlier = np.asarray(y_with_outlier)

        #down, up = np.percentile(y_with_outlier, self._percentile)
        # if down > up:
        #    down, up = up, down
        # ystd = y_with_outlier[
        #    np.logical_and(
        #        y_with_outlier > down,
        #        y_with_outlier < up)].std()
############################################################

        slope = pymc.Uniform('slope', -5.0, 5.0)
        intercept = pymc.Uniform('intercept', 0.0, 1.0)

        @pymc.deterministic
        def inlier(
                slope=slope,
                intercept=intercept,
                y_with_outlier=y_with_outlier,
                x=x):

            inlier = np.empty(len(y_with_outlier), dtype=bool)
            for i in range(len(y_with_outlier)):
                if y_with_outlier[i] < (x[i] * slope + intercept):
                    inlier[i] = True
                else:
                    inlier[i] = False

            return inlier

        def getLinalgSlope(y, x):
            if len(y) == 0:
                return float("Inf")
            A = np.vstack([x, np.ones(len(x))]).T
            slope, intercept = np.linalg.lstsq(A, y)[0]

            return slope

        def log_posterior_likelihood_of_outlier(
                y_with_outlier,
                inlier,
                slope,
                intercept):

            # Here requires multiple

            tempSlopes = getLinalgSlope(y_with_outlier[inlier], x[inlier])



            #inlier_posterior = np.log(
            #    1.0 / (np.sum((tempSlopes - slope) ** 2) + 1))

            inlier_posterior = pymc.distributions.normal_like(tempSlopes,
                                                                 slope,
                                                                 100)

            outlier_posterior = np.log(sum(inlier)) - np.log(len(inlier))

            return inlier_posterior + outlier_posterior

        outlier_distribution = pymc.stochastic_from_dist(
            'outlier_distribution',
            logp=log_posterior_likelihood_of_outlier,
            dtype=np.float,
            mv=True)

        outlier_dist = outlier_distribution('outlier_dist',
                                            inlier=inlier,
                                            slope=slope,
                                            intercept=intercept,
                                            observed=True,
                                            value=y_with_outlier)

        model = dict(outlier_dist=outlier_dist,
                     slope=slope,
                     inlier=inlier,
                     intercept=intercept
                     )

        mcmc = pymc.MCMC(model)
        mcmc.sample(100000, 20000)

        slope_trace = mcmc.trace('slope')[:]
        intercept_trace = mcmc.trace('intercept')[:]

        inlier_trace = mcmc.trace('inlier')[:]
        probability_of_points = inlier_trace.astype(float).mean(0)

        outlier_x = x[probability_of_points < self._prob_threshold]
        outlier_y = y_with_outlier[probability_of_points < self._prob_threshold]

        fix, ax = plt.subplots(figsize=(8, 8))

        ax = plt.subplot(211)
        plt.hist(slope_trace, histtype='stepfilled', bins=30, alpha=.7,
                 color="#A60628", normed=True)
        plt.legend(loc="upper left")
        plt.title("Posterior distributions")
        plt.ylim([0, 50])
        plt.ylabel("Slope")

        ax = plt.subplot(212)
        plt.hist(intercept_trace, histtype='stepfilled', bins=30,
                 alpha=.7, color="#7A68A6", normed=True)
        plt.legend(loc="upper left")
        plt.ylim([0, 2])
        plt.ylabel("Intercept")

        plt.show()

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(
            x,
            y_with_outlier,
            c='#7A68A6',
            lw=1,
            label='Original Signal')
        ax.scatter(
            outlier_x,
            outlier_y,
            facecolors='#A60628',
            edgecolors='#7A68A6',
            label='Outliers',
            lw=1,
            s=100,
            alpha=0.2)

        ax.set_xlim(np.min(x)-10, np.max(x) + 30)
        ax.set_ylim(np.min(y_with_outlier)-10, np.max(y_with_outlier) + 20)
        ax.set_title('Original Signal and Outliers')
        plt.legend()
        plt.show()

    def _getSampledData(self):
        if 0 == self._n:
            sampledSegs = self._data.segments
        else:
            sampledSegs = np.random.choice(self._data.segments, self._n)

        print "all sample: {}".format(len(sampledSegs))

        x0 = np.array(map(lambda seg: seg.gc * 100, sampledSegs))
        y0 = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                          np.log(seg.normal_reads_num + 1), sampledSegs))
        l = sorted(zip(y0, x0), reverse=True)
        y0, x0 = [list(t) for t in zip(*l)]

        return np.array(x0), np.array(y0)
