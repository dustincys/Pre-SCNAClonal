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

        def getSlopeMC(y, x):
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
        m, c = getSlopeMC(y_with_outlier, x)

        slope = pymc.Uniform('slope', m-1.0, m+1.0)
        #intercept = pymc.Uniform('intercept', c, c + 1.0)

        @pymc.deterministic
        def inlier(
                slope=slope,
                #intercept=intercept,
                y_with_outlier=y_with_outlier,
                x=x):

            intercept = 1
            #inlier_left_margin = constants.INLIER_LEFT_MARGIN * 100
            inlier = np.empty(len(y_with_outlier), dtype=bool)
            for i in range(len(y_with_outlier)):
                # if y_with_outlier[i] < (x[i] * slope + intercept + slope *
                # inlier_left_margin ):
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

        def get_inlier_posterior(y, x, slop, n):
            # divide y into n groups
            if n < 1:
                print "n less than 1"
                sys.exit(1)

            if len(y) < n:
                return -float("Inf")

            prob = 0

# tao = 1/ rho**2
            for i in range(1, n + 1):
                y_up = np.percentile(y, int(i * 100.0 / n))
                y_down = np.percentile(y, int((i - 1) * 100.0 / n))
                y_temp = y[np.logical_and(y > y_down, y < y_up)]
                x_temp = x[np.logical_and(y > y_down, y < y_up)]
                slope_temp = getLinalgSlope(y_temp, x_temp)
                prob = prob + 0.5 / n * pymc.distributions.normal_like(
                                                                      slope_temp,
                                                                      slop, 1 /
                                                                      (8 ** 2))

            slope_all = getLinalgSlope(y, x)

            prob = prob + 0.5 * pymc.distributions.normal_like(slope_all, slop,
                                                               1/(8 ** 2))

            return prob

        def get_outlier_posterior(inlier):
            alpha = 60
            beta = 10

            prob = pymc.distributions.beta_like(
                sum(inlier) * 1.0 / len(inlier), alpha, beta)

            return prob

        def correctY(y, x, slope, intercept):
            K = np.percentile(y, 50)
            A = slope * x + intercept
            return y - A + K

        def log_posterior_likelihood_of_outlier(
                y_with_outlier,
                inlier,
                slope
                ):

            # Here requires multiple

            # inlier_posterior = get_inlier_posterior(y_with_outlier[inlier],
            #                                        x[inlier],
            #                                        slope,
            #                                        2)

            #slope_mu, C = getSlopeMC(y_with_outlier, x)
            # intercept_posterior = pymc.distributions.normal_like(intercept, C,
            #                                                     1.0/ (0.05
            #                                                           **2))
            # slope_posterior = pymc.distributions.normal_like(slope,
            #                                                 slope_mu,
            #                                                 1.0 / (0.05 **2))

            # outlier_posterior = pymc.distributions.normal_like(
            #    sum(inlier) * 1.0 / len(inlier),
            #    0.6,
            #    1.0 / (0.002**2)
            #)

            y_corrected = correctY(y_with_outlier, x, slope, 0)
            y_density = gaussian_kde(y_corrected)

            y_30 = np.percentile(y_with_outlier, 30)
            y_70 = np.percentile(y_with_outlier, 70)
            y_xs = np.linspace(y_30, y_70, 20)

            y_ys = y_density(y_xs)
            index = np.argmax(y_ys)

            prob =  y_ys[index]

            return prob

        outlier_distribution = pymc.stochastic_from_dist(
            'outlier_distribution',
            logp=log_posterior_likelihood_of_outlier,
            dtype=np.float,
            mv=True)

        outlier_dist = outlier_distribution('outlier_dist',
                                            inlier=inlier,
                                            slope=slope,
                                            observed=True,
                                            value=y_with_outlier)

        model = dict(outlier_dist=outlier_dist,
                     slope=slope,
                     inlier=inlier
                     )

        mcmc = pymc.MCMC(model)
        mcmc.sample(10000, 2000)

        slope_trace = mcmc.trace('slope')[:]

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

        #ax = plt.subplot(212)
        #plt.hist(intercept_trace, histtype='stepfilled', bins=30,
        #         alpha=.7, color="#7A68A6", normed=True)
        #plt.legend(loc="upper left")
        #plt.ylim([0, 2])
        #plt.ylabel("Intercept")

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
