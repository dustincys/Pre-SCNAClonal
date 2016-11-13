#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: lmmodel.py
#          Desc: lm model for MCMC
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-28 15:29:18
#       History:
# =============================================================================
'''
import numpy as np
import scipy
import pymc
import random
from matplotlib import pyplot as plt

a = 2
x = np.arange(10,50,.5)
# Add some noise
x += np.random.normal(size=x.size)
y = a * x + 10

#y = a * np.sin(x) + 5
y_with_outlier = np.copy(y)

# Add some outliers
for ii in np.arange(len(x)/5, len(x), len(x)/5.):
    y_with_outlier[ii]= 60*(random.random()-.5) + y[ii]

y_with_outlier = np.asarray(y_with_outlier)
# assuming I do not know, just fill it with ones
spread = np.asarray([1 for _ in range(y_with_outlier.size)])

outlier_points = pymc.Uniform('outlier_points', 0.0, 1.0, value=0.1)
mean_outliers = pymc.Uniform('mean_outliers', -100.0, 100.0, value=0.0)
spread_outliers = pymc.Uniform('spread_outliers', -100.0, 100.0, value=0.0)

slope = pymc.Uniform('slope', -2.0, -1.0, value = -1.0)
intercept = pymc.Uniform('intercept', -0.2, 0.2, value = 0.0)
intercept2 = pymc.Uniform('intercept2', -5.0, 5.0, value = 0.0)

@pymc.deterministic
def model_(x=x, slope=slope, intercept=intercept):
    fit = slope * x + intercept
    return fit

@pymc.deterministic
def inlier(slope = slope, intercept2 = intercept2, y_with_outlier =
           y_with_outlier, x = x):

    print "length y = {}".format(len(y_with_outlier))
    #print "length x = {}".format(len(x))
    #print "shape y = {}".format(y_with_outlier.shape)

    inlier = np.empty(len(y_with_outlier), dtype = bool)
    for i in range(len(y_with_outlier)):
        if y_with_outlier[i] < (x[i] * slope + intercept2):
            inlier[i] = True
        else:
            inlier[i] = False
        #print "y_with_outlier[i] = {}".format(y_with_outlier[i])
        #print "x[i] = {}".format(x[i])

    #print "slope = {}".format(slope)
    #print "intercept2 = {}".format(intercept2)
    #print "inlier"
    #print inlier

    return inlier

def log_posterior_likelihood_of_outlier(
    y_with_outlier,
    mu,
    spread,
    inlier,
    mean_outliers,
     spread_outliers):

    inlier_posterior = np.sum(
        inlier *
        (np.log(2 * np.pi * spread ** 2) + (y_with_outlier - mu) ** 2 /
         (spread ** 2)))

    outlier_posterior = np.sum(
        (1 - inlier) *
        (np.log(2 * np.pi * (spread_outliers ** 2)) +
         (y_with_outlier - mean_outliers) ** 2 /
         (spread_outliers ** 2)
         )
    )

    return -0.5 * (inlier_posterior + outlier_posterior)

outlier_distribution = pymc.stochastic_from_dist(
    'outlier_distribution',
    logp=log_posterior_likelihood_of_outlier,
    dtype=np.float,
     mv=True)


outlier_dist = outlier_distribution('outlier_dist',
                                    mu=model_,
                                    spread=spread,
                                    mean_outliers=mean_outliers,
                                    spread_outliers=spread_outliers,
                                    inlier=inlier,
                                    observed=True,
                                    value=y_with_outlier)

model = dict(outlier_dist=outlier_dist,
             slope=slope,
             intercept=intercept,
             model_=model_,
             inlier=inlier,
             intercept2=intercept2,
             mean_outliers=mean_outliers,
             spread_outliers=spread_outliers
             )
