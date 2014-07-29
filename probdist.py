#!/opt/local/bin/python2.6

'''This file plots random numbers  '''

import numpy
import scipy
from scipy import stats
from matplotlib import *
import pylab
import sys
import math

def probdist1(l):
    '''This function creates an array of random numbers from Gaussian 
    distribution, and then plots a histogram of those random numbers and 
    the exptected probability distribution computed solely from the mean
    and variance of the data points.'''
    lstrand = []
    for i in range(l):
        lstrand.append(numpy.random.randn());
    mu = numpy.mean(lstrand)
    sigma = math.sqrt(numpy.var(lstrand))
    n, bins, patches = pylab.hist(lstrand, 100, normed=1)
    pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    y = pylab.normpdf(bins, mu, sigma)
    line = pylab.plot(bins, y, 'k--', linewidth=1.5)
    pylab.show()

def prob_threshold(thr):
    '''This function, using the distribution, calculates the probability
    for one of the data to be greater than a threshold, which is given 
    as an argument.'''
    return 1 - stats.norm.cdf(thr)

def probdist2(l, df):
    '''This function plots the histrogram of data sets which are squares of 
    data from normal distributions. '''
    lstrand = []
    for i in range(l):
        lstrand.append(numpy.random.randn() ** 2)

    n, bins, patches = pylab.hist(lstrand, 100)
    pylab.setp(patches, 'facecolor', 'b', 'alpha', '0.75')
    y = stats.chi2.pdf(bins, df)
    line = pylab.plot(bins, y, 'k--', linewidth=1.5)
    pylab.show()

def prob_threshold2(thr, df):
    '''This function, using the chi^2 distribution, calculates the probability
    for one of the data to be greater than a threshold, which is given as
    an argument.'''
    return 1 - stats.chi2.cdf(thr, df)

if __name__ == '__main__':
    if len(sys.argv) == 0:
        print("p for 1st, t for threshold, x for 2nd option")
    if sys.argv[1] == 'p':
        probdist1(int(sys.argv[2]))
    elif sys.argv[1] == 't':
        print prob_threshold(float(sys.argv[2]))
    elif sys.argv[1] == 'x':
        probdist2(float(sys.argv[2]), float(sys.argv[3]))
    elif sys.argv[1] == 'r':
        print prob_threshold2(float(sys.argv[2]), float(sys.argv[3]))
