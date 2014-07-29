#!/opt/local/bin/python2.6

'''This program creates streams of white Gaussian noise ni, computes FFT
nk, computes chi^2 variable xk in the frequency domain, plots
histogram of the chi^2 variable.'''

import numpy
import scipy
from scipy import stats
from matplotlib import *
import pylab
import sys
import math
import random

def createGaussian(l):
    '''This function creates stremas of white Gaussian noise
    of length l, given as an argument.'''
    return numpy.random.randn(l)

def computeFourier(ni):
    '''This function takes in the white Gaussian noise of length N, ni,
    and computes FFT of it.'''
    return scipy.fft(ni)

def computeChi(nk, N):
    '''This function computes chi^2 variable xk in the frequency domain with
    degrees of freedom N'''
    xk = []
    if N % 2 != 0:
        print('wrong degress of freedom: it must be even')    
    res = nk.real ** 2 + nk.imag ** 2
    for i in range(len(res)):
        s = 0
        for j in range(N / 2):
            s += random.choice(res)
        xk.append(s)
    return xk

def plotHistogram(xk, N):
    '''This function plots the historgram of the chi^2 variables.'''
    n, bins, patches = pylab.hist(xk, 100)
    pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
    y = stats.chi2.pdf(bins, N)
    line = pylab.plot(bins, y, 'k--', linewidth=1.5)
    pylab.xlabel('Frequency [Hz]',fontsize=10)
    pylab.ylabel('Count', fontsize=10)
    pylab.title("Frequency domain Chi^2 distribution",fontsize=10,weight="bold")
    pylab.show()

if __name__ == '__main__':
    l = int(sys.argv[1]) # Number of data
    N = int(sys.argv[2]) # Degrees of freedom
    ni = createGaussian(l)
    nk = computeFourier(ni)
    xk = computeChi(nk, N)
    plotHistogram(xk, N)
