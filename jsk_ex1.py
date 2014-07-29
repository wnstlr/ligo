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
import random
import math

def createGaussian(l):
    '''This function creates stremas of white Gaussian noise
    of length l, given as an argument.'''
    return numpy.random.randn(l)

def computeFourier(ni):
    '''This function takes in the white Gaussian noise of length N, ni,
    and computes FFT of it.'''
    return scipy.fft(ni)

def normalize(nk):
    '''This function takes in the fourier transform and computes power
    spectral density normalizing the FFT and averaging over many noise
    realizations.'''
    return (2 * (nk.real**2 + nk.imag**2) / len(nk))

def computeChi(pd, N):
    '''This function computes chi^2 variable xk in the frequency domain with
    degrees of freedom N'''
    xk = []
    if N % 2 != 0:
        print('wrong degress of freedom: it must be even')    
    for i in range(len(pd)):
        s = 0
        for j in range(N / 2):
            s += random.choice(pd)
        xk.append(s)
    return xk

def plotHistogram(xk, N):
    '''This function plots the historgram of the chi^2 variables along with the
    theoretical curve line.'''
    n, bins, patches = pylab.hist(xk, 100, normed=True)
    pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    y = stats.chi2.pdf(bins, N)
    line = pylab.plot(bins, y, 'k', linewidth=1.5)
    pylab.xlabel('$h(f)^2$ [ $1/Hz$ ]',fontsize=11)
    pylab.ylabel('Count', fontsize=11)
    pylab.title("$\chi^2$ Distribution of Gaussian Noise with " + str(N) + \
        " Degrees of Freedom",fontsize=13,weight="bold")
    pylab.show()

if __name__ == '__main__':
    l = int(sys.argv[1]) # Number of data
    N = int(sys.argv[2]) # Degrees of freedom
    f_s = int(sys.argv[3]) # Sampling Frequency
    # Create white noise
    ni = createGaussian(l)
    # Compute Fourier transform of the white noise
    nk = computeFourier(ni)
    # Compute power spectral density
    nki = normalize(nk)
    # Compute chi^2 variables
    xk = computeChi(nki, N)
    # Plot the distribution along with theoretical curveline
    plotHistogram(xk, N)
