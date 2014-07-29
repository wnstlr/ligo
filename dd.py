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
    '''This function creates streams of white Gaussian noise
    of length l, given as an argument.'''
    return numpy.random.randn(l)

def computeFourier(ni):
    '''This function takes in the white Gaussian noise of length N, ni,
    and computes FFT of it.'''
    return scipy.fft(ni)

def psd(nk, f_s):
    '''This function takes in the fourier transform and computes power
    spectral density normalizing the FFT and averaging over many noise
    realizations.'''
    return (2 * (nk.real**2 + nk.imag**2) / (f_s * len(ni)))

def plotHistogram(xk, N):
    '''This function plots the historgram of the chi^2 variables along with the
    theoretical curve line.'''
    n, bins, patches = pylab.hist(xk, 100, normed=1)
    pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    y = stats.chi2.pdf(bins, N)
    line = pylab.plot(bins, y, 'k', linewidth=1.5)
    pylab.xlabel('$h(f)^2$ [ $1/Hz$ ]',fontsize=11)
    pylab.ylabel('Count', fontsize=11)
    pylab.title("Time Domain $\chi^2$ Distribution with " + str(N) + \
        " Degrees of Freedom",fontsize=13,weight="bold")
    pylab.show()

if __name__ == '__main__':
    l = int(sys.argv[1]) # Number of data
    N = int(sys.argv[2]) # Degrees of freedom
    f_s = float(sys.argv[3]) # Sampling frequency
    # Create white noise
    ni = createGaussian(l)
    # Compute Fourier transform of the white noise
    nk = computeFourier(ni)
    # Compute Power Spectral Density
    xk = psd(nk, f_s)
    # Plot the distribution along with theoretical curveline
    plotHistogram(xk, N)