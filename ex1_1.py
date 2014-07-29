#!/opt/local/bin/python2.6

'''This program creates streams of white Gaussian noise ni, computes FFT
nk, computes chi^2 variable xk in the frequency domain, plots
histogram of the chi^2 variable.'''

import numpy
import scipy
from scipy import stats
from matplotlib import *
import pylab
import random

def createGaussian(l):
    '''This function creates streams of white Gaussian noise
    of length l, given as an argument.'''
    return numpy.random.randn(l)

def split(ni, l):
    '''Splits a list into segments with size l'''
    return [ni[i:i+l] for i in range(0, len(ni), l)]

def splitFourier(ni, n):
    '''This function takes in a data stream and splits into n segments
    and fourier transforms each segments and than averages them up on the
    time domain.'''
    table = []
    lst = []
    result = numpy.zeros(len(ni) / n)
    new_arr = ni[0:(len(ni) / n * n)]
    segments = split(new_arr, len(new_arr) / n)
    for i in range(len(segments)):
        table.append(scipy.fft(segments[i]).real ** 2 + \
            scipy.fft(segments[i]).imag ** 2)
    for i in range(len(segments[0])):
        for j in range(len(table)):
            lst.append(table[j][i])
        result[i] = numpy.mean(lst)
    return result          
        
def psd(nk, f_s):
    '''This function takes in the fourier transform and computes power
    spectral density normalizing the FFT and averaging over many noise
    realizations. The argument takes is already-split and averaged '''
    return 2 * nk / float(f_s * len(nk))

def plotHistogram(xk, N):
    '''This function plots the historgram of the chi^2 variables along with the
    theoretical curve line.'''
    n, bins, patches = pylab.hist(xk, 100, normed=1)
    pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    y = stats.chi2.pdf(bins, N)
    line = pylab.plot(bins, y, 'k', linewidth=5.0)
    pylab.xlabel('$h(f)^2$ [ $1/Hz$ ]',fontsize=11)
    pylab.ylabel('Count', fontsize=11)
    pylab.title("Time Domain $\chi^2$ Distribution with " + str(N) + \
        " Degrees of Freedom",fontsize=13,weight="bold")
    pylab.show()

if __name__ == '__main__':
    l = int(raw_input('Enter the Number of data: '))
    N = int(raw_input('Enter the desired degrees of freedom: '))
    while N % 2 != 0:
        print('Error: the degrees of freedom must be a multiple of 2')
        N = int(raw_input('Try again: '))
    f_s = float(raw_input('Enter the sampling frequency: '))
    n = int(raw_input('Enter the number of segments: '))

    # Create white noise
    ni = createGaussian(l)
    # Split the noise, compute FFT, and average over time domain
    nk = splitFourier(ni, n)
    # Compute Power Spectral Density
    xk = psd(nk, f_s)
    # Plot the distribution along with theoretical curveline
    plotHistogram(xk, N)