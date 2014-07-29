#!/opt/local/bin/python2.6

'''This program creates streams of white Gaussian noise ni, computes FFT
and averages up the segments to calculate the average norm of FFT,
plots the PSD for the noise.'''

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

if __name__ == '__main__':
    l = int(raw_input('Enter the Number of data in each segment: '))
    f_s = float(raw_input('Enter the sampling frequency: '))
    lst = [8, 16, 32, 64, 128]
    lines = []
    expl = ('8 segments', '16 segments', '32 segments', '64 segments',\
     '128 segments')
    t = numpy.linspace(0, 1, l)
    for i in range(len(lst)):
        ni = createGaussian(lst[i] * l)
        nk = splitFourier(ni, lst[i])
        xk = psd(nk, f_s)
        line = pylab.plot(t, xk, label=str(lst[i]) + " segments")
        lines.append(line)
    pylab.legend(loc="upper right", prop={'size':8})
    pylab.title('PSD of White Noise based on Number of Segments Averaged')
    pylab.xlabel('Frequency Bin')
    pylab.ylabel('Power')
    pylab.show()

