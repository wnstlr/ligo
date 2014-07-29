#!/opt/local/bin/python2.6

import numpy
import scipy
from matplotlib import *
import pylab
import math
import sys

def fourier(lst, delta):
    '''lst: input time-series data vector; delta: sampling rate. Plots 
    the magnitude of fft and hanning function.'''
    result = []
    result1 = []
    result2 = []
    for i in range(len(lst)):
        result.append(numpy.fft(lst[i]))
        result2.append(numpy.hanning(lst[i]))
    for i in range(len(result)):
        result1[i] = math.sqrt(result[i].real ** 2 + result[i].imag ** 2)

    t = numpy.linspace(0, len(lst), float(1 / delta) * len(lst))
    line1 = pylab.plot(t, result1, 'g')
    line2 = pylab.plot(t, result2, 'b')
    pylab.show()

def spectral_density(data):
    '''calculates spectral density by normalizing the FFT and averaging
    over many noise realizations.'''
    

if __name__ == '__main__':
    if sys.argv[1] == 'f':
        lst = 
        fourier(lst, sys.argv[2])
    elif sys.argv[2] == 'p':

    


