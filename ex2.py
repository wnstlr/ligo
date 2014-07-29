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

def createGaussian(l):
    '''This function creates stremas of white Gaussian noise
    of length l, given as an argument.'''
    return numpy.random.randn(l)

def computeFourier(ni):
    '''This function takes in the white Gaussian noise of length N, ni,
    and computes FFT of it.'''
    return scipy.fft(ni)

def computeChi(nk):
    '''This function computes chi^2 variable xk in the frequency domain with
    degrees of freedom N'''
    xk = []
    for i in range(len(nk)):
        xk.append(nk[i].real ** 2 + nk[i].imag ** 2)
    return xk

def chirp(phi0,f0,mchirp,distance,t):
    '''amplitude function of Newtonian chirp in time domain domain
    Equation (1),(2),(3),(4) and (5) in the iMDC primer'''
    tau = 5.0/(256.0*((numpy.pi*f0)**(8.0/3.0))*(mchirp**(5.0/3.0)))
    phival = phi0 + 2.0*numpy.pi*(tau**(3.0/8.0))*(f0)*\
    (integrate.quad(lambda x: (tau-x)**(-3.0/8.0), 0.0, t)[0])
    amp = 4*(3.6e-22)*((distance/100.0)**(-1.0))*((mchirp/1.22)**(5.0/3.0))*\
    (((f0/10.0)*(1.0-t/tau)**(-3.0/8.0))**(2.0/3.0))
    hval = amp*numpy.cos(phival)
    return hval

def addWave(ns, phi0, t):
    '''This function adds a Newtonian chirp signal to the noise created from 
    function white_noise(). Set M = 8 * solar_mass, D = 50Mpc, F_0 = 40Hz'''
    lst = []
    h0 = []
    for i in range(len(t)):
        h0.append(chirp(0, 40, 8 * 1.98892 * 10 ** 30, 50 * 3.08568025 * 10 ** 22,\
            t[i]))
    for i in range(len(ns)):
        ns[i] += h0[i]

    lst.append(ns)
    lst.append(h0)
    return lst

def plotHistogram(xk, N):
    '''This function plots the historgram of the chi^2 variables.'''
    n, bins, patches = pylab.hist(xk, len(xk), histtype='stepfilled')
    pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
    y = stats.chi2.pdf(bins, N)
    line = pylab.plot(bins, y, 'k--', linewidth=1.5)
    pylab.show()

if __name__ == '__main__':
    l = int(sys.argv[1]) # Number of data
    N = int(sys.argv[2]) # Degrees of freedom
    ni = createGaussian(l)
    nk = computeFourier(ni)
    xk = computeChi(nk)
    plotHistogram(xk, N)
