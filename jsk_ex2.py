#!/opt/local/bin/python2.6

'''This program creates streams of white Gaussian noise ni, adds chirp
   gravitational wave signal to the noise, computes FFT of ni, nk, 
   computes chi^2 variable xk in the frequency domain, plots
   histogram of the chi^2 variable.'''

import numpy
import scipy
from scipy import integrate
from scipy import stats
from matplotlib import *
import pylab
import random

def createGaussian(l):
    '''This function creates stremas of white Gaussian noise
    of length l, given as an argument.'''
    return numpy.random.randn(l) * 10 ** (-23)

def computeChi(nk, N):
    '''This function computes chi^2 variable xk in the frequency domain with
    degrees of freedom N'''
    xk = []
    if N % 2 != 0:
        print('Error: wrong degrees of freedom: it must be even')
        exit(1)    
    for i in range(len(nk)):
        s = 0
        for j in range(N / 2):
            s += random.choice(nk)
        xk.append(s)
    return xk

def normalize(nk):
    '''This function takes in the fourier transform and computes power
    spectral density normalizing the FFT.'''
    return (2 * (nk ** 2) / len(nk))

def amplitude(distance,mchirp,frequency):
    '''This function creates an amplitude function of Newtonian chirp in 
    frequency domain.'''
    amplitude = ((4.92549095e-9)/(1.0292712503e8))*(numpy.pi**(-2.0/3.0))*\
    (numpy.sqrt(5.0/24.0))*(distance**(-1.0))*(mchirp**(5.0/6.0))*\
    ((frequency)**(-7./6.))
    return amplitude

def addWave(ni, D, M, f):
    '''This function adds a Newtonian chirp signal to the noise created from 
    function createGaussian(). Takes in various physical properties of the 
    star as parameters. '''
    h0 = numpy.zeros(len(ni))
    for i in range(len(ni)):
        h0[i] = amplitude(D, M, f[i])
    return (ni + h0)

def plotHistogram(xk, N):
    '''This function plots the historgram of the chi^2 variables.'''
    n, bins, patches = pylab.hist(xk, 100, normed=1)
    pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    y = stats.chi2.pdf(bins, N)
    line = pylab.plot(bins, y, 'k', linewidth=1.5)
    pylab.xlabel('$h(f)^2$ [ $1/Hz$ ]',fontsize=11)
    pylab.ylabel('Count', fontsize=11)
    pylab.title("$\chi^2$ Distribution of White Noise with Chirp \nSignal \
Inserted with " + str(N) + " Degrees of Freedom",fontsize=13,weight="bold")
    pylab.show()

if __name__ == '__main__':
    # Number of data
    l = int(raw_input('Enter the desired number of data: '))
    # Degrees of freedom
    N = int(raw_input('Enter the desired degrees of freedom: ')) 
    # Sampling frequency
    # f_s = float(raw_input('Enter the desired sampling fequency: '))
    # Given in the exercise problem, we give the following values for 
    # the parameters of the chirp function:
    # chirp mass 8 * solarmass so M = 16 solar mass, D = 50Mpc, F0 = 40Hz
    f = numpy.linspace(5,1500,l)
    M = 16. # in solar mass
    D = 50000. # in kpc
    white = createGaussian(l)
    wn = numpy.sqrt(scipy.fft(white).real ** 2 + scipy.fft(white).imag ** 2)
    nk = addWave(wn, D, M, f)
    h0 = numpy.zeros(len(white))
    for i in range(len(nk)):
        h0[i] = amplitude(D, M, f[i])
    nk = (nk - h0) * 10 ** 23
    nki = normalize(nk)
    xk = computeChi(nki, N)
    plotHistogram(xk, N)
