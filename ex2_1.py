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
    res = nk ** 2
    for i in range(len(nk)):
        s = 0
        for j in range(N / 2):
            s += random.choice(res)
        xk.append(s)
    return xk

def normalize(nk):
    '''This function takes in the fourier transform and computes power
    spectral density normalizing the FFT and averaging over many noise
    realizations.'''
    return (2 * (nk.real**2 + nk.imag**2) / len(nk))

def chirpTime(phi0,f0,mchirp,distance,t):
    '''This function creates an amplitude function of Newtonian chirp in
    time domain.'''
    tau = 5.0/(256.0*((numpy.pi*f0)**(8.0/3.0))*(mchirp**(5.0/3.0)))
    phival = phi0 + 2.0*numpy.pi*(tau**(3.0/8.0))*(f0)*\
    (integrate.quad(lambda x: (tau-x)**(-3.0/8.0), 0.0, t)[0])
    amp = 4*(3.6e-22)*((distance/100.0)**(-1.0))*((mchirp/1.22)**(5.0/3.0))*\
    (((f0/10.0)*(1.0-t/tau)**(-3.0/8.0))**(2.0/3.0))
    hval = amp*numpy.cos(phival)
    return hval

def addWave(ni, phi0, f0, M, D, t):
    '''This function adds a Newtonian chirp signal to the noise created from 
    function createGaussian(). Takes in various physical properties of the 
    star as parameters. '''
    AmpT = numpy.empty(len(t))
    for i in range(len(t)):
        AmpT[i] = chirpTime(phi0,f0,M,D,t[i])
    return (ni + AmpT)

def plotHistogram(xk, N):
    '''This function plots the historgram of the chi^2 variables.'''
    n, bins, patches = pylab.hist(xk, 100, normed=1)
    pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    y = stats.chi2.pdf(bins, N)
    line = pylab.plot(bins, y, linewidth=1.5)
    pylab.xlabel('$h(f)^2$ [ $1/Hz$ ]',fontsize=11)
    pylab.ylabel('Count', fontsize=11)
    pylab.title("$\chi^2$ Distribution of White Noise with Chirp Signal Inserted"\
        ,fontsize=13,weight="bold")
    pylab.show()

if __name__ == '__main__':
    # Number of data
    l = int(raw_input('Enter the desired number of data: '))
    # Degrees of freedom
    N = int(raw_input('Enter the desired degrees of freedom: ')) 
    # Sampling frequency
    f_s = float(raw_input('Enter the desired sampling fequency: '))
    # Given in the exercise problem, we give the following values for 
    # the parameters of the chirp function:
    # M = 8 * solar_mass, D = 50Mpc, F0 = 40Hz
    f0 = 40
    phi0 = 0
    t = numpy.linspace(-100000,0,l)
    t = t * 4.92549095e-6
    M = 8 * 1.98892 * 10 ** 30
    D = 50 * 3.08568025 * 10 ** 22
    white = createGaussian(l)
    nk = addWave(white, phi0, f0, M, D, t)
    nki = scipy.fft(nk)
    nkii = normalize(nki)
    xk = computeChi(nkii, N)
    plotHistogram(xk, N)
