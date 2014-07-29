#!/opt/local/bin/python2.6

'''This program creates streams of white Gaussian noise ni, adds chirp
   gravitational wave signal to the noise, and considers small deviation 
   from the General Relativity Template.'''

import numpy
import scipy
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

def computeNewtonianChirpFD(distance,mchirp,frequency,phase,time):
    '''This function creates an h(f) of Newtonian chirp in 
    frequency domain.'''
    amplitude = ((4.92549095e-9)/(1.0292712503e8)) * (numpy.pi**(-2.0/3.0))*\
    (numpy.sqrt(5.0/24.0))*(distance**(-1.0))*(mchirp**(5.0/6.0))*\
    ((frequency)**(-7./6.))
    psi = 2 * numpy.pi * frequency * time - phase * numpy.ones(len(frequency))\
    - numpy.pi / 4. * numpy.ones(len(frequency)) + 3. / 128.\
    * (numpy.pi * mchirp * frequency) ** (-5./3.)
    return amplitude * numpy.exp(psi * complex(0,1))

def plotHistogram(xk, N, epsilon):
    '''This function plots the historgram of the chi^2 variables.'''
    n, bins, patches = pylab.hist(xk, 100, normed=1,histtype='step',\
        label="$\epsilon$=" + str(epsilon))
    pylab.gca().set_yscale("log")

if __name__ == '__main__':
    # Number of data
    l = int(raw_input('Enter the desired number of data: '))
    # Degrees of freedom
    N = int(raw_input('Enter the desired degrees of freedom: ')) 
    # Given in the exercise problem, we give the following values for 
    # the parameters of the chirp function:
    # chirp mass mc = 10 solar mass, d = 50Mpc
    phi0 = 0 # phase
    t0 = 0 # initial time
    mc = 10. # in solar mass
    d = 50000. # in kpc
    G = 6.67300 * 10 ** (-11)# Gravitational constant in SI
    c = 2.998 * 10 ** 8 # Speed of light in SI
    f = numpy.linspace(5,1500,l) # frequency bin    
    epsil = [0, 0.01, 0.1, 0.3, 0.5, 0.75, 1, 5]
    white = createGaussian(l)
    sigma = numpy.std(white)
    wfk = scipy.fft(white) / numpy.sqrt(len(white))
    h_GR = computeNewtonianChirpFD(d, mc, f, phi0, t0)
    for i in range(len(epsil)):
        nk = (wfk + epsil[i] * h_GR)
        nki = (nk.real ** 2 + nk.imag ** 2) / sigma ** 2
        #xk = computeChi(nki, N)
        plotHistogram(numpy.log10(nki), N, epsil[i])
    pylab.xlabel('$\log{(h(f)^2) }$ [ $1/Hz$ ]',fontsize=11)
    pylab.ylabel('Count', fontsize=11)
    pylab.title("Distribution of White Noise with\n" + str(N) + \
        " Degrees of Freedom (Deviation from the GR)",\
fontsize=13,weight="bold")
    pylab.legend(loc="upper right", prop={'size':10})
    pylab.savefig('deviation.pdf')
    pylab.close()
    
