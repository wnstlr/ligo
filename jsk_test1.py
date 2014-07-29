#!/opt/local/bin/python2.6

'''This module generates PSD, certain colored noise that follows such PSD,
computes chi^2 distribution of the data stream minus the template for General
relativity. Then, using False Alarm Probability, it plots the reliability of
the test using chi^2 method along the deviation from the GR template.'''

import numpy
import scipy
from scipy import stats
from matplotlib import *
import pylab

def generatePSD(f):
    '''This function calculates the Power Spectral Density using the fitted
    curve equation.'''
    x = f / 245.4
    S_h = 10 ** (-48) * (0.0152 * x ** (-4) + 0.2935 * x ** (9./4.) + 2.7951\
                             * x ** (3./2.) - 6.5080 * x ** (3./4.) + 17.7622)
    return S_h

def generateNoise(S_h):
    '''This function generates noise in the frequency domain that follows the
    PSD provided in the argument.'''	
    nreal = numpy.sqrt(S_h) * numpy.random.randn(len(S_h)) / 2.
    nimag = numpy.sqrt(S_h) * numpy.random.randn(len(S_h)) / 2.
    return (nreal + complex(0, 1) * nimag)

def computeNewtonianChirpFD(distance,mchirp,frequency,phase,time):
    '''This function creates an h(f) of Newtonian chirp in 
    frequency domain.'''
    amplitude = ((numpy.pi**(-2.0/3.0))*(numpy.sqrt(5.0/24.0))*(distance**(-1.0))\
                     *(mchirp**(5.0/6.0))*((frequency)**(-7./6.)))
    psi = 2 * numpy.pi * frequency * time - phase * numpy.ones(len(frequency))\
    - numpy.pi / 4. * numpy.ones(len(frequency)) + 3. / 128.\
    * (numpy.pi * mchirp * frequency) ** (-5./3.)
    return amplitude * numpy.exp(psi * complex(0,1))

def FAPthreshold(P, l, N):
    '''Given a false alarm probability, this function determines the threshold
    using the degrees of freedom N.'''
    lst = []
    for thr in range(l):
        if 1 - stats.chi2.cdf(thr, N) >= P:
            p = thr
    for i in numpy.linspace(p, p+1, 100001):
        lst.append(abs(P - (1 - stats.chi2.cdf(i, N))))
    return lst.index(min(lst)) * 0.00001 + p

if __name__ == '__main__':
    l = int(raw_input("Enter the number of data: "))
    N = 2
    fap = float(raw_input("Enter the desired FAP: "))
    MSOLAR_SI = 1.98892e30      # 1 solar mass in kg
    MPC_IN_SI = 3.08568025e22   # 1 mega parsec in meters
    G_NEWT = 6.67300 * 10 ** (-11) # gravitational constant in SI
    LIGHT_SPEED_SI = 2.998 * 10 ** 8 # Speed of light in SI
    MSOLAR_IN_SEC = G_NEWT * MSOLAR_SI / LIGHT_SPEED_SI ** 3
    m1 = 10. # mass1 in solar mass
    m2 = 10. # mass2 in solar mass
    d_mpc = 1000 # distance in megaparsecs
    d = d_mpc * MPC_IN_SI / LIGHT_SPEED_SI # distance in seconds
    
    m1 = m1 * MSOLAR_IN_SEC  # mass1 in seconds
    m2 = m2 * MSOLAR_IN_SEC  # mass2 in seconds
    m = m1 + m2  # total mass
    eta = m1 * m2 / m ** 2  # symmetric mass ratio   
    mc = m * eta ** (3./5.) # chirp mass
    
    phi0 = 0 # phase
    t0 = 0 # initial time
    f_isco = (1./6.) ** (3./2.) / (numpy.pi * m)
    f = numpy.linspace(20, f_isco, l)
    epsil = numpy.linspace(0, 0.5, 101)
    #epsil = [0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5]
    theta = []
    # thr = 9.2103 / 16. for fap = 0.01
    thr = FAPthreshold(fap, l, N) / 16.
    S_h = generatePSD(f)
    nk = generateNoise(S_h)
    h_GR = computeNewtonianChirpFD(d, mc, f, phi0, t0)
    sigmak = 2 * numpy.sqrt(S_h)
    for i in range(len(epsil)):
       	count = 0
        sum = 0
       	ck = (nk + epsil[i] * h_GR)
       	pk = (ck.real ** 2 + ck.imag ** 2) / (sigmak ** 2)
       	n, bins, patches = pylab.hist(pk, 500, normed=True)
       	pylab.close()
       	for j in range(len(bins)):
            if thr < bins[j]:
                count = j
                break
        for j in range(count,len(n)):
            sum = sum + (bins[j+1] - bins[j]) * n[j]
        theta.append(sum)
    pylab.plot(epsil, theta, 'go')	
    pylab.xlabel("Deviation from GR $\epsilon=h/h_{GR}-1$")
    pylab.ylabel("Efficiency of Detecting Deviations from GR")
    pylab.title("Efficiency of the Test based on Several Values of Deviations")
    pylab.show()
    #pylab.savefig("FAP" + str(fap) + "_" + str(d_mpc) + ".pdf")
    #pylab.close()
