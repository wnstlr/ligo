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

if __name__ == '__main__':
    l = int(raw_input("Enter the number of data: "))
    N = 2
    f_s = 4096. * 2
    # fap = float(raw_input("Enter the desired FAP: "))
    MSOLAR_SI = 1.98892e30 	# 1 solar mass in kg 
    MPC_IN_SI =	3.08568025e22	# 1 mega parsec in meters 
    G_NEWT = 6.67300 * 10 ** (-11)
    LIGHT_SPEED_SI = 2.998 * 10 ** 8 # Speed of light in SI

    MSOLAR_IN_SEC = G_NEWT*MSOLAR_SI/LIGHT_SPEED_SI**3 

    m1 = 10. 		# mass1 in solar masses 
    m2 = 10. 		# mass2 in solar masses 
    d_mpc = 1000	# distance in megaparsecs

    d = d_mpc * MPC_IN_SI / LIGHT_SPEED_SI # distance in seconds 
    m1 = m1 * MSOLAR_IN_SEC 	# mass1 in seconds 
    m2 = m2 * MSOLAR_IN_SEC 	# mass2 in seconds 
    m = m1 + m2 		# total mass
    eta = m1 * m2 / m ** 2 	# symmetric mass ratio 	
    mc = m * eta ** (3./5.)	# chirp mass 

    phi0 = 0 # phase
    t0 = 0 # initial time
    f0 = 20 # lower cut-off frequency
    f_isco = (1./6.) ** (3./2.) / (numpy.pi * m) 
    f = numpy.linspace(20, f_isco, l)
    epsil = [0, 0.01, 0.1, 0.5, 0.75]
    S_h = generatePSD(f)
    nk = generateNoise(S_h)
    h_GR = computeNewtonianChirpFD(d, mc, f, phi0, t0)
    SNR = numpy.sqrt(4 * numpy.sum((h_GR.real ** 2 + h_GR.imag ** 2) / S_h * \
                         ((f_isco - f0) / l)))
    sigmak = 2 * numpy.sqrt(S_h)
    for i in range(len(epsil)):
       	ck = (nk + epsil[i] * h_GR)
       	pk = (ck.real ** 2 + ck.imag ** 2) / (sigmak ** 2)
       	n, bins, patches = pylab.hist(pk, 500, normed=True, \
                                          histtype='step', \
                                          label="$\epsilon=$"+str(epsil[i]))
       	pylab.setp(patches, 'alpha', 0.5)
    y = stats.chi2.pdf(bins, N, scale=1./16.)
    pylab.loglog(bins, y, 'k--', linewidth=1.)
    pylab.legend(loc="lower left", prop={'size':10})
    pylab.xlabel("$h(f)^2$ [1/Hz]")
    pylab.ylabel("Count")
    pylab.title("Distribution of Noise with Chirp Inserted \n(Deviation from GR Template with Colored Noise)")
    pylab.figtext(0.6,0.8,'SNR = ' + str(SNR), fontsize=10)
    pylab.xlim(10 ** (-3), 5)
    pylab.ylim(10 ** (-2),10)
    pylab.show()
    #pylab.savefig("General_Noise_Deviation_loglog.pdf")
    #pylab.close

