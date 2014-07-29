#!/opt/local/bin/python2.6 

'''This module compares the wave form computed by MATLAB code MkIMRPhenomBWave.m
with that created by IMR_Psi_plot and IMR_PN_plot.'''

import numpy
import scipy
from matplotlib import *
import pylab

def computeIMRChirpAmplitudeFD(f,f1,f2,f3,sigma,c_constant,alpha2,alpha3,\
                          epsilon1,epsilon2,m,w_m,w_r):
    '''This function creates improved wave signal under GR using IMRPhenomB wave
    form. Takes in distance, frequency, time, phase, total mass, symmetric mass 
    ratio, a constant, alpha2, alpha3, epsilon1, epsilon2, total mass, and 
    normalization constant as arguments.'''
    nu = (numpy.pi * m * f) ** (1./3.)
    f_prime = f / f1
    lorentzian = 1. / (2 * numpy.pi) * sigma / ((f - f2) ** 2 + sigma ** 2 / 4.)
    amplitude = c_constant*f1**(-7./6.)
    if (f < f1):
        amplitude = amplitude*f_prime**(-7./6.)*(1+alpha2*nu**2+alpha3*nu**3)
    elif ((f1 <= f) & (f < f2)):
        amplitude = amplitude*w_m*f_prime**(-2./3.)*(1+epsilon1*nu+epsilon2*nu**2)
    elif ((f2 <= f) & (f < f3)):
        amplitude = amplitude*w_r*lorentzian
    else:
        amplitude = 0
    return amplitude

if __name__ == '__main__':
	### Importing MATLAB results
	imr_amp = numpy.loadtxt("./imr_amp_freq.txt")
	l = int(raw_input("Number of data: "))    
    ### Compute physical values of the binary system
    MSOLAR_SI = 1.98892e30
    MPC_IN_SI = 3.08568025e22   # 1 mega parsec in meters
    G_NEWT = 6.67300 * 10 ** (-11) # gravitational constant in SI
    LIGHT_SPEED_SI = 2.998 * 10 ** 8 # Speed of light in SI
    MSOLAR_IN_SEC = G_NEWT * MSOLAR_SI / LIGHT_SPEED_SI ** 3
    m1 = 10. # mass1 in solar mass
    m2 = 10. # mass2 in solar mass
    print ">> mass1 in solarmass=" + str(m1)
    print ">> mass2 in solarmass=" + str(m2)
    d_mpc = 1000. # distance in megaparsecs
    print ">> Distance in megaparsecs=" + str(d_mpc)
    d = d_mpc * MPC_IN_SI / LIGHT_SPEED_SI # distance in seconds
    m1 = m1 * MSOLAR_IN_SEC  # mass1 in seconds
    m2 = m2 * MSOLAR_IN_SEC  # mass2 in seconds
    m = m1 + m2  # total mass
    eta = m1 * m2 / m ** 2  # symmetric mass ratio   
    mc = m * eta ** (3./5.) # chirp mass
    print ">> Total mass in sec=" + str(m)
    print ">> Chirp mass in sec=" + str(mc)
    print ">> Symmetric mass ratio:=" + str(eta)
    print ">> Distance in sec=" + str(d)
    
    ### Compute phenomological phase parameters
    psi2 = 3715./756-920.9*eta+6742*eta**2-1.34e4*eta**3
    psi3 = -16*numpy.pi+1.702e4*eta-1.214e5*eta**2+2.386e5*eta**3
    psi4 = 15293365./508032.-1.254e5*eta+8.735e5*eta**2-1.694e6*eta**3
    psi5 = 0
    psi6 = -8.898e5*eta+5.981e6*eta**2-1.128e7*eta**3
    psi7 = 8.696e5*eta-5.838*eta**2+1.089e7*eta**3
    print ">> [psi2, psi3, psi4, psi5, psi6, psi7]="
    print psi2, psi3, psi4, psi5, psi6, psi7

    f1 = (1-4.455+3.521+0.6437*eta-0.05822*eta**2-7.092*eta**3)/(numpy.pi*m)
    f2 = ((1-0.63)/2.+0.1469*eta-0.0249*eta**2+2.325*eta**3)/(numpy.pi*m)
    f3 = (0.3236-0.1331*eta-0.2714*eta**2+4.922*eta**3)/(numpy.pi*m)
    sigma = ((1-0.63)/4.-0.4098*eta+1.829*eta**2-2.87*eta**3)/(numpy.pi*m)
    print ">> [f1, f2, f3, sigma]="
    print f1, f2, f3, sigma
    
    alpha2 = -323. / 224. + 451. * eta / 168.
    alpha3 = 0
    epsilon1 = -1.8897
    epsilon2 = 1.6557
    redshift = 0.21
    #d_lum = d * 5 * (1 + redshift) ** 2 / (1 + (2 + redshift) \
        * (1 + redshift + numpy.sqrt(1+redshift)))
	f0 = 20
    f_isco = (1./6.) ** (3./2.) / (numpy.pi * m)
    print ">> f_isco=" + str(f_isco)
    tau0 = 5.*m*MSOLAR_IN_SEC/(256.*eta*(numpy.pi*m*MSOLAR_IN_SEC*f0)^(8/3)); 
	n = 2^numpy.ceil(numpy.log2(1.2*tau0*4096))
    f = numpy.linspace(20, 1000, l)
    f_matlab = numpy.linspace(0, 2048, n/2. + 1)
    c_constant = (m ** (5./6.) / (d * numpy.pi ** (2./3.))) \
        * numpy.sqrt(5. * eta / 24.)
    print ">> C=" + str(c_constant)
    
    ### Compute normalization constants                                           
    vMerg = (numpy.pi * m * f1) ** (1./3.)
    vRing = (numpy.pi * m * f2) ** (1./3.)
    w_m =  1. + alpha2 * vMerg ** 2 + alpha3 * vMerg ** 3           
    w_m = w_m / (1. + epsilon1 * vMerg + epsilon2 * vMerg ** 2)
    w_r = w_m*(numpy.pi*sigma/2.)*(f2/f1)**(-2./3.)*(1.+epsilon1*vRing+epsilon2\
                                                        *vRing**2)
    ### Compute the IMR Amplitude
    IMR_amp = numpy.empty(len(f))
    for i in range(len(f)):
        IMR_amp[i] = computeIMRChirpAmplitudeFD(f[i],f1,f2,f3,sigma,c_constant,\
                                                    alpha2,alpha3,epsilon1,\
                                                    epsilon2,m,w_m,w_r)

    pylab.loglog(f, IMR_amp, 'r', label=From Python)
    pylab.loglog(f_matlab, imr_amp, 'b', label=From MATLAB)
    pylab.xlabel("Frequency (Hz)")
    pylab.ylabel("Amplitude $|A(f)|$")
    pylab.title("Insipiral Merger Ringdown Amplitude")
    pylab.legend(loc="upper right", prop={'size':10})
    pylab.grid(True)
    pylab.show()
    #pylab.savefig("matlab_python_amp_compare.pdf")
    #pylab.close()



