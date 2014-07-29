#!/opt/local/bin/python2.6 

'''Plots the Amplitude against frequency plot for Post-Newtonian Approach and 
Inspiral Merger Ringdown Template.'''

import numpy
import scipy
from matplotlib import *
import pylab

if __name__ == '__main__':
    imr_psi_matlab = numpy.loadtxt("./imr_psi_freq.txt")
    l = int(raw_input("Number of data: "))
    #m1 = float(raw_input("Mass1 in solar mass: "))
    #m2 = float(raw_input("Mass2 in solar mass: "))
    #d = float(raw_input("Distance in megaparsecs: "))
    
    ### Compute physical values of the binary system
    MSOLAR_SI = 1.98892e30      # 1 solar mass in kg
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
    psi7 = 8.696e5*eta-5.838e6*eta**2+1.089e7*eta**3
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
    lambda_g = 2.8e15 / LIGHT_SPEED_SI
    D = d * 5 * (1 + redshift) ** 2 / (1 + (2 + redshift) \
        * (1 + redshift + numpy.sqrt(1+redshift)))
    beta = numpy.pi * D / (lambda_g ** 2 * (1 + redshift))
    t0 = 0
    phi0 = 0
    phi_g = 0
    tau_g = 0
    f_isco = (1./6.) ** (3./2.) / (numpy.pi * m)
    print ">> f_isco=" + str(f_isco)
    f = numpy.linspace(10, 2048, l)
    f_matlab = numpy.linspace(0, 2048, len(imr_psi_matlab))
    c_constant = (m ** (5./6.) / (d * numpy.pi ** (2./3.))) \
        * numpy.sqrt(5. * eta / 24.)
    print ">> C=" + str(c_constant)

    #IMR_Psi = numpy.empty(len(f))
    #for i in range(len(f)):
    #    nu = (numpy.pi * m * f[i]) ** (1./3.)
    #    IMR_Psi[i] = 2*numpy.pi*f[i]*t0+phi0+3./(128.*eta*nu**5)*\
    #        (1+nu**2*psi2+nu**3*psi3+nu**4*psi4+nu**6*psi6+nu**7*psi7)

    nu = (numpy.pi * f * m) ** (1./3.)
    IMR_Psi = 2*numpy.pi*f*t0+phi0+3./(128.*eta*nu**5)*\
            (1+nu**2*psi2+nu**3*psi3+nu**4*psi4+nu**6*psi6+nu**7*psi7)

    Graviton_Psi = IMR_Psi - beta * 1./ f + phi_g + tau_g * f
    pylab.semilogx(f, IMR_Psi, 'b', label="Python Inspiral Merger Ringdown")
    pylab.semilogx(f, Graviton_Psi, 'r', label="Python Insipiral Merger Ringdown Graviton")
    #pylab.semilogx(f_matlab, imr_psi_matlab, 'g', alpha=0.5, label="Matlab IMR Graviton")
    pylab.xlabel("Frequency (Hz)")
    pylab.ylabel("Phase $\Psi(f)$")
    pylab.title("Insipiral Merger Ringdown Phase")
    pylab.legend(loc="upper right", prop={'size':10})
    pylab.grid(True)
    pylab.show()
    #pylab.savefig("IMR_psi_plot.pdf")
    #pylab.close()
