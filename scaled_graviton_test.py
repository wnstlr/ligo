#!/opt/local/bin/python2.6 

'''Tests with IMR template and graviton template. 
Improved the GR signal with Insipiral Merger Ringdown (IMR) PhenomB.
Improved the non-GR signal with massive graviton theory.
@ Modification type: scaledGR, graviton
@ GR signal type: PN (PostNewtonian), IMR (Inspiral Merger Ringdown)'''

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

def computeNewtonianChirpAmplitudeFD(distance,mchirp,frequency):
    '''This function creates an h(f) of Newtonian chirp in 
    frequency domain.'''
    amplitude = ((numpy.pi**(-2.0/3.0))*(numpy.sqrt(5.0/24.0))*(distance**(-1.0))\
                     *(mchirp**(5.0/6.0))*((frequency)**(-7./6.)))
    return amplitude

def computeIMRChirpAmplitudeFD(f,f1,f2,f3,sigma,c_constant,alpha2,alpha3,\
                          epsilon1,epsilon2,m,w_m,w_r):
    '''This function creates improved wave signal under GR using IMRPhenomB wave
    form. Takes in distance, frequency, time, phase, total mass, symmetric mass 
    ratio, a constant, alpha2, alpha3, epsilon1, epsilon2, nu, and normalization
    constatns as arguments'''
    f_prime = f / f1
    nu = (numpy.pi * m * f) ** (1./3.)
    lorentzian = 1. / (2 * numpy.pi) * sigma / ((f - f2) ** 2 + sigma ** 2 / 4.)
    amplitude = c_constant*f1**(-7./6.)
    if (f < f1):
        amplitude = amplitude*f_prime**(-7./6.)*(1+alpha2*nu**2+alpha2*nu**3)
    elif((f1 <= f) & (f < f2)):
        amplitude = amplitude*w_m*f_prime**(-2./3.)*(1+epsilon1*nu+epsilon2*nu**2)
    elif((f2 <= f) & (f < f3)):
        amplitude = amplitude*w_r*lorentzian
    else:
        amplitude = 0
    return amplitude

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
    ### List of possible options to select from
    modtype_options = ['scaledGR', 'graviton']
    grsignaltype_options = ['PN', "IMR"]
    ### Get input number of data
    l = int(raw_input("Number of data: "))
    print("Select " + str(modtype_options))
    modtype = raw_input("Modification type: ")
    while (modtype not in modtype_options):
        modtype = raw_input("Try again: ")
    print("Select " + str(grsignaltype_options))
    grsignaltype = raw_input("GR signal type: ")                        
    while (grsignaltype not in grsignaltype_options):
        grsignaltype = raw_input("Try again: ")
    ### Degrees of freedom
    N = 2
    #m1 = float(raw_input("Mass1 in solar mass: "))
    #m2 = float(raw_input("Mass2 in solar mass: "))
    #d = float(raw_input("Distance in megaparsecs: "))
    #fap = float(raw_input("Enter the desired FAP: "))
    
    ### Compute physical values of the binary system
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
    print ">> mass1 in seconds = " + str(m1)
    print ">> mass2 in seconds = " + str(m2)
    print ">> Symmetric mass ratio = " + str(eta)
    print ">> Chirp mass in seconds = " + str(mc) 
    print ">> Distance to the source in Mpc = " + str(d_mpc)
    
    ### Compute phenomological phase parameters
    psi2 = 3715./756-920.9*eta+6742*eta**2-1.34e4*eta**3
    psi3 = -16*numpy.pi+1.702e4*eta-1.214e5*eta**2+2.386e5*eta**3
    psi4 = 15293365./508032.-1.254e5*eta+8.735e5*eta**2-1.694e6*eta**3
    psi5 = 0
    psi6 = -8.898e5*eta+5.981e6*eta**2-1.128e7*eta**3
    psi7 = 8.696e5*eta-5.838e6*eta**2+1.089e7*eta**3
    f1 = (1-4.455+3.521+0.6437*eta-0.05822*eta**2-7.092*eta**3)/(numpy.pi*m)
    f2 = ((1-0.63)/2.+0.1469*eta-0.0249*eta**2+2.325*eta**3)/(numpy.pi*m)
    f3 = (0.3236-0.1331*eta-0.2714*eta**2+4.922*eta**3)/(numpy.pi*m)
    sigma = ((1-0.63)/4.-0.4098*eta+1.829*eta**2-2.87*eta**3)/(numpy.pi*m)
    print ">> [psi2, psi3, psi4, psi5, psi6, psi7] = "
    print psi2, psi3, psi4, psi5, psi6, psi7
    print ">> [f1, f2, f3, sigma] = "
    print f1, f2, f3, sigma

    ### Compute relevant constants
    alpha2 = -323. / 224. + 451. * eta / 168.
    alpha3 = 0
    epsilon1 = -1.8897
    epsilon2 = 1.6557
    redshift = 0.21
    phi0 = 0 # phase offset
    t0 = 0 # arrival time
    tau_g =  0 # redefinition of arrival time
    phi_g =  0 # redefinition of phase offset
    D = d * (1 + (2 + redshift) * (1 + redshift + numpy.sqrt(1+redshift))) /\
    (5 * (1 + redshift) ** 2)
    print ">> Distance parameter = "
    print D
    f_isco = (1./6.) ** (3./2.) / (numpy.pi * m)
    f_cut = 2048.
    f_pn = numpy.linspace(20, f_isco, l) # frequency domain
    f_imr = numpy.linspace(10, f_cut, l) # frequency domain for IMR and gravtion
    c_constant = (m ** (5./6.) / (d * numpy.pi ** (2./3.))) \
        * numpy.sqrt(5. * eta / 24.)
    
    ### Compute normalization constants
    vMerg = (numpy.pi * m * f1) ** (1./3.)
    vRing = (numpy.pi * m * f2) ** (1./3.)
    w_m =  1. + alpha2*vMerg**2+alpha3*vMerg**3 # normalization constant
    w_m = w_m / (1. + epsilon1 * vMerg + epsilon2 * vMerg ** 2)
    w_r = w_m*(numpy.pi*sigma/2.)*(f2/f1)**(-2./3.)*(1.+epsilon1*vRing+epsilon2\
                                                        *vRing**2)    
    
    ### Set up the procedure with necessary values
    #epsil = numpy.linspace(0, 0.5, 101)
    epsil = numpy.array([0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5])
    grav_list = numpy.logspace(16,16.74036,num=50) 
    grav_list_sec = grav_list / LIGHT_SPEED_SI
    #grav_list = numpy.array([2.8e15, 3e15, 5e15, 7.5e15, 1e16, 5e16, 1e17])\
    # / LIGHT_SPEED_SI 
    beta_list = numpy.pi * D / (grav_list_sec ** 2 * (1 + redshift))
    theta = []
    thr = 9.2103 / 16. 
    # thr = FAPthreshold(fap, l, N) / 16.

    ### Compute Psi of each cases
    PN_Psi = 2 * numpy.pi * f_pn * t0 - phi0 * numpy.ones(len(f_pn))\
    - numpy.pi / 4. * numpy.ones(len(f_pn)) + 3. / 128.\
    * (numpy.pi * mc * f_pn) ** (-5./3.)
    nu = (numpy.pi * m * f_imr) ** (1./3.)
    IMR_Psi = 2*numpy.pi*f_imr*t0+phi0+3./(128.*eta*nu**5)*\
            (1+nu**2*psi2+nu**3*psi3+nu**4*psi4+nu**6*psi6+nu**7*psi7)

    ### Compute GR signal type
    if (grsignaltype == 'PN'):
        S_h = generatePSD(f_pn)
        nk = generateNoise(S_h)
        sigmak = 2 * numpy.sqrt(S_h)
        signal_GR = computeNewtonianChirpAmplitudeFD(d, mc, f_pn) * \
            numpy.exp(PN_Psi * complex(0,1))
    elif (grsignaltype == 'IMR'):
        S_h = generatePSD(f_imr)
        nk = generateNoise(S_h)
        sigmak = 2 * numpy.sqrt(S_h)
        signal_GR = numpy.empty(len(f_imr))
        for i in range(len(f_imr)):
            signal_GR[i] = computeIMRChirpAmplitudeFD(f_imr[i],f1,f2,f3,sigma,\
                                               c_constant,\
                                               alpha2,alpha3,epsilon1,\
                                               epsilon2,m,w_m,w_r)
        signal_GR = signal_GR * numpy.exp(IMR_Psi * complex(0,1))
    ### Compute Non-GR signal amplitude
    ## Graviton template
    signal_mg_amp = numpy.empty(len(f_imr))
    for j in range(len(f_imr)):
        signal_mg_amp[j] = computeIMRChirpAmplitudeFD(f_imr[j],f1,f2,f3,sigma,\
                                                  c_constant,alpha2,alpha3,\
                                                  epsilon1,epsilon2,m,w_m,w_r)
    ### When signal buried has the form of (1 + epsil) * signal_GR
    if (modtype == 'scaledGR'):
        print "Scales for Deviation from GR: "
        print epsil
        for i in range(len(epsil)):
            count = 0
            s = 0
            ck = (nk + epsil[i] * signal_GR)
            pk = (ck.real ** 2 + ck.imag ** 2) / (sigmak ** 2)
            n, bins, patches = pylab.hist(pk, 500, normed=True)
            pylab.close()
            for j in range(len(bins)):
                if thr < bins[j]:
                    count = j
                    break
            for j in range(count,len(n)):
                s = s + ((bins[j+1] - bins[j]) * n[j])
            theta.append(s)
        pylab.plot(epsil, theta, 'go')
        pylab.xlabel("Deviation from GR $\epsilon=h/h_{GR}-1$")
        pylab.ylabel("Efficiency of Detecting Deviations from GR")
        pylab.title("Efficiency of the Test based on Values of Deviations")
        pylab.show()
        #pylab.savefig("FAP" + str(fap) + "_" + str(d_mpc) + ".pdf")
        #pylab.close()

    ### When signal buried follows massive graviton theory
    elif (modtype == 'graviton'):
        grav_list2_1 = numpy.logspace(12, 15.85715, 4)
        smooth = numpy.loadtxt("smooth3.txt")
        grav_list2_2 = smooth[:,0]
        grav_list2_3 = numpy.logspace(17, 30, 10)
        grav_list2 = numpy.append(grav_list2_1, grav_list2_2)
        grav_list2 = numpy.append(grav_list2, grav_list2_3)
        grav_list_sec2 = grav_list2 / LIGHT_SPEED_SI
        beta_list2 = numpy.pi * D / (grav_list_sec2 ** 2 * (1 + redshift))
        x_current = (numpy.array([2.8e15, 2.8e15]))
        y_current = numpy.array([0,0.03])
        for i in range(len(beta_list2)):
            count = 0
            s = 0
            Graviton_Psi = IMR_Psi - beta_list2[i] * 1./ f_imr + phi_g + tau_g * f_imr
            signal_mg = signal_mg_amp * numpy.exp(Graviton_Psi * complex(0,1))
            ck = nk + signal_mg - signal_GR
            pk = (ck.real ** 2 + ck.imag ** 2) / (sigmak ** 2)
            n, bins, patches = pylab.hist(pk, 500, normed = True)
            pylab.close()
            for j in range(len(bins)):
                if thr < bins[j]:
                    count = j
                    break
            for j in range(count,len(n)):
                s = s + ((bins[j+1] - bins[j]) * n[j])
            theta.append(s)
        pylab.semilogx(grav_list2, theta, "g")
        #numpy.savetxt("graviton_prob_.txt", theta)
        pylab.semilogx(x_current, y_current, 'k--', label="Current compton length boundry")
        pylab.xlabel("Compton Wavelength of the Graviton in Meters")
        pylab.ylabel("Efficiency of Detecting Deviations from GR")
        pylab.title("Efficeiceny of Detecting Deviation from GR with \n\
        Massive Graviton Template")
        pylab.legend(loc="upper right", prop={"size":8})
        pylab.show()
        #pylab.savefig("graviton_test.pdf")
        #pylab.close()
