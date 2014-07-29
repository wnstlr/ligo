#!/opt/local/bin/python2.6 

'''Plots the histogram of deviation from GR according to massive-graviton 
theory for various values of mass of graviton.'''

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
    ### Get input number of data
    l = int(raw_input("Number of data: "))
    ### Degrees of freedom
    N = 2
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
    print ">> Distance to the source in seconds = " + str(d)
    
    ### Compute phenomological phase parameters
    psi2 = 3715./756-920.9*eta+6742*eta**2-1.34e4*eta**3
    psi3 = -16*numpy.pi+1.702e4*eta-1.214e5*eta**2+2.386e5*eta**3
    psi4 = 15293365./508032.-1.2544e5*eta+8.7354e5*eta**2-1.6936e6*eta**3
    psi5 = 0
    psi6 = -8.898e5*eta+5.981e6*eta**2-1.128e7*eta**3
    psi7 = 8.696e5*eta-5.838e6*eta**2+1.089e7*eta**3
    f1 = (1-4.455+3.521+0.6437*eta-0.05822*eta**2-7.092*eta**3)/(numpy.pi*m)
    f2 = ((1-0.63)/2.+0.1469*eta-0.0249*eta**2+2.325*eta**3)/(numpy.pi*m)
    f3 = (0.3236-0.1331*eta-0.2714*eta**2+4.922*eta**3)/(numpy.pi*m)
    sigma = ((1-0.63)/4.-0.4098*eta+1.829*eta**2-2.87*eta**3)/(numpy.pi*m)
    print ">> [psi2, psi3, psi4, psi5, psi6, psi7] = "
    print [psi2, psi3, psi4, psi5, psi6, psi7]
    print ">> [f1, f2, f3, sigma] = "
    print [f1, f2, f3, sigma]

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
    print ">> D = "
    print D
    f_cut = 2048.
    f_imr = numpy.linspace(10, f_cut, l) # frequency domain for IMR and gravtion
    c_constant = (m ** (5./6.) / (d * numpy.pi ** (2./3.))) \
        * numpy.sqrt(5. * eta / 24.)
    print ">> alpha2, alpha3, epsilon1, epsilon2 =  "
    print [alpha2, alpha3, epsilon1, epsilon2]
    print ">> Frequency bins for IMR: "
    print [10, f_cut]
    ### Compute normalization constants
    vMerg = (numpy.pi * m * f1) ** (1./3.)
    vRing = (numpy.pi * m * f2) ** (1./3.)
    w_m =  1. + alpha2*vMerg**2+alpha3*vMerg**3 # normalization constant
    w_m = w_m / (1. + epsilon1 * vMerg + epsilon2 * vMerg ** 2)
    w_r = w_m*(numpy.pi*sigma/2.)*(f2/f1)**(-2./3.)*(1.+epsilon1*vRing+epsilon2\
                                                        *vRing**2)
    print ">> Normalization constant = "
    print w_m, w_r    
    
    ### Set up the procedure with necessary values
    grav_list = numpy.array([1,1.e14,2.8e15,1e16,1e30])
    grav_list_sec = grav_list / LIGHT_SPEED_SI 
    beta_list = numpy.pi * D / (grav_list_sec ** 2 * (1. + redshift))
    print ">> Beta List = "
    print beta_list

    ### Compute Psi of each cases
    nu = (numpy.pi * m * f_imr) ** (1./3.)
    IMR_Psi = 2*numpy.pi*f_imr*t0+phi0+3./(128.*eta*nu**5)*\
            (1+nu**2*psi2+nu**3*psi3+nu**4*psi4+nu**6*psi6+nu**7*psi7)
    print ">> IMR_Psi = "
    print IMR_Psi
    numpy.savetxt("IMR_Psi.txt", IMR_Psi)
    PN_Psi = 2 * numpy.pi * f_imr * t0 - phi0 * numpy.ones(len(f_imr))\
    - numpy.pi / 4. * numpy.ones(len(f_imr)) + 3. / 128.\
    * (numpy.pi * mc * f_imr) ** (-5./3.)
    print ">> PN_Psi = "
    print PN_Psi

    ### Compute PSD
    S_h = generatePSD(f_imr)
    ### Generate noise
    nk = generateNoise(S_h)
    ### Estimate std of noise
    sigmak = 2 * numpy.sqrt(S_h)

    amplitude = numpy.empty(len(f_imr))
    for i in range(len(f_imr)):
        amplitude[i] = computeIMRChirpAmplitudeFD(f_imr[i],f1,f2,f3,sigma,\
                                               c_constant,\
                                               alpha2,alpha3,epsilon1,\
                                               epsilon2,m,w_m,w_r)
    signal_GR = amplitude * numpy.exp(IMR_Psi * complex(0,1))
    print "Amplitude: "
    print amplitude
    numpy.savetxt("IMR_amp.txt", amplitude)
    print "GR signal"
    print signal_GR
    numpy.savetxt("signal_GR.txt", signal_GR)
    
    pylab.figure(0)
    ### When signal buried follows massive graviton theory
    for i in range(len(beta_list)):
        Graviton_Psi = IMR_Psi - beta_list[i] * 1./ f_imr + phi_g + tau_g * f_imr
        signal_mg = amplitude * numpy.exp(Graviton_Psi * complex(0,1))
        print ">> " + str(beta_list[i])
        print "Psi_Eff: "
        print Graviton_Psi
        numpy.savetxt("Graviton_Psi_" + str(grav_list[i]) + ".txt", Graviton_Psi)
        print "Graviton template: "
        print signal_mg
        numpy.savetxt("signal_MG_"+ str(grav_list[i]) + ".txt", signal_mg)
        ck = nk + signal_mg - signal_GR
        pk = (ck.real ** 2 + ck.imag ** 2) / (sigmak ** 2)
        print pk
        numpy.savetxt("pk_" + str(grav_list[i]) + ".txt", pk)
        if (grav_list[i] == 2.8e15):
            label_text = str(grav_list[i]) + "m (Current lower bound)" 
        else:
            label_text = str(grav_list[i]) + "m"
        n, bins, patches = pylab.hist(pk, 500, normed = True, alpha=0.5, \
            histtype='step', label=label_text)
    y = stats.chi2.pdf(bins, N, scale=1./16.)
    pylab.plot(bins, y, "k--", linewidth=2.0)
    pylab.xlabel("$|n_k+h_{mg}-h_{GR}|^2$")
    pylab.ylabel("Count")
    pylab.title("Massive Graviton Template's Deviation \n from $\chi^2$ distribution")
    pylab.legend(loc="upper right", title="Compton Wavelength",prop={"size": 10})

    ### When signal buried follows massive graviton theory 
    ## in loglog scale
    pylab.figure(1)
    for i in range(len(beta_list)):
        Graviton_Psi = IMR_Psi - beta_list[i] * 1./ f_imr + phi_g + tau_g * f_imr
        signal_mg = amplitude * numpy.exp(Graviton_Psi * complex(0,1))
        print ">> " + str(beta_list[i])
        print "Psi_Eff: "
        print Graviton_Psi
        print "Graviton template: "
        print signal_mg
        ck = nk + signal_mg - signal_GR
        pk = (ck.real ** 2 + ck.imag ** 2) / (sigmak ** 2)
        print pk
        if (grav_list[i] == 2.8e15):
            label_text = str(grav_list[i]) + "m (Current lower bound)" 
        else:
            label_text = str(grav_list[i]) + "m"
        n, bins, patches = pylab.hist(pk, 500, normed = True, alpha=0.5, \
            histtype='step', label=label_text)
    y = stats.chi2.pdf(bins, N, scale=1./16.)
    pylab.loglog(bins, y, "k--", linewidth=2.0)
    pylab.plot(bins, y, "k--", linewidth=2.0)
    pylab.xlabel("$|n_k+h_{mg}-h_{GR}|^2$")
    pylab.ylabel("Count")
    pylab.title("Massive Graviton Theory's Deviation \n from $\chi^2$ distribution")
    pylab.legend(loc="lower left", title="Compton Wavelength", prop={"size": 10})
    pylab.loglog(bins, y, "k--", linewidth=2.0)
    '''
    ### Plots h_mg - h_gr
    pylab.figure(2)
    pylab.subplot(1,2,1)
    Graviton_Psi = IMR_Psi - beta_list[2] / f_imr + phi_g + tau_g * f_imr
    signal_mg = amplitude * numpy.exp(Graviton_Psi * complex(0,1))
    residual = signal_mg - signal_GR 
    var1 = numpy.sqrt(residual.real ** 2 + residual.imag ** 2)
    var2 = numpy.sqrt(signal_GR.real ** 2 + signal_GR.imag ** 2)
    var3 = numpy.sqrt(signal_mg.real ** 2 + signal_mg.imag ** 2)
    pylab.loglog(f_imr, var1, 'b', label="Residual")
    pylab.loglog(f_imr, var2, 'r', alpha=0.5, label="GR signal (IMR)")
    pylab.loglog(f_imr, var3, 'g', alpha=0.5, label="Graviton Template")
    pylab.legend(loc="upper right", prop={"size":7})
    pylab.title("Residual between \nGraviton Template and GR Signal")
    pylab.xlabel("Frequency (Hz)")
    pylab.ylabel("Amplitude")
    pylab.subplot(1,2,2)
    pylab.semilogx(f_imr, abs(Graviton_Psi - IMR_Psi), 'b', label="Residual")
    pylab.semilogx(f_imr, IMR_Psi, 'r', label="IMR Phase")
    pylab.semilogx(f_imr, Graviton_Psi, 'g', label="Graviton Phase")
    pylab.title("Residual of Phase Between IMR and Graviton")
    pylab.xlabel("Frequency (Hz)")
    pylab.ylabel("Phase ($\Psi$)")
    pylab.legend(loc="upper right", prop={"size":7})

    ### Plots h_mg - h_gr
    pylab.figure(3)
    Graviton_Psi = IMR_Psi - beta_list[1] / f_imr + phi_g + tau_g * f_imr
    pylab.semilogx(f_imr, IMR_Psi, 'r', label="IMR Phase")
    pylab.semilogx(f_imr, Graviton_Psi, 'g', label="Graviton Phase")
    pylab.title("Phase of IMR and Graviton for $lambda=$" + str(grav_list[1]))
    pylab.xlabel("Frequency (Hz)")
    pylab.ylabel("Phase ($\Psi$)")
    pylab.legend(loc="upper right", prop={"size":7})

    pylab.figure(4)
    Graviton_Psi = IMR_Psi - beta_list[3] / f_imr + phi_g + tau_g * f_imr
    pylab.semilogx(f_imr, IMR_Psi, 'r', label="IMR Phase")
    pylab.semilogx(f_imr, Graviton_Psi, 'g', label="Graviton Phase")
    pylab.title("Phase of IMR and Graviton for $\lambda=$" + str(grav_list[3]))
    pylab.xlabel("Frequency (Hz)")
    pylab.ylabel("Phase ($\Psi$)")
    pylab.legend(loc="upper right", prop={"size":7})

    pylab.figure(5)
    Graviton_Psi = IMR_Psi - beta_list[4] / f_imr + phi_g + tau_g * f_imr
    pylab.semilogx(f_imr, IMR_Psi, 'r', label="IMR Phase")
    pylab.semilogx(f_imr, Graviton_Psi, 'g', label="Graviton Phase")
    pylab.title("Phase of IMR and Graviton for $\lambda=$" + str(grav_list[4]))
    pylab.xlabel("Frequency (Hz)")
    pylab.ylabel("Phase ($\Psi$)")
    pylab.legend(loc="upper right", prop={"size":7})

    pylab.figure(6)
    Graviton_Psi = IMR_Psi - beta_list[0] / f_imr + phi_g + tau_g * f_imr
    pylab.semilogx(f_imr, IMR_Psi, 'r', label="IMR Phase")
    pylab.semilogx(f_imr, Graviton_Psi, 'g', label="Graviton Phase")
    pylab.title("Phase of IMR and Graviton for $\lambda = $" + str(grav_list[0]))
    pylab.xlabel("Frequency (Hz)")
    pylab.ylabel("Phase ($\Psi$)")
    pylab.legend(loc="upper right", prop={"size":7})
'''
    pylab.show()