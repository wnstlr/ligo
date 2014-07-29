#!/opt/local/bin/python2.6

import numpy
import scipy
from scipy import integrate
from matplotlib import *
import pylab
import sys

def white_noise(n):
    '''This function creates an array of random l numbers that reflect 
    characteristics of Gaussian white noise with zero mean and unit 
    variance. It also creates a noise vector with variance 1e-46.'''
    wn = numpy.random.randn(n)
    ns = wn * 10 ** (-23)
    return ns

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

def add_chirp(ns, phi0, t):
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

def correlation(x, h):
    '''This function calculates the correlation function between the data(x) 
    and a normalized template(h).'''
    func = numpy.fft(x) * numpy.conjugate(numpy.fft(h))
    return numpy.ifft(func)

if __name__ == '__main__':  
    n = int(sys.argv[1]) # number of data in noise spectrum
    f_s = float(sys.argv[2]) # 2048 in this case
    phi0 = 0
    ns = white_noise(n * f_s)
    t = numpy.linspace(0, n, f_s * n)
    res = add_chirp(ns, phi0, t)
    corr = correlation(res[0], res[1])
    pylab.plot(t, res[0], "g", t, res[1], "r")
    pylab.show()
    print('%f', corr)

