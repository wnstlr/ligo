#!/opt/local/bin/python2.6

'''This module creates a chirp signal in both time and frequency domains.'''

import numpy
import matplotlib
matplotlib.use("Agg")
import pylab
from scipy import integrate
import scipy

## the 4km long detector's senstivity we need this just to show the 
##detector's sensitivity in contrast to the signal spectrum
## the design sensitvity file which is call in the line below
## can be downloaded from: http://gw-indigo.org/mdc-2011/tools/4k_design.txt
design4k = numpy.loadtxt("./4k_design.txt")

### amplitude function of Newtonian chirp in frequency domain
### Equation (6) in the iMDC primer
def A(distance,mchirp,frequency):
    amplitude = ((4.92549095e-9)/(1.0292712503e8))*(numpy.pi**(-2.0/3.0))*\
    (numpy.sqrt(5.0/24.0))*(distance**(-1.0))*(mchirp**(5.0/6.0))*\
    ((frequency)**(-7./6.))
    return amplitude


### amplitude function of Newtonian chirp in time domain domain
### Equation (1),(2),(3),(4) and (5) in the iMDC primer
def h(phi0,f0,mchirp,distance,t):
    tau = 5.0/(256.0*((numpy.pi*f0)**(8.0/3.0))*(mchirp**(5.0/3.0)))
    phival = phi0 + 2.0*numpy.pi*(tau**(3.0/8.0))*(f0)*\
    (integrate.quad(lambda x: (tau-x)**(-3.0/8.0), 0.0, t)[0])
    amp = 4*(3.6e-22)*((distance/100.0)**(-1.0))*((mchirp/1.22)**(5.0/3.0))*\
    (((f0/10.0)*(1.0-t/tau)**(-3.0/8.0))**(2.0/3.0))
    hval = amp*numpy.cos(phival)
    return hval

####### frequency domain plot
### choose starting frequency as 10 Hz
### choose phi0 = 0
f0 = 10.
phi0 = 0.
f = numpy.linspace(5,1500,1000)
AmpF = numpy.empty(len(f))
# choose distance to be 5 Mpc = 5000 kpc
d = 5000.
# choose chirp mass to be 5 solar mass
m = 10.
for i in range(len(f)):
    AmpF[i] = A(d,m,f[i])

pylab.loglog(design4k[:,0],design4k[:,1],"r-",label="Sensitivity design $\sqrt{S(f)}$ of the detector") 
pylab.loglog(f,AmpF,"b.",label="A(f) of Newtonian chirp") 
pylab.ylim(1e-24,1e-20) 
pylab.xlim(1e1,5e3) 
pylab.xlabel("Frequency [Hz]",fontsize=13,weight="bold") 
pylab.ylabel("h(f) [ $1/\sqrt{Hz} $ ] ",fontsize=13,weight="bold") 
pylab.title("Frequency domain plot",fontsize=13,weight="bold")
pylab.figtext(0.2,0.3,"Chirp mass = " + str(m) + " $M_{\odot}$",fontsize=13,weight="bold")
pylab.figtext(0.2,0.2,"Distance = " + str(d) + " kpc",fontsize=13,weight="bold")
pylab.xticks(size = 15,weight="bold")
pylab.yticks(size = 15,weight="bold")
pylab.legend(bbox_to_anchor=(0.95, 1), borderaxespad=0.)
pylab.savefig("amp_spectra.pdf") 
pylab.close()

####### time domain plot
t = numpy.linspace(-100000,0,10000)
AmpT = numpy.empty(len(t))
for i in range(len(t)):
    AmpT[i] = h(phi0,f0,m,d,t[i])

### in geometrized unit  1 solar mass = 4.92e-6 seconds

t = t*4.92549095e-6

pylab.plot(t,AmpT,"g.") 
pylab.xlabel("Time before merger [s]",fontsize=13,weight="bold") 
pylab.ylabel("Strain amplitude",fontsize=13,weight="bold") 
pylab.ylim([-6e-24,6e-24])
pylab.figtext(0.5,0.8,"Chirp mass = " + str(m) + " $M_{\odot}$",fontsize=13,weight="bold")
pylab.figtext(0.5,0.7,"Distance = " + str(d) + " kpc",fontsize=13,weight="bold")
pylab.xticks(size = 15,weight="bold")
pylab.yticks(size = 15,weight="bold")
pylab.title("Time domain wave plot",fontsize=13,weight="bold")
pylab.savefig("hoft.pdf") 
pylab.close()