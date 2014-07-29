#!/opt/local/bin/python2.6

import numpy
import scipy
from matplotlib import *
import pylab
from scipy import stats

thr = 9.2103
x1 = numpy.linspace(0, 10, 100000)
y1 = stats.chi2.pdf(x1, 2)
y2 = stats.chi2.pdf(x1, 4)
y3 = stats.chi2.pdf(x1, 6)
line1 = pylab.plot(x1, y1, 'r', alpha = 0.5)
line2 = pylab.plot(x1, y2, 'b', alpha = 0.5)
line3 = pylab.plot(x1, y3, 'g', alpha = 0.5)
X1 = numpy.linspace(thr, 10, 100000)
Y1 = stats.chi2.pdf(X1, 2) 
Y2 = stats.chi2.pdf(X1, 4)
Y3 = stats.chi2.pdf(X1, 6)
pylab.fill_between(X1,0,Y1, facecolor='r', alpha=0.5)
pylab.fill_between(X1,0,Y2, facecolor='b', alpha = 0.5)
pylab.fill_between(X1,0,Y3, facecolor='g', alpha=0.5)
pylab.show()