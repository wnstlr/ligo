#!/opt/local/bin/python2.6

import numpy
import scipy
from scipy import stats
from matplotlib import *
import pylab

def fap(P, n, N):
    '''Given a false alarm probability, this function determines the threshold
    using the degrees of freedom N.'''
    lst = []
    for thr in range(n):
        if 1 - stats.chi2.cdf(thr, N) >= P:
            p = thr
    for i in numpy.linspace(p, p+1, 10001):
        lst.append(abs(P - (1 - stats.chi2.cdf(i, N))))
    return lst.index(min(lst)) * 0.0001 + p

def plotChi2(thr, N):
    '''This function plots the chi^2 distribution with the shaded area 
    indicating the given false alarm probability, and computes the 
    threshold that we need to get the given FAP.'''
    x = numpy.linspace(0, 10, 100000)
    y = stats.chi2.pdf(x, N)
    line = pylab.plot(x, y, 'r')
    X = numpy.linspace(thr, 10, 100000)
    Y = stats.chi2.pdf(X, N) 
    pylab.fill_between(X,0,Y, facecolor='#C0C0CF')
    pylab.figtext(0.5,0.8,'The area of shaded region is ' + str(P), \
        fontsize=13)
    pylab.figtext(0.5,0.7,'The threshold is ' + str(ans), fontsize=13)
    pylab.title("The Threshold Indication with " + str(N) + \
        " degrees of freedom", fontsize=13, weight='bold')

if __name__ == '__main__':
    # Given false alarm probability
    P = float(raw_input('Enter the FAP: '))
    while (P > 1) or (P < 0):
        print('Error: The probability cannot exceed 1 or be negative.')
        P = float(raw_input('Enter again: '))
    # The length of the data
    n = int(raw_input('Enter the length of the data: '))
    while (n > 1000000):
        print('Error: Too many data input')
        n = int(raw_input('Enter again: '))
    # Degrees of freedom
    N = int(raw_input('Enter the desired degrees of freedom: '))
    while (N % 2 != 0):
        print("Error: degrees of freedom must be even. Enter again: ")
        N = int(raw_input())
    # Compute the threshold
    ans = 9.2103
    # Plot the graph
    plotChi2(ans, N)
    pylab.show()
