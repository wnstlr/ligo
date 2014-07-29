#!/opt/local/bin/python2.6 

'''Tests the given waveforms with Kolmogorov-Smirnov Test to quantify
the deviation from the known chi^2 distribution with 2 degrees of freedom.'''

import numpy
import scipy
from scipy import stats
from matplotlib import *
import pylab

def empricialDistFunc(data, threshold):
	'''Returns the empirical distributino function for n independent and
	identically distributed random variables using values from indicator
	function.'''
	arr = numpy.empty(len(data))
	for i in range(len(data)):
		if (data[i] <= threshold):
			arr[i] = 1
		else:
			arr[i] = 0
	return 1. / len(data) * numpy.sum(arr) 

def computeKSStatistic(empirical, cdffunc):
	'''This function takes in the empirical distribution function and a given 
	CDF fx to compute Kolmogorov-Smirnov Statistic.'''
	res = numpy.abs(empirical - cdffunc)
	return numpy.max(res)

def findKAlpha(alpha):
	'''Finds the value of K_alpha used for goodness-of-fit test 
	(Kolmogorov-Smirnov).'''
	if (alpha == 0.01):
		return 1.627624
	elif (alpha == 0.05):
		return 1.223848
	elif(alpha == 0.1):
		return 1.072749

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

def kstest(data, alpha, N, threshold):
	'''Performs KS test on the data input as an argument with alpha level of 
	conficence.'''
	#threshold = FAPthreshold(alpha, len(data), N) / 16.
	empirical = empiricalDistFunc(data, threshold)
	domain = numpy.linspace(0, threshold ,len(data))
	cdffunc = stats.chi2.cdf(domain,N,scale=1./16.)
	D_n = computeKSStatistic(empirical, cdffunc)
	K_alpha = findKAlpha(alpha)
	if (numpy.sqrt(len(data)) * D_n > K_alpha)
		print "The sample does not follow chi^2 distribution."
	else:
		print "The sample follows chi^2 distribution with the confidence \
		level of " + str(alpha)
