#!/opt/local/bin/python2.6

import matplotlib.numerix as nx
import pylab as p

def boltzman(x, xmid, tau):
    """
    evaluate the boltzman function with midpoint xmid and time constant tau
    over x
    """
    return 1. / (1. + nx.exp(-(x-xmid)/tau))

def fill_below_intersection(x, S, Z):
    """
    fill the region below the intersection of S and Z
    """
    #find the intersection point
    ind = nx.nonzero( nx.absolute(S-Z)==min(nx.absolute(S-Z)))[0]
    # compute a new curve which we will fill below
    Y = nx.zeros(S.shape, typecode=nx.Float)
    Y[:ind] = S[:ind]  # Y is S up to the intersection
    Y[ind:] = Z[ind:]  # and Z beyond it
    p.fill(x, Y, facecolor='blue', alpha=0.5)

x = nx.arange(-6, 6, .01)
S = boltzman(x, 0, 1)
Z = 1-boltzman(x, 0.5, 1)
p.plot(x, S, x, Z, color='red', lw=2)
fill_below_intersection(x, S, Z)
p.show()