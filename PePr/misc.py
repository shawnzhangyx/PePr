#!/usr/bin/env python
from math import exp
import logging
from scipy.stats import binom_test

debug = logging.debug


def binomial(n, p):
    #calculate the expected maximum number of replicated reads at a single position

    x = 1
    pvalue = 0
    while (binom_test(x,n,p) > 0.00001):
        x = x + 1
    if x >1:
        x= x - 1
    return x 


def median(list):
    '''Will return the median of the list of numbers '''
    list.sort()
    if len(list) % 2 == 0:
        med = (list[int(len(list)/2)]+list[int(len(list)/2-1)])/2
    else:
        med = list[int((len(list)-1)/2)]
    return med

def erfcc(x):
    """Complementary error function."""
    z = abs(x)
    t = 1. / (1. + 0.5*z)
    r = t * exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
    	t*(.09678418+t*(-.18628806+t*(.27886807+
    	t*(-1.13520398+t*(1.48851587+t*(-.82215223+
    	t*.17087277)))))))))
    if (x >= 0.):
    	return r
    else:
    	return 2. - r

def ncdf(x):
    return 1. - 0.5*erfcc(x/(2**0.5))

