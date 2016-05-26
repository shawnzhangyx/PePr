###############################################################
## script to calcualte the maximum number of duplicates to keep at single position
## usage: python cal_max_duplicates_to_keep.py genome_size library_size
## By Yanxiao Zhang <yanxiazh@umich.edu>
## Timestamp: 5/26/2016
#######################

import sys
from scipy.stats import binom_test

def binomial(n, p):
    #calculate the expected maximum number of replicated reads at a single position
    x = 1
    pvalue = 0
    while (binom_test(x,n,p) > 0.00001):
        x = x + 1
    if x >1:
        x= x - 1
    return x 


def cal_max_dup(argv):
    ''' The main function for calculating the max dup'''
    help = 'usage: python cal_max_duplicates_to_keep.py genome_size library_size'
    if len(argv) != 3:
        print help
        exit(1)
    genome_size = float(argv[1])
    read_total = float(argv[2])
    expected_dup = binomial(read_total, 1.0/genome_size)
    print "The expected maximum duplicate is {0}".format(expected_dup)
    return expected_dup

if __name__ == "__main__":
    cal_max_dup(sys.argv)    
