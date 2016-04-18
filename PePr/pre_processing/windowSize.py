#!/usr/bin/env python

import sys
import re
import math
import random
import numpy

from scipy.stats import rankdata
from logging import info 
from cal_normalization import parse_to_bin

def estimate_window_size(filename, parameter, iter=500):
    '''function to estimate the window size if not provided'''
    
    info ("begin estimating window size")
    bin_size = 20 
    array_window = parse_to_bin(filename, bin_size, parameter)  
    #print array_window
    #print numpy.sum(array_window)
    rowNum = len(array_window)
    peak_len_list = []
    for x in range(iter):
        peak = numpy.max(array_window)
        peak_idx = numpy.where(array_window == peak)[0][0]
        # set the peak to -1 so that it won't confuse the next iteration.
        array_window[peak_idx]= -1 
        if peak_idx == 0:  # Check if the window reaches the left boundary.
            left_boundary_not_reached = False
        else:
            left_boundary_not_reached = True    
        # Check if the window reaches of the right boundary.
        if peak_idx == rowNum-1: 
            right_boundary_not_reached = False
        else: 
            right_boundary_not_reached = True
        i_l = 1
        i_r = 1 
        # Be careful if anything equals to -1. 
        # Though it's unlikely that it will happen.
        while left_boundary_not_reached:

            while (array_window[peak_idx-i_l] >= 0.1*peak): 
                array_window[peak_idx-i_l] = -1
                if peak_idx-i_l==0:
                    left_boundary_not_reached = False
                    break
                i_l += 1
            if left_boundary_not_reached:
                if peak_idx-i_l==0:
                    array_window[peak_idx-i_l] = -1
                    left_boundary_not_reached = False
            
            if left_boundary_not_reached: 
            # Will continue counting if the next window(one window gap) 
            # has above the 10% mode reads. 
                if (array_window[peak_idx-i_l-1] >= 0.1*peak): 
                    array_window[peak_idx-i_l] = -1
                    i_l += 1
                else:
                    break

        while right_boundary_not_reached:
            while (array_window[peak_idx+i_r] >= 0.1*peak):
                array_window[peak_idx+i_r] = -1
                if peak_idx+i_r ==rowNum-1:
                    right_boundary_not_reached = False
                    break
                i_r += 1
            if right_boundary_not_reached:
                if peak_idx+i_r==rowNum-1:
                    array_window[peak_idx+i_r] = -1
                    right_boundary_not_reached = False

            if right_boundary_not_reached:
                if (array_window[peak_idx+i_r+1] >=0.1*peak): 
                    array_window[peak_idx+i_r] = -1
                    i_r += 1
                else:
                    break
        peak_len_list.append((1+i_l+i_r)*bin_size)
    # print peak_len_list
    window_size = int(numpy.median(peak_len_list))
    info ("Estimated window size is %d", window_size)
    if window_size < 100:
        window_size = 100
        info("Estimated window size too narrow, use 100 instead")
    elif window_size > 1000:
        window_size = 1000
        info("Estimated window size too large, use 1000 instead")
    parameter.window_size = window_size
    return 


