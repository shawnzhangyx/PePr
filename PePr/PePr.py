#!/usr/bin/env python
''' The main .py file for the PePr pipeline '''

# basic modules
import re
import os
import sys
import time
from logging import info

# local modules
import optParser
from pre_processing import initialize
import prepareData
import sigTests
from post_processing import post_process_PePr

def post_processing_module():
    '''setuptools entry_points'''
    post_process_PePr.main(sys.argv)

def pre_processing_module():
    opt = optParser.opt_parser(sys.argv)
    parameter = optParser.process_opt(opt)
    # read data and estimate the shiftsize, normalization constant
    initialize.preprocess(parameter)
    parameter.write_parameter_to_file()

def argless_main():
    '''setuptools entry_points take no arguments, wrap main instead'''
    main(sys.argv)

def main(argv):
    '''PePr main function. Including preprocing and peak calling. 
    No post-processing'''

    opt = optParser.opt_parser(argv)
    parameter = optParser.process_opt(opt)
    # read data and estimate the shiftsize, normalization constant
    initialize.preprocess(parameter)
    parameter.write_parameter_to_file()

    # read data again, begin to process. 
    
    prepareData.read_files_to_arrays(parameter)
    read_dict = prepareData.prepare_data(parameter)
    
    if parameter.difftest is False: 
        swap = False
        peakfilename = parameter.output_directory + parameter.name +"__PePr_peaks.bed"
        sigTests.negative_binomial(read_dict, peakfilename, swap, parameter)
    else: 
        
        up_peakfilename = parameter.output_directory + parameter.name+"__PePr_chip1_peaks.bed"
        swap = False
        sigTests.negative_binomial(read_dict, up_peakfilename, swap, parameter)
        down_peakfilename = parameter.output_directory + parameter.name+"__PePr_chip2_peaks.bed"
        swap = True
        sigTests.negative_binomial(read_dict, down_peakfilename,
                                   swap, parameter)

    info("PePr finished running, thanks for all the wait!")
	
if __name__ == "__main__":
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        print ("user interrupt me")
