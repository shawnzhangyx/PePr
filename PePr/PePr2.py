#!/usr/bin/env python
''' The main .py file for the PePr pipeline '''

# basic modules
import re
import os
import sys
import logging
import time

# local modules
import optParser
from pre_processing import preprocess
import prepareData
import sigTests
#import misc
from classDef import Parameters

debug = logging.debug
info = logging.info


def main(argv):

    opt = optParser.opt_parser(argv)
    parameter, readData = optParser.process_opt(opt)
    # read data and estimate the shiftsize, normalization constant
    preprocess(parameter)
    # read data again, begin to process. 
    
    prepareData.read_files_to_arrays(parameter)
    '''
    if parameter.difftest is False: 
        info ("peak-calling")
        swap = False
        peakfilename = parameter.name +"__PePr_peaks.bed"
        sigTests.negative_binomial(readData, peakfilename, swap, parameter)
    else: 
        up_peakfilename = parameter.name+"__PePr_chip1_peaks.bed"
        swap = False
        sigTests.negative_binomial(readData, up_peakfilename, swap, parameter)
        down_peakfilename = parameter.name+"__PePr_chip2_peaks.bed"
        swap = True
        sigTests.negative_binomial(readData, down_peakfilename,
                                   swap, parameter)
	parameter.write_parameter_to_file()
    '''
    info("PePr finished running, thanks for all the wait!")		
	
if __name__ == "__main__":
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        print ("user interrupt me")