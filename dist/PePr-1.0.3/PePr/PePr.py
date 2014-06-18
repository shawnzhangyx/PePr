#!/usr/bin/env python
''' The main .py file for the PePr pipeline '''

import re
import os
import sys
import logging

import logConfig
import optParser
import fileParser
import shiftSize
import windowSize
import sigTests
import misc
from classDef import Parameters

def main(argv):
    # initialize the logger 
    root_logger = logging.getLogger("")
    debug = root_logger.debug
    info = root_logger.info
    # performing the option parser
    opt = optParser.opt_parser(argv)
    parameter, readData = optParser.process_opt(opt)
    # 1. read and parse the data
    fileParser.parse(readData, parameter.file_format)
    # 2. remove the redundant reads
    if (parameter.remove_redundant):
        readData.remove_redundant_reads()    
    # 3. shiftSize estimation and shifting reads
    shiftSize.estimate_shift_size(readData, parameter)
    shiftSize.shift_reads(readData)
    # 4. calculating the normalization constant 
    windowSize.estimate_normalization_constant(readData, parameter)
    # 5. windowSize estimation and split reads into windows
    windowSize.estimate_window_size(readData, parameter)
    info (" The windowSize is %s", parameter.window_size)
    windowSize.separate_exact_by_window(readData, parameter) 
    # 6. calling peaks
    if parameter.difftest is False:
        swap = False
        peakfilename = parameter.name+"__PePr_peaks.bed"
        sigTests.negative_binomial(readData, peakfilename, swap, parameter)
    else: 
        up_peakfilename = parameter.name+"__PePr_up_peaks.bed"
        swap = False
        sigTests.negative_binomial(readData, up_peakfilename, swap, parameter)
        down_peakfilename = parameter.name+"__PePr_down_peaks.bed"
        swap = True
        sigTests.negative_binomial(readData, down_peakfilename,
                                   swap, parameter)
    # 7. Write to a file that record the command and parameters.     
    parameter.write_parameter_to_file() 
 
if __name__ == '__main__':
    try: main(sys.argv)
    except KeyboardInterrupt: 
        print "user interrupted me"
        
        

