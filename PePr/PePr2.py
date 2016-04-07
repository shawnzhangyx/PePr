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
#from pre_processing import shiftSize
#from pre_processing import windowSize
#import sigTests
#import misc
from classDef import Parameters

debug = logging.debug
info = logging.info

def main(argv):

    opt = optParser.opt_parser(argv)
    parameter, readData = optParser.process_opt(opt)
    preprocess(parameter)
    # read data and estimate the shiftsize, normalization constant
    

if __name__ == "__main__":
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        print ("user interrupt me")