### module for preprocessing the data to estimate the shift size and normalization constant. 
from collections import Counter
import logging 

from fileParser import parse
from shiftSize import estimate_shiftsize

info = logging.info 

def preprocess(parameter):
    ''' read file and estimate shift sizes and normazliation parameters. 
    '''
    
    # estimate the normalization constants for chip1 
    bin_dict = Counter()
    for idx,entry in enumerate(parameter.chip1):
        data = parse(parameter, entry[0])
        shift, data_bin = estimate_shiftsize(data)
        parameter.chip1[idx] = [entry[0],shift]
        del data 
        info ("shift size for %s is %d", entry[0],shift)
    
    # estimate the normalization constants for chip2 
    
    