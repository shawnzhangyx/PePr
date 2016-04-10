### module for preprocessing the data to estimate the shift size and normalization constant. 
from collections import Counter
from logging import info

from fileParser import parse
from fileParser import reads_to_bin
from fileParser import get_chromosome_info
from shiftSize import estimate_shiftsize
import cal_normalization


def preprocess(parameter):
    ''' read file and estimate shift sizes and normazliation parameters. 
    '''
    # initialize the dictionaries
    bin_dict = {}
    shift_dict = {}
    normalization_dict = {}
    
	# use one chip sample estimate the chromosome size
    chip_filename = parameter.chip1[0]
    get_chromosome_info(parameter, chip_filename)
	
	
    # estimate the shift size for chip1 and input1
    shift_chip1_list = []
    for sample_name in parameter.chip1:
        data = parse(parameter, sample_name)
        shift, data_bin = estimate_shiftsize(data)
        shift_dict[sample_name] = shift
        shift_chip1_list.append(shift)
        del data 
        info ("shift size for %s is %d", sample_name,shift)
        bin_dict[sample_name] = data_bin
    shift_chip1 = sum(shift_chip1_list)/len(shift_chip1_list)
        
    for sample_name in parameter.input1:
        data = parse(parameter, sample_name)
        data_bin = reads_to_bin(data)
        del data
        bin_dict[sample_name] = data_bin
        shift_dict[sample_name] = shift_chip1
        info ("shift size for %s is %d", sample_name,shift_chip1)
    
    normalization_dict = cal_normalization.cal_normalization_constant(parameter, bin_dict)
    parameter.shift_dict = shift_dict
    parameter.normalization_dict = normalization_dict
    del bin_dict  
    
    return 
    # estimate the normalization constants for chip2 
    
    