### module for preprocessing the data to estimate the shift size and normalization constant. 
from collections import Counter
from logging import info


from fileParser import parse
from fileParser import parse_to_bin
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
    
    
    # if shift size is not available, estimate the shift size for chip1 and input1
    if len(parameter.shift_dict) is 0:
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
    
        if len(parameter.input1) > 0:
            for sample_name in parameter.input1:
                data_bin = parse_to_bin(parameter, sample_name)
                bin_dict[sample_name] = data_bin
                shift_dict[sample_name] = shift_chip1
                info ("shift size for %s is %d", sample_name,shift_chip1)

        # estimate shift size for chip2 
        if len(parameter.chip2) > 0:
            shift_chip2_list = []
            for sample_name in parameter.chip2:
                data = parse(parameter, sample_name)
                shift, data_bin = estimate_shiftsize(data)
                shift_dict[sample_name] = shift 
                shift_chip2_list.append(shift)
                del data
                info ("shift size for %s is %d", sample_name,shift)
                bin_dict[sample_name] = data_bin
            shift_chip2 = sum(shift_chip2_list)/len(shift_chip2_list)
        # process the data for input 2
        if len(parameter.input2) > 0:
            for sample_name in parameter.input2:
                data_bin = parse_to_bin(parameter, sample_name)
                bin_dict[sample_name] = data_bin
                shift_dict[sample_name] = shift_chip2
                info ("shift size for %s is %d", sample_name,shift_chip2)
        
        parameter.shift_dict = shift_dict
         
    elif len(parameter.normalization_dict) is 0:# if not need to estimate shift size, directly process them into bins. 
        for sample_name in parameter.chip1 + parameter.chip2 \
                        + parameter.input1 + parameter.input2:
            data_bin = parse_to_bin(parameter, sample_name)
            bin_dict[sample_name] = data_bin
            
    if len(parameter.normalization_dict) is 0:
        parameter.normalization_dict = cal_normalization.cal_normalization_constant(parameter, bin_dict)

    return 
    # estimate the normalization constants for chip2 
    
    