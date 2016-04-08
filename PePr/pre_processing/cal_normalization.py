from PePr import misc
from logging import info

def cal_normalization_constant(parameter, bin_dict):
    info ("calculating normalization constants")
    if parameter.normalization == "scale":
        return scale(parameter,bin_dict)
    elif parameter.normalization == "compound":
        return compound(parameter,bin_dict)

def scale(parameter, bin_dict):
    sizes = {}
    norm_constants = {}
    for key in bin_dict:
        
        sizes[key] = sum(bin_dict[key].values())
    
    size_median = misc.median(sizes.values())
    
    for key in sizes:
        norm_constants[key] = float(size_median)/sizes[key]
        #print norm_constants
        info ("scale parameter for %s is %.2f",key, norm_constants[key])
    return norm_constants
        
        
### will need to work on the compound normalization later. 
def compound(parameter, bin_dict):
    return
    
    
def chip_tmm(parameter, bin_dict):
	return 