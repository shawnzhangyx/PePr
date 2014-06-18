import sys
import misc
import logging
import logConfig
import numpy


root_logger =  logging.getLogger("")
info = root_logger.info
debug = root_logger.debug
warning = root_logger.warning
shift_logger = logging.getLogger("shiftSizeEst")

def shift_size_per_chrom(forward, reverse, file=-1):
    #estimate the shift size for each chromosome separately. 

    shift_max = 200
    
    max_overlapping =  len(list(set(forward) & set(reverse)))
    optimum_shift = 0
    
    # iterate from 0 offset to shift_max, to find which offset produce the maximal overlap.  
    for offset in xrange(1, shift_max):
        
        overlapping = set([x+offset for x in forward]) & set([y-offset for y in reverse])
        overlapping = len(list(overlapping))

        if overlapping > max_overlapping:
            max_overlapping = overlapping
            optimum_shift = offset
        if file !=-1:
            file.write( str(offset)+'\t'+str(overlapping) + '\n')
    return optimum_shift

def estimate_shift_size(readData, parameter):
    # If shift size provided by the user, then skip estimating it.
    if parameter.shift_size != "-1": 
        shift_list = parameter.shift_size.split(',')
        if len(shift_list)== 1:
            for name in readData.filename_list: 
                readData.shift_size[name] = int(shift_list[0])
                info("%-10s %s", name, shift_list[0])
        else: 
            if parameter.difftest is True: 
                for idx, name in enumerate(readData.chip_filename_list):
                    readData.shift_size[name] = int(shift_list[idx])
                    info("%-10s %s", name, shift_list[idx])
                chip1_shift_list = [readData.shift_size[chip]
                                    for chip in readData.chip1_filename_list]
                input1_shift_size = sum(chip1_shift_list)/len(chip1_shift_list)
                for name in readData.input1_filename_list:
                    readData.shift_size[name] = input1_shift_size
                chip2_shift_list = [readData.shift_size[chip]
                                    for chip in readData.chip2_filename_list]
                input2_shift_size = sum(chip2_shift_list)/len(chip2_shift_list)
                for name in readData.input2_filename_list:
                    readData.shift_size[name] = input2_shift_size
            else: 
                for idx, chip in enumerate(readData.chip1_filename_list):
                    readData.shift_size[chip] = int(shift_list[idx])
                    info("%-10s %s", chip, shift_list[idx])
                for input in readData.input1_filename_list:
                    readData.shift_size[input] = \
                            sum([int(x) for x in shift_list])/len(shift_list)
                    info("%-10s %s", input, 
                          sum([int(x) for x in shift_list])/len(shift_list))
        return 
    info("begin estimating the shift size...")
    for chip_filename in readData.chip1_filename_list:
        info("estimating for %-10s", chip_filename)
        shift_list = []
        for count, chr in enumerate(readData.chr_list):
            # estimate the shift size for five chromosomes and take
            # the median of these as the estimator. 
            if count ==5:   
                break
            forward = readData.data_dict_by_strands[chr][chip_filename]['f']
            reverse = readData.data_dict_by_strands[chr][chip_filename]['r']
            shift = shift_size_per_chrom(forward, reverse)
            info("%-10s %d", chr, shift)
            shift_list.append(shift)
        shift_median = misc.median(shift_list[:])
        info("%-10s %d", chip_filename, shift_median)
        readData.shift_size[chip_filename] = shift_median
    chip1_shift_list = [readData.shift_size[chip] 
                        for chip in readData.chip1_filename_list]
    input1_shift_size = sum(chip1_shift_list)/len(chip1_shift_list)
    for input_filename in readData.input1_filename_list: 
        readData.shift_size[input_filename] = input1_shift_size
            
    if parameter.difftest is True:  # If we're calling differential binding. 
        for filename in readData.chip2_filename_list: 
            info("estimating for %-10s", filename)
            shift_list = []
            for count, chr in enumerate(readData.chr_list):
                if count==5: 
                    break
                forward = readData.data_dict_by_strands[chr][filename]['f']
                reverse = readData.data_dict_by_strands[chr][filename]['r']
                shift = shift_size_per_chrom(forward, reverse)
                info("%-10s %d", chr, shift)
                shift_list.append(shift)
            shift_median = misc.median(shift_list[:])
            info("%-10s %d", filename, shift_median)
            readData.shift_size[filename] = shift_median
        chip2_shift_list = [readData.shift_size[chip]
                            for chip in readData.chip2_filename_list]
        input2_shift_size = sum(chip2_shift_list)/len(chip2_shift_list)
        for input_filename in readData.input2_filename_list: 
            readData.shift_size[input_filename] = input2_shift_size
    return 

    
def shift_reads(readData): #shift the reads according the specific shift size estimated or provided by user
    readData.data_dict = {}
    for chr in readData.chr_list:
        debug("shifting the reads for %s...", chr)
        readData.data_dict[chr] ={}
        for file in readData.filename_list:
            forward=readData.data_dict_by_strands[chr][file]['f']
            reverse=readData.data_dict_by_strands[chr][file]['r']
            readData.data_dict_by_strands[chr][file]['f'] = numpy.array(readData.data_dict_by_strands[chr][file]['f'])
            readData.data_dict_by_strands[chr][file]['r'] = numpy.array(readData.data_dict_by_strands[chr][file]['r'])
            readData.data_dict[chr][file]=[x+readData.shift_size[file] for x in forward]+[y-readData.shift_size[file] for y in reverse]

#    return data_dict            

