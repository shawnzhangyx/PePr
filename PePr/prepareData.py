import numpy
from logging import info
import pysam
import multiprocessing
import itertools 

def read_bam(filename, data_dict, parameter):
    shift_size = parameter.shift_dict[filename]
    move_size = parameter.window_size/2
        
    infile = pysam.Samfile(filename, 'rb')
    for line in infile.fetch():
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            if line.is_reverse is False:
                pos = line.pos+ shift_size
            else:
                pos = line.pos - shift_size
            try: 
                data_dict[chr][int(pos/move_size)] += 1
            except IndexError: pass # index out of range at the end of chr. 
                
    return data_dict           
    
def read_sam(filename,data_dict, parameter):

    shift_size = parameter.shift_dict[filename]
    move_size = parameter.window_size/2
    
    infile = open(filename, 'r')
    # skip the header of the SAM file. 
    for line in infile:
        if not line.startswith("@"):
            break
    # start reading the real data
    for line in infile:
        words = line.strip().split()
        flag = int(words[1])
    
        if flag & 0x0004: #if not unmapped
            chr, pos =  words[2], int(words[3])-1
            if flag & 0x0010:
                pos -= shift_size
            else: 
                pos += shift_size
            try: 
                data_dict[chr][int(pos/move_size)] += 1
            except IndexError: pass # index out of range at end of chr.
                
                
    return data_dict 
    
    
def read_bed(filename, data_dict, parameter): 
    
    shift_size = parameter.shift_dict[filename]
    move_size = parameter.window_size/2
    
    infile = open(filename, 'r')
    num = 0
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        num += 1
        if num %1000000 == 0:
            print num
        if strand == "+":
            pos = int(start) + shift_size
        else: 
            pos = int(end) - shift_size
        try: data_dict[chr][int(pos/move_size)] += 1
        except IndexError:    pass
        
    return data_dict
        
def read_file_to_array(filename, parameter):
    ''' read file into arrays that can be used on significance testing. '''
    info ("processing %s", filename)
    move_size = parameter.window_size/2
    data_dict = {}
    for chr in parameter.chr_info:
        row_num = int(parameter.chr_info[chr]/move_size) - 1
        data_dict[chr] = numpy.zeros(row_num, dtype=numpy.float64)
    
    if parameter.file_format == "bed":
        data_dict = read_bed(filename, data_dict, parameter)
    elif parameter.file_format == "bam":
        data_dict = read_bam(filename, data_dict, parameter)
    elif paramter.file_format == "sam":
        data_dict = read_sam(filename, data_dict, parameter)
        
    for chr in parameter.chr_info: 
        data_dict[chr] = data_dict[chr] + numpy.roll(data_dict[chr],-1)
        data_dict[chr] = data_dict[chr] * parameter.normalization_dict[filename]
    return data_dict 
     
def read_file_to_array_wrapper(args):
    return read_file_to_array(*args)
     
def read_files_to_arrays(parameter):
    
    if parameter.num_procs <2:
        for file in parameter.get_filenames():
            parameter.array_dict[file] = read_file_to_array(file, parameter)
    else:
        pool = multiprocessing.Pool(parameter.num_procs)
        results = pool.map(read_file_to_array_wrapper, itertools.izip(parameter.get_filenames(), itertools.repeat(parameter)),1)
    
    return 
        
def prepare_data_peak_calling(readData, parameter):
    data_dict = {}
    for file in parameter.chip1:
        per_file_data_dict = read_file_to_array(file, parameter)
        for chr in parameter.chr_info:
            try: 
                data_dict[chr] = numpy.column_stack((data_dict[chr], per_file_data_dict[chr]))
            except KeyError:
                data_dict[chr] = per_file_data_dict[chr]
    
    for file in parameter.input1:
        per_file_data_dict = read_file_to_array(file, parameter)
        for chr in parameter.chr_info:
            data_dict[chr] = numpy.column_stack((data_dict[chr], per_file_data_dict[chr]))
        
    readData.reads_dict = data_dict
    return 
    
    
def prepare_data_diff_binding(readData, parameter):
    data_dict = {}
    chip1_array = {}
    for file in parameter.chip1:
        per_file_data_dict = read_file_to_array(file, parameter)
        for chr in parameter.chr_info:
            try: 
                chip1_array[chr] = numpy.column_stack((chip1_array[chr], per_file_data_dict[chr]))
            except KeyError:
                chip1_array[chr] = per_file_data_dict[chr]
                
    if len(parameter.input1)> 0:
        input1_array = {}
        if parameter.chip1_matched_input is True: 
            for file in parameter.input1:
                per_file_data_dict = read_file_to_array(file, parameter)
                for chr in parameter.chr_info:
                    try: 
                        input1_array[chr] = numpy.column_stack((input1_array[chr], per_file_data_dict[chr]))
                    except KeyError:
                        input1_array[chr] = per_file_data_dict[chr]                
        else: 
            for file in parameter.input1:
                per_file_data_dict = read_file_to_array(file, parameter)
                for chr in parameter.chr_info:
                    try: 
                        input1_array[chr] += per_file_data_dict[chr]
                    except KeyError:
                        input1_array[chr] = per_file_data_dict[chr] 
            for chr in parameter.chr_info:
                input1_array[chr] /= len(parameter.input1)
                input1_array[chr].resize(len(input1_array[chr]),1)
        for chr in parameter.chr_info:        
            chip1_array[chr] -= input1_array[chr]
            chip1_array[chr][chip1_array[chr]<0] = 0 
    
    # do the same thing for chip2 
    chip2_array = {}
    for file in parameter.chip2:
        per_file_data_dict = read_file_to_array(file, parameter)
        for chr in parameter.chr_info:
            try: 
                chip2_array[chr] = numpy.column_stack((chip2_array[chr], per_file_data_dict[chr]))
            except KeyError:
                chip2_array[chr] = per_file_data_dict[chr]
                
    if len(parameter.input2)> 0:
        input2_array = {}
        if parameter.chip2_matched_input is True: 
            for file in parameter.input2:
                per_file_data_dict = read_file_to_array(file, parameter)
                for chr in parameter.chr_info:
                    try: 
                        input2_array[chr] = numpy.column_stack((input2_array[chr], per_file_data_dict[chr]))
                    except KeyError:
                        input2_array[chr] = per_file_data_dict[chr]                
        else: 
            for file in parameter.input2:
                per_file_data_dict = read_file_to_array(file, parameter)
                for chr in parameter.chr_info:
                    try: 
                        input2_array[chr] += per_file_data_dict[chr]
                    except KeyError:
                        input2_array[chr] = per_file_data_dict[chr] 
            for chr in parameter.chr_info:
                input2_array[chr] /= len(parameter.input2)
                input2_array[chr].resize(len(input2_array[chr]),1)
        for chr in parameter.chr_info:
            chip2_array[chr] -= input2_array[chr]
            chip2_array[chr][chip2_array[chr]<0] = 0    
        
    for chr in parameter.chr_info:
        data_dict[chr] =numpy.column_stack((chip1_array[chr], chip2_array[chr]))
        
    readData.reads_dict = data_dict
    return  
    
    
def prepare_data(readData, parameter):
    ''' The wrapper function to prepare data. 
    Read and process data into arrays. '''
    
    data_by_window_dict = {}
    
    if parameter.difftest == False:
        prepare_data_peak_calling(readData, parameter)
    else: 
        prepare_data_diff_binding(readData, parameter)
        
        
        
        