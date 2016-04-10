import numpy
from logging import info

def read_bed(filename, data_dict, parameter): 
    
    shift_size = parameter.shift_dict[filename]
    move_size = parameter.window_size/2
    
    infile = open(filename, 'r')
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
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
        
    for chr in parameter.chr_info: 
        data_dict[chr] = data_dict[chr] + numpy.roll(data_dict[chr],-1)
        data_dict[chr] = data_dict[chr] * parameter.normalization_dict[filename]
    return data_dict 
     
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
    pass

    
    
def prepare_data(readData, parameter):
    ''' The wrapper function to prepare data. 
    Read and process data into arrays. '''
    
    data_by_window_dict = {}
    
    if parameter.difftest == "false":
        prepare_data_peak_calling(readData, parameter)
    else: 
        prepare_data_diff_binding(readData, parameter)
        
        
        
        