import numpy
from logging import info
import pysam
import multiprocessing
import itertools 
#from pre_processing.shiftSize import parse_bed_for_f_r,parse_bam_for_f_r,parse_sam_for_f_r
from pre_processing.fileParser import parse_file_by_strand

def read_bam(filename, data_dict, parameter):
    shift_size = parameter.shift_dict[filename]
    read_length = parameter.read_length_dict[filename]
    move_size = parameter.window_size/2
    num = 0 
    infile = pysam.Samfile(parameter.input_directory+filename, 'rb')
    for line in infile.fetch(until_eof = True):
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            if line.is_reverse is False:
                pos = line.pos + shift_size
            else:
                pos = line.pos + read_length- shift_size
            try: 
                data_dict[chr][int(pos/move_size)] += 1
            except (IndexError, KeyError) as e: pass # index out of range at the end of chr. 
                
    return data_dict           
    
def read_sam(filename,data_dict, parameter):

    shift_size = parameter.shift_dict[filename]
    read_length = parameter.read_length_dict[filename]
    move_size = parameter.window_size/2
    num = 0
    infile = open(parameter.input_directory+filename, 'r')
    # skip the header of the SAM file. 
    for line in infile:
        if not line.startswith("@"):
            break
    # start reading the real data
    for line in infile:
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        words = line.strip().split()
        flag = int(words[1])
    
        if not flag & 0x0004: #if not unmapped
            chr, pos =  words[2], int(words[3])-1
            if flag & 0x0010: #if reversed
                pos = pos + read_length - shift_size
            else: 
                pos += shift_size
            try: 
                data_dict[chr][int(pos/move_size)] += 1
            except (IndexError, KeyError) as e: pass # index out of range at end of chr.
                
                
    return data_dict 
    
    
def read_bed(filename, data_dict, parameter): 
    
    shift_size = parameter.shift_dict[filename]
    move_size = parameter.window_size/2
    
    infile = open(parameter.input_directory+filename, 'r')
    num = 0
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        if strand == "+":
            pos = int(start) + shift_size
        else: 
            pos = int(end) - shift_size
        try: data_dict[chr][int(pos/move_size)] += 1
        except (IndexError, KeyError) as e:    pass
        
    return data_dict

def remove_duplicate(forward, reverse, max):
    for chr in set(forward.keys())&set(reverse.keys()):
        forward_chr = forward[chr]
        forward_chr.sort()
        idx_keep = numpy.zeros(forward_chr.size)
        pos_pre = 0
        count = 1
        for idx,pos in enumerate(forward_chr):
            if pos == pos_pre:
                count += 1
                if count > max:
                    idx_keep[idx] = count 
            else:
                count = 1
                pos_pre = pos
        forward[chr] = forward_chr[numpy.where(idx_keep==0)]
        # do the same for the reverse strand. 
        reverse_chr = reverse[chr]
        reverse_chr.sort()
        idx_keep = numpy.zeros(reverse_chr.size)
        pos_pre = 0
        count = 1
        for idx,pos in enumerate(reverse_chr):
            if pos == pos_pre:
                count += 1
                if count > max:
                    idx_keep[idx] = count
            else:
                count = 1
                pos_pre = pos
        reverse[chr] = reverse_chr[numpy.where(idx_keep==0)]
        #print "before removing:",chr,str(len(forward_chr))
        #print "after:", chr, str(len(forward[chr]))
    return forward, reverse
        

def read_file_to_array(filename, parameter):
    ''' read file into arrays that can be used on significance testing. '''
    info ("processing %s", filename)
    move_size = parameter.window_size/2
    data_dict = {}
    for chr in parameter.chr_info:
        row_num = int(parameter.chr_info[chr]/move_size) - 1
        data_dict[chr] = numpy.zeros(row_num, dtype=numpy.float64)
    forward, reverse = parse_file_by_strand[parameter.file_format](filename,parameter)
    for chr in forward:
        forward[chr] = numpy.array(forward[chr])
    for chr in reverse:
        reverse[chr] = numpy.array(reverse[chr])
    if parameter.keep_max_dup > 0:
        forward, reverse = remove_duplicate(forward,reverse,parameter.keep_max_dup)
    for chr in forward:
        for pos in forward[chr]+parameter.shift_size:
            try:
                data_dict[chr][int(pos/move_size)] += 1
            # index out of range at the end of chr.
            except (IndexError, KeyError) as e: pass 
    for chr in reverse:
        for pos in reverse[chr]-parameter.shift_size+parameter.read_length_dict[filename]:
            try:
                data_dict[chr][int(pos/move_size)] += 1
            except (IndexError, KeyError) as e: pass 
        
    for chr in parameter.chr_info: 
        data_dict[chr] = data_dict[chr] + numpy.roll(data_dict[chr],-1)
        data_dict[chr] = data_dict[chr] * parameter.normalization_dict[filename]
    return data_dict 
     
def read_file_to_array_wrapper(args):
    try: 
        return read_file_to_array(*args)
    except KeyboardInterrupt, e:
        pass
     
def read_files_to_arrays(parameter):
    
    if parameter.num_procs <2:
        for filename in parameter.get_filenames():
            parameter.array_dict[filename] = read_file_to_array(filename, parameter)
    else:
        pool = multiprocessing.Pool(parameter.num_procs)
        p = pool.map_async(read_file_to_array_wrapper, itertools.izip(parameter.get_filenames(), itertools.repeat(parameter)),1)
        try: results = p.get()
        except KeyboardInterrupt:
            exit(1)
        for filename, result in itertools.izip(parameter.get_filenames(), results):
            parameter.array_dict[filename] = result
            
    return 
        
def prepare_data_peak_calling(parameter):
    info ("Begin peak-calling")
    data_dict = {}
    for filename in parameter.chip1:
        for chr in parameter.chr_info:
            try: 
                data_dict[chr] = numpy.column_stack((data_dict[chr], parameter.array_dict[filename][chr]))
            except KeyError:
                data_dict[chr] = parameter.array_dict[filename][chr]
    
    for filename in parameter.input1:
        for chr in parameter.chr_info:
            data_dict[chr] = numpy.column_stack((data_dict[chr], parameter.array_dict[filename][chr]))
        
    
    return data_dict
    
    
def prepare_data_diff_binding(parameter):
    info ("Begin differential binding analysis")
    data_dict = {}
    chip1_array = {}
    for filename in parameter.chip1:
        for chr in parameter.chr_info:
            try: 
                chip1_array[chr] = numpy.column_stack((chip1_array[chr], parameter.array_dict[filename][chr]))
            except KeyError:
                chip1_array[chr] = parameter.array_dict[filename][chr]
                
    if len(parameter.input1)> 0:
        input1_array = {}
        if parameter.chip1_matched_input is True: 
            for filename in parameter.input1:
                for chr in parameter.chr_info:
                    try: 
                        input1_array[chr] = numpy.column_stack((input1_array[chr], parameter.array_dict[filename][chr]))
                    except KeyError:
                        input1_array[chr] = parameter.array_dict[filename][chr]     
        else: 
            for filename in parameter.input1:
                for chr in parameter.chr_info:
                    try: 
                        input1_array[chr] += parameter.array_dict[filename][chr]
                    except KeyError:
                        input1_array[chr] = parameter.array_dict[filename][chr] 
            for chr in parameter.chr_info:
                input1_array[chr] /= len(parameter.input1)
                input1_array[chr].resize(len(input1_array[chr]),1)
        for chr in parameter.chr_info:        
            chip1_array[chr] -= input1_array[chr]
            chip1_array[chr][chip1_array[chr]<0] = 0 
    
    # do the same thing for chip2 
    chip2_array = {}
    for filename in parameter.chip2:
        for chr in parameter.chr_info:
            try: 
                chip2_array[chr] = numpy.column_stack((chip2_array[chr], parameter.array_dict[filename][chr]))
            except KeyError:
                chip2_array[chr] = parameter.array_dict[filename][chr]
                
    if len(parameter.input2)> 0:
        input2_array = {}
        if parameter.chip2_matched_input is True: 
            for filename in parameter.input2:
                for chr in parameter.chr_info:
                    try: 
                        input2_array[chr] = numpy.column_stack((input2_array[chr], parameter.array_dict[filename][chr]))
                    except KeyError:
                        input2_array[chr] = parameter.array_dict[filename][chr]                
        else: 
            for filename in parameter.input2:
                for chr in parameter.chr_info:
                    try: 
                        input2_array[chr] += parameter.array_dict[filename][chr]
                    except KeyError:
                        input2_array[chr] = parameter.array_dict[filename][chr] 
            for chr in parameter.chr_info:
                input2_array[chr] /= len(parameter.input2)
                input2_array[chr].resize(len(input2_array[chr]),1)
        for chr in parameter.chr_info:
            chip2_array[chr] -= input2_array[chr]
            chip2_array[chr][chip2_array[chr]<0] = 0    
        
    for chr in parameter.chr_info:
        data_dict[chr] =numpy.column_stack((chip1_array[chr], chip2_array[chr]))
        
    return data_dict
    
    
def prepare_data(parameter):
    ''' The wrapper function to prepare data. 
    Arrange the data into arrays for significance testing. '''
    if parameter.difftest == False:
        read_dict = prepare_data_peak_calling(parameter)
    else: 
        read_dict = prepare_data_diff_binding(parameter)
    # remove the array dict when it is read. 
    del parameter.array_dict 
    for chr in parameter.chr_info:
        read_array = read_dict[chr]
        read_array[numpy.where(read_array ==0)] = 1
    return read_dict
        
        
        
