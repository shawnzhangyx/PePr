### module for preprocessing the data to estimate the shift size and normalization constant. 
#from collections import Counter
from logging import info
import pysam

import shiftSize
import cal_normalization
import windowSize

def preprocess(parameter):
    ''' read file and estimate shift sizes and normazliation parameters. 
    '''
    # use one chip sample estimate the chromosome size
    chip_filename = parameter.chip1[0]
    get_chromosome_info(parameter, chip_filename)
    get_read_length_info(parameter)
    if parameter.window_size == -1:
        windowSize.estimate_window_size(chip_filename, parameter)
    # if shift size is not available, estimate the shift size for chip1 and input1
    if parameter.shift_size != -1:
        for filename in parameter.get_filenames():
                parameter.shift_dict[filename] = parameter.shift_size 
    elif len(parameter.shift_dict) is 0:
        shiftSize.estimate_shiftsizes(parameter)
    # if normalization constants are not available, estimate it.     
    if len(parameter.normalization_dict) is 0:
        cal_normalization.cal_normalization_constant(parameter)
    
    del parameter.bin_dict 
    return 
    

    
def get_chr_info_bam(parameter, filename):
    chr_info = {}
    num = 0
    infile = pysam.Samfile(parameter.input_directory+filename, 'rb')
    for line in infile.fetch(until_eof=True):
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            try: 
                chr_info[chr] = max(chr_info[chr], line.pos)
            except KeyError: 
                chr_info[chr] = line.pos
    parameter.chr_info = chr_info
    return 
    
def get_chr_info_sam(parameter, filename):
    chr_info = {}
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
            try: 
                chr_info[chr] = max(chr_info[chr], pos)
            except KeyError: 
                chr_info[chr] = pos
    parameter.chr_info = chr_info
    return  
    
def get_chr_info_bed(parameter, filename):
    chr_info = {}
    infile = open(parameter.input_directory+filename, 'r')
    num = 0
    for line in infile: 
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        chr,start,end,col3,col4,strand = line.strip().split()
        try: 
            chr_info[chr] = max(chr_info[chr], int(start))
        except KeyError:
            chr_info[chr] = int(start)
    parameter.chr_info = chr_info
    return


def get_chromosome_info(parameter, chip_filename):
    info ("getting chromosome info")
    if parameter.file_format == "bam":
         get_chr_info_bam(parameter, chip_filename)   
    if parameter.file_format == "sam":
         get_chr_info_sam(parameter, chip_filename)        
    if parameter.file_format == "bed":
         get_chr_info_bed(parameter, chip_filename)
    message = ' | '.join([str(key)+":"+str(value) for key,value in parameter.chr_info.iteritems()])
    info("Found "+message)
    return 

def get_read_length_info(parameter):
    if parameter.file_format == "bam":
        get_read_length_from_bam(parameter)
    if parameter.file_format == "sam":
        get_read_length_from_sam(parameter)
    if parameter.file_format == "bed":
        get_read_length_from_bed(parameter)
    
    message = '|'.join([ str(x) for x in parameter.read_length_dict.values()])
    info("Read length: "+message)

def get_read_length_from_bam(parameter):
    for filename in parameter.get_filenames():
        with pysam.Samfile(parameter.input_directory+filename, 'rb') as infile:
            length_list = []
            for idx in range(1000):
                line = infile.fetch(until_eof=True).__next__()
                length_list.append(line.alen)
            length = max(length_list)
        parameter.read_length_dict[filename] = length
    return

def get_read_length_from_sam(parameter):
    for filename in parameter.get_filenames():
        with open(parameter.input_directory+filename, 'r') as infile:
            idx = 0
            length_list = []
            for line in infile:
                words = line.strip().split()
                if not words[0].startswith("@"):
                    length_list.append(len(words[9]))
                    idx += 1
                    if idx == 1000:
                        break
            length = max(length_list)
        parameter.read_length_dict[filename] = length
    return 

def get_read_length_from_bed(parameter):
    for filename in parameter.get_filenames():
        with open(parameter.input_directory+filename, 'r') as infile:
            length_list = []
            for idx in range(1000):
                line = infile.readline()
                chr,start,end,col3,col4,strand = line.strip().split()
                length_list.append(int(end)-int(start))
            length = max(length_list)
        parameter.read_length_dict[filename] = length
    return


