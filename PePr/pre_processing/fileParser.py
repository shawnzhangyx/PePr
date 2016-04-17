from collections import Counter
from logging import info
import pysam
#from memory_profiler import profile

BIN = 1000

def bam_parse(filename):
    output = {}
        
    infile = pysam.Samfile(filename, 'rb')
    for line in infile.fetch():
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            if line.is_reverse is False:
                strand = "+"
            else:
                strand = "-"
            try: 
                output[(chr,line.pos)].append((line.pos%BIN,strand))
            except KeyError: 
                output[(chr,line.pos)] = [(line.pos%BIN,strand)]
                
    return output           
    
def sam_parse(filename):
    output = {}
    
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
                strand = "-"
            else: 
                strand = "+"
            try: 
                output[(chr,line.pos)].append((line.pos%BIN,strand))
            except KeyError: 
                output[(chr,line.pos)] = [(line.pos%BIN,strand)]
                
    return output 
            
def bed_parse(filename): 
    
    output = {}
    infile = open(filename, 'r')
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        if strand == "+":
            pos = int(start) 
        else:
            pos = int(end)
        try: 
            output[(chr, pos/BIN)].append((pos%BIN,strand))
        except KeyError:
            output[(chr,pos/BIN)] = [(pos%BIN,strand)]
        
        
    return output
    

def parse(parameter, filename):
    info ("start reading %s", filename)
    if parameter.file_format == "bed":
        data = bed_parse(filename)
    elif parameter.file_format == "bam":
        data = bam_parse(filename)
    elif parameter.file_format == "sam":
        data = sam_parse(filename)
    #info ("finished reading %s", filename)
    return data
    
def bam_parse_to_bin(filename):
    output = {}
        
    infile = pysam.Samfile(filename, 'rb')
    for line in infile.fetch():
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            try: 
                output[(chr,line.pos)] += 1
            except KeyError: 
                output[(chr,line.pos)] = 1
                
    return output           
    
def sam_parse_to_bin(filename):
    output = {}
    
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
            try: 
                output[(chr,line.pos)] += 1
            except KeyError: 
                output[(chr,line.pos)] = 1
                
    return output 
    
def bed_parse_to_bin(filename, bin_dict):
    # parse the bed files into bed, do not estimate the shift sizes. 

    with open(filename, 'r') as infile:
        #num = 0 
        for line in infile: 
            chr,start,end,col3,col4,strand = line.strip().split()
            if strand == "+":
                pos = int(start) 
            else:
                pos = int(end)
            try: 
                bin_dict[chr][int(pos/BIN)] += 1
            except IndexError:      pass
                
    return bin_dict
    
    
    
def parse_to_bin(parameter, filename):
    ''' read file into arrays that can be used on estimate the normalization.'''
    bin_dict = {}
    for chr in parameter.chr_info:
        row_num = int(parameter.chr_info[chr]/BIN) - 1
        bin_dict[chr] = numpy.zeros(row_num, dtype = numpy.float64)
    
    info ("start reading %s", filename)
    if parameter.file_format == "bam":
        bin_dict = bam_parse_to_bin(filename, bin_dict)
    if parameter.file_format == "sam":
        bin_dict = sam_parse_to_bin(filename, bin_dict)
    if parameter.file_format == "bed":
        bin_dict = bed_parse_to_bin(filename, bin_dict)
        
    return data
        
    
    
        