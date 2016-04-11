from collections import Counter
import logging
#from memory_profiler import profile

info = logging.info
BIN = 1000

#@profile
def bed_parse(filename):
    
    output = {}
    infile = open(filename, 'r')
    ##info ("reading sequences from file: %s", filename)
    #count = 1000000 will not report lines. 
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        pos = int(start)+int(end)/2
        try: 
            output[(chr, pos/BIN)].append((pos%BIN,strand))
        except KeyError:
            output[(chr,pos/BIN)] = [(pos%BIN,strand)]
        
        
    return output
    
def bed_parse_to_bin(filename):
    output = {}
    infile = open(filename, 'r')
    ##info ("reading sequences from file: %s", filename)
    #count = 1000000 will not report lines. 
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        pos = int(start)+int(end)/2
        try: 
            output[(chr, pos/BIN)] += 1
        except KeyError:
            output[(chr,pos/BIN)] = 1
        
        
    return output

    

def parse(parameter, filename):
    info ("start reading %s", filename)
    if parameter.file_format == "bed":
        data = bed_parse(filename)
    
    
    info ("finished reading %s", filename)
    return data
    

def reads_to_bin(data):
    counter = Counter()
    for key in data:
        counter[key] = len(data[key])
    return counter
    
    
def get_chr_info_bed(filename):
    chr_info = {}
    infile = open(filename, 'r')

    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        try: 
            chr_info[chr].append(int(end))
        except KeyError:
            chr_info[chr] = [int(end)]
            
    for chr in chr_info: 
        chr_info[chr] = max(chr_info[chr])
        info("length of %s is %d", chr, chr_info[chr] )
        
    return chr_info


def get_chromosome_info(parameter, chip_filename):
    info ("getting chromosome info")
    if parameter.file_format == "bed":
        parameter.chr_info = get_chr_info_bed(chip_filename)
        
    return 

    
        