from collections import Counter
import logging


info = logging.info
bin = 1000

def bed_parse(filename):
    
    output = {}
    infile = open(filename, 'r')
    ##info ("reading sequences from file: %s", filename)
    #count = 1000000 will not report lines. 
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        pos = int(start)+int(end)/2
        try: 
            output[(chr, pos/bin)].append((pos%bin,strand))
        except KeyError:
            output[(chr,pos/bin)] = [(pos%bin,strand)]
        
        
    return output
    


def parse(parameter, filename):
    info ("start reading %s", filename)
    if parameter.file_format == "bed":
        data = bed_parse(filename)
    
    
    info ("finished reading %s", filename)
    return data