#from logging import info
import multiprocessing
from scipy.stats import rankdata
from logging import info 
import numpy
import itertools 
import pysam 
import array

from fileParser import parse_file_by_strand

BIN = 10000

def estimate_shiftsizes(parameter):
    ''' the root function for estimating the shiftsizes'''
    
    ## estimate the shift size for chip samples. 
    if parameter.num_procs < 2:
        for chip in parameter.get_chip_filenames():
            shift_size, bin_array = estimate_shiftsize(chip, parameter)
            parameter.shift_dict[chip] = shift_size
            parameter.bin_dict[chip] = bin_array
    else: 
        pool = multiprocessing.Pool(parameter.num_procs)
        p = pool.map_async(estimate_shiftsize_wrapper, itertools.izip(parameter.get_chip_filenames(), itertools.repeat(parameter)),1)
        
        try: results = p.get()
        except KeyboardInterrupt:
            exit(1)
        for chip, result in itertools.izip(parameter.get_chip_filenames(), results):
            parameter.shift_dict[chip] = result[0]
            parameter.bin_dict[chip] = result[1]
           
    # put the shift size for input samples
    if len(parameter.input1) > 0:
        for input in parameter.input1: 
            shift_list = [parameter.shift_dict[file] for file in parameter.chip1]
            parameter.shift_dict[input] = sum(shift_list)/len(shift_list)
    if len(parameter.input2) > 0:
        for input in parameter.input2: 
            shift_list = [parameter.shift_dict[file] for file in parameter.chip2]
            parameter.shift_dict[input] = sum(shift_list)/len(shift_list) 
    return 
def estimate_shiftsize_wrapper(args):
    try:
        return estimate_shiftsize(*args)
    except KeyboardInterrupt, e:
        pass
        
def estimate_shiftsize(chip, parameter):
    ''' estimate the shiftsize for each file '''
    info("estimating the shift size for %s", chip)
    info_dict = {} # saving the info matrix for deriving the shift size data. 
    bin_dict = {}
    
    for chr in parameter.chr_info:
        row_num = int(parameter.chr_info[chr]/BIN) 
        bin_dict[chr] = numpy.zeros(row_num, dtype=numpy.float64)
        info_dict[chr] = numpy.zeros((row_num,4),dtype=numpy.int64)
    # parsing file into strand 
    forward, reverse = parse_file_by_strand[parameter.file_format](chip, parameter)
    shift_list =[]
    for chr in parameter.get_top3_chr():
        chr_f,chr_r = forward[chr],reverse[chr]
        shift_list.append(cross_cor(chr_f,chr_r))
    shift_list.sort()
    frag_size = shift_list[len(shift_list)/2]
    ### parse the reads into bins
    for chr in parameter.chr_info:
        if chr in forward:
            for pos in forward[chr]:
                try: bin_dict[chr][pos/BIN] += 1 
                except IndexError:
                    pass
        if chr in reverse:
            for pos in reverse[chr]:
                try: bin_dict[chr][pos/BIN] += 1
                except IndexError:
                    pass
        try: bin_array = numpy.append(bin_array, bin_dict[chr])
        except UnboundLocalError:
            bin_array = bin_dict[chr]

    frag_size += parameter.read_length_dict[chip]
    shift_size = frag_size/2
    info("The shift size for %s is %d", chip, shift_size)
    return (shift_size, bin_array)

def parse_bed_for_f_r(filename, parameter):
    infile = open(parameter.input_directory+filename, 'r')
    num = 0
    forward = {}
    reverse = {}
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        pos = int(start)
        if strand == "+":
            try: forward[chr].append(pos)
            except KeyError:
                forward[chr] = array.array('i',[pos])
        else: 
            try: reverse[chr].append(pos)
            except KeyError:
                reverse[chr] = array.array('i',[pos])
            
    return forward, reverse
        
def parse_bam_for_f_r(filename, parameter):
    num = 0
    forward = {}
    reverse = {}

    infile =pysam.Samfile(parameter.input_directory+filename, 'rb')
    for line in infile.fetch(until_eof = True):
        num += 1
        if num % 1000000 == 0 :
            print ("{0:,} lines processed in {1}".format(num, filename))
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            if line.is_reverse is False:
                try: forward[chr].append(line.pos)
                except KeyError:
                    forward[chr] = array.array('i',[line.pos])
            else:
                try: reverse[chr].append(line.pos)
                except KeyError:
                    reverse[chr] = array.array('i',[line.pos])
    return forward,reverse
    
def parse_sam_for_f_r(filename, parameter):
    infile = open(parameter.input_directory+filename, 'r')
    num = 0
    forward = {}
    reverse = {}
    # skip the header of the SAM file. 
    for line in infile:
        if not line.startswith("@"):
            break
    # start reading the real data
    for line in infile:
        num += 1
        if num % 10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        words = line.strip().split()
        flag = int(words[1])
    
        if not flag & 0x0004: #if not unmapped
            chr, pos =  words[2], int(words[3])-1
            if not flag & 0x0010: # if not reverse
                try:  forward[chr].apend(pos)
                except KeyError:
                    forward[chr] = array.array('i',[pos])
            else: 
                try:  reverse[chr].append(pos)
                except KeyError:
                    reverse[chr] = array.array('i',[pos])
    return forward, reverse


def cross_cor(f, r):
    npf = numpy.array(f)
    npr = numpy.array(r)
    cor_list = []
    for x in range(50,302,2):
        #print x
        y = len(numpy.intersect1d(npf+x,npr))
        cor_list.append(y)
    return range(50,302,2)[cor_list.index(max(cor_list))]
 
