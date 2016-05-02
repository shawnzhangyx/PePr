#from logging import info
import multiprocessing
from scipy.stats import rankdata
from logging import info 
import numpy
import itertools 
import pysam 
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
    
    if parameter.file_format == "bed":
        info_dict, bin_dict = parse_bed_for_shift_bin(chip, info_dict, bin_dict, parameter)
    elif parameter.file_format == "bam":
        info_dict, bin_dict = parse_bam_for_shift_bin(chip, info_dict, bin_dict, parameter)
    elif parameter.file_format == "sam":
        info_dict, bin_dict = parse_sam_for_shift_bin(chip, info_dict, bin_dict, parameter)
    
    #info_array = numpy.array([],dtype=numpy.float64)
    #bin_array = numpy.array([],dtype=numpy.float64)
    for idx,chr in enumerate(parameter.chr_info):
        try: 
            info_array = numpy.row_stack((info_array, info_dict[chr]))
            bin_array = numpy.append(bin_array, bin_dict[chr])
        except UnboundLocalError:
            info_array = info_dict[chr]
            bin_array = bin_dict[chr]
    
    #using the top 2000 peaks or top 10% of the windows to estimate the shift size, whichever is smaller. 
    SIZE = 20000
    top = min([SIZE, int(0.1*len(bin_array))])
    #print top
    rank = rankdata(-bin_array)
    order = numpy.argsort(rank)
    info_array_top = info_array[order][range(top)]
    # get rid of windows that have only reads coming from one strand
    info_array_top = info_array_top[numpy.where(info_array_top[:,1]>0)]
    info_array_top = info_array_top[numpy.where(info_array_top[:,3]>0)]
    shift_array = info_array_top[:,2]/info_array_top[:,3] - info_array_top[:,0]/info_array_top[:,1]
    
    ### output shift_size
    with open(chip+'.shift.txt', 'w') as fileout:
        for shift in shift_array:
            fileout.write(str(shift)+'\n')
           
    
    frag_size = int(numpy.median(shift_array))
    if frag_size < 0:
        frag_size = 0
    frag_size += parameter.read_length_dict[chip]
    shift_size = frag_size/2
    info("The shift size for %s is %d", chip, shift_size)
    return (shift_size, bin_array)

def parse_bed_for_shift_bin(filename, info_dict, bin_dict, parameter):
    infile = open(parameter.input_directory+filename, 'r')
    num = 0
    for line in infile: 
        chr,start,end,col3,col4,strand = line.strip().split()
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        pos = int(start)
        if strand == "+":
            try: 
                info_dict[chr][int(pos/BIN),0] += pos%BIN
                info_dict[chr][int(pos/BIN),1] += 1
            except (IndexError, KeyError) as e:    pass
            
        else: 
            try: 
                info_dict[chr][int(pos/BIN),2] += pos%BIN
                info_dict[chr][int(pos/BIN),3] += 1
            except (IndexError, KeyError) as e:    pass
            
        try: bin_dict[chr][int(pos/BIN)] += 1
        except (IndexError, KeyError):    pass
        
    return info_dict, bin_dict
        
def parse_bam_for_shift_bin(filename, info_dict, bin_dict, parameter):
    num = 0 
    infile = pysam.Samfile(parameter.input_directory+filename, 'rb')
    for line in infile.fetch(until_eof=True):
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            if line.is_reverse is False:
                try:
                    info_dict[chr][int(line.pos/BIN),0] += line.pos%BIN
                    info_dict[chr][int(line.pos/BIN),1] += 1
                except (IndexError, KeyError) as e: pass
            else:
                try:
                    info_dict[chr][int(line.pos/BIN),2] += line.pos%BIN
                    info_dict[chr][int(line.pos/BIN),3] += 1
                except (IndexError, KeyError) as e: pass

            try: 
                bin_dict[chr][int(line.pos/BIN)] += 1
            except (IndexError, KeyError) as e: pass # index out of range at the end of chr. 
                
    return info_dict, bin_dict           
    
def parse_sam_for_shift_bin(filename,info_dict, bin_dict, parameter):
    infile = open(parameter.input_directory+filename, 'r')
    num = 0
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
            if flag & 0x0010: # if reverse
                try:
                    info_dict[chr][int(pos/BIN),2] += pos%BIN
                    info_dict[chr][int(pos/BIN),3] += 1
                except (IndexError, KeyError) as e:    pass
            
            else: 
                try:
                    info_dict[chr][int(pos/BIN),0] += pos%BIN
                    info_dict[chr][int(pos/BIN),1] += 1
                except (IndexError, KeyError) as e: pass
            try: 
                bin_dict[chr][int(pos/BIN)] += 1
            except (IndexError, KeyError) as e: pass # index out of range at end of chr.
                
                
    return info_dict, bin_dict
            
