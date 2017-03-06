from logging import info
from logging import debug
import array
import pysam 
import numpy
from scipy.stats import rankdata
import multiprocessing
import itertools


BIN = 1000 # bin size for normalization purposes. 

def cal_normalization_constant(parameter):
    info ("calculating normalization constants")
    if parameter.normalization == "scale":
        scale(parameter)
    elif parameter.normalization == "intra-group":
        within_group_tmm(parameter)
    elif parameter.normalization == "inter-group":
        between_group_tmm(parameter)
    else: #no normalization
        no_normalization(parameter)

def no_normalization(parameter):
    for filename in parameter.get_filenames():
        parameter.normalization_dict[filename] = 1.0
    return

def scale(parameter):
    ''' scale the libaries so the total number of reads are the same'''

    if parameter.num_procs < 2:
        for filename in parameter.get_filenames():
            #print filename
            line_num = get_file_line_counts(filename, parameter)
            parameter.lib_size_dict[filename] = line_num
    else:
        pool =  multiprocessing.Pool(parameter.num_procs)
        p = pool.map_async(get_file_line_counts_wrapper, zip(parameter.get_filenames(), itertools.repeat(parameter)),1)
        try: results = p.get()
        except KeyboardInterrupt:
            exit(1)
        for filename, line_num in zip(parameter.get_filenames(), results):
            parameter.lib_size_dict[filename] = line_num

    lib_size_mean = sum(parameter.lib_size_dict.values())/len(parameter.lib_size_dict)
    for filename in parameter.get_filenames():
        parameter.normalization_dict[filename] = float(lib_size_mean)/parameter.lib_size_dict[filename]     
    return 
        
def get_file_line_counts_wrapper(args):
    try:
        return get_file_line_counts(*args)
    except KeyboardInterrupt:
        pass

def get_file_line_counts(filename, parameter):
    # keep the line count. 
    idx = 0
    #print filename
    if parameter.file_format == "bed":
        with open(parameter.input_directory + filename, 'r') as infile: 
            for idx, line in enumerate(infile):
                pass
            idx += 1    
    elif parameter.file_format in ["sam","sampe"]:
        with open(parameter.input_directory + filename, 'r') as infile:
            for line in infile:
                if not line.startswith('@'): #skip the header lines.
                    break
            for line in infile:
                words = line.strip().split()
                flag = int(words[1])
                if not flag & 0x004:
                    idx += 1
    elif parameter.file_format in ["bam","bampe"]:
        with pysam.Samfile(parameter.input_directory + filename, 'rb') as infile:
            for line in infile.fetch(until_eof = True):
                if line.is_unmapped is False:
                    idx += 1
                    
    return idx
        

### within-group normalization.  
def within_group_tmm(parameter):
    '''within_group normalization '''
    # process the reads into bins. 
    bin_dict = parameter.bin_dict 
    bin_size = BIN
    
    if parameter.num_procs < 2: 
        for filename in parameter.get_filenames_wo_bin_dict():
            bin_dict[filename] = parse_to_bin(filename, bin_size, parameter)
    
    else: 
        pool = multiprocessing.Pool(parameter.num_procs)
        p = pool.map_async(parse_to_bin_wrapper, zip(parameter.get_filenames_wo_bin_dict(), itertools.repeat(bin_size), itertools.repeat(parameter)),1)
        try: results = p.get()
        except KeyboardInterrupt:
            exit(1)
            
        for filename, array in zip(parameter.get_filenames_wo_bin_dict(), results):
            parameter.bin_dict[filename] = array
            
    # process chip1 array 
    for chip in parameter.chip1:
        try: 
            chip1_array_mixed = chip1_array_mixed + bin_dict[chip] 
            chip1_rep_rank_sum += rankdata(-bin_dict[chip])
        except UnboundLocalError:
            chip1_array_mixed = bin_dict[chip].copy()
            chip1_rep_rank_sum = rankdata(-bin_dict[chip])
            
    chip1_array_mixed /= len(parameter.chip1)
    
    # process chip2 array if there is any
    if len(parameter.chip2) > 0:
        for chip in parameter.chip2:
            try: 
                chip2_array_mixed = chip2_array_mixed + bin_dict[chip] 
                chip2_rep_rank_sum += rankdata(-bin_dict[chip])
            except UnboundLocalError:
                chip2_array_mixed = bin_dict[chip].copy()
                chip2_rep_rank_sum = rankdata(-bin_dict[chip])
                
        chip2_array_mixed /= len(parameter.chip2)
        # will scale chip1 and chip2 if their lib sizes are very different. 
        chip1_size = numpy.sum(chip1_array_mixed) 
        chip2_size = numpy.sum(chip2_array_mixed)
        chip1_array_mixed *= (chip1_size+chip2_size)/chip1_size/2
        chip2_array_mixed *= (chip1_size+chip2_size)/chip2_size/2
        
    for chip in parameter.chip1: 
        parameter.normalization_dict[chip] = chip_tmm(chip1_array_mixed, bin_dict[chip], chip1_rep_rank_sum)
    
    if len(parameter.chip2) > 0:
        for chip in parameter.chip2:
            parameter.normalization_dict[chip] = chip_tmm(chip2_array_mixed, bin_dict[chip], chip2_rep_rank_sum)
    if len(parameter.input1) > 0:
        for input in parameter.input1:
            parameter.normalization_dict[input] = input_ncis(chip1_array_mixed, bin_dict[input])
    if len(parameter.input2) > 0:
        for input in parameter.input2:
            parameter.normalization_dict[input] = input_ncis(chip2_array_mixed, bin_dict[input])

    return

def between_group_tmm(parameter):
    '''between group tmm method'''
    # process the reads into bins.
    bin_dict = parameter.bin_dict
    bin_size = BIN

    if parameter.num_procs < 2:
        for filename in parameter.get_filenames_wo_bin_dict():
            bin_dict[filename] = parse_to_bin(filename, bin_size, parameter)
    else:
        pool = multiprocessing.Pool(parameter.num_procs)
        p = pool.map_async(parse_to_bin_wrapper, zip(parameter.get_filenames_wo_bin_dict(), itertools.repeat(bin_size), itertools.repeat(parameter)),1)
        try: results = p.get()
        except KeyboardInterrupt:
            exit(1)

        for filename, array in zip(parameter.get_filenames_wo_bin_dict(), results):
            parameter.bin_dict[filename] = array

    # process chip array
    for chip in parameter.get_chip_filenames():
        try:
            chip_array_mixed = chip_array_mixed + bin_dict[chip]
            chip_rep_rank_sum += rankdata(-bin_dict[chip])
        except UnboundLocalError:
            chip_array_mixed = bin_dict[chip].copy()
            chip_rep_rank_sum = rankdata(-bin_dict[chip])

    chip_array_mixed /= len(parameter.get_chip_filenames())

    for chip in parameter.get_chip_filenames():
        parameter.normalization_dict[chip] = chip_tmm(chip_array_mixed, bin_dict[chip], chip_rep_rank_sum)

    if len(parameter.get_input_filenames()) > 0:
        for input in parameter.get_input_filenames():
            parameter.normalization_dict[input] = input_ncis(chip_array_mixed, bin_dict[input])

    return

    
    
def chip_tmm(ref, target, rep_rank_sum): 
    '''Estimate the ChIP normalization constant against the mixed ChIP
       using the modified TMM method'''
    N_PEAKS_GRID = numpy.array([1000, 5000, 10000, 20000, 30000, 40000, 50000])
    TRIM_M = 0.2
    TRIM_A = 0.05
    order = numpy.argsort(rep_rank_sum)
#    order = numpy.argsort(-ref)
    len_target_not_zero = len(numpy.where(target > 0)[0])
    #print len_target_not_zero
    ref = ref[order]
    target = target[order]
    tmm_array = numpy.array([])
    for idx,n in enumerate(N_PEAKS_GRID):
        if n > len_target_not_zero:
            if idx == 0:
                # if the non-zero list is smaller than 1000
                n = len_target_not_zero
            else:
                break
        ref_n = ref[range(n)]
        ref_n[ref_n==0] = 1
        target_n = target[range(n)]
        target_n[target_n==0] = 1
        Mg = numpy.log2(ref_n/target_n)
        Ag = 0.5*numpy.log2(ref_n*target_n)
        Mg_sorted = Mg.copy()
        Mg_sorted.sort()
        Mg_lower_bound = Mg_sorted[int(n*TRIM_M)]
        Mg_upper_bound = Mg_sorted[int(n*(1-TRIM_M))]
        Ag_sorted = Ag.copy()
        Ag_sorted.sort()
        Ag_lower_bound = Ag_sorted[int(n*TRIM_A)]
        Ag_upper_bound = Ag_sorted[int(n*(1-TRIM_A))]
        trim_index = numpy.where((Mg > Mg_lower_bound) &
            (Mg < Mg_upper_bound) & (Ag > Ag_lower_bound) &
            (Ag < Ag_upper_bound))
        ref_trim = ref_n[trim_index]
        target_trim = target_n[trim_index]
        Mgk = numpy.log2(ref_trim/target_trim)
        Wgk = 0.5*numpy.log2(ref_trim*target_trim)
        tmm = 2**(numpy.sum(Mgk*Wgk)/numpy.sum(Wgk))
        tmm_array = numpy.append(tmm_array, tmm)
        #debug("The TMM estiamted from top %s windows is %s", n, tmm)
    library_ratio = numpy.sum(ref)/numpy.sum(target)
    tmm_diff_array = numpy.abs(tmm_array - library_ratio)
    # print tmm_array, tmm_diff_array
    tmm_max = tmm_array[numpy.argmax(tmm_diff_array)]
    return tmm_max

    
def input_ncis(ref, target):
    '''Estimate the input normalization constant using the NCIS method'''
    combined = ref + target 
    READ_MAX = 200
    MIN_GENOME_COVERAGE = 0.75
    pre_ratio = 1
    #file = open(str(sum(target))+"input_norm.txt",'w')
    for r_cut in range(1, READ_MAX):
        index = numpy.where(combined <= r_cut)[0]
        percent_genome_covered = float(len(index))/len(combined)
        target_sum = numpy.sum(target[index]) or 1.0
        ratio = numpy.sum(ref[index])/target_sum
        #file.write(str(r_cut)+'\t'+str(numpy.sum(ref[index]))+ '\t' +
        #        str(target_sum) +'\t'+ str(ratio)+ '\n')
        if percent_genome_covered > MIN_GENOME_COVERAGE and ratio > pre_ratio:
            return ratio
        else: 
            pre_ratio = ratio
    #file.close()
    return ratio
    
def parse_to_bin_wrapper(args):
    try: 
        return parse_to_bin(*args)
    except KeyboardInterrupt:
        pass 
        
    
    
def parse_bed_to_bin(filename, bin_size, bin_dict, input_dir):
    ''' parse the bed files into bin '''
    infile = open(input_dir + filename, 'r')
    num = 0
    for line in infile: 
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        chr,start,end,col3,col4,strand = line.strip().split()
        pos = int(start)           
        try: bin_dict[chr][int(pos/bin_size)] += 1
        except (IndexError, KeyError) as e: pass # index out of range at the end of chr
    
            
    return bin_dict

def parse_bam_to_bin(filename, bin_size, bin_dict, input_dir):
    num = 0
    infile = pysam.Samfile(input_dir + filename, 'rb')
    for line in infile.fetch(until_eof=True):
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            try: 
                bin_dict[chr][int(line.pos/BIN)] += 1
            except (IndexError, KeyError) as e: pass # index out of range at the end of chr. 
    return bin_dict           
    
def parse_sam_to_bin(filename, bin_size, bin_dict, input_dir):
    infile = open(input_dir + filename, 'r')
    # skip the header of the SAM file. 
    for line in infile:
        if not line.startswith("@"):
            break
    # start reading the real data
    num = 0
    for line in infile:
        num += 1
        if num %10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        words = line.strip().split()
        flag = int(words[1])
    
        if not flag & 0x0004: #if not unmapped
            chr, pos =  words[2], int(words[3])-1
            try: 
                bin_dict[chr][int(pos/bin_size)] += 1
            except (IndexError, KeyError) as e: pass # index out of range at end of chr.
                
    return bin_dict 
    
def parse_sampe_to_bin(filename, bin_size, bin_dict, input_dir):
    infile = open(input_dir+filename, 'r')
    num = 0
    reads_dict = {}
    flen_dict = {}
    flen_list = array.array('i',[])
    # skip the header of the SAM file.
    for line in infile:
        if not line.startswith("@"):
            pre_name = line.strip().split()[0]
            break
    line_saved = False
    # start reading the real data
    for line in infile:
        num += 1
        if num % 10000000 == 0:
            print("{0:,} lines processed in {1}".format(num, filename))
        words = line.strip().split()
        name = words[0]
        # if the sequence name has already been processed
        if name == pre_name and line_saved == True:
                continue
        # else
        #initialize the condition
        pre_name = name
        line_saved = False
        # process the new sequence. 
        flag = int(words[1])
        if not flag & 0x0004: #if not unmapped
            chr, pos,flen = words[2], int(words[3])-1, int(words[8])
            if flen != 0: # if the other end mapped to the same chromosome
                try:
                    reads_dict[chr].append(pos+flen/2)
                    flen_dict[chr].append(flen)
                except KeyError:
                    reads_dict[chr] = array.array('i',[pos+flen/2])
                    flen_dict[chr] = array.array('i', [flen])
                flen_list.append(abs(flen))
                line_saved = True

    print(num)
    # will calculate the median fragment size and remove the reads that have larger fragment size than it. 
    flen_median = numpy.median(flen_list)
    # print flen_median

    for chr in reads_dict:
        flen_chr = numpy.array(flen_dict[chr])
        reads_chr = numpy.array(reads_dict[chr])
        reads_chr = reads_chr[numpy.where(numpy.abs(flen_chr) <= 2*flen_median)[0]]
        for pos in reads_chr:
            try:
                bin_dict[chr][int(pos/bin_size)] += 1
            except (IndexError, KeyError) as e: pass # index out of range
    
    return bin_dict

def parse_bampe_to_bin(filename, bin_size, bin_dict, input_dir):
    '''parse paired-end bam file'''

    infile =pysam.Samfile(input_dir+filename, 'rb')
    num = 0
    reads_dict = {}
    flen_dict = {}
    flen_list = array.array('i',[])
    # start proccessing the data.
    line_saved = False
    pre_name = ''
    for line in infile.fetch(until_eof = True):
        num += 1
        if num % 10000000 == 0 :
            print ("{0:,} lines processed in {1}".format(num, filename))
        name = line.query_name
        if name==pre_name and line_saved == True:
            continue
        # else 
        #initialize the condition
        pre_name = name
        line_saved = False
        if line.is_unmapped is False:
            chr = infile.getrname(line.tid)
            # pos = line.pos
            # flen = line.tlen
            if line.tlen !=0:
                try: 
                    reads_dict[chr].append(line.pos+line.tlen/2)
                    flen_dict[chr].append(line.tlen)
                except KeyError:
                    reads_dict[chr] = array.array('i',[line.pos+line.tlen/2])
                    flen_dict[chr] = array.array('i',[line.tlen])
                flen_list.append(abs(line.tlen))
                line_saved = True
    print(num)
    # calculate the median fragment size and remove reads that are two times larger.
    flen_median = numpy.median(flen_list)
    for chr in reads_dict:
        flen_chr = numpy.array(flen_dict[chr])
        reads_chr = numpy.array(reads_dict[chr])
        reads_chr = reads_chr[numpy.where(numpy.abs(flen_chr) <= 2*flen_median)[0]]
        for pos in reads_chr:
            try:
                bin_dict[chr][int(pos/bin_size)] += 1
            except (IndexError, KeyError) as e: pass # index out of range

    return bin_dict
    
def parse_to_bin(filename, bin_size, parameter):
    bin_dict = {}
    for chr in parameter.chr_info:
        row_num = int(parameter.chr_info[chr]/bin_size)
        bin_dict[chr] = numpy.zeros(row_num, dtype=numpy.float64)
        
    parse_to_bin_dict = {'bed':parse_bed_to_bin,'bam':parse_bam_to_bin,'sam':parse_sam_to_bin,'sampe':parse_sampe_to_bin, 'bampe':parse_bampe_to_bin}

    bin_dict = parse_to_bin_dict[parameter.file_format](filename, bin_size, bin_dict, parameter.input_directory)

    for chr in parameter.chr_info:
        try: 
            bin_array = numpy.append(bin_array, bin_dict[chr])
        except UnboundLocalError: # if bin_array does not exist.
            bin_array = bin_dict[chr]  
            
    return bin_array
    
