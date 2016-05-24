#!/usr/bin/env python

###########################################################
### script for removing artifactual peaks in the peak file. 
### Usage: python post_processing_PePr.py --peak peak_file --chip chip.files.sep.by.comma --input input.files.sep.by.comma --file-type type
### Yanxiao Zhang <yanxiazh@umich.edu> 
### Timestamp 5/24/2016
#######################

import re, os, sys,time
import logging
import numpy
from optparse import OptionParser

from PePr.pre_processing.initialize import get_read_length_info 
from PePr.pre_processing.fileParser import parse_file_by_strand

# optionParser
def opt_parser(argv):
    parser = OptionParser()
    parser.add_option(
        '--peak', action='store',type='string',
        dest='peak', default='',
        help='peak file')
    parser.add_option(
        '--chip', action='store',type='string',
        dest='chip',default='',
        help='chip files separated by comma')
    parser.add_option(
        '--input', action='store', type='string',
        dest='input', default='',
        help='input files separated by comma')
    parser.add_option(
        '--file-type', action='store',type='string',
        dest='type', default='',
        help='read file types. bed, sam, bam')
    (opt, args)=parser.parse_args(argv)
    if len(argv) == 1:
        parser.print_help()
        exit(1)
    return opt

def process_opt(opt):
    if opt.chip == '':
        raise Exception("Please specify chip samples.")
    if opt.input == '':
        raise Exception("Please specify input samples.")
    if opt.peak == '':
        raise Exception("No peak files.")
    if opt.type =='':
        raise Exception("File type not given.")
    ## start process chip file names.
    opt.chip = opt.chip.strip().split(',')
    opt.input = opt.input.strip().split(',')
    return opt

class ReadData:
    def __init__(self, opt):
        self.chip_filename_list = opt.chip
        self.input_filename_list = opt.input
        self.data_dict_by_strands = {}
        self.read_length_dict = {}
        self.file_format = opt.type
        self.input_directory = ''  

    def get_read_length(self):
        return self.read_length_dict.values()[0]

    def get_filenames(self):
        return self.chip_filename_list + self.input_filename_list


def main(argv):
    ## performing the option parser
    opt = opt_parser(argv)
    opt = process_opt(opt)
    readData = ReadData(opt)
    get_read_length_info(readData)
    readData.read_length = readData.get_read_length()
    print'read length is {0}'.format(readData.read_length)
    data_dict = {}

    for filename in opt.chip+opt.input:
        print "reading {0}".format(filename)
        data_dict[filename] = {}
        forward, reverse = parse_file_by_strand[readData.file_format](filename,readData)
        for chr in set(forward.keys())&set(reverse.keys()):
            data_dict[filename][chr] = {}
            data_dict[filename][chr]['f'] = numpy.array(forward[chr])
            data_dict[filename][chr]['r'] = numpy.array(reverse[chr])
    readData.data_dict_by_strands = data_dict     
    
    file_in = open(opt.peak, 'r')
    file_passed = open(opt.peak+'_passed', 'w')
    file_failed = open(opt.peak+'_failed', 'w')
    post_processing(readData,file_in,file_passed,file_failed)

def post_processing(readData, peak, passed, filtered, remove_artefacts=True, narrow_peak=False):
    print(" Begin post-processing.")
    chip_list =  readData.chip_filename_list
    input_list = readData.input_filename_list
    strands_dict = readData.data_dict_by_strands
    removed_count = 0
    for idx,line in enumerate(peak):
        if idx%100 == 0:
            print '{0} peaks processed'.format(idx)
        words = line.strip().split('\t')
        chr = words[0]
        start = int(words[1])
        end = int(words[2])
        rest = words[3:]
        (start, end, chip_input_ratio, overlap_orig,
        overlap_roll) = post_processing_per_peak(
            strands_dict, chip_list, input_list, chr, start, end,
            readData.read_length,
            narrow_peak, remove_artefacts)
        # print start,end, chip_input_ratio, overlap_orig, overlap_roll
        if remove_artefacts and (chip_input_ratio > 0.5 or
                (overlap_orig > 0.2 and
                overlap_roll/overlap_orig < 0.5)):
            filtered.write(chr+"\t" + str(start) + '\t' + str(end) + '\t' +
                        str(end-start)+'\t')
            for item in rest: 
                filtered.write(item +'\t')
            filtered.write(str(chip_input_ratio) + '\t' +
                        '\t' + str(overlap_orig) +
                        '\t' + str(overlap_roll) + '\n')
            removed_count += 1
            continue
        passed.write(chr+"\t" + str(start) + '\t' + str(end) + '\t' +
                     str(end-start)+'\t')
        for item in rest:
            passed.write(item +'\t')
        passed.write('\n')
        ## will not write additional information
        #passed.write(str(chip_input_ratio) + '\t' +
        #                  str(overlap_orig) +
        #                 '\t' + str(overlap_roll) + '\n')
    print 'Done. {0} out of {1} peaks failed the test and removed'.format(removed_count, idx+1)
    return 


def post_processing_per_peak(strands_dict, chip_list, input_list, chr,
                             start, end, readLength, narrow_peak,
                             remove_artefacts):
    ''' Remove artefacts that may be caused by PCR duplicates'''

    chip_forward = numpy.zeros(end-start)
    chip_reverse = numpy.zeros(end-start)
    input_forward = numpy.zeros(end-start)
    input_reverse = numpy.zeros(end-start)
    for chip in chip_list:
        forward = strands_dict[chip][chr]['f']
        reverse = strands_dict[chip][chr]['r']
        forward_read = forward[numpy.where( (forward >= start) &
                (forward < end) )]
        reverse_read = reverse[numpy.where( (reverse >= start+readLength) &
                (reverse < end+readLength) )]
        for read in forward_read:
            try: chip_forward[(read-start)] +=1
            except IndexError:
                debug("index error ignored. end-start: %d. read-start: %d.", end-start, read-start)
        for read in reverse_read:
            try: chip_reverse[(read-start-readLength)] +=1
            except IndexError:
                debug("index error ignored. end-start: %d. read-start-readLength: %d.", end-start, read-start-readLength)
    # using input reads to remove artefacts
    if remove_artefacts is True:
        for input in input_list:
            forward = strands_dict[input][chr]['f']
            reverse = strands_dict[input][chr]['r']
            forward_read = forward[numpy.where( (forward >=
                    start) & (forward < end) )]
            reverse_read = reverse[numpy.where( (reverse >=
                    start+readLength) & (reverse < end+readLength) )]
            for read in forward_read:
                try: input_forward[(read-start)] +=1
                except IndexError:
                    debug("index error ignored. end-start: %d. read-start: %d.", end-start, read-start)
            for read in reverse_read:
                try: input_reverse[(read-start-readLength)] +=1
                except IndexError:
                    debug("index error ignored. end-start: %d. read-start-readLength: %d.", end-start, read-start-readLength)

        chip_both = chip_forward + chip_reverse
        input_both = input_forward + input_reverse
        if sum(chip_both) == 0:
            pass
        else:
            chip_both = chip_both/sum(chip_both)
        if sum(input_both) ==0:
            pass
        else:
            input_both = input_both/sum(input_both)
        overlap_chip_input = numpy.sum(numpy.minimum(chip_both, input_both))
        if sum(chip_reverse) != 0:
            chip_reverse = chip_reverse/sum(chip_reverse)
        chip_forward_roll = numpy.roll(chip_forward,readLength)
        if sum(chip_forward) != 0:
            chip_forward = chip_forward/sum(chip_forward)
            chip_forward_roll = chip_forward_roll/sum(chip_forward)
        overlap_orig = numpy.sum(numpy.minimum(chip_forward, chip_reverse))
        overlap_orig = numpy.max([overlap_orig, 1e-5])
        overlap_roll = numpy.sum(numpy.minimum(chip_forward_roll, chip_reverse))
    else:
        overlap_chip_input = 0
        overlap_orig = 1e-5
        overlap_roll = 0

    '''if narrow_peak is True:
        sum_forward = 0
        sum_reverse = 0
        if sum(chip_forward) > 0:
            chip_forward = chip_forward/sum(chip_forward)
            for i in xrange(end-start):
                sum_forward += chip_forward[i]
                if sum_forward > 0.2:
                    new_start = start + i
                    break
        if sum(chip_reverse) > 0:
            chip_reverse = chip_reverse/sum(chip_reverse)
            for i in xrange(end-start-1, -1, -1):
                sum_reverse += chip_reverse[i]
                if sum_reverse > 0.2:
                    new_end = start + i
                    break
        if sum(chip_forward) == 0:
            new_start = new_end - 2*shiftSize
        if sum(chip_reverse) == 0:
            new_end = new_start + 2*shiftSize
        start = new_start
        end = new_end'''

    return (start, end, overlap_chip_input,
            overlap_orig, overlap_roll)



if __name__ == '__main__':
    try: main(sys.argv)
    except KeyboardInterrupt:
        print "user interrupted me"



