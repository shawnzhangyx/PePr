#!/usr/bin/env python
import re, os, sys,time
import logging
import numpy

import logConfig
import optParser
import fileParser
import shiftSize
import windowSize
import sigTests
import misc

# initialize the logger
root_logger = logging.getLogger("")
debug = root_logger.debug
info = root_logger.info


def main(argv):
        ## performing the option parser
    opt = optParser.opt_parser(argv)
    parameter, readData = optParser.process_opt(opt)
    ## read in the data
    fileParser.parse(readData, parameter.file_format)

    ## remove the redundant reads
    if (parameter.remove_redundant):
        readData.remove_redundant_reads()

    ## shiftSize estimation and shifting reads
    shiftSize.estimate_shift_size(readData,parameter)
    shiftSize.shift_reads(readData)

    file_in = open(argv[1], 'r')
    file_passed = open(argv[1]+'_passed', 'w')
    file_filtered = open(argv[1]+'_filtered', 'w')
    post_processing(readData,file_in,file_passed,file_filtered)

def post_processing(readData, peak, passed, filtered, remove_artefacts=True, narrow_peak=False):
    info(" Begin post-processing.")
    chip_list =  readData.chip_filename_list
    input_list = readData.input_filename_list
    strands_dict = readData.data_dict_by_strands

    for line in peak: 
        words = line.strip().split('\t')
        chr = words[0]
        start = int(words[1])
        end = int(words[2])
        rest = words[3:]
        (start, end, chip_input_ratio, tag_monopoly, overlap_orig,
        overlap_roll) = post_processing_per_peak(
            strands_dict, chip_list, input_list, chr, start, end,
            readData.shift_size, readData.read_length,
            narrow_peak, remove_artefacts)
        if remove_artefacts and (chip_input_ratio > 0.5 or
                # tag_monopoly > 0.5 or 
                (overlap_orig > 0.2 and
                overlap_roll/overlap_orig < 0.5)):
            filtered.write(chr+"\t" + str(start) + '\t' + str(end) + '\t' +
                        str(end-start)+'\t')
            for item in rest: 
                filtered.write(item +'\t')
            filtered.write(str(chip_input_ratio) + '\t' +
                        str(tag_monopoly) + '\t' + str(overlap_orig) +
                        '\t' + str(overlap_roll) + '\n')
            continue
        passed.write(chr+"\t" + str(start) + '\t' + str(end) + '\t' +
                     str(end-start)+'\t')
        for item in rest:
            passed.write(item +'\t')
        passed.write(str(chip_input_ratio) + '\t' +
                         str(tag_monopoly) + '\t' + str(overlap_orig) +
                         '\t' + str(overlap_roll) + '\n')
    return 


def post_processing_per_peak(strands_dict, chip_list, input_list, chr,
                             start, end, shiftSize, readLength, narrow_peak,
                             remove_artefacts):
    ''' Remove artefacts and refine peak width.'''
    chip_forward = numpy.zeros(end-start)
    chip_reverse = numpy.zeros(end-start)
    input_forward = numpy.zeros(end-start)
    input_reverse = numpy.zeros(end-start)
    for chip in chip_list:
        forward = strands_dict[chr][chip]['f']
        reverse = strands_dict[chr][chip]['r']
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
            forward = strands_dict[chr][input]['f']
            reverse = strands_dict[chr][input]['r']
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
        chip_both.sort()
        chip_3_maximum = numpy.sum(chip_both[-3:])
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
        chip_3_maximum = 0
        overlap_orig = 1e-5
        overlap_roll = 0

    if narrow_peak is True:
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
        end = new_end

    return (start, end, overlap_chip_input, chip_3_maximum,
            overlap_orig, overlap_roll)



if __name__ == '__main__':
    try: main(sys.argv)
    except KeyboardInterrupt:
        print "user interrupted me"


