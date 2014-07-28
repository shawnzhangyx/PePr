#!/usr/bin/env python

import sys
import re
import math
import random
import numpy
import numpy.lib.stride_tricks as std
from scipy.stats import rankdata
import misc
import logging

root_logger = logging.getLogger("")
debug = root_logger.debug
info = root_logger.info
error = root_logger.error

window_logger = logging.getLogger("windowSizeEst")

def get_window_size2(array, bin=20, iter=100):
    # this is the function that estimate the window size. 
    chr_len = array.size
    rowNum = chr_len/bin
    array_window = std.as_strided(array, (rowNum, bin), 
            (bin*array.itemsize, 1*array.itemsize))
    array_window = numpy.sum(array_window, 1)
    peak_len_list = []
    for x in range(iter):
        peak = numpy.max(array_window)
        peak_idx = numpy.where(array_window == peak)[0][0]
        # set the peak to -1 so that it won't confuse the next iteration.
        array_window[peak_idx]= -1 
        if peak_idx == 0:  # Check if the window reaches the left boundary.
            left_boundary_not_reached = False
        else:
            left_boundary_not_reached = True    
        # Check if the window reaches of the right boundary.
        if peak_idx == rowNum-1: 
            right_boundary_not_reached = False
        else: 
            right_boundary_not_reached = True
        i_l = 1
        i_r = 1 
        # Be careful if anything equals to -1. 
        # Though it's unlikely that it will happen.
        while left_boundary_not_reached:

            while (array_window[peak_idx-i_l] >= 0.1*peak): 
                array_window[peak_idx-i_l] = -1
                if peak_idx-i_l==0:
                    left_boundary_not_reached = False
                    break
                i_l += 1
            if left_boundary_not_reached:
                if peak_idx-i_l==0:
                    array_window[peak_idx-i_l] = -1
                    left_boundary_not_reached = False
            
            if left_boundary_not_reached: 
            # Will continue counting if the next window(one window gap) 
            # has above the 10% mode reads. 
                if (array_window[peak_idx-i_l-1] >= 0.1*peak): 
                    array_window[peak_idx-i_l] = -1
                    i_l += 1
                else:
                    break

        while right_boundary_not_reached:
            while (array_window[peak_idx+i_r] >= 0.1*peak):
                array_window[peak_idx+i_r] = -1
                if peak_idx+i_r ==rowNum-1:
                    right_boundary_not_reached = False
                    break
                i_r += 1
            if right_boundary_not_reached:
                if peak_idx+i_r==rowNum-1:
                    array_window[peak_idx+i_r] = -1
                    right_boundary_not_reached = False

            if right_boundary_not_reached:
                if (array_window[peak_idx+i_r+1] >=0.1*peak): 
                    array_window[peak_idx+i_r] = -1
                    i_r += 1
                else:
                    break
        peak_len_list.append((1+i_l+i_r)*bin)

    return misc.median(peak_len_list)

def estimate_window_size(readData, parameter):
    if (parameter.window_size != -1):
        # window size should be even number 
        if parameter.window_size % 2 ==1:
            info("warning: window size should be even number; adding 1 to it.")
            parameter.window_size = parameter.window_size + 1
        return parameter.window_size

    info( "begin window_size calculation...")
    window_size_list = []
    for chr in readData.chr_list:
        reads = []
        for chip_filename in readData.chip_filename_list:
            reads.extend(readData.data_dict[chr][chip_filename])
        coord_array = numpy.zeros(readData.chr_length_dict[chr],
                dtype=numpy.int64)
        for x in reads:
            try: coord_array[x] += 1
            except IndexError: pass #debug("coordinates out of range")

        window = get_window_size2(coord_array)
#        if window >1000:  # limit the window size to be less than 1000 bases.
#            window = 1000
        window_size_list.append(window)
        window_logger.debug("%-10s %s", chr, window)
    debug("length of windowsize list is "+str(len(window_size_list)))
    window_median = misc.median(window_size_list)
    # Window size cannot be odd number since we need 
    # to use overlapping half windows. 
    if window_median % 2 == 1:
        window_median = window_median+1
    window_logger.debug("window_size: %+10s", window_median)
    info ( "finishing window size calculation...")
    parameter.window_size = window_median
    return  


def count_reads_by_file(read_data, bin, row_num):
    array = numpy.zeros(row_num, dtype=numpy.float64)
    for x in read_data:
        try: array[x/bin] += 1
        except IndexError:    pass  # out of range
    if any(array<0):
        error("array has reads less than 0")
    array = array + numpy.roll(array, -1)
    return array


def separate_exact_by_window(readData, parameter, normalize = "Large"):
    ''' split genome into windows'''
    window_size = parameter.window_size
    data_by_window_dict = {}
    index_dict = {}
    move_size = window_size/2

    for chr in readData.chr_list:
        info( "partitioning by window on \t%s...", chr)
        chr_len =readData.chr_length_dict[chr]
        row_num = chr_len/move_size - 1
        if row_num <= 0:
            continue
        if parameter.difftest is False: 
            for idx, filename in enumerate(readData.filename_list):
                array = count_reads_by_file(readData.data_dict[chr][filename],
                                            move_size, row_num)
                array = (
                        array*readData.normalization_constant[filename])
                if idx ==0:
                    chr_array = array
                else:
                    chr_array = numpy.column_stack((chr_array, array))   
        else: 
            for idx, filename in enumerate(readData.chip1_filename_list):
                array = count_reads_by_file(readData.data_dict[chr][filename],
                                            move_size, row_num)
                array = (
                        array*readData.normalization_constant[filename])
                if idx ==0:
                    chip1_array = array
                else:
                    chip1_array = numpy.column_stack((chip1_array, array))
            if parameter.chip1_matched_input is True: 
                for idx, filename in enumerate(readData.input1_filename_list):
                    array = count_reads_by_file(
                            readData.data_dict[chr][filename],
                            move_size, row_num)
                    array = (
                            array*readData.normalization_constant[filename])
                    if idx ==0:
                        input1_array = array
                    else:
                        input1_array = numpy.column_stack((input1_array, array))
            else:
                for idx, filename in enumerate(readData.input1_filename_list):
                    array = count_reads_by_file(
                            readData.data_dict[chr][filename],
                            move_size, row_num)
                    array = (
                            array*readData.normalization_constant[filename])
                    if idx ==0:
                        input1_array = array
                    else:
                        input1_array += array
                input1_array /= len(readData.input1_filename_list)
                input1_array.resize(row_num, 1)
            chip1_array -= input1_array
            chip1_array[chip1_array<0] = 0
            for idx, filename in enumerate(readData.chip2_filename_list):
                array = count_reads_by_file(readData.data_dict[chr][filename],
                                            move_size, row_num)
                array = (
                        array*readData.normalization_constant[filename])
                if idx ==0:
                    chip2_array = array
                else:
                    chip2_array = numpy.column_stack((chip2_array, array))
            if parameter.chip2_matched_input is True:
                for idx, filename in enumerate(readData.input2_filename_list):
                    array = count_reads_by_file(readData.data_dict[chr][filename],
                                                move_size, row_num)
                    array = (
                            array*readData.normalization_constant[filename])
                    if idx ==0:
                        input2_array = array
                    else:
                        input2_array = numpy.column_stack((input2_array, array))
            else:
                for idx, filename in enumerate(readData.input2_filename_list):
                    array = count_reads_by_file(readData.data_dict[chr][filename],
                                                move_size, row_num)
                    array = (
                            array*readData.normalization_constant[filename])
                    if idx ==0:
                        input2_array = array
                    else:
                        input2_array += array
                input2_array /= len(readData.input2_filename_list)
                input2_array.resize(row_num, 1)
            chip2_array -= input2_array
            chip2_array[chip2_array<0] = 0
            chr_array = numpy.column_stack((chip1_array, chip2_array))
        data_by_window_dict[chr] = chr_array
    readData.reads_dict = data_by_window_dict
#    for chr in readData.chr_list:
#        fileopen = open(chr+".reads.window", 'w')
#        for row in data_by_window_dict[chr]:
#            for reads in row: 
#                fileopen.write(str(reads)+'\t')
#            fileopen.write('\n')
#        fileopen.close()
    del readData.data_dict # this is not useful anymore
    return


def estimate_normalization_constant(readData, parameter): 
    '''Estimate the normalization constant for all samples'''
    # Split the genome into 1kb windows.
    bin = 1000
    array_dict = {}
    for file in readData.filename_list: 
        array = numpy.array([], dtype=numpy.float64)
        for chr in readData.chr_list:
            row_num = readData.chr_length_dict[chr]/bin
            array_by_chr = numpy.zeros(row_num, dtype=numpy.float64)
            for x in readData.data_dict[chr][file]:
                try: array_by_chr[x/bin] += 1
                except IndexError: pass 
            array = numpy.append(array, array_by_chr)
        array_dict[file] = array
    # Create a mixed chip sample and use it as the reference 
    mixed_chip_array = numpy.array([], dtype=numpy.float64)
    for idx, chip in enumerate(readData.chip_filename_list): 
        if idx == 0: 
            mixed_chip_array = array_dict[chip].copy()
            rep_rank_sum = rankdata(-array_dict[chip])
        else: 
            mixed_chip_array += array_dict[chip]
            rep_rank_sum += rankdata(-array_dict[chip])
    mixed_chip_array /= len(readData.chip_filename_list)
    # Estimate the input normalization constant using NCIS 
    for input in readData.input_filename_list: 
        norm_constant = estimate_input_normalization(
                mixed_chip_array, array_dict[input])
        readData.normalization_constant[input] = norm_constant
        debug("The scaling factor for %s is %s", input, norm_constant)
    # Estiamte the chip normalization cosntant using modified TMM method
    for chip in readData.chip_filename_list: 
        norm_constant = estimate_chip_normalization(
                mixed_chip_array, array_dict[chip], rep_rank_sum)
        readData.normalization_constant[chip] = norm_constant
        debug("The scaling factor for %s is %s", chip, norm_constant)
    return
def estimate_input_normalization(ref, target):
    '''Estimate the input normalization constant using the NCIS method'''
    combined = ref + target 
    READ_MAX = 200
    MIN_GENOME_COVERAGE = 0.75
    pre_ratio = 1
    file = open("input_norm.txt",'w')
    for r_cut in xrange(1, READ_MAX):
        index = numpy.where(combined <= r_cut)[0]
        percent_genome_covered = float(len(index))/len(combined)
        ratio = numpy.sum(ref[index])/numpy.sum(target[index])
        file.write(str(percent_genome_covered)+'\t'+str(ratio)+'\n')
        if percent_genome_covered > MIN_GENOME_COVERAGE and ratio > pre_ratio:
            return ratio
        else: 
            pre_ratio = ratio
    file.close()
    return ratio

def estimate_chip_normalization(ref, target, rep_rank_sum): 
    '''Estimate the ChIP normalization constant against the mixed ChIP
       using the modified TMM method'''
    N_PEAKS_GRID = numpy.array([1000, 5000, 10000, 20000, 30000, 40000, 50000])
    TRIM_M = 0.2
    TRIM_A = 0.05
    order = numpy.argsort(rep_rank_sum)
#    order = numpy.argsort(-ref)
    len_target_not_zero = len(numpy.where(target > 0)[0])
    ref = ref[order]
    target = target[order]
    tmm_array = numpy.array([])
    for n in N_PEAKS_GRID: 
        if n > len_target_not_zero: 
            break
        ref_n = ref[range(n)]
        ref_n[ref_n==0] = 1
        target_n = target[range(n)]
        target_n[target_n==0] = 1
        Mg = numpy.log2(ref_n/target_n)
        Ag = 0.5*numpy.log2(ref_n*target_n)
        Mg_sorted = Mg.copy()
        Mg_sorted.sort()
        Mg_lower_bound = Mg_sorted[n*TRIM_M]
        Mg_upper_bound = Mg_sorted[n*(1-TRIM_M)]
        Ag_sorted = Ag.copy()
        Ag_sorted.sort()
        Ag_lower_bound = Ag_sorted[n*TRIM_A]
        Ag_upper_bound = Ag_sorted[n*(1-TRIM_A)]
        trim_index = numpy.where((Mg > Mg_lower_bound) &
            (Mg < Mg_upper_bound) & (Ag > Ag_lower_bound) &
            (Ag < Ag_upper_bound))
        ref_trim = ref_n[trim_index]
        target_trim = target_n[trim_index]
        Mgk = numpy.log2(ref_trim/target_trim)
        Wgk = 0.5*numpy.log2(ref_trim*target_trim)
        tmm = 2**(numpy.sum(Mgk*Wgk)/numpy.sum(Wgk))
        tmm_array = numpy.append(tmm_array, tmm)
        debug("The TMM estiamted from top %s windows is %s", n, tmm)
    library_ratio = numpy.sum(ref)/numpy.sum(target)
    tmm_diff_array = numpy.abs(tmm_array - library_ratio)
    tmm_max = tmm_array[numpy.argmax(tmm_diff_array)]
    return tmm_max
