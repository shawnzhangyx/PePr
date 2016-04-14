#!/usr/bin/env python26

import sys
import re
from logging import info
import multiprocessing
import numpy
import pysam
import itertools



class FileParser:
    ''' have the functions to prcess different types of files'''
    ## BIN_size used for normalization. 
    BIN_SIZE = 1000
    
    def __init__(self, name, file_format):
        self.name = name
        self.file_format = file_format
        #self.window_size = 90
        #self.shift_size 
        #self.window_size 
        #self.normalization_constant = parameter.normalization_dict[name]
        
    
    def get_chr_info(self):
        ''' get the chromosome information from file '''
        
        info ("getting chromosome info")
        if self.file_format == "bam":
            chr_info = self.get_chr_info_bam()   
        if self.file_format == "sam":
            chr_info = self.get_chr_info_sam()        
        if self.file_format == "bed":
            chr_info = self.get_chr_info_bed()
        # change the class chromosome info. 
        FileParser.chr_info = chr_info 
        return chr_info
        
        
    def get_chr_info_bed(self):
        chr_info = {}
        with open(self.name , 'r') as infile:
            #num = 0
            for line in infile: 
                #num += 1
                #if num % 1000000 == 0:
                #    print num
                chr,start,end,col3,col4,strand = line.strip().split()
                try: 
                    chr_info[chr] = max(chr_info[chr],int(end))
                except KeyError:
                    chr_info[chr] = int(end)
        return chr_info    

    def get_chr_info_bam(self):
        chr_info = {}
        with pysam.Samfile(self.name, 'rb') as infile:
            ### you can get info from the header of a bam file actually. 
            ### will check for it. Some bam files do not have header. 
            for line in infile.fetch():
                if line.is_unmapped is False:
                    chr = infile.getrname(line.tid)
                    try: 
                        chr_info[chr] = max(chr_info[chr], line.pos)
                    except KeyError: 
                        chr_info[chr] = line.pos
        return chr_info
        
    def get_chr_info_sam(self):
        chr_info = {}
        with open(self.name, 'r') as infile:
        # skip the header of the SAM file. 
            for line in infile:
                # can actually get the info from the sam file header. 
                if not line.startswith("@"):
                    break
            # start reading the real data
            for line in infile:
                words = line.strip().split()
                flag = int(words[1])
            
                if flag & 0x0004: #if not unmapped
                    chr, pos =  words[2], int(words[3])-1
                    try: 
                        chr_info[chr] = max(chr_info[chr],pos)
                    except KeyError: 
                        chr_info[chr] = pos
        return chr_info
       
    def parse_data_bin(self, windowsize, shift_size, num_procs):
        info ("start reading %s", self.name)
        slide_size = windowsize/2
        self.data_dict = {}
        for chr in FileParser.chr_info:
            row_num = int(FileParser.chr_info[chr]/slide_size) - 1
            self.data_dict[chr] = numpy.zeros(row_num, dtype=numpy.float64)
        # 
        if self.file_format == "bed":
            self.parse_bed_bin(slide_size, shift_size)
        elif self.file_format == "bam":
            self.parse_bam_bin(slide_size, shift_size)
        elif self.file_format == "sam":
            self.parse_sam_bin(slide_size, shift_size)
            
        for chr in parameter.chr_info: 
            data_dict[chr] = data_dict[chr] + numpy.roll(data_dict[chr],-1)
            data_dict[chr] = data_dict[chr] * parameter.normalization_dict[filename]
        # info ("finished reading %s", self.name)
        #self.data_dict = data
        
    def parse_bed_bin(self, slide_size, shift_size):
        
        with open(self.name, 'r') as infile:
            for line in infile: 
                chr,start,end,col3,col4,strand = line.strip().split()
                if strand == "+":
                    pos = int(start) + shift_size
                else: 
                    pos = int(end) - shift_size
                try: self.data_dict[chr][int(pos/slide_size)] += 1
                except IndexError:    pass    

        return 
            

    