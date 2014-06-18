#!/usr/bin/env python26

import sys
import re
import logging

import numpy
import pysam


root_logger = logging.getLogger("")
info = root_logger.info


def bam_parse(filename_list):
    '''Parsing bam format files.'''
    data_dict = {}
    infile = pysam.Samfile(filename_list[0], 'rb')
    line = infile.fetch().next()
    length = line.alen
    infile.close()
    info("the lenghth of the reads is %s", length)

    for filename in filename_list: 
        infile = pysam.Samfile(filename, 'rb')
        info("retrieving reads from file : %s", filename)
        for line in infile.fetch():
            if line.is_unmapped is False: 
                chr = infile.getrname(line.tid)
                pos = line.pos
                if line.is_reverse is False: 
                    try:
                        data_dict[chr][filename]['f'].append(pos)
                    except KeyError:
                        if chr not in data_dict:
                            data_dict[chr] = {}
                        data_dict[chr][filename] = {}
                        data_dict[chr][filename]['f'] =[pos]
                        data_dict[chr][filename]['r'] =[]
                else:
                    try:
                        data_dict[chr][filename]['r'].append(pos + length +1)
                    except KeyError:
                        if chr not in data_dict:
                            data_dict[chr] = {}
                        data_dict[chr][filename] = {}
                        data_dict[chr][filename]['f'] =[]
                        data_dict[chr][filename]['r'] =[pos + length + 1]
        infile.close()
    return data_dict, length



def sam_parse(filename_list):
    '''Parsing sam format files.'''
    # Store the position and strand data for each file,
    # seperated by chromosomes.
    data_dict = {}
    infile = open(filename_list[0], 'r')

    for line in infile:
        words = line.strip().split()
        if not words[0].startswith("@"):
            length = len(words[9])
            break
    infile.close()
    info("the lenghth of the reads is %s", length)

    for filename in filename_list:
        infile = open(filename, 'r')
        info("retrieving reads from file : %s", filename)
        for line in infile: 
            if not line.startswith("@"):
                break
        for line in infile:
            words = line.strip().split()
            words[1] = int(words[1])

            if words[1] is not 0x0004:
                chr = words[2]
                pos = int(words[3])-1

                if words[1] is not 16:
                    try: 
                        data_dict[chr][filename]['f'].append(pos)
                    except KeyError:
                        if chr not in data_dict:
                            data_dict[chr] = {}
                        data_dict[chr][filename] = {}
                        data_dict[chr][filename]['f'] =[pos]
                        data_dict[chr][filename]['r'] =[]
                else:
                    try:
                        data_dict[chr][filename]['r'].append(pos + length +1)
                    except KeyError:
                        if chr not in data_dict:
                            data_dict[chr] = {}
                        data_dict[chr][filename] = {}
                        data_dict[chr][filename]['f'] =[]
                        data_dict[chr][filename]['r'] =[pos + length + 1]

                    except KeyError: pass
        infile.close()
    return data_dict, length

def bowtie_parse(filename_list, chr_list):
    #parsing bowtie format files
    data_dict = {} 
    infile = open(filename_list[0], 'r')
    line = infile.readline()
    words = line.strip().split()
    length = len(words[4])
    info("the lenghth of the reads is %s", length)
    infile.close()

    for filename in filename_list:
        infile = open(filename, 'r')
        info("retrieving reads from file : %s", filename)

        for line in infile:
            line = line.strip().split()
            chr = line[2]
            strand = line[1]
            pos = int(line[3])

            if strand == '+':
                try: data_dict[chr][filename]['f'].append(pos)
                except KeyError:
                    if chr not in data_dict:
                        data_dict[chr] = {}
                    data_dict[chr][filename] = {}
                    data_dict[chr][filename]['f'] =[pos]
                    data_dict[chr][filename]['r'] =[]
            elif strand == '-':
                try: data_dict[chr][filename]['r'].append(pos + length)
                except KeyError:
                    if chr not in data_dict:
                        data_dict[chr] = {}
                    data_dict[chr][filename] = {}
                    data_dict[chr][filename]['f'] =[]
                    data_dict[chr][filename]['r'] =[pos + length + 1]
            else:
                print("strand error")
        infile.close()
    return data_dict, length

def bed_parse(filename_list):
    '''Parsing BED format files'''
    data_dict={} 
    for filename in filename_list:
        infile = open(filename, 'r')
        info("retrieving reads from file: %s", filename)
        for line in infile:
            line = line.strip().split()
            chr=line[0]
            strand = line[5]
            if strand == '+':
                pos = int(line[1])  # genomic position
                try: data_dict[chr][filename]['f'].append(pos)
                except KeyError:
                    if chr not in data_dict:
                        data_dict[chr] = {}
                    data_dict[chr][filename] = {}
                    data_dict[chr][filename]['f'] =[pos]
                    data_dict[chr][filename]['r'] =[]
            elif strand == '-':
                pos = int(line[2])-1
                try: data_dict[chr][filename]['r'].append(pos)
                except KeyError:
                    if chr not in data_dict:
                        data_dict[chr] = {}
                    data_dict[chr][filename] = {}
                    data_dict[chr][filename]['f'] =[]
                    data_dict[chr][filename]['r'] =[pos]
            else:
                print("strand error")
        infile.close()
    return data_dict, int(line[2])-int(line[1])  # read_length

def process_illumina_match(align, mismatch, length, format):
    words = align.split(',')
    if format == "multi":
        for match in words:
            if match.startswith('c'):
                str, match =match.split(':')
                chr = re.match(r'(chr(\d\d?|\w))', str)
                chr = chr.group(0)
            if match.endswith(mismatch):
                strand = match[-2]
                pos = int(match[0:-2])
                if strand =='R':
                    pos = pos+length
                return chr, strand, pos
    elif format == "extended":
        for match in words:
            if match.startswith('c'):
                chr_part, match = match.split(':')
                chr = re.match(r'(chr(\d\d?|\w))', chr_part)
                chr = chr.group(0)
                m1 = re.search(r'([F|R])', match)
                strand = m1.group(1)
                pos = int(re.sub(r'[F|R].*', '', match))

                if strand == 'R':
                    pos = pos + length

                match = re.sub(r'(\^\w+\$)', '', match)
                m2 = re.findall(r'[A-Z]', match)

                if len(m2) - 1 == int(mismatch):
                    return chr, strand, pos


# Parsing eland format files
def eland_parse(filename_list, chr_list, format="default"):
    #store the position and strand data for each file, seperated by chromosomes.
    data_dict={}
    infile = open(filename_list[0], 'r')
    line = infile.readline()
    words=line.split()
    length=len(words[1])
    info ("the length of the reads is %s", length)
    infile.close()

    for filename in filename_list:
        infile = open(filename, 'r')
        info("retrieving reads from file:%s", filename)
        for line in infile:
            words = line.strip().split()
            mismatch = -1
            
            if len(words) == 4 and ':' in words[2]:
                num_mismatch = words[2].split(':')
                
                if num_mismatch[0]== '1':
                    mismatch = '0'
                elif (num_mismatch[0]=='0') and (num_mismatch[1] == '1'):
                    mismatch = '1'
                elif (num_mismatch[0]=='0' and num_mismatch[1]=='0' and
                        num_mismatch[2] == '1'):
                    mismatch = '2'
                if mismatch !=-1:
                    try:
                        chr, strand, pos = process_illumina_match(
                                words[3], mismatch, length, format)
                        if strand == 'F':
                            try: data_dict[chr][filename]['f'].append(pos)
                            except KeyError:
                                if chr not in data_dict:
                                    data_dict[chr] = {}
                                data_dict[chr][filename] = {}
                                data_dict[chr][filename]['f'] =[pos]
                                data_dict[chr][filename]['r'] =[]
                        elif strand == 'R':
                            try: data_dict[chr][filename]['r'].append(pos)
                            except KeyError:
                                if chr not in data_dict:
                                    data_dict[chr] = {}
                                data_dict[chr][filename] = {}
                                data_dict[chr][filename]['f'] =[]
                                data_dict[chr][filename]['r'] =[pos + length + 1]
                    except TypeError:
                        pass
        infile.close()
    return data_dict, length


def parse(readData, file_format):
    '''retrieve the data form files and make it into appropriate data 
        format'''
    info("begin parsing files...")
    data_dict = {}
    if file_format == "bed":
        data_dict, read_length = bed_parse(readData.filename_list)
    elif file_format == "eland_multi":
        data_dict, read_length = eland_parse(
                readData.filename_list, readData.chr_list, format="multi")
    elif file_format == "eland_extended":
        data_dict, read_length = eland_parse(
                readData.filename_list, readData.chr_list, format="extended")
    elif file_format == "bowtie":
        data_dict, read_length = bowtie_parse(
                readData.filename_list, readData.chr_list)
    elif file_format == "sam":
        data_dict, read_length = sam_parse(readData.filename_list)
    elif file_format == "bam":
        data_dict, read_length = bam_parse(readData.filename_list)
    else: 
        pass

    #--- here can add other file format support ---#

    # fill in the data structures. 
    readData.data_dict_by_strands = data_dict
    readData.read_length = read_length
    readData.refine_chr_list()
    readData.cal_read_total()
    for file in readData.filename_list: 
        info("Total #reads for %s is %s", file, 
                readData.read_total_per_file[file])
    readData.cal_genome_size()
    info( "leaving file_parser.py")

