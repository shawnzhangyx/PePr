import pysam
import array
import numpy

def parse_bam_for_f_r(filename, input_dir):
    num = 0
    forward = {}
    reverse = {}

    infile =pysam.Samfile(input_dir+filename, 'rb')
    for line in infile.fetch(until_eof = True):
        num += 1
        if num % 10000000 == 0 :
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

def parse_sam_for_f_r(filename, input_dir):
    infile = open(input_dir+filename, 'r')
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
                try:  forward[chr].append(pos)
                except KeyError:
                    forward[chr] = array.array('i',[pos])
            else:
                try:  reverse[chr].append(pos)
                except KeyError:
                    reverse[chr] = array.array('i',[pos])
    return forward, reverse

def parse_bed_for_f_r(filename, input_dir):
    infile = open(input_dir+filename, 'r')
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

def parse_sampe(filename,input_dir):
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

    # will calculate the median fragment size and remove the reads that have larger fragment size than it. 
    flen_median = numpy.median(flen_list)
    # print flen_median

    for chr in reads_dict:
        flen_chr = numpy.array(flen_dict[chr])
        reads_chr = numpy.array(reads_dict[chr])
        reads_dict[chr] = reads_chr[numpy.where(numpy.abs(flen_chr) <= 2*flen_median)[0]]
        # print chr, len(reads_chr),len(reads_dict[chr])

    return reads_dict 

def parse_bampe(filename, input_dir):
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
        if num % 1000000 == 0 :
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
    # calculate the median fragment size and remove reads that are two times larger.
    flen_median = numpy.median(flen_list)
    for chr in reads_dict:
        flen_chr = numpy.array(flen_dict[chr])
        reads_chr = numpy.array(reads_dict[chr])
        reads_dict[chr] = reads_chr[numpy.where(numpy.abs(flen_chr) <= 2*flen_median)[0]]
    return reads_dict 

### function dictionary. 
parse_file_by_strand = {'bed':parse_bed_for_f_r,'bam':parse_bam_for_f_r,'sam':parse_sam_for_f_r}
parse_file_pe = {'sampe':parse_sampe, 'bampe':parse_bampe}
