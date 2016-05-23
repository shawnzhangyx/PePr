#!/user/bin/env python
# use this script to find the mode of a peak. 
# Input: peak file, read file and format(bed,bam,sam)
# Output: peak file with additional column denoting the peak_mode
### usage: python find_peak_mode.py Peaks_file read_file file_format fragment_size
# Created By Yanxiao Zhang, 5/19/2016
############################################
import sys
import array
import numpy
import pysam 

def parse_sam_for_f_r(filename):
    infile = open(filename, 'r')
    num = 0
    forward = {}
    reverse = {}
    # skip the header of the SAM file.
    for line in infile:
        if not line.startswith("@"):
            break

    words = infile.readline().strip().split()
    length = len(words[9])
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
    return forward, reverse, readlength

def parse_bam_for_f_r(filename):
    num = 0
    forward = {}
    reverse = {}

    infile =pysam.Samfile(filename, 'rb')
    line = infile.fetch(until_eof=True).__next__()
    length = line.alen
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
    return forward,reverse, length

def parse_bed_for_f_r(filename):
    infile = open(filename, 'r')
    num = 0
    forward = {}
    reverse = {}
    words = infile.readline().strip().split()
    length = int(words[2]) - int(words[1])
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

    return forward, reverse, length

def cross_cor(f, r):
    npf = numpy.array(f)
    npr = numpy.array(r)
    cor_list = []
    for x in range(50,302,2):
        #print x
        y = len(numpy.intersect1d(npf+x,npr))
        cor_list.append(y)
    return range(50,302,2)[cor_list.index(max(cor_list))]


def estimate_shift_size(forward,reverse):
    shift_list = []
    print 'estimating shift size:'
    for chr in set(forward.keys())&set(reverse.keys()):
        print chr
        chr_f, chr_r = forward[chr], reverse[chr]
        shift_list.append(cross_cor(chr_f, chr_r))
    shift_size = shift_list[len(shift_list)/2]
    return shift_size/2
    
def shift_reads(forward, reverse,shift):
    print 'shifting reads'
    reads_by_chr = {}
    for chr in set(forward.keys())&set(reverse.keys()):
        chr_f = numpy.array(forward[chr])+shift
        chr_r = numpy.array(reverse[chr])-shift
        reads_by_chr[chr] = numpy.append(chr_f, chr_r)
    return reads_by_chr

def find_peak_mode(peak_filename, reads_by_chr,frag_size):
    out_filename = peak_filename + '.w_mode.txt'
    file_in = open(peak_filename, 'r')
    file_out = open(out_filename, 'w')
    for idx,line in enumerate(file_in):
        if idx > 0 and idx % 1000 == 0:
            print idx, 'peaks processed'
        words = line.strip().split('\t')
        chr = words[0]
        start = int(words[1])
        end = int(words[2])
        others = words[3:]
        window = numpy.zeros(end-start)
        read_chr = reads_by_chr[chr]
        read = read_chr[numpy.where((read_chr > start) & (read_chr < end) )]
        for item in read:
            if item-start-frag_size<0:
                window[0:(item-start+frag_size)]+=1
            else: 
                window[(item-start-frag_size):(item-start+frag_size)] +=1
        mode = numpy.max(window)
        where_mode = numpy.where(window == mode)[0][0]
        where_mode += start
        out_str = '\t'.join(words)+'\t'+str(where_mode)+'\n'
        file_out.write(out_str)

def main(argv):
    ''' usage: python find_peak_mode.py Peaks_file read_file file_format'''
    if len(argv) < 4:
        print '''usage: python find_peak_mode.py Peaks_file read_file file_format'''
        exit(1)

    peak_filename = argv[1]
    read_filename = argv[2]
    file_format = argv[3]
    process_file_for_f_r_dict = {'bam':parse_bam_for_f_r,'sam':parse_sam_for_f_r,'bed':parse_bed_for_f_r}
    forward, reverse,read_length = process_file_for_f_r_dict[file_format](read_filename)
    try: 
        frag_size = int(argv[4])
        shift_size = (frag_size -read_length)/2
    except IndexError:
        shift_size = estimate_shift_size(forward, reverse)
        frag_size = shift_size*2 + read_length
    reads_by_chr = shift_reads(forward,reverse,shift_size)
    find_peak_mode(peak_filename, reads_by_chr,frag_size/2)

if __name__ == "__main__":
    main(sys.argv)
