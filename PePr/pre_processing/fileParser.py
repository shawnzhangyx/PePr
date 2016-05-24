import pysam
import array


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

### function dictionary. 
parse_file_by_strand = {'bed':parse_bed_for_f_r,'bam':parse_bam_for_f_r,'sam':parse_sam_for_f_r}
