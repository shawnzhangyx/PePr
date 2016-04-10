import logging

import logConfig
import misc

info = logging.info
debug = logging.debug


class ReadData:
    ''' A data structure for the read data'''
    def __init__(self, chip1, input1, chip2, input2, diff_test):
        self.data_dict_by_strands = {}
        self.data_dict = {}
        self.reads_dict = {}
        self.chip1_filename_list = chip1
        self.input1_filename_list = input1
        if '' in self.input1_filename_list: 
            self.input1_filename_list.remove('')
        self.chip2_filename_list = chip2
        self.input2_filename_list = input2
        if '' in self.input2_filename_list:
            self.input2_filename_list.remove('')
        self.chip_filename_list = chip1[:]
        self.input_filename_list = input1[:]
        if diff_test is True:
            self.chip_filename_list += chip2
            self.input_filename_list += input2
        self.filename_list = self.chip_filename_list + self.input_filename_list
        #print(self.filename_list)
        
        # remove duplicated file names. 
        self.filename_list = list(set(self.filename_list))
        self.chr_list = []
        self.chr_length_dict = {}
        self.read_total_per_file = {}
        self.normalization_constant = {}
        self.genome_size = 0
        self.shift_size = {}
        self.read_length = 0

    def cal_genome_size(self):  # calculate the genome size.
        for chr in self.chr_list:
            self.genome_size += self.chr_length_dict[chr]
        debug("The genome size is %d", self.genome_size)

    def cal_read_total(self):
        # calculate the total number of reads for each sample.
        for chr in self.chr_list:
            for filename in self.filename_list:
                for strand in ["f", "r"]:
                    try:
                        self.read_total_per_file[filename] += \
                        len(self.data_dict_by_strands[chr][filename][strand])
                    except KeyError:
                        self.read_total_per_file[filename] = \
                        len(self.data_dict_by_strands[chr][filename][strand])

    def remove_redundant_reads(self):
        ''' remove additional redundant reads that are not warranted
        by a binomial test.'''
        def rm_max(list, max):
            # remove reads that are more than the maximum specified.
            list_new = []
            pre = 0
            count = 0

            for x in list:
                if x == pre:
                    count += 1
                else:
                    pre = x
                    count = 1
                if count > max:
                    continue
                list_new.append(x)
            return list_new

        for filename in self.filename_list:
            max = misc.binomial(
                self.read_total_per_file[filename], 1.0/self.genome_size)
            # calculate maximum duplicates for each genomic location.
            info("The maximum at one position for %s is %d", filename, max)
            total_before = self.read_total_per_file[filename]
            total_after = 0
            for chr in self.chr_list:
                for strand in ["f", "r"]:
                    reads = self.data_dict_by_strands[chr][filename][strand]
                    reads.sort()
                    reads = rm_max(reads, max)
                    total_after += len(reads)
                    self.data_dict_by_strands[chr][filename][strand] = reads
            self.read_total_per_file[filename] = total_after
            debug("Total # of reads before removing redundancies is %d",
                  total_before)
            debug("Total # of reads after removing redundancies is %d",
                  total_after)
            info("The percentage of redundant reads for %s is %f ",
                 filename, 1-float(total_after)/total_before)

    def refine_chr_list(self):
        '''Delete the chromosomes that having 0 reads. Find the maximum
           genomic location in each chromosome'''
        chr_list = []
        chr_length_dict = {}
        for chr in self.data_dict_by_strands:
            temp_max_chr = 0
            all_file_has_this_chr = True
            for file in self.filename_list:
                if file not in self.data_dict_by_strands[chr]:
                    all_file_has_this_chr = False
                    break
                max_read_f = max(self.data_dict_by_strands[chr][file]['f']+[0])
                max_read_r = max(self.data_dict_by_strands[chr][file]['r']+[0])
                temp_max_chr = max(temp_max_chr, max_read_f, max_read_r)
            if all_file_has_this_chr is True: 
                chr_list.append(chr)
                chr_length_dict[chr] = temp_max_chr
        chr_list.sort()
        self.chr_list = chr_list
        self.chr_length_dict = chr_length_dict


class Parameters:
    "store the parameters that will be used in many processes."
    # ---- input parameters ---- #
    # ---- function specific parameters ---- #

    def __init__(self, opt):
        self.chip1 = []
        self.chip2 = []
        self.input1 = []
        self.input2 = []
        self.file_format = None
        self.shift_size = -1
        self.window_size = -1
        self.difftest = None
        self.name = None
        self.remove_redundant = None
        self.threshold = None
        self.peaktype = None
        self.normalization = None
        self.shift_dict = {}
        self.normalization_dict = {}
        # deprecated parameters. 
        self.remove_artefacts = False
        self.narrow_peak_width = False 
        
        if opt.parameter is not "":
            self.process_parameter_file(opt.parameter)
        else:
            self.process_command_line_option(opt)
        self.validate_parameters()
        logConfig.startLog(self.name)
        self.print_parameters()
        # --- initialize logging --- #

    def process_parameter_file(self, parameter_file):
        # using a case/switch commands
        file = open(parameter_file, 'r')
        for line in file:
            items = line.split('#')[0].split()
            # skip empty lines and lines starting with #
            if len(items) == 0:
                continue
            if len(items) > 4:
                raise Exception("There are more than 4 columns in your parameter file.")
            key = items[0].strip().lower()
            value = items[1:]
            if key == "chip1":
                self.chip1.append(value[0])
            if key == "chip2":
                self.chip2.append(value[0])
            if key == "input1":
                self.input1.append(value[0])
            if key == "input2":
                self.input2.append(value[0])
            # read the shift size and normalization constants. 
            if key in ['chip1','chip2','input1','input2']:
                if len(value) > 1: 
                    self.shift_dict[value[0]] = int(value[1])
                if len(value) > 2:
                    self.normalization_dict[value[0]] = int(value[2])
                    
            if key == "file_format":
                self.file_format = value[0].lower()
            if key == "shiftsize":
                self.shift_size = int(value[0])
            if key == "windowsize":
                self.window_size = int(value[0])
            if key == "peaktype":
                self.peaktype = value[0].lower()
            if key == "difftest":
                self.difftest = value[0].lower()
            if key == "name":
                self.name = value[0]
            if key == "remove_redup":
                self.redundant = value[0]
            if key == "threshold":
                self.threshold = float(value[0])
            if key == "normalization":
                self.normalization = value[0].lower()
                
    def process_command_line_option(self, opt):
        self.chip1 = opt.chip1
        self.input1 = opt.input1
        self.chip2 = opt.chip2
        self.input2 = opt.input2
        self.file_format = opt.file_format.lower()
        self.shift_size = opt.shift_size
        self.window_size = opt.window_size
        self.difftest = opt.difftest
        self.name = opt.name
        self.remove_redundant = opt.remove_redundant
        self.threshold = opt.threshold
        self.peaktype = opt.peaktype.lower()
        self.normalization = opt.normalization

    def validate_parameters(self):
        # if there are files with the same name, raise an exception. 
        
        # if any of the required files are missing, raise an exception.
        if self.chip1 == ['']:
            raise Exception("Please specify ChIP-1 sample names")
        if self.input1 == [''] and self.difftest is False:
            raise Exception("Please specify input-1 sample names")
        if self.difftest is True:  # also check group 2 files if difftest
            if self.chip2 == ['']:
                raise Exception("Please specify ChIP-2 sample names")
        # will not test if input2 is available. 
        #    if self.input2 == ['']:
        #        raise Exception("Please specify input-2 sample names")
        if len(self.chip1) < 2 and len(self.input1) <2: 
            raise Exception('''Only 1 replicates detected. To use PePr, at least
                two replicates are required''')
        if self.chip2 != [''] and self.input2 != [''] and  self.difftest is False: 
            raise Exception('''chip2 and input2 detected, but --diff parameter
                is missing. Please add --diff if you want to run differential
                analysis''')
        if not self.file_format:
            raise Exception('''Please specify a file format:  bed,
                eland_multi, eland_extended, bowtie, sam, or bam''')
        if self.file_format not in [
                'bed', 'eland_multi', 'eland_extended',
                'bowtie', 'sam', 'bam'
                ]:
            raise Exception('''Please specify a valid file format: bed,
            eland_multi, eland_extended, bowtie, sam, or bam''')
        if self.peaktype not in ['sharp', 'broad']:
            raise Exception('''please specify a peak type: sharp or broad.
            Typically, sharp works for TF better and broad
            for histone modifications.''')
        
        if self.difftest is True:
            if len(self.chip1) == \
                    len(self.input1):
                self.chip1_matched_input = True
            else:
                self.chip1_matched_input = False
            if len(self.chip2) == \
                    len(self.input2):
                self.chip2_matched_input = True
            else:
                self.chip2_matched_input = False

    def print_parameters(self):
        info('printing parameters')
        info('the chip1 files are...')
        print self.chip1
        
    def write_parameter_to_file(self):
        fileout = open(self.name+"__PePr_parameters.txt", 'w')
        attributes = vars(self)
        keys_list = list(attributes.keys())
        keys_list.sort()
        for keys in keys_list:
            fileout.write(keys+'\t'+str(attributes[keys])+'\n')
        fileout.close()


class Peak:
    "data structure that contains the significant peaks"

    def __init__(self, chr, index, pvalue, qvalue):
        self.chr = chr
        self.index = index
        self.pvalue = pvalue
        self.qvalue = qvalue

    def __repr__(self):
        return repr((self.chr, self.index, self.pvalue, self.qvalue))
