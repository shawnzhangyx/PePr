import os
from logging import info, debug, warning
import collections
import pysam

from . import logConfig
## global variable read_dict
read_dict = {}
def init_dict():
    global read_dict
    read_dict = {}

       
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
        self.difftest = False 
        self.name = None
        self.threshold = None
        self.peaktype = None
        self.normalization = '' 
        self.input_directory = ''
        self.output_directory = ''
        self.keep_max_dup = -1
 
        ## dictionaries storing file related structures. 
        self.read_length_dict = {} # store the read length
        self.shift_dict = {} #store the shift sizes
        self.bin_dict = {} #store the binned reads for normalization.
        self.lib_size_dict = {} # store the library size for normalization.
        self.normalization_dict = {} # store the normalization constants.
        self.array_dict = {} # store the processed read array for significance testing. 
        
        self.num_procs = 1
        # reads specific parameters. 
        self.read_length = -1
        # deprecated parameters. 
        self.remove_artefacts = False
        self.narrow_peak_width = False 
        
        if opt.parameter is not "":
            self.process_parameter_file(opt.parameter)
        else:
            self.process_command_line_option(opt)
        self.validate_parameters()
        self.validate_files()
        logConfig.startLog(self.output_directory + self.name)
        self.check_file_formats()
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
            elif key == "chip2":
                self.chip2.append(value[0])
            elif key == "input1":
                self.input1.append(value[0])
            elif key == "input2":
                self.input2.append(value[0])
            # read the shift size and normalization constants. 
            if key in ['chip1','chip2','input1','input2']:
                if len(value) > 1: 
                    self.shift_dict[value[0]] = int(value[1])
                if len(value) > 2:
                    self.normalization_dict[value[0]] = float(value[2])
                    
            elif key == "file-format":
                self.file_format = value[0].lower()
            elif key == "shiftsize":
                self.shift_size = int(value[0])
            elif key == "windowsize":
                self.window_size = int(value[0])
            elif key == "peaktype":
                self.peaktype = value[0].lower()
            elif key == "difftest":
                self.difftest = value[0].lower() == "true"
            elif key == "name":
                self.name = value[0]
            elif key == "threshold":
                self.threshold = float(value[0])
            elif key == "keep-max-dup":
                self.keep_max_dup = int(value[0])
            elif key == "normalization":
                self.normalization = value[0].lower()
            elif key == "num-processors":
                self.num_procs = int(value[0])
            elif key == "input-directory":
                self.input_directory = value[0]
            elif key == "output-directory":
                self.output_directory = value[0]                
            else:
                raise Exception("Incorrect or unknown parameter: {0}".format(line))
        return 
        
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
        self.threshold = opt.threshold
        self.peaktype = opt.peaktype.lower()
        self.normalization = opt.normalization
        self.keep_max_dup = opt.keep_max_dup
        self.input_directory = opt.input_directory
        self.output_directory = opt.output_directory
        self.num_procs = opt.num_procs
        return 
        
    def validate_parameters(self):
        # if there are files with the same name, raise an exception. 
        
        # if any of the required files are missing, raise an exception.
        if self.chip1 == []:
            raise Exception("Please specify ChIP-1 sample names")
        if self.input1 == [] and self.difftest is False:
            raise Exception("Please specify input-1 sample names")
        if self.difftest is True:  # also check group 2 files if difftest
            if self.chip2 == []:
                raise Exception("Please specify ChIP-2 sample names")
        # will not test if input2 is available. 
        #    if self.input2 == ['']:
        #        raise Exception("Please specify input-2 sample names")
        if len(self.chip1) < 2 and len(self.input1) <2: 
            raise Exception('''Only 1 replicates detected. To use PePr, at least
                two replicates are required''')
        if self.chip2 != [] and self.input2 != [] and  self.difftest is False: 
            raise Exception('''chip2 and input2 detected, but --diff parameter
                is missing. Please set difftest to True if you want to run differential
                analysis''')
        if not self.file_format:
            raise Exception('''Please specify a file format:  bed,
                sam, bam, sampe or bampe''')
        if self.file_format not in [
                'bed', 'sam', 'bam', 'sampe', 'bampe'
                ]:
            raise Exception('''Please specify a valid file format: bed,
            sam, bam, sampe or bampe''')
        if self.peaktype not in ['sharp', 'broad']:
            raise Exception('''please specify a peak type: sharp or broad.
            Typically, sharp works for TF better and broad
            for histone modifications.''')
       
        ## check normalization
        # if normalization constants not provided
        if  len(self.normalization_dict) == 0:
            if  self.normalization =='':
                if self.difftest is False:
                    self.normalization = "intra-group"
                else: 
                    self.normalization = "inter-group"
            if self.normalization not in ['inter-group','intra-group','scale','no']:
                raise Exception('''Please specify a normalization method: inter-group, intra-group, or scale. put 'no' if you don't want to normalize''')
            
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
        
        # add '/' to the end of directories if not. 
        # no '/' at beginning of filenames if input directory is present.
        if self.input_directory != '':
            for filename in self.get_filenames():
                if filename.startswith('/'):
                    raise Exception('If input directory is specified, absolute path("/") in file name is invalid')
            if not self.input_directory.endswith('/'):
                self.input_directory += '/'
        if self.output_directory != '' and \
                not self.output_directory.endswith('/'):
            self.output_directory += '/'

        # can change relative directory to absolute directory. 

    def validate_files(self):
        if self.input_directory != '':
            if os.path.isdir(self.input_directory) is False:
                raise Exception('Input directory: {0} does not exist'.format(self.input_directory))
        if self.output_directory != '':
            if os.path.isdir(self.output_directory) is False:
                raise Exception('Output directory: {0} does not exist'.format(self.output_directory))
        for filename in self.get_filenames():
            if os.path.isfile(self.input_directory+filename) is False:
                print("File:",self.input_directory+filename, " not found")
                exit(1)

    @staticmethod
    def check_sampe_sorted(filename, input_dir):
        with open(input_dir + filename,'r') as infile:
            for line in infile:
                if not line.startswith("@"):
                    break
            count_list = []
            count = 1
            pre_name = ''
            for idx,line in enumerate(infile): 
                name = line.strip().split()[0]
                if name == pre_name:
                    count += 1
                else: 
                    count_list.append(count)
                    count = 1
                    pre_name = name
                if idx == 999:
                    break

        count1 = len([i for i in count_list if i==1])
        ratio = float(count1)/sum(count_list)
        #print filename, ratio
        if ratio > 0.8:
            warning("%s may not be sorted by read name. Please check.",filename)

    @staticmethod
    def check_bampe_sorted(filename, input_dir):
        infile = pysam.Samfile(input_dir + filename, 'rb')
        count_list = []
        count = 1
        pre_name = ''
        for idx,line in enumerate(infile.fetch(until_eof = True)):
            if line.query_name == pre_name:
                count += 1
            else:
                count_list.append(count)
                count = 1
                pre_name = line.query_name
            if idx == 999:
                break

        count1 = len([i for i in count_list if i==1])
        ratio = float(count1)/sum(count_list)
        #print filename, ratio
        if ratio > 0.8:
            warning("%s may not be sorted by read name. Please check.",filename)

    def check_file_formats(self):
        check_name_sorted = {'sampe':self.check_sampe_sorted, 'bampe':self.check_bampe_sorted}
        if self.file_format in ['sampe', 'bampe']:
            for filename in self.get_filenames():
                check_name_sorted[self.file_format](filename, self.input_directory)
        

    def write_parameter_to_file(self):
        '''write the current parameters to the files so user can repeat the analysis'''
        # check if the file name has already be taken. 
        # never mind, will only implement it if really needed.
        output_str = '\t'.join(['#filetype', 'filename', 'shift_size','normalization_factor'])+'\n' 
        
        for filename,filetype in zip(self.get_filenames(), self.get_filetypes()):
            output_str += '\t'.join([filetype,filename,str(self.shift_dict[filename]), str(self.normalization_dict[filename])]) +'\n'
                
        output_str += 'file-format\t'+self.file_format+'\n'
        output_str += 'peaktype\t'+self.peaktype+'\n'
        output_str += 'windowsize\t'+str(self.window_size)+'\n'
        output_str += 'difftest\t'+str(self.difftest)+'\n'
        output_str += 'threshold\t'+str(self.threshold)+'\n'
        output_str += 'normalization\t'+self.normalization+'\n'
        if self.keep_max_dup > 0:
            output_str +='keep-max-dup\t'+str(self.keep_max_dup)+'\n'
        if self.input_directory != '':
            output_str += 'input-directory\t'+self.input_directory+'\n'
        if self.output_directory !='':
            output_str += 'output-directory\t'+self.output_directory+'\n'
        output_str += 'name\t'+self.name+'\n'
        output_str += 'num-processors\t'+str(self.num_procs)+'\n'
        #print and save the parameters into a file. 
        self.print_parameters(output_str)
        with open(self.output_directory + self.name+"__PePr_parameters.txt", 'w') as outfile:
            outfile.write(output_str)
        return 
            
    def print_parameters(self, output_str):
        info("printing running parameters:")
        print('\n\n'+ output_str + '\n')    
        return 
        
    def get_filenames(self):
        return self.chip1+self.chip2+self.input1+self.input2
    
    def get_uniq_filenames(self):
       return list(set(self.get_filenames()))
    
    def get_filetypes(self):
        return ['chip1']*len(self.chip1) + ['chip2'] * len(self.chip2) + ['input1'] * len(self.input1) + ['input2'] * len(self.input2)
    
    def get_filenames_wo_bin_dict(self):
        return [filename for filename in self.get_filenames() if filename not in self.bin_dict ]
    
    def get_chip_filenames(self):
        return self.chip1+self.chip2

    def get_input_filenames(self):
        return self.input1+self.input2
        
    def get_genome_size(self):
        return sum(self.chr_info.values())
    
    def get_top3_chr(self):
        values = [i for i in self.chr_info.values()]
        values.sort()
        values.reverse()
        values_top3 = values[0:3]
        chr_top3 = []
        for value in values_top3:
            for chr in self.chr_info:
                if self.chr_info[chr] == value:
                    chr_top3.append(chr)
        return chr_top3
        
class Peak:
    "data structure that contains the significant peaks"

    def __init__(self, chr, index, group1_count, group2_count, pvalue, qvalue):
        self.chr = chr
        self.index = index
        self.g1_count = group1_count
        self.g2_count = group2_count
        self.pvalue = pvalue
        self.qvalue = qvalue

    def __repr__(self):
        return repr((self.chr, self.index, self.pvalue, self.qvalue))
