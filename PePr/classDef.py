import os
import logConfig
from logging import info, debug

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
        self.difftest = None
        self.name = None
        self.threshold = None
        self.peaktype = None
        self.normalization = None
        self.input_directory = './'
        self.output_directory = './'
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
                    self.normalization_dict[value[0]] = float(value[2])
                    
            if key == "file_format":
                self.file_format = value[0].lower()
            if key == "shiftsize":
                self.shift_size = int(value[0])
            if key == "windowsize":
                self.window_size = int(value[0])
            if key == "peaktype":
                self.peaktype = value[0].lower()
            if key == "difftest":
                self.difftest = value[0].lower() == "true"
            if key == "name":
                self.name = value[0]
            if key == "remove_redup":
                self.redundant = value[0]
            if key == "threshold":
                self.threshold = float(value[0])
            if key == "keep_max_dup":
                self.keep_max_dup = int(value[0])
            if key == "normalization":
                self.normalization = value[0].lower()
            if key == "num_processors":
                self.num_procs = int(value[0])
            if key == "input_directory":
                self.input_directory = value[0]
            if key == "output_directory":
                self.output_directory = value[0]                
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
                sam, or bam''')
        if self.file_format not in [
                'bed', 'sam', 'bam'
                ]:
            raise Exception('''Please specify a valid file format: bed,
            sam, or bam''')
        if self.peaktype not in ['sharp', 'broad']:
            raise Exception('''please specify a peak type: sharp or broad.
            Typically, sharp works for TF better and broad
            for histone modifications.''')
        
        if  len(self.normalization_dict) == 0 and self.normalization not in ['scale','compound','none']:
            raise Exception('''Please specify a normalization method: scale, compound, none. 
            Or give normalization constants.''')
            
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
        if not self.input_directory.endswith('/'):
            self.input_directory += '/'
        if not self.output_directory.endswith('/'):
            self.output_directory += '/'
        # can change relative directory to absolute directory. 

    def validate_files(self):
        for filename in self.get_filenames():
            if os.path.isfile(self.input_directory+filename) is False:
                print "File:",self.input_directory+filename, " not found"
                exit(1)
    def write_parameter_to_file(self):
        '''write the current parameters to the files so user can repeat the analysis'''
        # check if the file name has already be taken. 
        # never mind, will only implement it if really needed.
        output_str = '\t'.join(['#filetype', 'filename', 'shift_size','normalization_factor'])+'\n' 
        
        for filename,filetype in zip(self.get_filenames(), self.get_filetypes()):
            output_str += '\t'.join([filetype,filename,str(self.shift_dict[filename]), str(self.normalization_dict[filename])]) +'\n'
                
        output_str += 'file_format\t'+self.file_format+'\n'
        output_str += 'peaktype\t'+self.peaktype+'\n'
        output_str += 'windowsize\t'+str(self.window_size)+'\n'
        output_str += 'difftest\t'+str(self.difftest)+'\n'
        output_str += 'threshold\t'+str(self.threshold)+'\n'
        output_str += 'normalization\t'+self.normalization+'\n'
        if self.keep_max_dup > 0:
            output_str +='keep_max_dup\t'+str(self.keep_max_dup)+'\n'
        output_str += 'input_directory\t'+self.input_directory+'\n'
        output_str += 'output_directory\t'+self.output_directory+'\n'
        output_str += 'name\t'+self.name+'\n'
        output_str += 'num_processors\t'+str(self.num_procs)+'\n'
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
    
    def get_filetypes(self):
        return ['chip1']*len(self.chip1) + ['chip2'] * len(self.chip2) + ['input1'] * len(self.input1) + ['input2'] * len(self.input2)
    
    def get_filenames_wo_bin_dict(self):
        return [filename for filename in self.get_filenames() if filename not in self.bin_dict ]
    
    def get_chip_filenames(self):
        return self.chip1+self.chip2
        
    def get_genome_size(self):
        return sum(self.chr_info.values())
    
    def get_top3_chr(self):
        values = self.chr_info.values()
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
