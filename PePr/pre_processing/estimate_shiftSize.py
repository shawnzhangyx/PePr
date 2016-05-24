###### you can use this script to estimate shift size for each sample.
# example: python estimate_shiftSize.py file_format filename
import sys
from shiftSize import estimate_shiftsize
from initialize import get_chromosome_info
from initialize import get_read_length_info


class TestParameter:
    def __init__(self,argv):
        self.file_format = argv[1]
        self.filename = [argv[2]]
        self.chr_info = {} 
        self.input_directory = ''
        self.read_length_dict = {}

    def get_filenames(self):
        return self.filename

if __name__ == "__main__":
    parameter = TestParameter(sys.argv)
    print sys.argv
    get_chromosome_info(parameter, sys.argv[2])
    get_read_length_info(parameter)
    print parameter.chr_info
    estimate_shiftsize(sys.argv[2], parameter)

