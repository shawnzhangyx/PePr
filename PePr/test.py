from fileParser import FileParser
import logging 
logging.basicConfig(level=logging.DEBUG, format = "%(levelname)-6s %(name)s %(asctime)s | %(message)s", datefmt = "%m/%d/%Y %I:%M:%S %p")
if __name__ == "__main__":
	file = FileParser("/data/yanxiao/test/pepr/new_pepr/chr7/Sample_16028_R1.chr7.bed", "bed")
	file.get_chr_info()
	file.parse_data_bin(200, 90, 2)
