#logConfig.py


import logging
import sys
import time
import os

def startLog(filename):
    #configure root settings.
    logging.basicConfig(level=logging.DEBUG, format = "%(levelname)-6s %(name)s %(asctime)s | %(message)s", datefmt = "%m/%d/%Y %I:%M:%S %p")
    format_basic = logging.Formatter("%(levelname)-7s %(asctime)s | %(message)s", datefmt = "%m/%d/%Y %I:%M:%S %p")

    #Create handlers that prints DEBUG messages to log files
    log_filename = '-'.join([filename, time.strftime("%m-%d-%Y-%I.%M"), str(os.getpid()), "debug.log"])
    # check if the same file name existk, if yes, add an appendix to the file name and then create the file. 
    # although not very likely
   
    file_hand = logging.FileHandler(log_filename, 'w')
    file_hand.setLevel(logging.DEBUG)
    file_hand.setFormatter(format_basic)
    logging.getLogger().addHandler(file_hand)


