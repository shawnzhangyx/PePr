#logConfig.py


import logging
import sys
import time

#configure root settings.
#logging.basicConfig(level=logging.DEBUG, format = "%(levelname)s %(name)s %(asctime)s %(message)s", datefmt = "%m/%d/%Y %I:%M:%S %p")


#Create logging format for the handlers
#format_basic = logging.Formatter("%(levelname)-10s %(name)s %(asctime)s # %(message)s", datefmt = "%m/%d/%Y %I:%M:%S")
format_basic = logging.Formatter("%(levelname)-7s %(asctime)s # %(message)s", datefmt = " %I:%M:%S")
format_log =logging.Formatter("%(message)s")

#Create handler that prints INFO messages to stderr
info_hand = logging.StreamHandler(sys.stderr)
info_hand.setLevel(logging.DEBUG)
info_hand.setFormatter(format_basic)


#Create handlers that prints DEBUG messages to log files
log_filename = time.strftime("%m-%d-%Y-%I.%M")+"-debug.log"
debug_hand = logging.FileHandler(log_filename, 'w')
debug_hand.setLevel(logging.DEBUG)
debug_hand.setFormatter(format_basic)

'''
shiftSize_hand = logging.FileHandler("shift.log", 'w')
shiftSize_hand.setLevel(logging.DEBUG)
shiftSize_hand.setFormatter(format_log)

windowSize_hand = logging.FileHandler("window.log", 'w')
windowSize_hand.setLevel(logging.DEBUG)
windowSize_hand.setFormatter(format_log)
'''


#Create a root logger
root_log = logging.getLogger("")
root_log.setLevel(logging.DEBUG)
root_log.addHandler(info_hand)
root_log.addHandler(debug_hand)

'''
#Create a logger that handles the shift sizes
shiftSizeEst_log = logging.getLogger("shiftSizeEst")
shiftSizeEst_log.addHandler(shiftSize_hand)

#Create a logger that handles the window sizes
windowSizeEst_log = logging.getLogger("windowSizeEst")
windowSizeEst_log.addHandler(windowSize_hand)
'''


