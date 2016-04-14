#from logging import info
from collections import Counter




def estimate_shiftsize(data):
    counter = Counter()
    for key in data:
        counter[key] = len(data[key])

    top = counter.most_common(2000)
    shift_list = []
    for window in top:
        key = window[0]
        reads = data[key]
        pos = [item[0] for item in reads if item[1] is '+']
        neg = [item[0] for item in reads if item[1] is '-']
        if len(neg) == 0 or len(pos) == 0: 
            continue
        shift = sum(neg)/len(neg) - sum(pos)/len(pos)
        shift_list.append(shift)
    
    return sum(shift_list)/len(shift_list), counter

