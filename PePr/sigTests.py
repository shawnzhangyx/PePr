
# initialize libraries.
import os
import numpy
import logging
from classDef import Peak
from scipy.special import psi
from scipy.stats.distributions import norm
from operator import attrgetter, itemgetter
import multiprocessing
from multiprocessing import sharedctypes
import itertools
import time
import gc 

debug = logging.debug
info = logging.info

def get_candidate_window2( x, y, repx, repy, threshold):
    # using PHI = 1e6 to prescreen the genome 
    PHI = 1e6
    GAMMA_HAT = 1.0
    tau = numpy.sqrt( y*( repx*x*(PHI+y) + repy*y*(PHI+x))/repx/repy/PHI/x**3)
    gamma = y/x
    z = (numpy.log(gamma)-numpy.log(GAMMA_HAT))*gamma/tau
    pvalue = norm.cdf(-z)
    pre_idx_list = numpy.where(pvalue[10:-10]<threshold)[0]+10
    return numpy.array(pre_idx_list)


def weighted_log_likelihood(v_hat, m, n, reads, diff_test):
    '''This function returns the derivative of the log likelihood
       of current and adjacent windows using a triangular weight.'''
    equation = 0
    n_window = len(reads)
    baseline = numpy.mean(reads[n_window/2,0:m])
    for idx in range(n_window):
        x = reads[idx, 0:m]  # x/m refer the test sample
        y = reads[idx, m:(m+n)]  # y/n refer to the control sample
        weight_x = numpy.mean(x)/baseline
        weight_y = numpy.mean(y)/baseline
        if n == 1: 
            weight_y = 0
        log_likelihood = (-(m*weight_x+n*weight_y)*psi(v_hat) +
                numpy.sum(psi(v_hat+x))*weight_x + 
                numpy.sum(psi(v_hat+y))*weight_y + 
                m*numpy.log(v_hat/(v_hat+numpy.mean(x)))*weight_x + 
                n*numpy.log(v_hat/(v_hat+numpy.mean(y)))*weight_y)
        equation = (equation + 
                log_likelihood*(1-(abs(n_window/2-idx)/(n_window/2+1))))
    return equation


def estimate_area_dispersion_factor(read, m, n, idx, peaktype, diff_test):
    '''This function implements a simple method to find the root
       of the derivative of the log likelihood.'''
    BASE = 10
    N_ITERATION = 15
    if peaktype == "sharp":
        N_WINDOWS = 1 
    elif peaktype  == "broad":
        N_WINDOWS = 5
    left = -1.0
    right = 6.0
    if weighted_log_likelihood(BASE**right, m, n, 
            read[(idx-N_WINDOWS):(idx+N_WINDOWS+1)], diff_test) > 0:
        return BASE**right
    elif weighted_log_likelihood(BASE**left, m, n,
            read[(idx-N_WINDOWS):(idx+N_WINDOWS+1)], diff_test) < 0:
        return BASE**left
    for times in range(N_ITERATION):
        temp = (right+left)/2
        if weighted_log_likelihood(10**temp, m, n,
                read[(idx-N_WINDOWS):(idx+N_WINDOWS+1)], diff_test) > 0:
            left = temp
        else:
            right = temp
    return BASE**right


def cal_FDR(peak_list, num_tests):
    '''Calculate the Benjamini-Hochberg FDR'''
    if len(peak_list)==0:
        return peak_list
    peak_list = sorted(peak_list, key=attrgetter('pvalue'))
    #calculate BH q-values
    q_list = [item.pvalue*num_tests/(idx+1) 
            for idx, item in enumerate(peak_list)]
    q_min = min(q_list[-1],1) # q-value should be less than 1. 
    for i in range(len(q_list)-1, -1, -1):
        if q_list[i] >= q_min:
            q_list[i] = q_min
        else:
            q_min = q_list[i]
    for i in range(len(q_list)):
        peak_list[i] = Peak(peak_list[i].chr, peak_list[i].index,
                            peak_list[i].g1_count, peak_list[i].g2_count,
                            peak_list[i].pvalue, q_list[i])
    peak_list = sorted(peak_list, key=attrgetter('chr', 'index'))
    return peak_list


    
def per_chr_nbtest(read_array, chr, swap,threshold, peaktype,difftest, start1,end1,start2,end2,test_rep,control_rep):
    t1 = time.time()
    sig_peaks_list = []
    y_bar_array = numpy.mean(read_array[:, start1:end1], 1)
    x_bar_array = numpy.mean(read_array[:, start2:end2], 1)
    if swap: #swap the chip and control reads.
        x_bar_array, y_bar_array = y_bar_array, x_bar_array
    cand_index = get_candidate_window2( x_bar_array,
                    y_bar_array, control_rep, test_rep, threshold)
    debug("There are %d candidate windows for %s (PID:%d)", len(cand_index), chr, os.getpid())
    if not swap:
        disp_list = numpy.array([estimate_area_dispersion_factor(read_array,
                test_rep, control_rep, idx, peaktype, difftest)
                for idx in cand_index])
    else:
        disp_list = numpy.array([estimate_area_dispersion_factor(read_array,
                control_rep, test_rep, idx, peaktype, difftest)
                for idx in cand_index])
    #debug("finished estimating dispersion for %s", chr)
   # return []
    cand_x_bar_array = x_bar_array[cand_index]
    cand_y_bar_array = y_bar_array[cand_index]
    gamma_array = cand_y_bar_array / cand_x_bar_array
    tau_hat_array = numpy.sqrt(cand_y_bar_array*
            ((control_rep*cand_x_bar_array*(disp_list+cand_y_bar_array)) +
            (test_rep*cand_y_bar_array*(disp_list+cand_x_bar_array)))/
            (test_rep*control_rep*disp_list*(cand_x_bar_array**3)))

    gamma_hat = 1.0 #Null hypothesis
    z_score_array = ((numpy.log(gamma_array)-numpy.log(gamma_hat))*
            gamma_array/tau_hat_array)
    pval_array = norm.cdf(-z_score_array)
    test_index = numpy.where(pval_array<threshold)
    test_index = test_index[0]
    sig_index = cand_index[test_index]
    sig_pval = pval_array[test_index]
    sig_group1_count = cand_y_bar_array[test_index]
    sig_group2_count = cand_x_bar_array[test_index]
    #sig_disp = disp_list[test_index]
    for i, a in enumerate(test_index):
        sig_peaks_list.append(Peak(chr, sig_index[i], sig_group1_count[i], sig_group2_count[i], sig_pval[i], 0))
    t2 = time.time()
    debug ("Analysis finished for %s, used %f sec CPU time", chr, t2-t1)
    return sig_peaks_list
	
def per_chr_nbtest_helper(args):
	### help the workers to take multiple arguments 
    try: 
        return per_chr_nbtest(*args)
    except KeyboardInterrupt, e:
        pass 

def mock_helper(args):
    return mock(*args)

def mock(d,p):
    return []        
 
def negative_binomial(read_dict,peakfilename, swap, parameter):
    '''the main function that test for significant windows.'''
    print len(read_dict)
    # Initialize the parameters
    peaktype = parameter.peaktype
    threshold = parameter.threshold
    windowsize = parameter.window_size
    # Indicate the data 
    if parameter.difftest is True: 
        test_list = parameter.chip1
        control_list = parameter.chip2
    else: 
        test_list = parameter.chip1
        control_list = parameter.input1
    num_tests = parameter.get_genome_size()/windowsize
    

    #compute number of replicates
    test_rep = len(test_list)
    control_rep = len(control_list)
    start1 = 0
    end1 = start2 = test_rep
    end2 = test_rep+control_rep
    # if swap 
    if swap is True: 
        test_rep, control_rep = control_rep, test_rep
        #test_list, control_list = control_list, test_list

    # initialize basic array structures
    sig_peaks_list = []    
    
    # single-core version. 
    if parameter.num_procs <2:
        for chr in parameter.chr_info:
            read_array = read_dict[chr]
            sig_peaks_list.extend(per_chr_nbtest(read_array, chr, swap,threshold, peaktype,parameter, start1,end1,start2,end2,test_rep,control_rep))
    # multi-core version
    else: 
        result_list = []
        def log_result(result):
            result_list.append(result)
        try: 
            import sharedmem
            for chr in parameter.chr_info:
                read_array = read_dict[chr]
                read_dict[chr] = sharedmem.copy(read_array)
        except ImportError:
            print "Import sharedmem package failed"
            
        pool = multiprocessing.Pool(processes=parameter.num_procs)#,maxtasksperchild=1)
        for chr in parameter.chr_info:
            read_array = read_dict[chr]
            pool.apply_async(per_chr_nbtest, (read_array, chr, swap,threshold, peaktype, parameter.difftest, 
                       start1,end1,start2,end2,test_rep,control_rep),callback=log_result)
        pool.close()
        pool.join()
        sig_peaks_list = list(itertools.chain(*result_list))
            
    #calculate the BH FDR. 
    debug("begin estimating FDR")
    sig_peaks_list = cal_FDR(sig_peaks_list, num_tests)
    debug("finished estimating FDR")

    # merge adjacent significant peaks. 
    info ("Merging adjacent significant windows...")
    final_peak_list = []
    for chr in read_dict:
        sig_peak_list_by_chr = \
                [item for item in sig_peaks_list if item.chr==chr]
        if len(sig_peak_list_by_chr) == 0:
            continue  # if there is no significant peak in this chromosome, skip it.
        sig_index = [item.index for item in sig_peak_list_by_chr]
        sig_pval = [item.pvalue for item in sig_peak_list_by_chr]
        sig_qval = [item.qvalue for item in sig_peak_list_by_chr]
        sig_g1_count = [item.g1_count for item in sig_peak_list_by_chr]
        #print sig_g1_count, len(sig_g1_count)
        sig_g2_count = [item.g2_count for item in sig_peak_list_by_chr]
        sig_start, sig_end, sig_fc, sig_pval, sig_qval = merge_sig_window(sig_index,
                sig_g1_count, sig_g2_count, sig_pval, sig_qval, peaktype)
        for idx in range(len(sig_start)):
            final_peak = [chr, sig_start[idx]*windowsize/2, 
                          sig_end[idx]*windowsize/2+windowsize, sig_fc[idx],
                          sig_pval[idx], sig_qval[idx]]
            final_peak_list.append(final_peak) 
    # sort the peak list
    final_peak_list = sorted(final_peak_list, key=itemgetter(4))
    info("%d peaks called.", len(final_peak_list))
    if len(final_peak_list) == 0:
        return 
    #start output peaks. 
    all_fc = [peak[3] for peak in final_peak_list] 
    #print all_fc
    max_fc = max(all_fc)
    # write results to peak file.
    peakfile = open(peakfilename, 'w')
    for idx, final_peak in enumerate(final_peak_list): 
        chr = final_peak[0]
        start = final_peak[1]
        end = final_peak[2]
        fc = final_peak[3]
        pval = final_peak[4]
        qval = final_peak[5]
        # tentatively, assign the normalized fold change as the score in 5th column
        score = fc/max_fc*1000 # range from 0 to 1000
        peakfile.write( '\t'.join([chr, str(start), str(end), 
                        ("chip2" if swap else "chip1") + "_peak_" +str(idx+1), 
                        str(score), '.', str(fc), str(pval), str(qval)]) + '\n')
    return 

def merge_sig_window(index_list, g1_count_list, g2_count_list, pval_list, qval_list, peaktype):
    ''' Merge significant windows that are nearby. '''
    if peaktype=="sharp":
        # Mininal number of windows required. 
        MIN_WINDOW = 1 
        # The maximal gap (max-1) of significant windows allowed to merge.
        MAX_WINDOW = 2 
    elif peaktype =="broad":
        MIN_WINDOW = 1 
        MAX_WINDOW = 5
    sig_peak_start = []
    sig_peak_end = []
    sig_fc = []
    sig_pval = []
    sig_qval = []
    for idx, pos in enumerate(index_list):
        if idx == 0:
            sig_peak_start = []
            sig_peak_end = []
            start = pos
            pre_start = pos
            peak_pval_list = [pval_list[idx]]
            peak_qval_list = [qval_list[idx]]
            g1_count = g1_count_list[idx]
            g2_count = g2_count_list[idx]
        else: 
            if pos - pre_start <= MAX_WINDOW:
                pre_start = pos
                peak_pval_list.append(pval_list[idx])
                peak_qval_list.append(qval_list[idx])
                g1_count += g1_count_list[idx]
                g2_count += g2_count_list[idx]
            elif pos-pre_start > MAX_WINDOW:
                end = pre_start
                if (end-start>=MIN_WINDOW-1):
                    sig_peak_start.append(start)
                    sig_peak_end.append(pre_start)
                    sig_pval.append(min(peak_pval_list))
                    sig_qval.append(min(peak_qval_list))
                    sig_fc.append(g1_count/g2_count)
                start = pos
                pre_start = pos
                peak_pval_list = [pval_list[idx]]
                peak_qval_list = [qval_list[idx]]
                g1_count = g1_count_list[idx]
                g2_count = g2_count_list[idx]
                
    if (pre_start-start>MIN_WINDOW-1):
        sig_peak_start.append(start)
        sig_peak_end.append(pre_start)
        sig_pval.append(min(peak_pval_list))
        sig_qval.append(min(peak_qval_list))
        sig_fc.append(g1_count/g2_count)
        
    return sig_peak_start, sig_peak_end, sig_fc, sig_pval, sig_qval

