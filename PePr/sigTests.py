
# initialize libraries.
import numpy
import logging
from classDef import Peak
from scipy.special import psi
from scipy.stats.distributions import norm
from operator import attrgetter

root_logger = logging.getLogger("")
debug = root_logger.debug
info = root_logger.info

def get_candidate_window2(read, x, y, repx, repy, threshold):
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
        N_WINDOWS = 10
    left = -1.0
    right = 6.0
    if weighted_log_likelihood(BASE**right, m, n, 
            read[(idx-N_WINDOWS):(idx+N_WINDOWS+1)], diff_test) > 0:
        return BASE**right
    elif weighted_log_likelihood(BASE**left, m, n,
            read[(idx-N_WINDOWS):(idx+N_WINDOWS+1)], diff_test) < 0:
        return BASE**left
    for times in xrange(N_ITERATION):
        temp = (right+left)/2
        if weighted_log_likelihood(10**temp, m, n,
                read[(idx-N_WINDOWS):(idx+N_WINDOWS+1)], diff_test) > 0:
            left = temp
        else:
            right = temp
    return BASE**right


def cal_FDR(peak_list, num_tests):
    '''Calculate the Benjamini-Hochberg FDR'''
    peak_list = sorted(peak_list, key=attrgetter('pvalue'))
    #calculate BH q-values
    q_list = [item.pvalue*num_tests/(idx+1) 
            for idx, item in enumerate(peak_list)]
    q_min = q_list[-1]
    for i in xrange(len(q_list)-1, -1, -1):
        if q_list[i] >= q_min:
            q_list[i] = q_min
        else:
            q_min = q_list[i]
    for i in range(len(q_list)):
        peak_list[i] = Peak(peak_list[i].chr, peak_list[i].index,
                            peak_list[i].pvalue, q_list[i])
    peak_list = sorted(peak_list, key=attrgetter('chr', 'index'))
    return peak_list


def negative_binomial(readData, peakfilename, swap, parameter):
    '''the main function that test for significant windows.'''
    # Initialize the parameters
    peaktype = parameter.peaktype
    remove_artefacts = parameter.remove_artefacts
    narrow_peak = parameter.narrow_peak_width
    threshold = parameter.threshold
    windowsize = parameter.window_size
    # Indicate the data 
    read_dict = readData.reads_dict 
    strands_dict = readData.data_dict_by_strands
    chr_list = readData.chr_list
    if parameter.difftest is True: 
        test_list = readData.chip1_filename_list
        control_list = readData.chip2_filename_list
    else: 
        test_list = readData.chip1_filename_list
        control_list = readData.input1_filename_list
    num_tests = readData.genome_size/windowsize
    peakfile = open(peakfilename, 'w')

    #compute number of replicates
    test_rep = len(test_list)
    control_rep = len(control_list)
    start1 = 0
    end1 = start2 = test_rep
    end2 = test_rep+control_rep
    # if swap 
    if swap is True: 
        test_rep, control_rep = control_rep, test_rep
        test_list, control_list = control_list, test_list

    # initialize basic array structures
    #sig_index_dict = {}
    sig_peaks_list = []    

    for chr in read_dict:

        read_array = read_dict[chr]
        read_array[numpy.where(read_array ==0)] = 1 
        y_bar_array = numpy.mean(read_array[:, start1:end1], 1)
        x_bar_array = numpy.mean(read_array[:, start2:end2], 1)
        if swap: #swap the chip and control reads. 
            x_bar_array, y_bar_array = y_bar_array, x_bar_array
        # setting the minimum # of reads in each window to 1 
        # so that they won't have arithmetic errors. 
    
        cand_index = get_candidate_window2(read_array, x_bar_array,
                y_bar_array, control_rep, test_rep, threshold)
        debug("there are %d candidate windows for %s", len(cand_index), chr)
        if not swap: 
            disp_list = numpy.array([estimate_area_dispersion_factor(read_array, 
                    test_rep, control_rep, idx, peaktype, parameter.difftest)
                    for idx in cand_index])
        else: 
            disp_list = numpy.array([estimate_area_dispersion_factor(read_array,
                    control_rep, test_rep, idx, peaktype, parameter.difftest)
                    for idx in cand_index])
        debug("finished estimating dispersion for %s", chr)
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
        # record results for potential windows
#        window_out = open(chr+"candidate.out",'w')
#        for i in range(len(cand_index)):
#            window_out.write(chr+'\t'+str(cand_index[i])+'\t'+str(disp_list[i])+'\t'+str(pval_array[i])+'\n')
#        window_out.close()
        # Record the indices of the windows that have p-value smaller
        # than the threshold. 
        test_index = numpy.where(pval_array<threshold)  
        test_index = test_index[0]
        sig_index = cand_index[test_index]
        sig_pval = pval_array[test_index]
        sig_disp = disp_list[test_index]
        for i, a in enumerate(test_index):
            sig_peaks_list.append(Peak(chr, sig_index[i], sig_pval[i], 0))

    #calculate the BH FDR. 
    debug("begin estimating fdr")
    sig_peaks_list = cal_FDR(sig_peaks_list, num_tests)
    debug("finished estimating fdr")

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
        sig_start, sig_end, sig_pval, sig_qval = merge_sig_window(sig_index,
                sig_pval, sig_qval, peaktype)
        for idx in xrange(len(sig_start)):
            final_peak = [chr, sig_start[idx]*windowsize/2, 
                          sig_end[idx]*windowsize/2+windowsize,
                          sig_pval[idx], sig_qval[idx]]
            final_peak_list.append(final_peak) 

    # post-processing 
    if remove_artefacts is True: 
         dump_file = open(peakfilename[:-10]+"_filtered_peaks.txt","w")
         dump_file.write("chr\tstart\tend\twidth\tp-value\tq-value\t")
         dump_file.write("chip_input_ratio\ttag_monopoly\t")
         dump_file.write("strand_overlap_original\tstrand_overlap_shifted\n")
    info("begin post-processing...")
    for final_peak in final_peak_list: 
        chr = final_peak[0]
        start = final_peak[1]
        end = final_peak[2]
        pval = final_peak[3]
        qval = final_peak[4]

        if swap is False: 
            chip_list = readData.chip1_filename_list
            input_list = readData.input1_filename_list
        else: 
            chip_list = readData.chip2_filename_list
            input_list = readData.input2_filename_list
        if narrow_peak or remove_artefacts: 
            (start, end, chip_input_ratio, tag_monopoly, overlap_orig, 
             overlap_roll) = post_processing_per_peak(
                    strands_dict, chip_list, input_list, chr, start, end,
                    readData.shift_size, readData.read_length,
                    narrow_peak, remove_artefacts)
            if remove_artefacts and (chip_input_ratio > 0.5 or 
                    tag_monopoly > 0.5 or (overlap_orig > 0.2 and
                overlap_roll/overlap_orig < 0.5)):
                dump_file.write(chr+"\t" + str(start) + '\t' + str(end) + '\t' +
                        str(end-start) + '\t' + str(pval) + '\t' + str(qval) + 
                        '\t' + str(chip_input_ratio) + '\t' + 
                        str(tag_monopoly) + '\t' + str(overlap_orig) +
                        '\t' + str(overlap_roll) + '\n')
                continue
        peakfile.write(chr + "\t" + str(start) + '\t' + str(end) + '\t' + 
                str(end-start) + '\t' + str(pval) + '\t' + str(qval) + '\n')
    return 

def merge_sig_window(index_list, pval_list, qval_list, peaktype):
    ''' Merge significant windows that are nearby. '''
    if peaktype=="sharp":
        # Mininal number of windows required. 
        MIN_WINDOW = 1 
        # The maximal gap (max-1) of significant windows allowed to merge.
        MAX_WINDOW = 2 
    elif peaktype =="broad":
        MIN_WINDOW = 2 
        MAX_WINDOW = 5
    sig_peak_start = []
    sig_peak_end = []
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
        else: 
            if pos - pre_start <= MAX_WINDOW:
                pre_start = pos
                peak_pval_list.append(pval_list[idx])
                peak_qval_list.append(qval_list[idx])

            elif pos-pre_start > MAX_WINDOW:
                end = pre_start
                if (end-start>=MIN_WINDOW-1):
                    sig_peak_start.append(start)
                    sig_peak_end.append(pre_start)
                    sig_pval.append(min(peak_pval_list))
                    sig_qval.append(min(peak_qval_list))
                start = pos
                pre_start = pos
                peak_pval_list = [pval_list[idx]]
                peak_qval_list = [qval_list[idx]]

    if (pre_start-start>MIN_WINDOW-1):
        sig_peak_start.append(start)
        sig_peak_end.append(pre_start)
        sig_pval.append(min(peak_pval_list))
        sig_qval.append(min(peak_qval_list))

    return sig_peak_start, sig_peak_end, sig_pval, sig_qval


def post_processing_per_peak(strands_dict, chip_list, input_list, chr,
                             start, end, shiftSize, readLength, narrow_peak,
                             remove_artefacts):
    ''' Remove artefacts and refine peak width.'''
    chip_forward = numpy.zeros(end-start)
    chip_reverse = numpy.zeros(end-start)
    input_forward = numpy.zeros(end-start)
    input_reverse = numpy.zeros(end-start)

    for chip in chip_list:
        forward = strands_dict[chr][chip]['f']
        reverse = strands_dict[chr][chip]['r']
        forward_read = forward[numpy.where( (forward > start-shiftSize[chip]) &
                (forward < end-shiftSize[chip]) )]
        reverse_read = reverse[numpy.where( (reverse > start+shiftSize[chip]) &
                (reverse < end+shiftSize[chip]) )]
        for read in forward_read:
            try: chip_forward[(read-start)] +=1
            except IndexError:
                debug("index error ignored. end-start: %d. read-start: %d.", end-start, read-start)
        for read in reverse_read:
            try: chip_reverse[(read-start-readLength)] +=1
            except IndexError:
                debug("index error ignored. end-start: %d. read-start-readLength: %d.", end-start, read-start-readLength)
    # using input reads to remove artefacts
    if remove_artefacts is True: 
        for input in input_list:
            forward = strands_dict[chr][input]['f']
            reverse = strands_dict[chr][input]['r']
            forward_read = forward[numpy.where( (forward > 
                    start-shiftSize[chip]) & (forward < end-shiftSize[chip]) )]
            reverse_read = reverse[numpy.where( (reverse > 
                    start+shiftSize[chip]) & (reverse < end+shiftSize[chip]) )]
            for read in forward_read:
                try: input_forward[(read-start)] +=1
                except IndexError: 
                    debug("index error ignored. end-start: %d. read-start: %d.", end-start, read-start)
            for read in reverse_read:
                try: input_reverse[(read-start-readLength)] +=1
                except IndexError:
                    debug("index error ignored. end-start: %d. read-start-readLength: %d.", end-start, read-start-readLength)

        chip_both = chip_forward + chip_reverse
        input_both = input_forward + input_reverse
        if sum(chip_both) == 0:
            pass
        else:
            chip_both = chip_both/sum(chip_both)
        if sum(input_both) ==0:
            pass
        else:
            input_both = input_both/sum(input_both)
        overlap_chip_input = numpy.sum(numpy.minimum(chip_both, input_both))
        chip_both.sort()
        chip_3_maximum = numpy.sum(chip_both[-3:])
        if sum(chip_reverse) != 0:
            chip_reverse = chip_reverse/sum(chip_reverse) 
        chip_forward_roll = numpy.roll(chip_forward,readLength)
        if sum(chip_forward) != 0:
            chip_forward = chip_forward/sum(chip_forward)
            chip_forward_roll = chip_forward_roll/sum(chip_forward)
        overlap_orig = numpy.sum(numpy.minimum(chip_forward, chip_reverse))
        overlap_orig = numpy.max([overlap_orig, 1e-5])
        overlap_roll = numpy.sum(numpy.minimum(chip_forward_roll, chip_reverse))
    else: 
        overlap_chip_input = 0
        chip_3_maximum = 0 
        overlap_orig = 1e-5
        overlap_roll = 0

    if narrow_peak is True:
        sum_forward = 0
        sum_reverse = 0
        if sum(chip_forward) > 0: 
            chip_forward = chip_forward/sum(chip_forward)
            for i in xrange(end-start): 
                sum_forward += chip_forward[i]
                if sum_forward > 0.2: 
                    new_start = start + i
                    break
        if sum(chip_reverse) > 0:
            chip_reverse = chip_reverse/sum(chip_reverse)
            for i in xrange(end-start-1, -1, -1):
                sum_reverse += chip_reverse[i]
                if sum_reverse > 0.2:
                    new_end = start + i
                    break
        if sum(chip_forward) == 0:
            new_start = new_end - 2*shiftSize
        if sum(chip_reverse) == 0:
            new_end = new_start + 2*shiftSize
        start = new_start
        end = new_end

    return (start, end, overlap_chip_input, chip_3_maximum, 
            overlap_orig, overlap_roll)


def shift_size_per_peak(strands_dict, chip_list, chr, start, end,
                        shiftSize, readLength, narrow_peak):
    '''Deprecated function. Used to estimate the shift size of 
       peak'''
    shift_list = []
    start_list = []
    end_list = []

    for chip in chip_list:
        forward = strands_dict[chr][chip]['f']
        reverse = strands_dict[chr][chip]['r']

        forward_read = forward[numpy.where( (forward > start-shiftSize[chip]) &
                (forward < end-shiftSize[chip]) )]
        reverse_read = reverse[numpy.where( (reverse > start+shiftSize[chip]) &
                (reverse < end+shiftSize[chip]) )]
        if len(reverse_read) ==0: 
            fold = 100  
        # If there is no reads in the reverse strand, define 
        # the fold (forward over reverse) to be 100 to avoid 
        # arithmetic errors.
        else:
            fold = float(len(forward_read))/len(reverse_read)
        # estimate the shift size
        forward = numpy.zeros(end-start)
        reverse = numpy.zeros(end-start)
        for read in forward_read:
            try: forward[(read-start):(read-start+readLength)] +=1
            except IndexError:
                debug("index error ignored. end-start: %d. read-start+readLength: %d.", end-start, read-start+readLength)
        for read in reverse_read:
            try: reverse[(read-start-readLength):(read-start)] +=1
            except IndexError: 
                debug("index error ignored. end-start: %d. read-start: %d.", end-start, read-start)
        shade_max = 0
        shift_max = 0
        # Iterate from 0 to 250 bp, find the optimum shift that produce 
        # the most overlap between the two strands. 
        for shift in xrange(-readLength, 250-readLength): 
            forward_temp = numpy.roll(forward, shift)
            shade = numpy.sum(numpy.min(numpy.vstack(
                    [forward_temp, reverse]), 0))
            if shade > shade_max:
                shift_max = shift
                shade_max = shade
        shift_list.append(shift_max+readLength)
        if narrow_peak:
            forward_read.sort()
            reverse_read.sort()
            forward_len = len(forward_read)
            reverse_len = len(reverse_read)
            if forward_len>0:
                start_list.append(forward_read[numpy.floor(forward_len*0.2)])
            if reverse_len>0: 
                end_list.append(reverse_read[-numpy.floor(reverse_len*0.2)])
    if narrow_peak:
        if (len(start_list) ==0) | (len(end_list)==0):
            return -1
        else:
            return min(start_list), max(end_list), shift_list, fold
    else:
        return start, end, shift_list, fold
