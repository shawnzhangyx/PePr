from optparse import OptionParser
from classDef import Parameters


def opt_parser(argv):
    parser = OptionParser()
    parser.add_option(
            "-p", "--parameter-file", action="store", type="string",
            dest="parameter", default="",
            help="provide a file that contain the parameters",
            metavar="PARAMETER-FILE")
    parser.add_option(
            "-c", "--chip1", action="store", type="string", dest="chip1",
            default="", help="chip1 file names separated by comma",
            metavar="CHIP1")
    parser.add_option(
            "-i", "--input1", action="store",
            type="string", dest="input1", default="",
            help="input1 file names separated by comma",
            metavar="INPUT1")
    parser.add_option(
            "--chip2", action="store", type="string", dest="chip2",
            default="", help="chip2 file names separated by comma",
            metavar="CHIP2")
    parser.add_option(
            "--input2", action="store", type="string", dest="input2",
            default="", help="input2 file names separated by comma",
            metavar="INPUT2")
    parser.add_option(
            "-f", "--file-format", action="store", type="string",
            dest="file_format",
            help="bed, sam, bam,...",
            metavar="FORMAT")
    parser.add_option(
            "-s", "--shiftsize", action="store",
            type="int", dest="shift_size", default=-1,
            help="Half the fragment size.", metavar="SHIFTSIZE")
    parser.add_option(
            "-w", "--windowsize", action="store",
            type="int", dest="window_size", default=-1,
            help="Window sizes",
            metavar="WINDOWSIZE")
    parser.add_option(
            "--diff", action="store_true",
            dest = "difftest", default=False,
            help="Perform differential binding instead of peak-calling")
    parser.add_option(
            "-n", "--name", action = "store",
            type="string", dest="name", default = "NA",
            help = "the experimental name. NA if none provided",
            metavar="NAME")
    parser.add_option(
            "--threshold", action ="store",
            type='float', dest="threshold", default=1e-5,
            help="p-value threshold. Default 1e-5.")
    parser.add_option(
            "--peaktype", action="store",
            type="string", dest="peaktype", default="broad",
            help="sharp or broad. Default broad.")
    parser.add_option(
            "--num-processors", action="store",
            type="int", dest="num_procs", default =1,
            help="number of CPUs use on it.")
    parser.add_option(
            "--input-directory", action="store",
            type="string", dest="input_directory", default='./',
            help="where the data files are")
    parser.add_option(
            "--output-directory", action="store",
            type="string", dest="output_directory", default='./',
            help="where you want the output files to be")
    parser.add_option(
            "--custom_normalization", action="store",
            type="string", dest="normalization", default="scale", 
            help='''Normalization method. scale or compound''')
    '''parser.add_option(
            "--no_log", action="store_true",
            dest = "unsave_log", default=False,
            help = "Disable saving the log files")
    '''
    parser.add_option(
            "--version", action="store_true",
            dest="version",default=False,
            help="Show version information and exit")

    (opt, args)=parser.parse_args(argv)
    if len(argv)==1:
        parser.print_help()
        exit(1)
    if opt.version == True:
        from  __init__ import __version__
        print "PePr "+ __version__
    exit(1)

    return opt

def process_opt(opt):
    ''' validate the parameters that the user specified'''
    # initial process the filenames.
    opt.chip1 = [] if opt.chip1 is '' else opt.chip1.strip().split(',') 
    opt.chip2 = [] if opt.chip2 is '' else opt.chip2.strip().split(',')
    opt.input1 = [] if opt.input1 is '' else opt.input1.strip().split(',')
    opt.input2 = [] if opt.input2 is '' else opt.input2.strip().split(',')
    opt.normalization = opt.normalization.lower()
    parameter = Parameters(opt)

    #add shift size validations
    return parameter
