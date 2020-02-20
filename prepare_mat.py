import argparse
from os.path import dirname, isdir
from os import mkdir, getcwd
import shutil
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--maf", dest="maf", help="MAF file from which to extract matrix.", required=True)
    parser.add_argument("-i", "--indels", action="store_true", dest="indels", help="Extract indel matrix, as well as SBS", default=False)
    parser.add_argument("-e", "--exome", dest="exome", help="Exome data - restrict genome to exome regions", default=False)
    parser.add_argument("-p", "--project", dest="project", default="PROJECT", help="Project name for output.")
    parser.add_argument("-d", "--directory", dest="directory", default=None, help="Input/Output directory")
    parser.add_argument("-P", "--plot", dest="plot", action="store_true", help="Output plots of input data.")
    parser.add_argument("-c", "--cushion", dest="cushion", type=int, default=100, help="Adds an Xbp cushion to the exome/bed_file ranges for downsampling the mutations.")

    return parser.parse_args()


if __name__ == "__main__":
    
    args = parse_args()
    
    if args.directory is None:
        current = getcwd()
        args.directory = current + "/" + "sigprof_input"
    try:
        if not isdir(args.directory):
            mkdir(args.directory)
    except:
        print("ERROR: creation of directory", args.directory, "failed. Please use the -d option to create a valid directory.")
    try:
        shutil.copyfile(args.maf, args.directory + "/" + args.maf)
    except:
        print("File copy failed.", args.maf, args.directory + "/" + args.maf)

    matrices = matGen.SigProfilerMatrixGeneratorFunc(args.project, "GRCh37", args.directory, plot=args.plot, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=args.cushion)
