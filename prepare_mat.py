import argparse
from os.path import dirname
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--maf", dest="maf", help="MAF file from which to extract matrix.", required=True)
    parser.add_argument("-i", "--indels", action="store_true", dest="indels", help="Extract indel matrix, as well as SBS", default=False)
    parser.add_argument("-e", "--exome", dest="exome", help="Exome data - restrict genome to exome regions", default=False)
    parser.add_argument("-p", "--project", dest="project", default="PROJECT", help="Project name for output.")
    parser.add_argument("-d", "--directory", dest="directory", default="input", help="Input/Output directory")
    parser.add_argument("-P", "--plot", dest="plot", action="store_true", help="Output plots of input data.")

    return parser.parse_args()


if __name__ == "__main__":
    
    args = parse_args()

    matrices = matGen.SigProfilerMatrixGeneratorFunc(args.project, "GRCh37", dirname(args.maf), plot=args.plot, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
