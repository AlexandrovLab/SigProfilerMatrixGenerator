import argparse
import sys
from collections import defaultdict
from pyfaidx import Fasta

"""
A debugging function.
Given a list of values with any type,
write those values to stderr separated by spaces.
"""
def write_err(text):
    sys.stderr.write(" ".join([str(i) for i in text]) + "\n")

"""
Creates a new header dictionary from a tab-separated header.
This can then be used to get the index of a given column name.
"""
def reheader(headerline):
    header = defaultdict(int)
    headertokens = headerline.strip("#").strip().split("\t")

    for i in range(0, len(headertokens)):
        header[headertokens[i]] = i

    return header 
"""
Undoes alignment trimming of indels,
replacing bases represented with a "-" with
the correct base in the reference
and prepending the reference base to the non-dash allele.
"""
def untrim_indel(chrom, start, end, ref, alt, fa_ref):
    ## Get the REF allele
    real_ref = str(fa_ref[chrom][int(start)-1-1:int(start)-1])
    if alt == "-":
        ref = real_ref + ref
        alt = real_ref
        start = int(start) - 1
    elif ref == "-":
        ref = real_ref
        alt = real_ref + alt
        start = int(start) - 1
    else:
        write_err(["Invalid allele: ", ref, "->", alt, "."])
        raise Exception("Exiting.")
    return chrom, start, end, ref, alt

def write_matlab_csv(d, features):
    return

def write_python_tsv(dict, features):
    return

def sbs96(tokens, sbs_d, ref):
    return
def sbs192(tokens, sbs_d, ref):
    return
def sbs1536(tokens, sbs_d, ref):
    return

def indel83(tokens, id_d, ref):
    return
def indel28(tokens, id_d, ref):
    return


"""
Creates a minimal representation of a variant that is compatible with the SigProfiler text format.
"""
def make_minimal_record(cancer_type, sample, assay, genome, variant_type, chrom, pos, end, ref, alt, mutation_type):
    vals = [cancer_type, sample, assay, genome, variant_type, chrom, pos, end, ref, alt, mutation_type]
    return "\t".join(vals)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--maf", required=True, dest="maf", help="MAF file from which to extract matrix.")
    parser.add_argument("-i", "--indels", action="store_true", dest="indels", help="Extract indel matrix, as well as SBS", default=False)
    parser.add_argument("-s", "--sigprofiler", help = "Output in sigprofiler format", action = "store_true")
    parser.add_argument("-e", "--exome", dest="exome", help="Exome data - restrict genome to exome regions", default=False)
    parser.add_argument("-p", "--project", dest="project", default="PROJECT", help="Project name for output.")
    parser.add_argument("-d", "--directory", dest="directory", default="input", help="Input/Output directory")
    parser.add_argument("-f", "--fasta", required=True, dest="ref", help="A fasta file reference.")
    parser.add_argument("-C", "--custom-id", dest="custom_id", help="Column name of custom, preferred ID field", default=None, type = str)

    return parser.parse_args()


if __name__ == "__main__":
    ## Holds the sample->mutational counts vectors
    ## Each Key is a sample name
    ## and each value is an abstracted list of lists,
    ## where sample_d["SAMP"]["SBS96"] holds the mutational count vector for sample SAMP
    ## and mutational type (e.g. ID83, SBS96, etc.) "SBS96".
    sample_d = defaultdict(lambda: defaultdict(list))   


    header_d = {
    "Chromosome" : 4,
    "Start_position" : 5,
    "End_position" : 6,
    "Strand" : 7,
    "Variant_Type" : 9,
    "Reference_Allele" : 10,
    "Tumor_Seq_Allele1" : 11,
    "Tumor_Seq_Allele2" : 12,
    "Tumor_Sample_Barcode": 15,
    "Matched_Norm_Sample_Barcode" : 16
    }

    args = parse_args()

    header_line = None
    ref = None
    id_field = None
    if args.ref is not None:
        ref = Fasta(args.ref)
    if args.custom_id is not None:
        id_field = args.custom_id
    else:
        id_field = "Tumor_Sample_Barcode"
    

    with open(args.maf, "r") as mfi:
        for line in mfi:
            line = line.strip()
            tokens = line.split("\t")
            if "Hugo_Symbol" not in line:
                vtype = tokens[header_d["Variant_Type"]]

                if vtype == "SNP" or vtype == "DNP" or vtype == "TNP" or vtype == "QNP" or vtype == "MNP":
                    chrom = tokens[header_d["Chromosome"]]
                    sample = tokens[header_d[id_field]]
                    start_pos = tokens[header_d["Start_position"]]
                    end_pos = tokens[header_d["End_position"]]
                    ref_allele = tokens[header_d["Reference_Allele"]]
                    alt_allele = tokens[header_d["Tumor_Seq_Allele2"]]
                    fasta_allele = ref[chrom][int(start_pos)-1:int(end_pos)]
                    if str(fasta_allele) == str(ref_allele) and args.sigprofiler:
                        ## make_minimal_record(cancer_type, sample, assay, genome, variant_type, chrom, pos, ref, alt, mutation_type = "Somatic"):
                        print(make_minimal_record(args.project, sample, "WGS", "GRCh37", vtype, chrom, start_pos, end_pos, ref_allele, alt_allele, "SOMATIC"))
                    elif str(fasta_allele) == str(ref_allele):
                        write_err(["Error: non-sigprofiler output not yet implemented."])
                    else:
                        write_err(["Error: ref allele [", 
                        ref_allele,
                         "] and fasta allele [", fasta_allele,
                          "] don't match." ])
                        raise Exception("Error: mismatched fasta and ref alleles.")
                    
                elif vtype == "DEL" or vtype == "INS":
                    chrom = tokens[header_d["Chromosome"]]
                    sample = tokens[header_d[id_field]]
                    start_pos = tokens[header_d["Start_position"]]
                    end_pos = tokens[header_d["End_position"]]
                    ref_allele = tokens[header_d["Reference_Allele"]]
                    alt_allele = tokens[header_d["Tumor_Seq_Allele2"]]

                    fasta_allele = ref[chrom][int(start_pos)-1:int(end_pos)]
                    #print(start_pos, fasta_allele, ref[chrom][int(start_pos) - 1: int(end_pos)], ref_allele, alt_allele)
                    #assert fasta_allele == ref_allele

                    chrom, start_pos, end_pos, ref_allele, alt_allele = untrim_indel(chrom, start_pos, end_pos, ref_allele, alt_allele, ref)
                    if args.sigprofiler:
                        print(make_minimal_record(args.project, sample, "WGS", "GRCh37", vtype, chrom, str(start_pos), str(end_pos), ref_allele, alt_allele, "SOMATIC"))

                    #print(chrom, start_pos, end_pos, ref_allele, alt_allele)

                else:
                    write_err(["Invalid mutation type", vtype])
            else:
                header_line = line
                if (args.custom_id):
                    header_d = reheader(header_line)
            

    




    
