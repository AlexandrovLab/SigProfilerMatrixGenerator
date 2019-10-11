import sys
from collections import defaultdict


"""
Usage: python python_to_matlab_csv.py <mutationalCountsFile> <ProjectName>
"""
if __name__ == "__main__":

    header = None
    project = sys.argv[2]
    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            if "MutationType" in line:
                header = line.strip().split("\t")
                header[0] = "Mutation Type"
                header.insert(1, "Trinucleotide")
                header = [project + ":" + header[i] if i > 1 else header[i] for i in range(0, len(header))]
                header = ["\"" + i + "\"" for i in header]
                print(",".join(header))
            else:
                line = line.strip()
                tokens = line.split("\t")
                mtype = tokens[0].strip("ACTG").strip("[]")
                trinuc = tokens[0][0] + tokens[0][2] + tokens[0][6]
                tokens[0] = mtype
                tokens.insert(1, trinuc)
                tokens = ["\"" + tokens[i] + "\"" if i < 2 else i for i in range(0, len(tokens))]
                print(",".join([str(i) for i in tokens]))
