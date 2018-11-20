#!/usr/bin/env python3
# This source code file is a part of SigProfilerMatrixGenerator

# SigProfilerMatrixGenerator is a tool included as part of the SigProfiler

# computational framework for comprehensive analysis of mutational

# signatures from next-generation sequencing of cancer genomes.

# SigProfilerMatrixGenerator provides a standard tool for displaying all 

# types of mutational signatures as well as all types of mutational 

# patterns in cancer genomes. The tool seamlessly integrates with 

#other SigProfiler tools.

# Copyright (C) 2018 Erik Bergstrom

#

# SigProfilerMatrixGenerator is free software: you can redistribute it and/or modify

# it under the terms of the GNU General Public License as published by

# the Free Software Foundation, either version 3 of the License, or

# (at your option) any later version.

#

# SigProfilerMatrixGenerator is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of

# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

# GNU General Public License for more details.

#

# You should have received a copy of the GNU General Public License

# along with SigProfilerMatrixGenerator.  If not, see http://www.gnu.org/licenses/
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

import pickle
import os
import argparse
import re

def save_chrom_tsb_separate (genome, ref_dir):
	chromosome_path = ref_dir + '/references/chromosomes/tsb/' + genome + '/'
	chromosome_BED_path = ref_dir + '/references/chromosomes/tsb_BED/' + genome + '/'
	
	chromosomes = ['X', 'Y', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
				   '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

	tsb_ref = {0:['N','A'], 1:['N','C'], 2:['N','G'], 3:['N','T'],
			   4:['T','A'], 5:['T','C'], 6:['T','G'], 7:['T','T'],
			   8:['U','A'], 9:['U','C'], 10:['U','G'], 11:['U','T'],
			   12:['B','A'], 13:['B','C'], 14:['B','G'], 15:['B','T'],
			   16:['N','N'], 17:['T','N'], 18:['U','N'], 19:['B','N']}
			   
	tsb_ref_original = {'N':0, 'T':1, 'U':2, 'B':3}


	if genome == 'mm10' or genome == 'mm9':
		chromosomes = chromosomes[:21]

	if not os.path.exists(chromosome_BED_path):
		os.makedirs(chromosome_BED_path)

	for chrom in chromosomes:
		with open (chromosome_path+chrom+".txt", "rb") as f, open(chromosome_BED_path + chrom + "_BED_TSB.txt", 'w') as out:
			print("<CHROM>\t<START>\t<END>\t<TSB>", file=out)
			chrom_tsb = f.read()
			first_tsb = tsb_ref_original[tsb_ref[chrom_tsb[0]][0]]
			current_range = [0]
			for i in range (1, len(chrom_tsb), 1):
				if tsb_ref_original[tsb_ref[chrom_tsb[i]][0]] != first_tsb:
					current_range.append(i-1)
					print(chrom + "\t" + str(current_range[0]) + "\t" + str(current_range[1]) + "\t" + str(tsb_ref_original[tsb_ref[chrom_tsb[i-1]][0]]), file=out)
					first_tsb = tsb_ref_original[tsb_ref[chrom_tsb[i]][0]]
					current_range = [i]
				else:
					continue
		print("chromosome ", chrom, "done")


def main ():
	parser = argparse.ArgumentParser(description="Provide the necessary arguments to install the reference files.")
	parser.add_argument("-g", "--genome", nargs='?', help="Optional parameter instructs script to install the custom genome.")
	args = parser.parse_args()
	genome = args.genome

	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)

	save_chrom_tsb_separate(genome, ref_dir)


if __name__ == '__main__':
	main()