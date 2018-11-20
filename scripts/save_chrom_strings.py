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

import os
import re
import sys
import argparse


def save_chrom_strings (genome):
	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)
	chrom_fasta_path = ref_dir + '/references/chromosomes/fasta/' + genome + "/"
	chrom_string_path = ref_dir + '/references/chromosomes/chrom_string/' + genome + '/'
	chrom_fasta_files = os.listdir(chrom_fasta_path)

	if not os.path.exists(chrom_string_path):
		os.makedirs(chrom_string_path)

	for files in chrom_fasta_files:
		if not files.startswith('.'):
			file_name = files.split(".")
			try:
				chromosome = file_name[-2]
			except:
				continue
			if file_name[-4] == 'dna':
				with open(chrom_fasta_path + files) as chrom, open(chrom_string_path + chromosome + ".txt", 'w') as out:
					next(chrom)
					chromosome_final = ''
					for lines in chrom:
						line = lines.strip()
						chromosome_final += line


					print(chromosome_final, file=out)
					print("The string file for Chromosome " + chromosome + " has been created.")

def main():
	parser = argparse.ArgumentParser(description="Provide the necessary arguments to save the chromosomes as strings for later reference.")
	parser.add_argument("--genome", "-g",help="Provide a reference genome. (ex: GRCh37, GRCh38, mm10)")
	args=parser.parse_args()
	genome = args.genome

	save_chrom_strings(genome)


if __name__ == '__main__':
	main()