#!/usr/bin/env python3
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

from __future__ import print_function
import os
import re
import sys
import argparse


def save_chrom_strings (genome):
	'''
	Saves the fasta files for each chromosome as single strings. This
	is required to generate the reference files.

	Parameters:
		genome  -> reference genome

	Returns:
		None

	Ouput:
		Text files that contain single strings for each chromosome. Saved
		into the "/references/chromosomes/chrom_string/[genome]/" path
	'''

	# Instantiates all of the relevant paths/creates them if not present
	ref_dir, tail = os.path.split(os.path.dirname(os.path.abspath(__file__)))
	chrom_fasta_path = ref_dir + '/references/chromosomes/fasta/' + genome + "/"
	chrom_string_path = ref_dir + '/references/chromosomes/chrom_string/' + genome + '/'
	chrom_fasta_files = os.listdir(chrom_fasta_path)
	if not os.path.exists(chrom_string_path):
		os.makedirs(chrom_string_path)

	# Iterates through all of the fasta files for the given genome and saves
	# each into a single string.
	first = True
	i = 1
	if ".DS_Store" in chrom_fasta_files:
		chrom_fasta_files.remove('.DS_Store')
	chrom_total = len(chrom_fasta_files)
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
					for lines in chrom:
						line = lines.strip()
						print(line, file=out, end='')

					if first:
						print("The string file has been created for Chromosome: " + chromosome + " (" + str(i) + "/" + str(chrom_total) + ")")
						first = False
					else:
						print("                                                 " + chromosome + " (" + str(i) + "/" + str(chrom_total) + ")")
		i += 1

def main():
	parser = argparse.ArgumentParser(description="Provide the necessary arguments to save the chromosomes as strings for later reference.")
	parser.add_argument("--genome", "-g",help="Provide a reference genome. (ex: GRCh37, GRCh38, mm10)")
	args=parser.parse_args()
	genome = args.genome

	save_chrom_strings(genome)


if __name__ == '__main__':
	main()