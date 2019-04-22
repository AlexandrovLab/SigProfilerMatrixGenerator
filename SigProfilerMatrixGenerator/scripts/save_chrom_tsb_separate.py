#!/usr/bin/env python3
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

from __future__ import print_function
import pickle
import os
import argparse
import re

def save_chrom_tsb_separate (genome, ref_dir):
	'''
	Saves the transcriptional strand bias information for a given genome 
	into a BED file.

	Parameters:
		 genome  -> reference genome
		ref_dir  -> directory from which the parent script is being called.

	Returns:
		  None
	
	Outputs:
		BED files that contain the TSB information and their associated ranges.
	'''

	# Instantiates all of the relevant path and reference varibales
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


	# Truncates the chromosome variable if a mouse genome is provided
	if genome == 'mm10' or genome == 'mm9':
		chromosomes = chromosomes[:21]

	# Creates the output path if it does not exist
	if not os.path.exists(chromosome_BED_path):
		os.makedirs(chromosome_BED_path)

	# Iterates through each chromosome, saving a new BED file for each
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