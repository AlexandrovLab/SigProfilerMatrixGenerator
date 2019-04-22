#!/usr/bin/env python3

# This source code file is a part of SigProfilerSimulator
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

from __future__ import print_function
import sys
import os
import re
import argparse

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])


def context_distribution (context_input, output_file, chromosome_path, chromosomes, tsb_ref):
	'''
	Creates a csv file for the distribution of nucleotides given a specific context. 
	This csv file needs to be created before simulating mutationalsigantures for the 
	given context.

	Requires: 
		Chromosomes saved in individual text files ('X.txt','Y.txt', '1.txt', etc)
		Transcriptional data saved in binary files for each chromosome. These files
			can be created using the script: "save_tsb_192.py"

	Parameters:
			  context_input  -> simulation context of interest (ex: 96, 192, 1536, 3072, DINUC, INDEL)
				output_file  -> file where the distribution for the given nucleotide context is saved (csv file)
			chromosome_path  -> path to the reference chromosomes
				chromosomes  -> list of chromosomes for the species of interest

	Returns:
		None

	Outputs:
		CSV file with each chromosome represented as a column and reach row 
		represented as a nucleotide. Under each chromosome for a given nucleotide
		is the proportion of the total length of that chromosome associated with that 
		nucleotide. 
	'''


	dinuc_types = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'GA', 'GC', 'TA']


	# Set the context parameter based upon the user input
	if context_input == '96' or context_input == '192':
		context = 3
	elif context_input == '1536' or context_input == '3072':
		context = 5
	elif context_input == 'DINUC':
		context = 2
	else:
		print('Not a valid context')
		sys.exit()

	count = 0 
	probs = {}
	chromosome_lengths = []

	# Iterate through each chromosome and open the associated file
	for chrom in chromosomes:
		with open (chromosome_path + chrom + ".txt", "rb") as f:

			chromosome = f.read().strip()
			chromosome_lengths.append(len(chromosome))
			print(chrom, len(chromosome))
			# If desired context is for TSB, open the transcriptional file, too
			# if context_input == '192' or context_input == '3072':
			# 	with open (chromosome_TSB_path + chrom + "_192.txt", 'rb') as tsb:
			# 		chromosome_tsb = tsb.read()

			# Iterate through the chromosome base by base
			for i in range (0, (len(chromosome)-context), 1):
				nuc = ''
				for l in range (i, i+context, 1):
					nuc += tsb_ref[chromosome[l]][1]
				base = nuc[int(context/2)]

				count += 1
				if count == 100000:
					print(i)
					count = 0
				# Skip the base if unknown 
				if "N" in nuc:
					pass

				else:
					if context_input != "DINUC":
						# Only save the pyrimidine context (canonical)
						if base == 'A' or base == 'G':
							nuc = revcompl(nuc)


						# Adjust the nucleotide representaiton if TSB is desired
						if context_input == '192' or context_input == '3072':
							bias = tsb_ref[chromosome[i+int(context/2)]][0]
							nuc = bias + ":" + nuc
							# if bias == 0:
							# 	nuc = 'N:' + nuc
							# elif bias == 1:
							# 	nuc = 'T:' + nuc
							# elif bias == 2:
							# 	nuc = 'U:' + nuc
							# else:
							# 	nuc = 'B:' + nuc 

						# Update the dictionary for the current nucleotide
						if nuc not in probs.keys():
							probs[nuc] = {chrom:1}
						else:
							if chrom not in probs[nuc].keys():
								probs[nuc][chrom] = 1
							else:
								probs[nuc][chrom] += 1

					else:
						# Populate the dictionary if desired context is DINUC
						for dinuc in dinuc_types:
							if dinuc not in probs.keys():
								probs[dinuc] = {}

						# Update the dictionary for the current nucleotide
						if nuc not in probs.keys():
							nuc = revcompl(nuc)
							
						if chrom not in probs[nuc]:
							probs[nuc][chrom] = 1
						else:
							probs[nuc][chrom] += 1


		print("chrom ", chrom, "done")
	print(probs, chromosome_lengths)
	# Write the resulting dictionary to the csv file
	with open (output_file, 'w') as out:
		print(' ,', end='', file=out)
		for chrom in chromosomes[:-1]:
			print(chrom + ',', end='',file=out)
		print(chromosomes[-1],file=out)
		for nuc in probs.keys():
			nuc_sum = sum(probs[nuc].values())
			print (nuc + ',', end='', file=out)
			for i in range (0, len(chromosomes[:-1]), 1):
				try:
					print(str(probs[nuc][chromosomes[i]]/nuc_sum) + ',', end='', file=out)
					out.flush()
				except:
					print(str(0) + ',', end='', file=out)
			print(probs[nuc][chromosomes[-1]]/nuc_sum, file=out)


	# Sort the file so that the nucleotides are in alphabetical order
	sort_command_1 = "sort -t ',' -k 1,1 "
	sort_command_2 = " -o "
	os.system (sort_command_1 + output_file + sort_command_2 + output_file)


def context_distribution_BED (context_input, output_file, chromosome_path, chromosomes, bed, bed_file, exome, exome_file, genome, ref_dir, tsb_ref):
	'''
	Creates a csv file for the distribution of nucleotides given a specific context and BED file. 
	This csv file needs to be created before simulating mutationalsigantures for the given 
	context.

	Requires: 
		Chromosomes saved in individual text files ('X.txt','Y.txt', '1.txt', etc)
		Transcriptional data saved in binary files for each chromosome. These files
			can be created using the script: "save_tsb_192.py"

	Parameters:
			  context_input  -> simulation context of interest (ex: 96, 192, 1536, 3072, DINUC, INDEL)
				output_file  -> file where the distribution for the given nucleotide context is saved (csv file)
			chromosome_path  -> path to the reference chromosomes
				chromosomes  -> list of chromosomes for the species of interest
						bed  -> flag that determines if the user has provided a BED file with specific ranges to simulate
				   bed_file  -> BED file that contains the ranges of interest

	Returns:
		None

	Outputs:
		CSV file with each chromosome represented as a column and reach row 
		represented as a nucleotide. Under each chromosome for a given nucleotide
		is the proportion of the total length of that chromosome associated with that 
		nucleotide. 
	'''


	dinuc_types = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'GA', 'GC', 'TA']

	# Set the context parameter based upon the user input
	if context_input == '96' or context_input == '192':
		context = 3
	elif context_input == '1536' or context_input == '3072':
		context = 5
	elif context_input == 'DINUC':
		context = 2
	else:
		print('Not a valid context')
		sys.exit()
	count = 0 
	probs = {}
	chromosome_lengths = {}
	first_line = True
	chrom_length = 0
	if exome:
		file_to_open = ref_dir + "/references/chromosomes/exome/" + genome + "/" + exome_file
	else:
		file_to_open = ref_dir + "/references/vcf_files/BED/" + bed_file
	with open(file_to_open) as b_file:
		next(b_file)
		for lines in b_file:
			line = lines.strip().split()
			chrom = line[0]
			if len(chrom) > 1 and chrom[0:3].upper() == 'CHR':
				chrom = chrom[3:]
			start = int(line[1])
			end = int(line[2])

			if first_line:
				chrom_initial = chrom
				first_line = False 
				f = open (chromosome_path + chrom + ".txt", "rb")
				chromosome = f.read()
				# if context_input == '192' or context_input == '3072':
				# 	tsb = open (chromosome_TSB_path + chrom + "_192.txt", 'rb')
				# 	chromosome_tsb = tsb.read()

			if chrom == chrom_initial:
				chrom_length += end-start
				for i in range(start, end+1-context, 1):
					nuc = ''
					for l in range(i, i+context,1):
						nuc += tsb_ref[chromosome[l]][1]
					base = nuc[int(context/2)]

					count += 1
					if count == 100000:
						#print(i)
						count = 0
					# Skip the base if unknown 
					if "N" in nuc:
						pass

					else:
						if context_input != "DINUC":
							# Only save the pyrimidine context (canonical)
							if base == 'A' or base == 'G':
								nuc = revcompl(nuc)

							# Adjust the nucleotide representaiton if TSB is desired
							if context_input == '192' or context_input == '3072':
								bias = tsb_ref[chromosome[i+int(context/2)]][0]
								nuc = bias + ":" + nuc
								# if bias == 0:
								# 	nuc = 'N:' + nuc
								# elif bias == 1:
								# 	nuc = 'T:' + nuc
								# elif bias == 2:
								# 	nuc = 'U:' + nuc
								# else:
								# 	nuc = 'B:' + nuc 

							# Update the dictionary for the current nucleotide
							if nuc not in probs.keys():
								probs[nuc] = {chrom:1}
							else:
								if chrom not in probs[nuc].keys():
									probs[nuc][chrom] = 1
								else:
									probs[nuc][chrom] += 1

						else:
							# Populate the dictionary if desired context is DINUC
							for dinuc in dinuc_types:
								if dinuc not in probs.keys():
									probs[dinuc] = {}

							# Update the dictionary for the current nucleotide
							if nuc not in probs.keys():
								nuc = revcompl(nuc)
								
							if chrom not in probs[nuc]:
								probs[nuc][chrom] = 1
							else:
								probs[nuc][chrom] += 1
			else:
				f.close()
				print("yes")
				print("Chromosome ", chrom_initial, "done")
				chromosome_lengths[chrom_initial] = chrom_length
				chrom_length = end-start
				chrom_initial = chrom
				print(chrom_initial)
				f = open (chromosome_path + chrom + ".txt", "rb")
				chromosome = f.read()
				#chromosome_lengths.append(len(chromosome))
				# if context_input == '192' or context_input == '3072':
				# 	tsb = open (chromosome_TSB_path + chrom + "_192.txt", 'rb')
				# 	chromosome_tsb = tsb.read()

				for i in range(start, end+1-context, 1):
					nuc = ''
					for l in range(i,i+context,1):
						nuc += tsb_ref[chromosome[l]][1]
					base = nuc[int(context/2)]

					count += 1
					if count == 100000:
						#print(i)
						count = 0
					# Skip the base if unknown 
					if "N" in nuc:
						pass

					else:
						if context_input != "DINUC":
							# Only save the pyrimidine context (canonical)
							if base == 'A' or base == 'G':
								nuc = revcompl(nuc)

							# Adjust the nucleotide representaiton if TSB is desired
							if context_input == '192' or context_input == '3072':
								bias = tsb_ref[chromosome[i+int(context/2)]][0]
								nuc = bias + ":" + nuc
									# if bias == 0:
									# 	nuc = 'N:' + nuc
									# elif bias == 1:
									# 	nuc = 'T:' + nuc
									# elif bias == 2:
									# 	nuc = 'U:' + nuc
									# else:
									# 	nuc = 'B:' + nuc 

							# Update the dictionary for the current nucleotide
							if nuc not in probs.keys():
								if chrom == '22':
									print(chrom)
								probs[nuc] = {chrom:1}
							else:
								if chrom not in probs[nuc].keys():
									if chrom == '22':
										print(chrom, "first")
									probs[nuc][chrom] = 1
								else:
									if chrom == '22':
										print(chrom, "multiple")
									probs[nuc][chrom] += 1

						else:
							# Populate the dictionary if desired context is DINUC
							for dinuc in dinuc_types:
								if dinuc not in probs.keys():
									probs[dinuc] = {}

							# Update the dictionary for the current nucleotide
							if nuc not in probs.keys():
								nuc = revcompl(nuc)
								
							if chrom not in probs[nuc]:
								probs[nuc][chrom] = 1
							else:
								probs[nuc][chrom] += 1

		chromosome_lengths[chrom_initial] = chrom_length

	# Write the resulting dictionary to the csv file
	with open (output_file, 'w') as out:
		print(' ,', end='', file=out)
		for chrom in chromosomes[:-1]:
			print(chrom + ',', end='',file=out)
		print(chromosomes[-1],file=out)
		for nuc in probs.keys():
			nuc_sum = sum(probs[nuc].values())
			print (nuc + ',', end='', file=out)
			for chroms in chromosomes[:-1]:
				try:
					print(str(probs[nuc][chroms]/nuc_sum) + ',', end='', file=out)
					out.flush()
				except:
					print(str(0) + ',', end='', file=out)
			try:
				print(probs[nuc][chromosomes[-1]]/nuc_sum, file=out)
			except:
				print(str(0),file=out)


	# Sort the file so that the nucleotides are in alphabetical order
	sort_command_1 = "sort -t ',' -k 1,1 "
	sort_command_2 = " -o "
	os.system (sort_command_1 + output_file + sort_command_2 + output_file)




def main():
	bed = False
	bed_file = None
	exome = False
	exome_file = None
	gender = 'male'

	parser = argparse.ArgumentParser(description="Provide the necessary arguments to save the nucleotide distributions for each chromosome.")
	parser.add_argument("--genome", "-g",help="Provide a reference genome. (ex: GRCh37, GRCh38, mm10)")
	parser.add_argument("--context", "-c",help="Whole genome context by default")
	parser.add_argument("-b", "--bed", nargs='?', help="Optional parameter instructs script to simulate on a given set of ranges (ex: exome). Whole genome context by default")
	parser.add_argument("-e", "--exome", nargs='?', help="Optional parameter instructs script to simulate only on exome). Whole genome context by default")
	parser.add_argument("-gD", "--gender", help="Optional parameter instructs script to create the context files based on female (two x chromosomes.", action='store_true')

	tsb_ref = {0:['N','A'], 1:['N','C'], 2:['N','G'], 3:['N','T'],
			   4:['T','A'], 5:['T','C'], 6:['T','G'], 7:['T','T'],
			   8:['U','A'], 9:['U','C'], 10:['U','G'], 11:['U','T'],
			   12:['B','A'], 13:['B','C'], 14:['B','G'], 15:['B','T'],
			   16:['N','N'], 17:['T','N'], 18:['U','N'], 19:['B','N']}

	args=parser.parse_args()
	genome = args.genome
	context = args.context
	if args.bed:
		bed = True
		bed_file = args.bed

	if args.exome:
		exome = True
		exome_file = args.exome

	chromosomes = ['X', 'Y', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
				   '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

	if args.gender:
		gender = 'female'
		chromosomes.remove('Y')


	if genome.upper() == 'MM10' or genome.upper() == 'MM9':
		chromosomes = chromosomes[:21]

	script_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', script_dir)
	#chromosome_path = ref_dir + "/references/chromosomes/chrom_string/" + genome + "/"
	#chromosome_path = "/Users/ebergstr/Desktop/test_bi/references/chromosomes/tsb/GRCh37/"
	chromosome_path = ref_dir + "/references/chromosomes/tsb/" + genome + "/"

	output_path = ref_dir + '/references/chromosomes/context_distributions/'
	if not os.path.exists(output_path):
		os.makedirs(output_path)

	if bed:
		output_file = ref_dir + '/references/chromosomes/context_distributions/context_distribution_' + genome + "_" + context + "_" + gender + '_BED.csv'
	else:
		if exome:
			output_file = ref_dir + '/references/chromosomes/context_distributions/context_distribution_' + genome + "_" + context + "_" + gender +'_exome.csv'
		else:
			output_file = ref_dir + '/references/chromosomes/context_distributions/context_distribution_' + genome + "_" + context + "_" + gender +'.csv'

	if bed or exome:
		context_distribution_BED(context, output_file, chromosome_path, chromosomes, bed, bed_file, exome, exome_file, genome, ref_dir, tsb_ref)
	else:
		context_distribution(context, output_file, chromosome_path, chromosomes, tsb_ref)
	
if __name__ == '__main__':
	main()

