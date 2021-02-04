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
revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1','U':'T','T':'U','B':'B','N':'N', 'Q':'Q'}[B] for B in x][::-1])


def context_distribution (context_input, output_file, chromosome_path, chromosomes, tsb_ref, genome):
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
	dinuc_tsb = ['AA', 'AG', 'CC', 'GA',]
	dinuc_non_tsb = ['Q:AC', 'Q:AT', 'Q:CA', 'Q:CG', 'Q:GC', 'Q:TA']
	tsb_bias = ['T','U','B','N']

	# Set the context parameter based upon the user input
	if context_input == '96' or context_input == '192' or context_input == "384":
		context = 3
	elif context_input == '1536' or context_input == '3072' or context_input == '6144':
		context = 5
	elif context_input == 'DINUC' or context_input == 'DBS186' or context_input == "DBS":
		context = 2
	elif context_input == '6' or context_input == '24':
		context = 1
	else:
		print('Not a valid context')
		sys.exit()

	count = 0 
	probs = {}
	chromosome_lengths = []


	# Populate the dictionary if desired context is DINUC
	if context_input == 'DBS' or context_input == 'DINUC': 
		for dinuc in dinuc_types:
			probs[dinuc] = {}

	if context_input == 'DBS186':
		for dinuc in dinuc_tsb:
			for bias in tsb_bias:
				dinuc_t = bias + ":" + dinuc
				probs[dinuc_t] = {}
		for dinuc in dinuc_non_tsb:
			probs[dinuc] = {}


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
				if count == 1000000:
					print(i)
					count = 0
				# Skip the base if unknown 
				if "N" in nuc:
					pass

				else:
					if context_input != "DINUC" and context_input != "DBS186":
						# Only save the pyrimidine context (canonical)
						if base == 'A' or base == 'G':
							nuc = revcompl(nuc)


						# Adjust the nucleotide representaiton if TSB is desired
						if context_input == '192' or context_input == '3072' or context_input == '384' or context_input == '6144' or context_input == '24':
							bias = tsb_ref[chromosome[i+int(context/2)]][0]
							nuc = bias + ":" + nuc

						# Update the dictionary for the current nucleotide
						if nuc not in probs:
							probs[nuc] = {chrom:1}
						else:
							if chrom not in probs[nuc]:
								probs[nuc][chrom] = 1
							else:
								probs[nuc][chrom] += 1

					else:
						if context_input == 'DBS186':
							bias = tsb_ref[chromosome[i+int(context/2)]][0]
							if nuc not in dinuc_types:
								nuc = revcompl(nuc)
								bias = revbias(bias)
							if nuc not in dinuc_tsb:
								bias = 'Q'
							nuc = bias + ":" + nuc

						else:
							if nuc not in dinuc_types:
								nuc = revcompl(nuc)
						# Update the dictionary for the current nucleotide
						if chrom not in probs[nuc]:
							probs[nuc][chrom] = 1
						else:
							probs[nuc][chrom] += 1


						# if context_input == 'DBS186':
						# 	bias = tsb_ref[chromosome[i+int(context/2)]][0]
						# 	nuc = bias + ":" + nuc
						# # Populate the dictionary if desired context is DINUC
						# for dinuc in dinuc_types:
						# 	if dinuc not in probs:
						# 		probs[dinuc] = {}

						# # Update the dictionary for the current nucleotide
						# if nuc not in probs:
						# 	nuc = revcompl(nuc)
							
						# if chrom not in probs[nuc]:
						# 	probs[nuc][chrom] = 1
						# else:
						# 	probs[nuc][chrom] += 1


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
			try:
				print(probs[nuc][chromosomes[-1]]/nuc_sum, file=out)
			except:
				print(str(0), file=out)

	counts_file = os.path.dirname(output_file)	
	with open (counts_file + "/context_counts_" + genome + "_" + context_input + ".csv", 'w') as out:
		print(' ,', end='', file=out)
		for chrom in chromosomes[:-1]:
			print(chrom + ',', end='',file=out)
		print(chromosomes[-1],file=out)
		for nuc in probs.keys():
			nuc_sum = sum(probs[nuc].values())
			print (nuc + ',', end='', file=out)
			for chroms in chromosomes[:-1]:
				try:
					print(str(probs[nuc][chroms]) + ',', end='', file=out)
					out.flush()
				except:
					print(str(0) + ',', end='', file=out)
			try:
				print(str(probs[nuc][chroms]), file=out)
			except:
				print(str(0), file=out)
			# try:
			# 	print(probs[nuc][chromosomes[-1]]/nuc_sum, file=out)
			# except:
			# 	print(str(0),file=out)
	# chromosomes.remove("Y")
	# chromosomes.remove("MT")
	# chromosomes.remove("M")

	
	# Sort the file so that the nucleotides are in alphabetical order
	sort_command_1 = "sort -t ',' -k 1,1 "
	sort_command_2 = " -o "
	os.system (sort_command_1 + output_file + sort_command_2 + output_file)


def context_distribution_BED (context_input, output_file, chromosome_path, chromosomes, bed, bed_file, exome, exome_file, genome, ref_dir, tsb_ref, gender):
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
	dinuc_tsb = ['AA', 'AG', 'CC', 'GA',]
	dinuc_non_tsb = ['Q:AC', 'Q:AT', 'Q:CA', 'Q:CG', 'Q:GC', 'Q:TA']
	tsb_bias = ['T','U','B','N']

	# Set the context parameter based upon the user input
	if context_input == '96' or context_input == '192' or context_input == '384':
		context = 3
	elif context_input == '1536' or context_input == '3072' or context_input == '6144':
		context = 5
	elif context_input == 'DINUC' or context_input == 'DBS' or context_input == 'DBS186':
		context = 2
	elif context_input == '6' or context_input == '24':
		context = 1
	else:
		print('Not a valid context')
		sys.exit()
	count = 0 
	probs = {}
	chromosome_lengths = {}
	first_line = True
	chrom_length = 0
	if exome:
		# file_to_open = ref_dir + "/references/chromosomes/exome/" + genome + "/" + exome_file
		file_to_open = exome_file
	else:
		# file_to_open = ref_dir + "/references/vcf_files/BED/" + bed_file
		file_to_open = bed_file

	chromosomes_sort = chromosomes
	if "Y" not in chromosomes:
		chromosomes_sort.append("Y")
	if "MT" not in chromosomes:
		chromosomes_sort.append("MT")
	if "M" not in chromosomes:
		chromosomes_sort.append("M")

	# Populate the dictionary if desired context is DINUC
	if context_input == 'DBS' or context_input == 'DINUC': 
		for dinuc in dinuc_types:
			probs[dinuc] = {}

	if context_input == 'DBS186':
		for dinuc in dinuc_tsb:
			for bias in tsb_bias:
				dinuc_t = bias + ":" + dinuc
				probs[dinuc_t] = {}
		for dinuc in dinuc_non_tsb:
			probs[dinuc] = {}


	with open(file_to_open) as f:
		lines = [line.strip().split() for line in f]
	output = open(file_to_open, 'w')
	print('\t'.join(lines[0]), file=output)
	if len(lines) > 1:
		if len(lines[1][0]) > 2:
			for line in sorted(lines[1:], key = lambda x: (chromosomes_sort.index(x[0][3:]), int(x[1]), int(x[2]))):
				print('\t'.join(line), file=output)
		else:
			for line in sorted(lines[1:], key = lambda x: (chromosomes_sort.index(x[0]), int(x[1]), int(x[2]))):
				print('\t'.join(line), file=output)
	output.close()

	with open(file_to_open) as b_file:
		next(b_file)
		for lines in b_file:
			line = lines.strip().split()
			chrom = line[0]
			if len(chrom) > 1 and chrom[0:3].upper() == 'CHR':
				chrom = chrom[3:]
			if gender == 'female' and chrom == 'Y':
				continue
			start = int(line[1])
			end = int(line[2])

			if first_line:
				chrom_initial = chrom
				first_line = False 
				f = open (chromosome_path + chrom + ".txt", "rb")
				chromosome = f.read()

			if chrom == chrom_initial:
				chrom_length += end-start
				for i in range(start, end+1-context, 1):
					nuc = ''
					for l in range(i, i+context,1):
						nuc += tsb_ref[chromosome[l]][1]
					base = nuc[int(context/2)]

					# Skip the base if unknown 
					if "N" in nuc:
						pass

					else:
						if context_input != "DINUC" and context_input != 'DBS186':
							# Only save the pyrimidine context (canonical)
							if base == 'A' or base == 'G':
								nuc = revcompl(nuc)

							# Adjust the nucleotide representaiton if TSB is desired
							if context_input == '192' or context_input == '3072' or context_input == '384' or context_input == '6144' or context_input == '24':
								bias = tsb_ref[chromosome[i+int(context/2)]][0]
								nuc = bias + ":" + nuc

							# Update the dictionary for the current nucleotide
							if nuc not in probs:
								probs[nuc] = {chrom:1}
							else:
								if chrom not in probs[nuc]:
									probs[nuc][chrom] = 1
								else:
									probs[nuc][chrom] += 1

						else:
							if context_input == 'DBS186':
								bias = tsb_ref[chromosome[i+int(context/2)]][0]
								if nuc not in dinuc_types:
									nuc = revcompl(nuc)
									bias = revbias(bias)
								if nuc not in dinuc_tsb:
									bias = 'Q'
								nuc = bias + ":" + nuc

							else:
								if nuc not in dinuc_types:
									nuc = revcompl(nuc)
							# Update the dictionary for the current nucleotide
							if chrom not in probs[nuc]:
								probs[nuc][chrom] = 1
							else:
								probs[nuc][chrom] += 1
			else:
				f.close()
				print("          Chromosome ", chrom_initial, "done")
				chromosome_lengths[chrom_initial] = chrom_length
				chrom_length = end-start
				chrom_initial = chrom
				try:
					f = open (chromosome_path + chrom + ".txt", "rb")
					chromosome = f.read()
				except:
					continue

				for i in range(start, end+1-context, 1):
					nuc = ''
					for l in range(i,i+context,1):
						nuc += tsb_ref[chromosome[l]][1]
					base = nuc[int(context/2)]

					# Skip the base if unknown 
					if "N" in nuc:
						pass

					else:
						if context_input != "DINUC" and context_input != 'DBS186':
							# Only save the pyrimidine context (canonical)
							if base == 'A' or base == 'G':
								nuc = revcompl(nuc)

							# Adjust the nucleotide representaiton if TSB is desired
							if context_input == '192' or context_input == '3072' or context_input == '384' or context_input == '6144' or context_input == '24':
								bias = tsb_ref[chromosome[i+int(context/2)]][0]
								nuc = bias + ":" + nuc

							# Update the dictionary for the current nucleotide
							if nuc not in probs:
								probs[nuc] = {chrom:1}
							else:
								if chrom not in probs[nuc]:
									probs[nuc][chrom] = 1
								else:
									probs[nuc][chrom] += 1
						else:
							if context_input == 'DBS186':
								bias = tsb_ref[chromosome[i+int(context/2)]][0]
								if nuc not in dinuc_types:
									nuc = revcompl(nuc)
									bias = revbias(bias)
								if nuc not in dinuc_tsb:
									bias = 'Q'
								nuc = bias + ":" + nuc

							else:
								if nuc not in dinuc_types:
									nuc = revcompl(nuc)
							# Update the dictionary for the current nucleotide
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

	counts_file = os.path.dirname(output_file)	
	with open (counts_file + "/context_counts_" + genome + "_" + context_input + "_exome.csv", 'w') as out:
		print(' ,', end='', file=out)
		for chrom in chromosomes[:-1]:
			print(chrom + ',', end='',file=out)
		print(chromosomes[-1],file=out)
		for nuc in probs.keys():
			nuc_sum = sum(probs[nuc].values())
			print (nuc + ',', end='', file=out)
			for chroms in chromosomes[:-1]:
				try:
					print(str(probs[nuc][chroms]) + ',', end='', file=out)
					out.flush()
				except:
					print(str(0) + ',', end='', file=out)
			# try:
			# 	print(probs[nuc][chromosomes[-1]]/nuc_sum, file=out)
			# except:
			# 	print(str(0),file=out)
			try:
				print(str(probs[nuc][chroms]), file=out)
			except:
				print(str(0), file=out)
	# chromosomes.remove("Y")
	# chromosomes.remove("MT")
	# chromosomes.remove("M")

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
	if genome.upper() == 'RN6':
		chromosomes = chromosomes[:22]

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
		context_distribution(context, output_file, chromosome_path, chromosomes, tsb_ref, genome)
	
if __name__ == '__main__':
	main()

