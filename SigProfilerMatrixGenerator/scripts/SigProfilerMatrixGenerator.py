#!/usr/bin/env python3
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

import os
import sys
import re
import argparse
import itertools
import pandas as pd
from itertools import chain
import time
import datetime
from scipy import stats
import statsmodels.stats.multitest as sm
from SigProfilerMatrixGenerator.scripts import convert_input_to_simple_files as convertIn
import sigProfilerPlotting as sigPlt
import uuid
from collections import defaultdict
from collections import OrderedDict
import numpy as np
#np.set_printoptions(threshold=np.nan)

################# Functions and references ###############################################
def perm(n, seq):
	'''
	Generates a list of all available permutations of n-mers.

	Parameters:
			   n  -> length of the desired permutation string
			 seq  -> list of all possible string values

	Returns:
		  permus  -> list of all available permutations
	'''
	permus = []
	for p in itertools.product(seq, repeat=n):
		permus.append("".join(p))
	return(permus)

def reference_paths (genome):
	'''
	Returns the path to the reference genomes installed with SigProfilerMatrixGenerator

	Parameters:
		genome  -> reference genome

	Returns:
		chrom_path  -> path to the reference genome's chromosome files
	'''
	current_dir = os.path.realpath(__file__)
	ref_dir = re.sub('\/scripts/SigProfilerMatrixGenerator.py$', '', current_dir)
	chrom_path =ref_dir + '/references/chromosomes/tsb/' + genome + "/"

	return(chrom_path)

def BED_filtering (bed_file_path):
	'''
	Creates ranges from a bed file for generating the matrix.

	Parameters:
		 bed_file_path  -> path to the desired bed file

	Returns:
		  ranges_final  -> dictionary of all ranges for each chromosome.
	'''
	ranges = {}
	ranges_final = {}
	with open(bed_file_path) as f:
		#next(f)
		for lines in f:
			if lines[0] == "#" or lines[0] == '@':
				next(f)
			else:
				line = lines.strip().split()
				chrom = line[0]
				if len(chrom) > 2:
					chrom = chrom[3:]
				start = int(line[1])
				end = int(line[2])
				if chrom not in ranges.keys():
					ranges[chrom] = []
				ranges[chrom].append((start, end))

	for chroms in ranges.keys():
		ranges_final[chroms] = set(chain(*(range(start, end+1) for start, end in ranges[chroms])))

	return(ranges_final)

def gene_range (files_path, indel=False):
	'''
	Creates a dictionary of gene ranges and gene names across the 
	given reference genome. 

	Parameters:
							files_path  -> path to the transcript files
								 indel  -> flag that will construct the data structures for indels

	Returns:
						   gene_ranges  -> dictionary that contains the gene ranges on a chromosome basis.
						   gene_counts  -> dictionary that contains the number of mutations found for a given gene.
										   This value broken into a list of two integers ([transcribed, untranscribed])
							gene_names  -> dictionary that contains all of the gene names on a chromosome basis
			sample_mut_counts_per_gene  -> dictionary that contains all of the genes. It will store the number
										   of mutations associated with each gene per sample.
		sample_mut_counts_per_mut_type  -> dictionary that contains the total mutation count for each gene per mutation type
	'''
	gene_ranges = {}
	gene_counts = {}
	gene_names = {}
	sample_mut_counts_per_gene = {}
	sample_mut_counts_per_mut_type = {}

	for file in os.listdir(files_path):
		name = file.split("_")
		chrom = name[0]
		gene_ranges[chrom] = []
		gene_names[chrom] = []
		if file == '.DS_Store':
			continue
		else:
			with open(files_path + file) as f:
				next(f)
				for lines in f:
					line = lines.strip().split("\t")
					gene, start, end, strand, chrom = line[6], line[4], line[5], line[3], line[2]
					start, end = int(start), int(end)
					if gene not in gene_names[chrom]:
						gene_counts[gene] = OrderedDict()
						gene_ranges[chrom].append((start, end, strand))
						gene_names[chrom].append(gene)
						if indel:
							gene_counts[gene] = {'T':0, 'U':0, 'samples':[]}
						else:
							gene_counts[gene] = {'T:C>A':0, 'T:C>G':0,'T:C>T':0,'T:T>A':0,'T:T>C':0,'T:T>G':0,
												 'U:C>A':0, 'U:C>G':0,'U:C>T':0,'U:T>A':0,'U:T>C':0,'U:T>G':0,
												 'samples':[]}
						sample_mut_counts_per_gene[gene] = {}
						sample_mut_counts_per_mut_type[gene] = {}
					else:
						lst = list(gene_ranges[chrom][-1])
						if start < lst[0]:
							lst[0] = start
						if end > lst[1]:
							lst[1] = end
						gene_ranges[chrom][-1] = tuple(lst) 

	return(gene_ranges, gene_counts, gene_names, sample_mut_counts_per_gene, sample_mut_counts_per_mut_type)




def catalogue_generator_single (vcf_path, vcf_path_original, vcf_files, bed_file_path, chrom_path, project, output_matrix, context, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, transcript_path, tsb_stat, gs, log_file):
	'''
	Generates the mutational matrix for 96, 1536, 384, and 6144 context using a single vcf file with all samples of interest.

	Parameters:
					vcf_path  -> path to vcf file of interest
		   vcf_path_original  -> path to the original vcf files
				   vcf_files  -> actual vcf file
			   bed_file_path  -> path to the user-provided bed files
				  chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
								file name: '1.txt', '2.txt', etc.
					 project  -> unique name given to the set of samples (ex. 'BRCA') 
			   output_matrix  -> path where the final mutational matrix is stored
					 context  -> desired context (ex. 96, 1536, 384, 6144)
					   exome  -> flag that generates the matrix in only the exome 
					  genome  -> the desired reference genome 
				  ncbi_chrom  -> dictionary that allows for the converstion of ncbi chromosome names to standard format
								for the mm10 assembly.
				functionFlag  -> flag that is used when calling this function from an alternate script
						 bed  -> parameter used to filter the mutations on a user-provided BED file
				  bed_ranges  -> dictionary that contains all of the ranges for each chromosome dictated by the user's input BED file
				 chrom_based  -> flag that generates the matrices on a chromosome basis
						plot  -> flag that plots the matrices after they are generated
					tsb_stat  -> performs a transcriptional strand bias test for the 24, 384, and 6144 contexts. The output is
								 saved into the output/TSB directory

					 tsb_ref  -> dictionary that allows for switching between binary and biologically relevant strings
			 transcript_path  -> path to the transcript files
						  gs  -> flag that generates a file for the strand bias on a gene basis.
					log_file  -> path to the output log file

	Returns:
		If called as a function, returns a nested dictionary of panda data frames for each matrix

	Outputs:
		Outputs a mutational matrix for the desired context within [user-provided input path]/output/[mut_type]/

	'''
	out = open(log_file, 'a')

	# Small functions to provide reverse complements of TSB and sequence info:
	revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])
	revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1','U':'T','T':'U','B':'B','N':'N'}[B] for B in x][::-1])
	
	# Provides the sorting order for the TSB matrices
	bias_sort = {'T':0,'U':1,'N':3,'B':2}
	tsb = ['T','U','N','B']
	bases = ['A','C','G','T']
	
	# Pre-fills the mutation types variable
	size = 5
	mut_types_initial = perm(size, "ACGT")
	mut_types = []
	for tsbs in tsb:
		for mut in mut_types_initial:
			current_base = mut[int(size/2)]
			if current_base == 'C' or current_base == 'T':
				for base in bases:
					if base != current_base:
						mut_types.append(tsbs+":"+mut[0:int(size/2)] + "[" + current_base+">"+ base+"]"+mut[int(size/2)+1:])



	# Instantiates all relevant variables
	types = []
	samples = []
	mutation_dict = {}
	flag = True
	i = 0
	#file = vcf_files[0]
	sample_start = None
	gene_counts = {}
	sample_count = {}
	skipped_count = 0
	total_analyzed = 0
	sequence = ''

	if exome:
		exome_temp_file = "exome_temp.txt"
		exome_file = open(vcf_path + exome_temp_file, 'w')

	if bed:
		bed_temp_file = "bed_temp.txt"
		bed_file = open(vcf_path + bed_temp_file, 'w')


	# Opens the input vcf file
	for file in vcf_files:
		chrom_start = file.split("_")[0]
		chrom = chrom_start
		with open (vcf_path + file) as f, open(chrom_path + chrom_start + ".txt", "rb") as f2:
				range_index = 0
				chrom_string = f2.read()
				# first_line = f.readline()

				# line = first_line.strip().split('\t')
				# chrom = line[1]
				# chrom_start = chrom
				# f.seek(0)

				# # Opens each chromosome string path and bias string path
				# try:
				# 	with open(chrom_path + chrom_start + ".txt", "rb") as f2:
				# 		chrom_string = f2.read()

				# except:
				# 	print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...", file=out)
				# 	out.flush()

				# Creates a gene strand bias output file
				if gs:
					out = open(output_matrix + "gene_strand_bias_counts_SNV.txt", "w")
					out_hot = open(output_matrix + "gene_strand_bias_counts_hotspots_SNV.txt", "w")
					gene_ranges, gene_counts, gene_names, sample_mut_counts_per_gene, sample_mut_counts_per_mut_type = gene_range(transcript_path)
				
				# Skips any header lines
				for lines in f:

					# Saves all relevant data from the file
					try:
						line = lines.strip().split('\t')
						# sample = line[1]
						# chrom = line[5]
						# # if chrom in ncbi_chrom.keys():
						# # 	chrom = ncbi_chrom[chrom]

						# start = int(line[6])
						# ref = line[8][0].upper()
						# mut = line[9][0].upper()

						sample = line[0]
						#chrom = line[1]
						# if chrom in ncbi_chrom.keys():
						# 	chrom = ncbi_chrom[chrom]

						start = int(line[2])
						ref = line[3][0].upper()
						mut = line[4][0].upper()

						# Saves the sample name for later reference
						# if sample not in samples:
						# 	samples.append(sample)
						# 	mutation_dict[sample] = {}
						# 	sample_count[sample] = [0,0]

						# if sample not in mutation_dict.keys():
						# 	mutation_dict[sample] = {}


						if sample not in mutation_dict:
							mutation_dict[sample] = {}
							samples.append(sample)
							mutation_dict[sample] = {}
							sample_count[sample] = [0,0]

						# # Opens the next chromosome data once a new chromosome is reached
						# if chrom != chrom_start:
						# 	range_index = 0
						# 	print("Chromosome " + chrom_start + " done",file=out)
						# 	out.flush()
						# 	if chrom_based:
						# 		matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start, tsb)
						# 		mutation_dict = {}
						# 		mutation_dict[sample] = {}
						# 	chrom_start = chrom
						# 	try:
						# 		with open(chrom_path + chrom_start + ".txt", "rb") as f:
						# 			chrom_string = f.read()
						# 		i += 1
						# 	except:
						# 		print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...", file=out)
						# 		out.flush()
						# 		continue

						# Pulls out the relevant sequence depending on the context
						try:
							sequence = ''.join([tsb_ref[chrom_string[start-3]][1],tsb_ref[chrom_string[start-2]][1], tsb_ref[chrom_string[start-1]][1], tsb_ref[chrom_string[start]][1], tsb_ref[chrom_string[start+1]][1]])
							# for i in range (start-3, start+2, 1):
							# 	sequence += tsb_ref[chrom_string[i]][1]
						except:
							print("The position is out of range. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
							out.flush()
							skipped_count += 1
							continue						

						bias = tsb_ref[chrom_string[start]][0]
						char = sequence[int(len(sequence)/2)]

						# Prints the sequence and position if the pulled sequence doesn't match
						# the variant from the file
						if char != ref:# and revcompl(char) != ref:
							print("The reference base does not match the reference genome. Skipping this mutation: "+ chrom + " "+ str(start) + " "+ ref + " "+ mut,file=out)
							out.flush()
							skipped_count += 1
							continue				
						
						# Saves the sequence/mutation type if it matched the reference/reverse strand  
						else:
							strand = '1'
							#bias_before = bias
							#ref_before = ref
							if ref == 'A' or ref == 'G':
								strand = '-1'
								bias = revbias(bias)
								ref = revcompl(ref)
								mut = revcompl(mut)
								sequence = revcompl(sequence)



							# Performs the gene strand bias test if desired
							if gs:
								for ranges in gene_ranges[chrom_start][range_index:]:
									dict_key = ref + ">" + mut
									if start < ranges[0]:
										break
									if ranges[0] <= start <= ranges[1]:
										gene_index = gene_ranges[chrom_start].index(ranges)
										gene = gene_names[chrom_start][gene_index]
										if int(strand) + int(ranges[2]) == 0:
											dict_key = 'T:' + dict_key
											gene_counts[gene][dict_key] += 1
											if sample not in gene_counts[gene]['samples']:
												gene_counts[gene]['samples'].append(sample)
												sample_mut_counts_per_gene[gene][sample] = 1
												sample_mut_counts_per_mut_type[gene][sample] = {'T:C>A':0, 'T:C>G':0,'T:C>T':0,'T:T>A':0,'T:T>C':0,'T:T>G':0,
												 'U:C>A':0, 'U:C>G':0,'U:C>T':0,'U:T>A':0,'U:T>C':0,'U:T>G':0}
												sample_mut_counts_per_mut_type[gene][sample][dict_key] += 1

											else:
												sample_mut_counts_per_gene[gene][sample] += 1
												sample_mut_counts_per_mut_type[gene][sample][dict_key] += 1
										elif strand == ranges[2]:
											dict_key = 'U:' + dict_key
											gene_counts[gene][dict_key] += 1
											if sample not in gene_counts[gene]['samples']:
												gene_counts[gene]['samples'].append(sample)
												sample_mut_counts_per_gene[gene][sample] = 1
												sample_mut_counts_per_mut_type[gene][sample] = {'T:C>A':0, 'T:C>G':0,'T:C>T':0,'T:T>A':0,'T:T>C':0,'T:T>G':0,
												 'U:C>A':0, 'U:C>G':0,'U:C>T':0,'U:T>A':0,'U:T>C':0,'U:T>G':0}
												sample_mut_counts_per_mut_type[gene][sample][dict_key] += 1
											else:
												sample_mut_counts_per_gene[gene][sample] += 1
												sample_mut_counts_per_mut_type[gene][sample][dict_key] += 1
									

							# Saves the mutation key for the current variant
							#mut_key = bias + ":" + sequence[0:int(len(sequence)/2)] + '[' + ref + '>' + mut + ']' + sequence[int(len(sequence)/2+1):]
							mut_key = ''.join([bias,":",sequence[0:int(len(sequence)/2)],'[',ref,'>',mut,']',sequence[int(len(sequence)/2+1):]])
							if mut_key not in mutation_dict[sample]:
								mutation_dict[sample][mut_key] = 1
							else:
								mutation_dict[sample][mut_key] += 1
							total_analyzed += 1


							# If exome is specified, it will write the variant to a temporary exome file.
							if exome:
								exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + mut_key + "\t" + ref + "\t" + mut + "\n")
							if bed:
								bed_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + mut_key + "\t" + ref + "\t" + mut + "\n")

					except:
						print("There appears to be an error in this line. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut,file=out)
						out.flush()
						skipped_count += 1

				# Once all variants are accounted for, complete the gene strand bias test/output to the final file. 
				if gs:
					pvals = []
					qvals = []
					pvals_hot = []
					qvals_hot = []
					hotspots = {}
					gene_bias = []
					gene_bias_hotspots = []
					for gene in gene_counts:
						if gene not in gene_bias:
							gene_bias.append(gene)
						total_count = sum(sample_mut_counts_per_gene[gene].values())
						for sample in sample_mut_counts_per_gene[gene]:
							mut_count = sample_mut_counts_per_gene[gene][sample]
							if mut_count > 3: # and mut_count/total_count > 0.5:
								if gene not in gene_bias_hotspots:
									gene_bias_hotspots.append(gene)
								if gene not in hotspots:
									hotspots[gene] = {}
									for mut, count in sample_mut_counts_per_mut_type[gene][sample].items():
										hotspots[gene][mut] = count
									hotspots[gene]['samples'] = [sample]
									for mut, count in sample_mut_counts_per_mut_type[gene][sample].items():
										gene_counts[gene][mut] -= count
									gene_counts[gene]['samples'].remove(sample)
								else:
									for mut, count in sample_mut_counts_per_mut_type[gene][sample].items():
										hotspots[gene][mut] += count
										gene_counts[gene][mut] -= count
									gene_counts[gene]['samples'].remove(sample)
									hotspots[gene]['samples'].append(sample)

						sum_tran = 0
						sum_untran = 0
						for mut, counts in gene_counts[gene].items():
							if mut[0] == 'T':
								sum_tran += counts
							elif mut[0] == 'U':
								sum_untran += counts
						pvals.append(stats.binom_test([sum_tran, sum_untran]))

						sum_tran_hot = 0
						sum_untran_hot = 0
						if gene in hotspots:
							for mut, counts in hotspots[gene].items():
								if mut[0] == 'T':
									sum_tran_hot += counts
								elif mut[0] == 'U':
									sum_untran_hot += counts
						pvals_hot.append(stats.binom_test([sum_tran_hot, sum_untran_hot]))

					qvals = sm.fdrcorrection(pvals)[1]
					qvals_hot = sm.fdrcorrection(pvals_hot)[1]
					ind = pvals.index('BMP7')
					ind2 = pvals_hot.index('BMP7')

					gene_ind = 0
					for gene in gene_bias:
						gene_counts[gene]['samples'] = len(gene_counts[gene]['samples'])
						print(gene, end='',file=out, flush=False)
						sum_tran = 0
						sum_untran = 0
						for mut, counts in gene_counts[gene].items():
							if mut[0] == 'T':
								sum_tran += counts
							elif mut[0] == 'U':
								sum_untran += counts
							print("\t" + str(counts), end='', file=out, flush=False)
						print("\t" + str(sum_tran) + "\t" + str(sum_untran) + "\t" + str(qvals[gene_ind]), flush=False, file=out)
						gene_ind += 1
					out.close()
					with open(output_matrix + "gene_strand_bias_counts_SNV.txt") as f2:
						lines = [line.strip().split() for line in f2]
					output = open(output_matrix + "gene_strand_bias_counts_SNV.txt", 'w')
					print('GENE\tT:C>A\tT:C>G\tT:C>T\tT:T>A\tT:T>C\tT:T>G\tU:C>A\tU:C>G\tU:C>T\tU:T>A\tU:T>C\tU:T>G\tSampleCount\tTranscribed_total\tUntranscribedTotal\tq_value', file=output)
					for line in sorted(lines, key = lambda x: (float(x[-1])), reverse=False):
						print('\t'.join(line), file=output)
					output.close()

					# Gene strand bias test for hot spot samples.
					gene_ind = 0
					for gene in gene_bias_hotspots:
						hotspots[gene]['samples'] = len(hotspots[gene]['samples'])
						print(gene, end='',file=out_hot, flush=False)
						sum_tran_hot = 0
						sum_untran_hot = 0
						for mut, counts in hotspots[gene].items():
							if mut[0] == 'T':
								sum_tran_hot += counts
							elif mut[0] == 'U':
								sum_untran_hot += counts
							print("\t" + str(counts), end='', file=out_hot, flush=False)
						print("\t" + str(sum_tran_hot) + "\t" + str(sum_untran_hot) + "\t" + str(qvals_hot[gene_ind]), flush=False, file=out_hot)
						gene_ind += 1
					out_hot.close()
					with open(output_matrix + "gene_strand_bias_counts_hotspots_SNV.txt") as f2:
						lines = [line.strip().split() for line in f2]
					output = open(output_matrix + "gene_strand_bias_counts_hotspots_SNV.txt", 'w')
					print('GENE\tT:C>A\tT:C>G\tT:C>T\tT:T>A\tT:T>C\tT:T>G\tU:C>A\tU:C>G\tU:C>T\tU:T>A\tU:T>C\tU:T>G\tSampleCount\tTranscribed_total\tUntranscribedTotal\tq_value', file=output)
					for line in sorted(lines, key = lambda x: (float(x[-1])), reverse=False):
						print('\t'.join(line), file=output)
					output.close()

		print("Chromosome " + chrom_start + " done",file=out)
		out.flush()
		if chrom_based:
			matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start, tsb)
			mutation_dict = {}
			mutation_dict[sample] = {}


	print("Chromosome " + chrom_start + " done", file=out)
	out.flush()
	# Generate the matrix for the final chromosome if specified by user.
	if chrom_based:
		matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start, tsb_stat)
		mutation_dict = {}

	# Organizes the variants based upon the exome
	if exome:
		exome_file.close()
		sort_file = exome_temp_file
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()
		mutation_dict, samples2 = exome_check(genome, vcf_path + exome_temp_file, output_matrix, project)

	# Organizes the variants based upon the user-provided bed file
	if bed:
		bed_file.close()
		sort_file = bed_temp_file
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()
		mutation_dict, samples2 = panel_check(genome, vcf_path + bed_temp_file, output_matrix, bed_file_path, project)

	# Calls the function to generate the final matrix
	if functionFlag:
		chrom_start = None
		chrom_start=None
		matrices = matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start, functionFlag, plot, tsb_stat)
		out.close()

		return(matrices, skipped_count, total_analyzed, len(samples))
	else:
		if not chrom_based:
			chrom_start=None
			matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start, functionFlag, plot, tsb_stat)
	out.close()

def catalogue_generator_DINUC_single (vcf_path, vcf_path_original, vcf_files, bed_file_path, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, log_file):
	'''
	Generates the mutational matrix for the dinucleotide context.

	Parameters:
						vcf_path  -> path to vcf file of interest
		   vcf_path_original  -> path to the original vcf files
				   vcf_files  -> actual vcf file
			   bed_file_path  -> path to the user-provided bed files
				  chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
								file name: '1.txt', '2.txt', etc.
					 project  -> unique name given to the set of samples (ex. 'BRCA') 
			   output_matrix  -> path where the final mutational matrix is stored
					   exome  -> flag that generates the matrix in only the exome 
					  genome  -> the desired reference genome 
				  ncbi_chrom  -> dictionary that allows for the converstion of ncbi chromosome names to standard format
								for the mm10 assembly.
				functionFlag  -> flag that is used when calling this function from an alternate script
						 bed  -> parameter used to filter the mutations on a user-provided BED file
				  bed_ranges  -> dictionary that contains all of the ranges for each chromosome dictated by the user's input BED file
				 chrom_based  -> flag that generates the matrices on a chromosome basis
						plot  -> flag that plots the matrices after they are generated
					 tsb_ref  -> dictionary that allows for switching between binary and biologically relevant strings
						  gs  -> flag that generates a file for the strand bias on a gene basis.
					log_file  -> path to the output log file
	
	Returns:
		If called as a function, returns a nested dictionary of panda data frames for each matrix

	Outputs:
		Outputs a mutational matrix for the desired context within [user-provided input path]/output/[mut_type]/

	'''

	out = open(log_file, 'a')
	# Small functions to provide reverse complements of TSB and sequence info:
	revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','[':'[',']':']','>':'>'}[B] for B in x][::-1])
	revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1','U':'T','T':'U','B':'B','N':'N'}[B] for B in x][::-1])
	bases = ['A','C','G','T']
	bias_sort = {'T':0,'U':1,'N':3,'B':2}
	tsb = ['T','U','N','B']


	# Instantiates the necessary variables/data structures
	dinucs = {}
	dinucs_context = {}
	dinucs_context_tsb = {}
	dinucs_tsb = {}
	samples = []

	mutation_types = ['AC>CA','AC>CG','AC>CT','AC>GA','AC>GG','AC>GT','AC>TA','AC>TG','AC>TT',
					  'AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA','CC>AA','CC>AG','CC>AT',
					  'CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT','CG>AT','CG>GC','CG>GT',
					  'CG>TA','CG>TC','CG>TT','CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG',
					  'CT>TA','CT>TC','CT>TG','GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA',
					  'TA>AT','TA>CG','TA>CT','TA>GC','TA>GG','TA>GT','TC>AA','TC>AG','TC>AT',
					  'TC>CA','TC>CG','TC>CT','TC>GA','TC>GG','TC>GT','TG>AA','TG>AC','TG>AT',
					  'TG>CA','TG>CC','TG>CT','TG>GA','TG>GC','TG>GT','TT>AA','TT>AC','TT>AG',
					  'TT>CA','TT>CC','TT>CG','TT>GA','TT>GC','TT>GG']



	# Organizes all of the mutation types for DINUCs
	mutation_types_context = []
	mutation_types_tsb_context = []
	mutation_types_tsb = []
	for base in bases:
		for mut in mutation_types:
			for base2 in bases:
				mutation_types_context.append(base + "[" + mut + "]" + base2)
				for base3 in tsb:
					mutation_types_tsb_context.append(''.join([base3,":",base,"[",mut,"]",base2]))
					mut_tsb = base3 + ":" + mut
					if mut_tsb not in mutation_types_tsb:
						mutation_types_tsb.append(base3 + ":" + mut)

	# Instantiates the required variables and paths
	types = []
	samples = []
	mutation_dict = {}
	flag = True
	i = 0
	#file = vcf_files[0]
	sample_start = None
	skipped_count = 0
	total_analyzed = 0
	if exome:
		exome_temp_file = "exome_temp.txt"
		exome_temp_file_context = "exome_temp_context.txt"
		exome_temp_file_tsb = "exome_temp_tsb.txt"
		exome_temp_file_context_tsb = "exome_temp_context_tsb.txt"

		exome_file = open(vcf_path + exome_temp_file, 'w')
		exome_file_context = open(vcf_path + exome_temp_file_context, 'w')
		exome_file_tsb = open(vcf_path + exome_temp_file_tsb, 'w')
		exome_file_context_tsb = open(vcf_path + exome_temp_file_context_tsb, 'w')

	if bed:
		bed_temp_file = "bed_temp.txt"
		bed_temp_file_context = "bed_temp_context.txt"
		bed_temp_file_tsb = "bed_temp_tsb.txt"
		bed_temp_file_context_tsb = "bed_temp_context_tsb.txt"

		bed_file = open(vcf_path + bed_temp_file, 'w')
		bed_file_context = open(vcf_path + bed_temp_file_context, 'w')
		bed_file_tsb = open(vcf_path + bed_temp_file_tsb, 'w')
		bed_file_context_tsb = open(vcf_path + bed_temp_file_context_tsb, 'w')

	# Iterates throught the provided vcf files
	for file in vcf_files:
		chrom_start = file.split("_")[0]
		chrom = chrom_start
		with open (vcf_path + file) as data, open(chrom_path + chrom_start + ".txt", "rb") as f2:
				range_index = 0
				chrom_string = f2.read()
		#with open (vcf_path + file) as data:
				# prev_line = None
				initial_line = data.readline()
				initial_line_data = initial_line.strip().split('\t')
				# previous_sample = initial_line_data[1]
				# previous_chrom = initial_line_data[5] 
				# previous_start = int(initial_line_data[6])
				# previous_ref = initial_line_data[8]
				# previous_mut = initial_line_data[9]

				previous_sample = initial_line_data[0]
				# previous_chrom = initial_line_data[1] 
				previous_start = int(initial_line_data[2])
				previous_ref = initial_line_data[3]
				previous_mut = initial_line_data[4]

				# opens the first chromosome file
				sample_start = previous_sample
				# chrom_start = previous_chrom
				# try:
				# 	with open(chrom_path + chrom_start + ".txt", "rb") as f:
				# 		chrom_string = f.read()
				# except:
				# 	print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...", file=out)

				# flag = False
						
				
				for lines in data:

					# Saves the relevant data from each line
					try:
						line = lines.strip().split('\t')
						# sample = line[1]
						# chrom = line[5]
						# # if chrom in ncbi_chrom.keys():
						# # 	chrom = ncbi_chrom[chrom]

						# start = int(line[6])
						# ref = line[8][0].upper()
						# mut = line[9][0].upper()
						sample = line[0]
						#chrom = line[1]
						# if chrom in ncbi_chrom.keys():
						# 	chrom = ncbi_chrom[chrom]

						start = int(line[2])
						ref = line[3][0].upper()
						mut = line[4][0].upper()

						# # Exception handling for incorrect input format
						# if ref not in 'ACGT-':
						# 	print("The ref base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						# 	skipped_count += 1
						# 	continue

						# if mut not in 'ACGT-':
						# 	print("The mutation base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						# 	skipped_count += 1
						# 	continue

						# if ref == mut:
						# 	print("The ref base appears to match the mutated base. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						# 	skipped_count += 1
						# 	continue

						# if line == prev_line:
						# 	print("There appears to be a duplicate single base substitution. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						# 	skipped_count += 1
						# 	continue

						# prev_line = line
						if sample not in samples:
							samples.append(sample)

						# # opens the first chromosome file
						# if flag:
						# 	sample_start = sample
						# 	chrom_start = chrom
						# 	try:
						# 		with open(chrom_path + chrom_start + ".txt", "rb") as f:
						# 			chrom_string = f.read()
						# 	except:
						# 		print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...", file=out)
						# 		continue

						# 	flag = False

						# Breaks the loop once a new chromosome is reached
						# if chrom != chrom_start:
						# 	print("Chromosome " + chrom_start + " done", file=out)
						# 	if chrom_based:
						# 		all_dinucs = {'78':dinucs, '312':dinucs_tsb, '1248':dinucs_context, '4992':dinucs_context_tsb}
						# 		all_mut_types = {'78':mutation_types, '312':mutation_types_tsb, '1248':mutation_types_context, '4992':mutation_types_tsb_context}
						# 		matrix_generator_DINUC (output_matrix, samples, bias_sort, all_dinucs, all_mut_types, dinucs, project, exome, bed, chrom_start, plot)
						# 		dinucs = {}
						# 	chrom_start = chrom
						# 	try:
						# 		with open(chrom_path + chrom_start + ".txt", "rb") as f:
						# 			chrom_string = f.read()
						# 	except:
						# 		print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...", file=out)
						# 		continue


						# Grabs the slice of interest and saves the relevant DINUC with context and TSB information
						else:
							if start == previous_start + 1:
								dinuc = previous_ref + ref + ">" + previous_mut + mut

								try:
									dinuc_seq = "".join([tsb_ref[chrom_string[start-3]][1],"[",dinuc,"]",tsb_ref[chrom_string[start]][1]])
									bias = tsb_ref[chrom_string[start-1]][0]
								except:
									print("The position is out of range. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
									skipped_count += 1
									continue
								dinuc_seq_tsb = bias + ":" + dinuc_seq
								dinuc_tsb = bias + ":" + dinuc

								if sample not in dinucs.keys():
									dinucs[sample] = {}
									dinucs_context[sample] = {}
									dinucs_context_tsb[sample] = {}
									dinucs_tsb[sample] = {}

									for dinucl in mutation_types:
										dinucs[sample][dinucl]=0
									for dinucl in mutation_types_context:
										dinucs_context[sample][dinucl]=0
									for dinucl in  mutation_types_tsb_context:
										dinucs_context_tsb[sample][dinucl]=0
									for dinucl in mutation_types_tsb:
										dinucs_tsb[sample][dinucl]=0

								# Saves the respective DINUC variable into the appropriate dictionary 
								if dinuc in mutation_types:
									dinucs[sample][dinuc] += 1 
								else:
									dinuc = revcompl(previous_ref + ref) + '>' + revcompl(previous_mut + mut)
									dinucs[sample][dinuc] += 1

								if dinuc_seq in mutation_types_context:
									dinucs_context[sample][dinuc_seq] += 1
								else:
									dinuc_seq = "".join([revcompl(dinuc_seq[0]),"[",revcompl(previous_ref + ref),'>',revcompl(previous_mut + mut),"]",revcompl(dinuc_seq[8])])
									dinucs_context[sample][dinuc_seq] += 1

								if dinuc_seq_tsb in mutation_types_tsb_context:
									dinucs_context_tsb[sample][dinuc_seq_tsb] += 1
								else:
									dinuc_seq_tsb = "".join([revbias(dinuc_seq_tsb[0]),":",revcompl(dinuc_seq_tsb[-1]),"[",revcompl(dinuc_seq_tsb[4:6]),">",revcompl(dinuc_seq_tsb[7:9]),"]",revcompl(dinuc_seq_tsb[2])])
									dinucs_context_tsb[sample][dinuc_seq_tsb] += 1

								if dinuc_tsb in mutation_types_tsb:
									dinucs_tsb[sample][dinuc_tsb] += 1
								else:
									dinuc_tsb = "".join([revbias(dinuc_tsb[0]),":",revcompl(dinuc_tsb[2:4]),">",revcompl(dinuc_tsb[5:])])
									dinucs_tsb[sample][dinuc_tsb] += 1

								# Saves the DINUC into temporary files for exome sorting
								if exome:
									exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc + "\t" + ref + "\t" + mut + "\n")
									exome_file_context.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc_seq + "\t" + ref + "\t" + mut + "\n")
									exome_file_tsb.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc_tsb + "\t" + ref + "\t" + mut + "\n")
									exome_file_context_tsb.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc_seq_tsb + "\t" + ref + "\t" + mut + "\n")

								# Saves the DINUC into temporary files for region sorting
								if bed:
									bed_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc + "\t" + ref + "\t" + mut + "\n")
									bed_file_context.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc_seq + "\t" + ref + "\t" + mut + "\n")
									bed_file_tsb.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc_tsb + "\t" + ref + "\t" + mut + "\n")
									bed_file_context_tsb.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc_seq_tsb + "\t" + ref + "\t" + mut + "\n")

								total_analyzed += 1

						previous_sample = sample
						#previous_chrom = chrom 
						previous_start = start
						previous_ref = ref
						previous_mut = mut

					except:
						print("There appears to be an error in this line. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						skipped_count += 1

		print("Chromosome " + chrom_start + " done", file=out)
		if chrom_based:
			all_dinucs = {'78':dinucs, '312':dinucs_tsb, '1248':dinucs_context, '4992':dinucs_context_tsb}
			all_mut_types = {'78':mutation_types, '312':mutation_types_tsb, '1248':mutation_types_context, '4992':mutation_types_tsb_context}
			matrix_generator_DINUC (output_matrix, samples, bias_sort, all_dinucs, all_mut_types, dinucs, project, exome, bed, chrom_start, plot)
			dinucs = {}


	# Organizes the required dictionaries for the final matrix generation.
	all_dinucs = {'78':dinucs, '312':dinucs_tsb, '1248':dinucs_context, '4992':dinucs_context_tsb}
	all_mut_types = {'78':mutation_types, '312':mutation_types_tsb, '1248':mutation_types_context, '4992':mutation_types_tsb_context}
	print("Chromosome " + chrom_start + " done", file=out)
	if chrom_based:
		matrix_generator_DINUC (output_matrix, samples, bias_sort, all_dinucs, all_mut_types, dinucs, project, exome, bed, chrom_start, plot)
		dinucs = {}

	# Sorts all of the temporary DINUC files used for the exome filtering
	if exome:
		exome_file.close()
		sort_file = exome_temp_file
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs['78'], samples2 = exome_check(genome, vcf_path + exome_temp_file, output_matrix, project)

		exome_file_context.close()
		sort_file = exome_temp_file_context
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs['1248'], samples2 = exome_check(genome, vcf_path + exome_temp_file_context, output_matrix, project)

		exome_file_tsb.close()
		sort_file = exome_temp_file_tsb
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs['312'], samples2 = exome_check(genome, vcf_path + exome_temp_file_tsb, output_matrix, project)

		exome_file_context_tsb.close()
		sort_file = exome_temp_file_context_tsb
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs['4992'], samples2 = exome_check(genome, vcf_path + exome_temp_file_context_tsb, output_matrix, project)

	# Sorts all of the temporary DINUC files used for the region filtering	
	if bed:
		bed_file.close()
		sort_file = bed_temp_file
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs['78'], samples2 = panel_check(genome, vcf_path + bed_temp_file, output_matrix, bed_file_path, project)

		bed_file_context.close()
		sort_file = bed_temp_file_context
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs['1248'], samples2 = panel_check(genome, vcf_path + bed_temp_file_context, output_matrix, bed_file_path, project)

		bed_file_tsb.close()
		sort_file = bed_temp_file_tsb
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs['312'], samples2 = panel_check(genome, vcf_path + bed_temp_file_tsb, output_matrix, bed_file_path, project)

		bed_file_context_tsb.close()
		sort_file = bed_temp_file_context_tsb
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs['4992'], samples2 = panel_check(genome, vcf_path + bed_temp_file_context_tsb, output_matrix, bed_file_path, project)

	# Calls the function to generate the final mutational matrix
	if not chrom_based:
		chrom_start=None
		matrix_generator_DINUC (output_matrix, samples, bias_sort, all_dinucs, all_mut_types, dinucs, project, exome, bed, chrom_start, plot)
		if functionFlag:
			return(pd.DataFrame.from_dict(dinucs), skipped_count, total_analyzed, len(samples))
		

		
	

def catalogue_generator_INDEL_single (vcf_path, vcf_path_original, vcf_files, bed_file_path, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, limited_indel, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, transcript_path, gs, log_file):
	'''
	Generates the mutational matrix for the INDEL context.

	Parameters:
					vcf_path  -> path to vcf file of interest
		   vcf_path_original  -> path to the original vcf files
				   vcf_files  -> actual vcf file
			   bed_file_path  -> path to the user-provided bed files
				  chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
								file name: '1.txt', '2.txt', etc.
					 project  -> unique name given to the set of samples (ex. 'BRCA') 
			   output_matrix  -> path where the final mutational matrix is stored
					   exome  -> flag that generates the matrix in only the exome 
					  genome  -> the desired reference genome 
				  ncbi_chrom  -> dictionary that allows for the converstion of ncbi chromosome names to standard format
								for the mm10 assembly.
			   limited_indel  -> flag that creates the matrix based on limited indels
				functionFlag  -> flag that is used when calling this function from an alternate script
						 bed  -> parameter used to filter the mutations on a user-provided BED file
				  bed_ranges  -> dictionary that contains all of the ranges for each chromosome dictated by the user's input BED file
				 chrom_based  -> flag that generates the matrices on a chromosome basis
						plot  -> flag that plots the matrices after they are generated
					 tsb_ref  -> dictionary that allows for switching between binary and biologically relevant strings
			 transcript_path  -> path to the transcript files
						  gs  -> flag that generates a file for the strand bias on a gene basis.
					log_file  -> path to the output log file

	Returns:
		If called as a function, returns a nested dictionary of panda data frames for each matrix

	Outputs:
		Outputs a mutational matrix for the desired context within [user-provided input path]/output/[mut_type]/

	'''

	out = open(log_file, 'a')
	# Instantiates the necessary variables/data structures
	indel_dict = {}
	indel_tsb_dict = {}
	indel_simple_dict = {}
	samples = []
	indel_types_tsb = []
	indel_types_simple = []
	range_index = 0
	tsb_abrev = ['T','U','B','N']

	revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])

						# Single point mutations
	indel_types = ['1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5',
				   '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5',
				   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
				   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', 
						# >1bp INDELS
				   '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
				   '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
				   '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
				   '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
				   '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', 
				   '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', 
				   '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
				   '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
						#MicroHomology INDELS
				   '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
				   '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5', '2:Ins:M:1', 
				   '3:Ins:M:1', '3:Ins:M:2', '4:Ins:M:1', '4:Ins:M:2', '4:Ins:M:3', '5:Ins:M:1', 
				   '5:Ins:M:2', '5:Ins:M:3', '5:Ins:M:4', '5:Ins:M:5', 'complex', 'non_matching']

	for indels in indel_types[:24]:
		for tsbs in tsb_abrev:
			indel_types_tsb.append(tsbs + ":" + indels)

	if limited_indel:
		indel_types = indel_types[:-13]

	indel_types_simple = indel_types[:24]
	indel_types_simple.append('long_Del')
	indel_types_simple.append('long_Ins')
	indel_types_simple.append('MH')
	indel_types_simple.append('complex')

	# Instantiates the remaining varibales/data structures
	i = 0
	chrom_string = None
	count = 0
	non_matching = 0
	complex_muts = 0
	skipped_count = 0
	total_analyzed = 0

	# Creates files for exome and region sorting
	if exome:
		exome_temp_file = "exome_temp.txt"
		exome_file = open(vcf_path + exome_temp_file, 'w')

	if bed:
		bed_temp_file = "bed_temp.txt"
		bed_file = open(vcf_path + bed_temp_file, "w")

	# Opens the input vcf files
	with open (vcf_path + vcf_files[0]) as data:
		prev_line = None

		# Creates files for the genes strand bias test
		if gs:
			out = open(output_matrix + "gene_strand_bias_counts_indel.txt", "w")
			out_hot = open(output_matrix + "gene_strand_bias_counts_hotspots_indel.txt", "w")
			gene_ranges, gene_counts, gene_names, sample_mut_counts_per_gene, sample_mut_counts_per_mut_type = gene_range(transcript_path, True)

		first_flag = True

		# Saves the relevant data from each line
		for lines in data:
			try:
				line = lines.strip().split('\t')
				sample = line[1]
				chrom = line[5]
				# if chrom in ncbi_chrom.keys():
				# 	chrom = ncbi_chrom[chrom]

				start = int(line[6])
				ref = line[8].upper()
				mut = line[9].upper()
				if ref == '-':
					mut = '-' + mut
				if mut == '-':
					ref = '-' + ref
				sub_type_length = 0

				mut_type = None

				# Exception handling for incorrect input formats
				skip = False
				for r in ref:
					if r not in 'ACGT-':
						print("The ref contains bases that are not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						skipped_count += 1
						skip = True
						break
				if skip:
					continue

				skip = False
				for m in mut:
					if m not in 'ACGT-':
						print("The mutation contains bases that are not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						skipped_count += 1
						skip = True
						break
				if skip:
					continue

				if ref == mut:
					print("The ref base appears to match the mutated base. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
					skipped_count += 1
					continue

				if line == prev_line:
					print("There appears to be a duplicate single base substitution. Skipping this mutation: " + chrom + str(start) + ref + mut, file=out)
					skipped_count += 1
					continue

				prev_line = line

				# Opens the first chromosome reference file
				if first_flag:
					initial_chrom = chrom
					try:
						with open (chrom_path + initial_chrom + '.txt', "rb") as f:
							chrom_string = f.read().strip()
						first_flag = False
					except:
						print(initial_chrom + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...", file=out)
						continue

				# Opens the next chromosome when a new chromosome is reached in the file
				if chrom != initial_chrom:
					if chrom_based:
						matrix_generator_INDEL(output_matrix, samples, indel_types, indel_types_tsb, indel_types_simple, indel_dict, indel_tsb_dict, indel_simple_dict, project, exome, limited_indel, bed, initial_chrom)
						indel_dict = {}
						indel_tsb_dict = {}
						indel_simple_dict = {}
					initial_chrom = chrom
					try:
						with open (chrom_path + initial_chrom + '.txt', "rb") as f:
							chrom_string = f.read().strip()
						i += 1
						print("Chromosome "+ chrom + " complete", file=out)
					except:
						print(initial_chrom + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...", file=out)
						continue

				# Saves the relevant chromosome information from the reference file
				try:
					base = tsb_ref[chrom_string[start-1]][1]
				except:
					print("The start position is out of range. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
					skipped_count += 1
					continue

				if ref[0] != base and ref[0] != '-':
					print("The reference base does not match the reference chromosome position. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
					skipped_count += 1
					continue

				if ref[0] == base or ref[0] == '-':
					bias = tsb_ref[chrom_string[start-1]][0]

					# Saves the mutation type for the given variant
					if len(ref) - len(mut) == len(ref)-1:
						mut_type = 'Del'
						ref_base = ref[1]
						if ref_base == 'G' or ref_base == 'A':
							ref_base = revcompl(ref_base)
						if ref_base == tsb_ref[chrom_string[start-1]][1]:
							strand = '1'
						else:
							strand = '-1'
					elif len(mut) - len(ref) == len(mut)-1:
						mut_type = 'Ins'
						mut_base = mut[1]
						if mut_base == 'G' or mut_base == 'A':
							mut_base = revcompl(mut_base)
						if mut_base == tsb_ref[chrom_string[start-1]][1]:
							strand = '1'
						else:
							strand = '-1'
					else:
						# Add code to count complex mutations
						if sample not in indel_dict.keys():
							indel_dict[sample] = {}
							indel_tsb_dict[sample] = {}
							indel_simple_dict[sample] = {}

						if not functionFlag:
							if 'complex' not in indel_dict[sample].keys():
								indel_dict[sample]['complex'] = 1
								indel_simple_dict[sample]['complex'] = 1
							else:
								indel_dict[sample]['complex'] += 1
								indel_simple_dict[sample]['complex'] += 1
						continue
					

					type_sequence = ''
					
					# Pulls out the mutation subtype for deletions
					if mut_type == 'Del': 
						type_sequence = ref[1:]   

						type_length = len(type_sequence)
						sequence = type_sequence
						pos = start + type_length 
						pos_rev = start 
						actual_seq = ''
						for i in range (pos_rev-type_length, pos_rev, 1):
							actual_seq += tsb_ref[chrom_string[i]][1]                    
						while pos_rev - type_length > 0 and actual_seq == type_sequence:
							sequence = actual_seq + sequence
							pos_rev -= type_length
							actual_seq = ''
							for i in range (pos_rev-type_length, pos_rev, 1):
								actual_seq += tsb_ref[chrom_string[i]][1]
						
						new_seq = ''
						for i in range(pos, pos+type_length, 1):
							new_seq += tsb_ref[chrom_string[i]][1]
						while pos + type_length < len(chrom_string) and new_seq == type_sequence:
							sequence += new_seq
							pos += type_length
							new_seq = ''
							for i in range(pos, pos+type_length, 1):
								new_seq += tsb_ref[chrom_string[i]][1]
						
						# Pulls out possible microhomology deletions
						if type_length > 1 and len(sequence) == type_length:
							forward_homology = ref[1:-1]
							reverse_homology = ref[2:]

							
							for_hom = False
							pos = start + type_length
							for i in range (len(forward_homology), 0, -1):
								seq = ''
								for l in range (pos, pos + i, 1):
									seq += tsb_ref[chrom_string[l]][1]
								if seq == forward_homology[:i]:
									sequence += forward_homology[:i]
									mut_type += '_Micro_for'
									for_hom = True
									break

							if for_hom != True:
								pos = start
								for i in range (len(reverse_homology), 0, -1):
									seq = ''
									for l in range (pos-i, pos, 1):
										seq += tsb_ref[chrom_string[l]][1]

									if seq == reverse_homology[-i:]:
										sequence = reverse_homology[-i:] + sequence
										mut_type += '_Micro_rev'
										break
					
					# Pulls out the mutation subtype for insertions
					elif mut_type == 'Ins':        
						type_sequence = mut[1:]
						type_length = len(type_sequence)
						sequence = type_sequence

						pos = start
						pos_rev = start
						seq = ''
						for i in range(pos_rev-type_length, pos_rev, 1):
							seq += tsb_ref[chrom_string[i]][1]
						while pos_rev - type_length > 0 and seq == type_sequence:
							sequence = seq + sequence
							pos_rev -= type_length
							seq = ''
							for i in range(pos_rev-type_length, pos_rev, 1):
								seq += tsb_ref[chrom_string[i]][1]

						seq = ''
						for i in range(pos, pos + type_length, 1):
							seq += tsb_ref[chrom_string[i]][1]
						while pos + type_length < len(chrom_string) and seq == type_sequence:
							sequence += seq
							pos += type_length
							seq = ''
							for i in range(pos, pos + type_length, 1):
								seq += tsb_ref[chrom_string[i]][1]  

						# Pulls possible microhomology for insertions
						if type_length > 1 and len(sequence) == type_length:
							forward_homology = mut[1:-1]
							reverse_homology = mut[2:]
							
							for_hom = False
							pos = start
							for i in range (len(forward_homology), 0, -1):
								seq = ''
								for i in range (pos, pos + i, 1):
									seq += tsb_ref[chrom_string[i]][1]
								if seq == forward_homology[:i]:
									sequence += forward_homology[:i]
									mut_type += '_Micro_for'
									for_hom = True
									break

							if for_hom != True:
								pos = start
								for i in range (len(reverse_homology), 0, -1):
									seq = ''
									for i in range (pos-i, pos, 1):
										seq += tsb_ref[chrom_string[i]][1]
									if seq == reverse_homology[-i:]:
										sequence = reverse_homology[-i:] + sequence
										mut_type += '_Micro_rev'
										break

					# Saves the sample name for later reference
					if sample not in indel_dict.keys():
						indel_dict[sample] = {}
						indel_tsb_dict[sample] = {}
						indel_simple_dict[sample] = {}
					if sample not in samples:
						samples.append(sample)

					# Instantiates variables used to create the unique INDEL keys
					indel_key_1 = None
					indel_key_2 = None
					indel_key_3 = None
					indel_key_4 = None
					indel_key = 'blah'

					output_sequence = None

					# Creates the INDEL key for all deletions
					if mut_type[0:3] == 'Del': 
						indel_key_2 = 'Del'

						# Includes deletions of >1 bp
						if len(ref)-1 > 1: 
							key_1 = len(ref)-1
							if key_1 < 5:
								indel_key_1 = key_1
							else:
								indel_key_1 = 5

							# Only regular deleletions
							if mut_type == 'Del': 
								indel_key_3 = 'R'
								key_4 = int(len(sequence)/key_1 - 1)
								if key_4 < 5:
									indel_key_4 = key_4
								else: 
									indel_key_4 = 5

							# Only for microhomologies
							else:
								indel_key_3 = 'M'
								key_4 = len(sequence) - (len(ref)-1) 
								if key_4 > 5:
									indel_key_4 = 5
								elif key_4 < 0:
									print(lines)
								else:
									indel_key_4 = key_4
					
						# For deletions of 1bp
						else:
							indel_key_1 = 1
							key_4 = len(sequence) -1 
							if key_4 > 5:
								indel_key_4 = 5
							else:
								indel_key_4 = key_4
							
							if ref[1] == 'C' or ref[1] == 'G':
								indel_key_3 = 'C'
							
							else:
								indel_key_3 = 'T'

								
					# Creates the INDEL key for all insertions
					elif mut_type[0:3] == 'Ins':
						indel_key_2 = 'Ins'

						#Includes insertions of >1bp
						if len(mut)-1 > 1:
							key_1 = len(mut)-1
							if key_1<5:
								indel_key_1 = key_1
							else:
								indel_key_1 = 5
								
							# Only regular insertions
							if mut_type == 'Ins':
								indel_key_3 = 'R'
								key_4 = int(len(sequence)/key_1 - 1)
								if key_4 < 5:
									indel_key_4 = key_4
								else:
									indel_key_4 = 5
							# Only for microhomologies
							else:
								indel_key_3 = 'M'
								key_4 = len(sequence) - (len(mut)-1) 
								if key_4 >= 5:
									indel_key_4 = 5
								elif key_4 < 0:
									print(lines)
								else:
									indel_key_4 = key_4
								
						# Includes insertions of 1bp
						else:
							indel_key_1 = 1
							key_4 = len(sequence)-1
							if key_4 >= 5:
								indel_key_4 = 5
							else:
								indel_key_4 = key_4
								
							if mut[1] == 'C' or mut[1] == 'G':
								indel_key_3 = 'C'
							else:
								indel_key_3 = 'T'

					# Counts the number of "complex" mutations
					else:
						non_matching += 1

					# Creates the final INDEl key and saves it into the data structure
					indel_key_simple = ''
					if limited_indel and indel_key_2 == 'Ins' and indel_key_3 == 'M':
							indel_key = str(indel_key_1) + ':' + indel_key_2 + ':' + 'R' + ':' + '0'
					else:        
						indel_key = str(indel_key_1) +':'+indel_key_2+':'+indel_key_3+':'+str(indel_key_4)

					if int(indel_key_1) > 1:
						if indel_key_3 == 'M':
							indel_key_simple = 'MH'
						else:
							indel_key_simple = 'long_' + indel_key_2
					else:
						indel_key_simple = indel_key


					indel_key_tsb = bias + ":" + indel_key


					# Performs the gene strand bias test if desired
					if gs:
						if indel_key[0] == '1':	
							continue_flag = False
							for ranges in gene_ranges[initial_chrom][range_index:]:
								if start < ranges[0]:
									break
								if ranges[0] <= start <= ranges[1]:
									continue_flag = True
									gene_index = gene_ranges[initial_chrom].index(ranges)
									gene = gene_names[initial_chrom][gene_index]
									if strand == ranges[2]:
										dict_key = 'T'
										gene_counts[gene][dict_key] += 1
										if sample not in gene_counts[gene]['samples']:
											gene_counts[gene]['samples'].append(sample)
											sample_mut_counts_per_gene[gene][sample] = 1
											sample_mut_counts_per_mut_type[gene][sample] = {'T':0, 'U':0}
											sample_mut_counts_per_mut_type[gene][sample][dict_key] += 1

										else:
											sample_mut_counts_per_gene[gene][sample] += 1
											sample_mut_counts_per_mut_type[gene][sample][dict_key] += 1
									elif int(strand) + int(ranges[2]) == 0:
										dict_key = 'U'
										gene_counts[gene][dict_key] += 1
										if sample not in gene_counts[gene]['samples']:
											gene_counts[gene]['samples'].append(sample)
											sample_mut_counts_per_gene[gene][sample] = 1
											sample_mut_counts_per_mut_type[gene][sample] = {'T':0, 'U':0}
											sample_mut_counts_per_mut_type[gene][sample][dict_key] += 1
										else:
											sample_mut_counts_per_gene[gene][sample] += 1
											sample_mut_counts_per_mut_type[gene][sample][dict_key] += 1




					# Saves the mutation type into the correct dictionary key
					if indel_key not in indel_dict[sample].keys():
						indel_dict[sample][indel_key] = 1
					else:
						indel_dict[sample][indel_key] += 1

					if indel_key_tsb not in indel_tsb_dict[sample].keys():
						indel_tsb_dict[sample][indel_key_tsb] = 1
					else:
						indel_tsb_dict[sample][indel_key_tsb] += 1

					if indel_key_simple not in indel_simple_dict[sample].keys():
						indel_simple_dict[sample][indel_key_simple] = 1
					else:
						indel_simple_dict[sample][indel_key_simple] += 1


					total_analyzed += 1

					# Writes the INDEL to a temporary file for exome/region sorting
					if exome:
						exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + indel_key + "\t" + ref + "\t" + mut + "\n")
					if bed:
						bed_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + indel_key + "\t" + ref + "\t" + mut + "\n")

				else:
					if not functionFlag:
						if 'non_matching' not in indel_dict[sample].keys():
							indel_dict[sample]['non_matching'] = 1
						else:
							indel_dict[sample]['non_matching'] += 1
			except:
				print("There appears to be an error in this line. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
				skipped_count += 1


		# Once all of the variants have been account for, complete the gene strand bias test
		if gs:
			pvals = []
			qvals = []
			pvals_hot = []
			qvals_hot = []
			hotspots = {}
			for gene in gene_counts:
				total_count = sum(sample_mut_counts_per_gene[gene].values())
				for sample in sample_mut_counts_per_gene[gene]:
					mut_count = sample_mut_counts_per_gene[gene][sample]
					if mut_count > 10 and mut_count/total_count > 0.5:
						if gene not in hotspots:
							hotspots[gene] = {}
							for mut, count in sample_mut_counts_per_mut_type[gene][sample].items():
								hotspots[gene][mut] = count
							hotspots[gene]['samples'] = [sample]
							for mut, count in sample_mut_counts_per_mut_type[gene][sample].items():
								gene_counts[gene][mut] -= count
							gene_counts[gene]['samples'].remove(sample)
						else:
							for mut, count in sample_mut_counts_per_mut_type[gene][sample].items():
								hotspots[gene][mut] += count
								gene_counts[gene][mut] -= count
							gene_counts[gene]['samples'].remove(sample)
							hotspots[gene]['samples'].append(sample)

				sum_tran = 0
				sum_untran = 0
				for mut, counts in gene_counts[gene].items():
					if mut[0] == 'T':
						sum_tran += counts
					elif mut[0] == 'U':
						sum_untran += counts
				pvals.append(stats.binom_test([sum_tran, sum_untran]))

				sum_tran_hot = 0
				sum_untran_hot = 0
				if gene in hotspots:
					for mut, counts in hotspots[gene].items():
						if mut[0] == 'T':
							sum_tran_hot += counts
						elif mut[0] == 'U':
							sum_untran_hot += counts
				pvals_hot.append(stats.binom_test([sum_tran_hot, sum_untran_hot]))

			qvals = sm.fdrcorrection(pvals)[1]
			qvals_hot = sm.fdrcorrection(pvals_hot)[1]
			gene_ind = 0
			for gene in gene_counts:
				gene_counts[gene]['samples'] = len(gene_counts[gene]['samples'])
				print(gene, end='',file=out, flush=False)
				sum_tran = 0
				sum_untran = 0
				for mut, counts in gene_counts[gene].items():
					if mut[0] == 'T':
						sum_tran += counts
					elif mut[0] == 'U':
						sum_untran += counts
					print("\t" + str(counts), end='', file=out, flush=False)
				print("\t" + str(qvals[gene_ind]), flush=False, file=out)
				gene_ind += 1
			out.close()
			with open(output_matrix + "gene_strand_bias_counts_indel.txt") as f2:
				lines = [line.strip().split() for line in f2]
			output = open(output_matrix + "gene_strand_bias_counts_indel.txt", 'w')
			print('GENE\tTranscribed\tUntranscribed\tSampleCount\tq_value', file=output)
			for line in sorted(lines, key = lambda x: (float(x[-1])), reverse=False):
				print('\t'.join(line), file=output)
			output.close()

			# Performs gene strand bias test hot spot samples
			gene_ind = 0
			for gene in hotspots:
				hotspots[gene]['samples'] = len(hotspots[gene]['samples'])
				print(gene, end='',file=out_hot, flush=False)
				sum_tran_hot = 0
				sum_untran_hot = 0
				for mut, counts in hotspots[gene].items():
					if mut[0] == 'T':
						sum_tran_hot += counts
					elif mut[0] == 'U':
						sum_untran_hot += counts
					print("\t" + str(counts), end='', file=out_hot, flush=False)
				print("\t" + str(qvals_hot[gene_ind]), flush=False, file=out_hot)
				gene_ind += 1
			out_hot.close()
			with open(output_matrix + "gene_strand_bias_counts_hotspots_indel.txt") as f2:
				lines = [line.strip().split() for line in f2]
			output = open(output_matrix + "gene_strand_bias_counts_hotspots_indel.txt", 'w')
			print('GENE\tTranscribed\tUntranscribed\ttSampleCount\tq_value', file=output)
			for line in sorted(lines, key = lambda x: (float(x[-1])), reverse=False):
				print('\t'.join(line), file=output)
			output.close()

	# Prints the total number of complex mutations
	print("Non-matching mutations: " + str(non_matching), file=out)
	if chrom_based:
		matrix_generator_INDEL(output_matrix, samples, indel_types, indel_types_tsb, indel_dict, indel_tsb_dict, indel_simple_dict, project, exome, limited_indel, bed, initial_chrom, plot)

	# Performs the final filter on the variants base upon the exome if desired
	if exome:
		exome_file.close()
		sort_file = exome_temp_file
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		indel_dict, samples2 = exome_check(genome, vcf_path + exome_temp_file, output_matrix, project)

	# Performs the final filter on the variants base upon the region if desired
	if bed:
		bed_file.close()
		sort_file = bed_temp_file
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()
		indel_dict, samples2 = panel_check(genome, vcf_path + bed_temp_file, output_matrix, bed_file_path, project)

	# Calls the function to generate the final mutational matrix
	if not chrom_based:
		initial_chrom=None
		matrix_generator_INDEL(output_matrix, samples, indel_types, indel_types_tsb, indel_types_simple, indel_dict, indel_tsb_dict, indel_simple_dict, project, exome, limited_indel, bed, initial_chrom, plot)
		if functionFlag:
			return(pd.DataFrame.from_dict(indel_dict), skipped_count, total_analyzed, len(samples))

def exome_check (genome, exome_temp_file, output_matrix, project):
	'''
	Filters the variants for those present within the exome. 

	Parameters:
				 genome  -> name of the genome of interest (ex: GRCh37)
		exome_temp_file  -> The temporary file that contains all of the variants used for filtering
		  output_matrix  -> path to the final output matrix folder
				project  -> name of the unique project

	Returns:
		  mutation_dict  -> updated mutation dictionary for each sample for each mutation type post filtering
				samples  -> updated list of samples that still contain mutations post filtering

	''' 

	# Instantiates the relevant variables/data structures
	base_cushion = 200
	mutation_dict = {}
	samples = []

	initial = True
	udpate_chrom = False
	current_dir = os.path.realpath(__file__)
	#current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts/SigProfilerMatrixGenerator.py$', '', current_dir)

	exome_file = ref_dir + "/references/chromosomes/exome/" + genome + "/" + genome + "_exome.interval_list"

	#exome_output_path = ref_dir + "/references/vcf_files/" + project + "_exome/SNV/"
	exome_output_path = output_matrix + "vcf_files/SNV/"
	exome_output = exome_output_path + project + "_exome.vcf"
	#exome_output = ref_dir + "/references/vcf_files/" + project + "_exome/SNV/" + project + "_exome.vcf"
	if not os.path.exists(exome_output_path):
		os.makedirs(exome_output_path)

	with open(exome_temp_file) as f, open(exome_file) as exome, open(exome_output, "w") as out:
		previous_chrom_ref = None
		previous_chrom_start = None
		previous_chrom_end = None

		chrom_ref = None
		start_ref = None
		end_ref = None

		read = True

		for lines in f:
			# Saves the relevant data for the current variant for later reference
			line = lines.strip().split()
			sample = line[0]
			chrom = line[1]
			start = int(line[2])
			mut_type = line[3]
			ref = line[4]
			mut = line[5]

			# Saves a value for the x and y chromosomes as a numeric reference
			if chrom == 'X':
				chrom_value = -1
			elif chrom == 'Y':
				chrom_value = 0
			elif chrom == 'MT' or chrom == 'M':
				chrom_value = 100
			else:
				chrom_value = int(chrom)


			if initial:
				chrom_start = chrom
				initial = False

			stop = False
			while not stop:
				if chrom == previous_chrom_ref:
					if start >= previous_chrom_start - base_cushion and start <= previous_chrom_end + base_cushion:
						if sample not in mutation_dict.keys():
							samples.append(sample)
							mutation_dict[sample] = {}
						if mut_type not in mutation_dict[sample].keys():
							mutation_dict[sample][mut_type] = 1
						else:
							mutation_dict[sample][mut_type] += 1
						read = True
						print('\t'.join([chrom, str(start), ".", ref, mut ]), file=out)
						break

				if read:
					lines2 = exome.readline()
				try:
					if lines2[0] == "@":
						continue
				except:
					break
				else:
					line2 = lines2.strip().split('\t')
					chrom_ref = line2[0]
					if len(chrom_ref) > 2:
						chrom_ref = chrom_ref[3:]
					start_ref = int(line2[1])
					end_ref = int(line2[2])

					if chrom_ref == 'X':
						ref_chrom_value = -1
					elif chrom_ref == 'Y':
						ref_chrom_value = 0
					else:
						ref_chrom_value = int(chrom_ref)

					if chrom == chrom_ref:

						if start > (start_ref - base_cushion and end_ref + base_cushion):
							read = True
							continue
						elif start >= start_ref - base_cushion and start <= end_ref + base_cushion: 
							if sample not in mutation_dict.keys():
								samples.append(sample)
								mutation_dict[sample] = {}
							if mut_type not in mutation_dict[sample].keys():
								mutation_dict[sample][mut_type] = 1
							else:
								mutation_dict[sample][mut_type] += 1
							read = True
							print('\t'.join([chrom, str(start), ".", ref, mut ]), file=out)
							break
						elif start < (start_ref - base_cushion):
							read = False
							break


					else:
						if chrom_value < ref_chrom_value:
							read = False
							break
						elif chrom_value > ref_chrom_value:
							read = True
							continue



			chrom_start = chrom
			previous_chrom_ref = chrom_ref
			previous_chrom_start = start_ref
			previous_chrom_end = end_ref

	# logging.info("Exome filtering is complete. Proceeding with the final catalogue generation...")
	# print("Exome filtering is complete. Proceeding with the final catalogue generation...")
	return(mutation_dict, samples)

def panel_check (genome, bed_temp_file, output_matrix, bed_file_path, project):
	'''
	Filters the variants for those present within the exome. 

	Parameters:
				 genome  -> name of the genome of interest (ex: GRCh37)
		exome_temp_file  -> The temporary file that contains all of the variants used for filtering
		  output_matrix  -> path to the final output matrix folder
		  bed_file_path  -> path to the bed file 
				project  -> unique project name


	Returns:
		  mutation_dict  -> updated mutation dictionary for each sample for each mutation type post filtering
				samples  -> updated list of samples that still contain mutations post filtering

	''' 

	# Instantiates the relevant variables/data structures
	base_cushion = 200
	mutation_dict = {}
	samples = []

	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)


	initial = True
	udpate_chrom = False
	current_dir = os.path.realpath(__file__)
	panel_file = bed_file_path
	panel_output_path = output_matrix + "vcf_files/SNV/"
	panel_output = panel_output_path + project + "_panel.vcf"
	if not os.path.exists(panel_output_path):
		os.makedirs(panel_output_path)


	with open(bed_temp_file) as f, open(panel_file) as exome, open(panel_output, "w") as out:
		previous_chrom_ref = None
		previous_chrom_start = None
		previous_chrom_end = None

		chrom_ref = None
		start_ref = None
		end_ref = None

		read = True

		for lines in f:
			# Saves the relevant data for the current variant for later reference
			line = lines.strip().split()
			sample = line[0]
			chrom = line[1]
			start = int(line[2])
			mut_type = line[3]
			ref = line[4]
			mut = line[5]

			# Saves a value for the x and y chromosomes as a numeric reference
			if chrom == 'X':
				chrom_value = -1
			elif chrom == 'Y':
				chrom_value = 0
			elif chrom == 'MT' or chrom == 'M':
				chrom_value = 100
			else:
				chrom_value = int(chrom)


			if initial:
				chrom_start = chrom
				initial = False

			stop = False
			while not stop:
				if chrom == previous_chrom_ref:
					if start >= previous_chrom_start - base_cushion and start <= previous_chrom_end + base_cushion:
						if sample not in mutation_dict.keys():
							samples.append(sample)
							mutation_dict[sample] = {}
						if mut_type not in mutation_dict[sample].keys():
							mutation_dict[sample][mut_type] = 1
						else:
							mutation_dict[sample][mut_type] += 1
						read = True
						print('\t'.join([chrom, str(start), ".", ref, mut ]), file=out)
						break

				if read:
					lines2 = exome.readline()
				try:
					if lines2[0] == "@":
						continue
				except:
					break
				else:
					line2 = lines2.strip().split('\t')
					chrom_ref = line2[0]
					if len(chrom_ref) > 2:
						chrom_ref = chrom_ref[3:]
					start_ref = int(line2[1])
					end_ref = int(line2[2])

					if chrom_ref == 'X':
						ref_chrom_value = -1
					elif chrom_ref == 'Y':
						ref_chrom_value = 0
					else:
						ref_chrom_value = int(chrom_ref)

					if chrom == chrom_ref:

						if start > (start_ref - base_cushion and end_ref + base_cushion):
							read = True
							continue
						elif start >= start_ref - base_cushion and start <= end_ref + base_cushion: 
							if sample not in mutation_dict.keys():
								samples.append(sample)
								mutation_dict[sample] = {}
							if mut_type not in mutation_dict[sample].keys():
								mutation_dict[sample][mut_type] = 1
							else:
								mutation_dict[sample][mut_type] += 1
							read = True
							print('\t'.join([chrom, str(start), ".", ref, mut ]), file=out)
							break
						elif start < (start_ref - base_cushion):
							read = False
							break


					else:
						if chrom_value < ref_chrom_value:
							read = False
							break
						elif chrom_value > ref_chrom_value:
							read = True
							continue



			chrom_start = chrom
			previous_chrom_ref = chrom_ref
			previous_chrom_start = start_ref
			previous_chrom_end = end_ref

	# logging.info("Panel filtering is complete. Proceeding with the final catalogue generation...")
	# print("Panel filtering is complete. Proceeding with the final catalogue generation...")
	return(mutation_dict, samples)



def matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start=None, functionFlag=False, plot=False, tsb_stat=False):
	'''
	Writes the final mutational matrix given a dictionary of samples, mutation types, and counts

	Parameters:
					context  -> desired context (ex. 96, 1536, 384, 6144)
			  output_matrix  -> path where the final mutational matrix is stored
					project  -> unique name given to the set of samples (ex. 'BRCA') 
					 samples -> a list of all sample names
				   bias_sort -> dictionary that provides the sorting order for the TSB matrices
			 muatation_dict  -> dictionary with the counts for each mutation type for each sample
					  types  -> list of the mutation types for the given context 
					  exome  -> Boolean for whether the catalogue should be generated across the whole
								genome or just the exome
				  mut_types  -> list with all possible mutation types for the given context
						bed  -> parameter used to filter the mutations on a user-provided BED file
				chrom_start  -> current chromosome when generating the matrix on a chromosome basis
			   functionFlag  -> flag that will return the matrix in memory.
					   plot  -> flag that will generate the plots for each context
				   tsb_stat  -> performs a transcriptional strand bias test for the 24, 384, and 6144 contexts. The output is
								saved into the output/TSB directory


	Returns:
		None

	Output:
		Write the final mutational matrix for 96, 384, 1536, 6144 contexts

	'''

	# Prepares all of the required data structures and files
	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)

	contexts = ['96', '384', '1536', '6', '24']
	mut_types_all = {'96':[], '384':[], '1536':[], '6':[], '24':[], '6_pvalue':[], '7_pvalue':[]}
	mut_count_all = {'96':{}, '384':{}, '1536':{}, '6':{}, '24':{}, '6_pvalue':{}, '7_pvalue':{}, '6_pvalue_temp':{}, '7_pvalue_temp':{}}

	mut_6144 = pd.DataFrame.from_dict(mutation_dict)
	mut_6144 = mut_6144.fillna(0)
	mut_6144.index.name = 'MutationType'

	zeros_data = np.zeros(len(mut_6144.columns.values))
	available_types = list(mut_6144.index)
	for mut in mut_types:
		if mut not in available_types:
			mut_6144.loc[mut] = zeros_data
	mut_6144 = mut_6144.astype(int)
	mut_count_all['6144'] = mut_6144
	mut_count_all['1536'] = mut_6144.groupby(mut_6144.index.str[2:]).sum()
	mut_count_all['96'] = mut_6144.groupby(mut_6144.index.str[3:10]).sum()
	mut_count_all['6'] = mut_6144.groupby(mut_6144.index.str[5:8]).sum()
	mut_count_all['384'] = mut_6144.groupby(mut_6144.index.str[0:2] + mut_6144.index.str[3:10]).sum()
	mut_count_all['24'] = mut_6144.groupby(mut_6144.index.str[0:2] + mut_6144.index.str[5:8]).sum()

	mut_count_all['1536'].index.name = 'MutationType'
	mut_count_all['96'].index.name = 'MutationType'
	mut_count_all['6'].index.name = 'MutationType'
	mut_count_all['384'].index.name = 'MutationType'
	mut_count_all['24'].index.name = 'MutationType'

	strandBiasOut = output_matrix + "TSB/"
	output_matrix_SBS = output_matrix + "SBS/"

	if not os.path.exists(output_matrix_SBS):
		os.mkdir(output_matrix_SBS)

	file_prefix = project + ".SBS" + context
	if exome:
		output_file_matrix = output_matrix_SBS + file_prefix + ".exome"
	else:
		if bed:
			output_file_matrix = output_matrix_SBS + file_prefix + ".region"
		else:
			output_file_matrix = output_matrix_SBS + file_prefix + ".all"

	if chrom_start != None:
		output_file_matrix += ".chr" + chrom_start

	types = list(mut_6144.index)
	types = sorted(types, key=lambda val: (bias_sort[val[0]], val[2:]))
	mut_6144 = mut_6144.reindex(types)

	types = list(mut_count_all['384'].index)
	types = sorted(types, key=lambda val: (bias_sort[val[0]], val[2:]))
	mut_count_all['384'] = mut_count_all['384'].reindex(types)

	types = list(mut_count_all['24'].index)
	types = sorted(types, key=lambda val: (bias_sort[val[0]], val[2:]))
	mut_count_all['24'] = mut_count_all['24'].reindex(types)

	mut_6144.to_csv(output_file_matrix, header=True, sep='\t')

	# TSB test:
	if tsb_stat:

		types_7pN = list(mut_count_all['24'].index)
		types.append('T:T>CpN')
		types.append('U:T>CpN')
		types_7pN = list(mut_count_all['384'].index)


		mut_count_all['6_pvalue_temp'] = mut_count_all['24']
		mut_count_all['7_pvalue_temp'] = mut_count_all['384']
		mut_count_all['6_pvalue'] = mut_count_all['6']
		mut_count_all['7_pvalue'] = mut_count_all['6']

		mut_count_all['6_pvalue'] = mut_count_all['6'].to_dict('dict')
		mut_count_all['7_pvalue'] = mut_count_all['6'].to_dict('dict')

		drop_muts = ['T:A[T>C]A', 'T:A[T>C]C', 'T:A[T>C]G', 'T:A[T>C]T',
					 'U:A[T>C]A', 'U:A[T>C]C', 'U:A[T>C]G', 'U:A[T>C]T',
					 'B:A[T>C]A', 'B:A[T>C]C', 'B:A[T>C]G', 'B:A[T>C]T',
					 'N:A[T>C]A', 'N:A[T>C]C', 'N:A[T>C]G', 'N:A[T>C]T']
		add_counts_T = []
		add_counts_U = []
		add_counts_B = []
		add_counts_N = []

		for mut in drop_muts:
			mut_count_all['7_pvalue_temp'] = mut_count_all['7_pvalue_temp'].drop(mut)
			if mut[0] == 'T':
				add_counts_T.append(mut_count_all['384'].loc[mut])
			elif mut[0] == 'U':
				add_counts_U.append(mut_count_all['384'].loc[mut])
			elif mut[0] == 'B':
				add_counts_B.append(mut_count_all['384'].loc[mut])
			elif mut[0] == 'N':
				add_counts_N.append(mut_count_all['384'].loc[mut])


		# don't combine the categories like that. Need to keep TSB categories
		mut_count_all['7_pvalue_temp'] = mut_count_all['7_pvalue_temp'].groupby(mut_count_all['7_pvalue_temp'].index.str[0:2] + mut_count_all['7_pvalue_temp'].index.str[4:7]).sum()

		mut_count_all['7_pvalue_temp'].loc['T:T>CpN'] = sum(add_counts_T)
		mut_count_all['7_pvalue_temp'].loc['U:T>CpN'] = sum(add_counts_U)
		mut_count_all['7_pvalue_temp'].loc['B:T>CpN'] = sum(add_counts_B)
		mut_count_all['7_pvalue_temp'].loc['N:T>CpN'] = sum(add_counts_N)

		function_tsb_test = ['6_pvalue_temp', '7_pvalue_temp']
		for cont in function_tsb_test:
			cont_save = cont[:8]
			if cont == '6_pvalue_temp':
				types = list(mut_count_all['24'].index)
			else:
				types = list(mut_count_all['24'].index)
				types.append('T:T>CpN')
				types.append('U:T>CpN')
			types = list(set(types))
			current_tsb = mut_count_all[cont]
			for sample in samples:
				pvals = []
				for mut_type in types:
					if mut_type[0] == 'T':
						pval = stats.binom_test([current_tsb.loc[mut_type][sample], current_tsb.loc['U:'+mut_type[2:]][sample]])
						if current_tsb.loc[mut_type][sample] >= current_tsb.loc['U:'+mut_type[2:]][sample]:
							strand_test = 1
						else:
							strand_test = -1
						mut_count_all[cont_save][sample][mut_type[2:]] = [pval, strand_test]


		mut_count_all['6_pvalue'] = pd.DataFrame.from_dict(mut_count_all['6_pvalue'])
		mut_count_all['7_pvalue'] = pd.DataFrame.from_dict(mut_count_all['7_pvalue'])


		if not os.path.exists(strandBiasOut):
			os.mkdir(strandBiasOut)
		significant_tsb = open(strandBiasOut + "significantResults_strandBiasTest.txt", 'w')
		with open (strandBiasOut + "strandBiasTest_6144.txt", 'w') as out2, open (strandBiasOut + "strandBiasTest_384.txt", 'w') as out384, open (strandBiasOut + "strandBiasTest_24.txt", 'w') as out24:
			print("Sample\tMutationType\tEnrichment[Trans/UnTrans]\tp.value\tFDR_q.value",file=out2)
			tsb_6144 = mut_6144[mut_6144.index.str[0] == 'T']
			tsb_6144_U = mut_6144[mut_6144.index.str[0] == 'U']
			tsb_index = [x[2:] for x in list(tsb_6144.index)]
			tsb_6144.index = tsb_index
			tsb_6144_U.index = tsb_index

			tsb_384 = mut_count_all['384'][mut_count_all['384'].index.str[0] == 'T']
			tsb_384_U = mut_count_all['384'][mut_count_all['384'].index.str[0] == 'U']
			tsb_index384 = [x[2:] for x in list(tsb_384.index)]
			tsb_384.index = tsb_index384
			tsb_384_U.index = tsb_index384

			tsb_24 = mut_count_all['24'][mut_count_all['24'].index.str[0] == 'T']
			tsb_24_U = mut_count_all['24'][mut_count_all['24'].index.str[0] == 'U']
			tsb_index24 = [x[2:] for x in list(tsb_24.index)]
			tsb_24.index = tsb_index24
			tsb_24_U.index = tsb_index24


			enr = tsb_6144/tsb_6144_U
			enr_384 = tsb_384/tsb_384_U
			enr_24 = tsb_24/tsb_24_U

			for sample in samples:
				pvals = []
				pvals_384 = []
				pvals_24 = []
				enrichment = list(enr[sample])
				enrichment_384 = list(enr[sample])
				enrichment_24 = list(enr[sample])

				for i in range(0, len(tsb_6144), 1):
					pvals.append(stats.binom_test([tsb_6144[sample][i], tsb_6144_U[sample][i]]))
				for i in range(0, len(tsb_384), 1):
					pvals_384.append(stats.binom_test([tsb_384[sample][i], tsb_384_U[sample][i]]))
				for i in range(0, len(tsb_24), 1):
					pvals_24.append(stats.binom_test([tsb_24[sample][i], tsb_24_U[sample][i]]))

				qvals = sm.fdrcorrection(pvals)[1]
				qvals_384 = sm.fdrcorrection(pvals)[1]
				qvals_24 = sm.fdrcorrection(pvals)[1]
				p_index = 0
				for mut in tsb_index:
					print(sample + "\t" + mut + "\t" + str(enrichment[p_index]) + "\t" + str(pvals[p_index]) + "\t" + str(qvals[p_index]), file=out2)
					if qvals[p_index] < 0.01:
						 print(sample + "\t" + mut + "\t" + str(enrichment[p_index]) + "\t" + str(pvals[p_index]) + "\t" + str(qvals[p_index]), file=significant_tsb)
				p_index += 1

				p_index = 0
				for mut in tsb_index384:
					print(sample + "\t" + mut + "\t" + str(enrichment_384[p_index]) + "\t" + str(pvals_384[p_index]) + "\t" + str(qvals_384[p_index]), file=out384)
					if qvals_384[p_index] < 0.01:
						 print(sample + "\t" + mut + "\t" + str(enrichment_384[p_index]) + "\t" + str(pvals_384[p_index]) + "\t" + str(qvals_384[p_index]), file=significant_tsb)
				p_index += 1

				p_index = 0
				for mut in tsb_index24:
					print(sample + "\t" + mut + "\t" + str(enrichment_24[p_index]) + "\t" + str(pvals_24[p_index]) + "\t" + str(qvals_24[p_index]), file=out24)
					if qvals_24[p_index] < 0.01:
						 print(sample + "\t" + mut + "\t" + str(enrichment_24[p_index]) + "\t" + str(pvals_24[p_index]) + "\t" + str(qvals_24[p_index]), file=significant_tsb)
				p_index += 1
		significant_tsb.close()


	# Generates the matrices for the remaining matrices (1536, 384, 96, 24, 6) by 
	# summing the counts from the 6144 matrix
	for cont in contexts:
		mutation_dict = mut_count_all[cont]

		file_prefix = project + ".SBS" + cont
		if exome:
			output_file_matrix = output_matrix_SBS + file_prefix + ".exome"
		else:
			if bed:
				output_file_matrix = output_matrix_SBS + file_prefix + ".region"
			else:
				output_file_matrix = output_matrix_SBS + file_prefix + ".all"

		if chrom_start != None:
			output_file_matrix += ".chr" + chrom_start

		mutation_dict.to_csv(output_file_matrix, header=True, sep='\t')
						


		if plot:
			output_path = output_matrix + "plots/"
			if not os.path.exists(output_path):
				os.makedirs(output_path)
			if cont == '96':
				try:
					sigPlt.plotSBS(output_file_matrix, output_path, project, '96', False)
				except:
					pass
			elif cont == '384':
				try:
					sigPlt.plotSBS(output_file_matrix, output_path, project, '384', False)
				except:
					pass
			elif cont == '6':
				try:
					sigPlt.plotSBS(output_file_matrix, output_path, project, '6', False)
				except:
					pass
			elif cont == '24':
				try:
					sigPlt.plotSBS(output_file_matrix, output_path, project, '24', False)
				except:
					pass


	# If this code is run as an imported function, delete the physcial matrix.
	if functionFlag:		
		return(mut_count_all)

def matrix_generator_INDEL (output_matrix, samples, indel_types, indel_types_tsb, indel_types_simple, indel_dict, indel_tsb_dict, indel_simple_dict, project, exome, limited_indel, bed, initial_chrom=None, plot=False):
	'''
	Writes the final mutational matrix for INDELS given a dictionary of samples, INDEL types, and counts

	Parameters:
			  output_matrix  -> path where the final mutational matrix is stored
					samples  -> a list of all sample names
				indel_types  -> list of the INDEL types 
			indel_types_tsb  -> list of the INDEL types for the TSB matrix
		 indel_types_simple  -> list of the INDEL types for the simple categories matrix
				 indel_dict  -> dictionary with the counts for each INDEL type for each sample
			 indel_tsb_dict  -> dictionary with the TSB counts for each INDEL type (only 1bp INDELs)
		  indel_simple_dict  -> dictionary with the counts for each simple INDEL type
					project  -> unique name given to the set of samples (ex. 'BRCA') 
					  exome  -> Boolean for whether the catalogue should be generated across the whole
								genome or just the exome
			  limited_indel  -> flag that instructs the function to create the a limited indel matrix
						bed  -> parameter used to filter the mutations on a user-provided BED file
			  initial_chrom  -> the current chromosome to generate the matrix for.
					   plot  -> flag that will generate the INDEL plots for the provided samples.

	Returns:
		None

	Output:
		Write the final mutational matrix for INDELS

	'''

	# Instantiates all of the required data structures and output files
	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)

	bias_sort = {'T':0,'U':1,'N':3,'B':2}

	output_matrix_INDEL = output_matrix + "INDEL/"
	if not os.path.exists(output_matrix_INDEL):
		os.mkdir(output_matrix_INDEL)

	if limited_indel:
		file_prefix = project + ".INDEL83"
		file_prefix_tsb = project + ".INDEL96"
		file_prefix_simple = project + ".INDEL28"
	else:
		file_prefix = project + ".INDEL94"
		file_prefix_tsb = project + ".INDEL96"
		file_prefix_simple = project + ".INDEL28"
	
	if exome:
		output_file_matrix = output_matrix_INDEL + file_prefix + ".exome"
		output_file_matrix_tsb = output_matrix_INDEL + file_prefix_tsb + ".exome"
		output_file_matrix_simple = output_matrix_INDEL + file_prefix_simple + ".exome"
	else:
		if bed:
			output_file_matrix = output_matrix_INDEL + file_prefix + ".region"
			output_file_matrix_tsb = output_matrix_INDEL + file_prefix_tsb + ".region"
			output_file_matrix_simple = output_matrix_INDEL + file_prefix_simple + ".region"
		else:
			output_file_matrix = output_matrix_INDEL + file_prefix + ".all"
			output_file_matrix_tsb = output_matrix_INDEL + file_prefix_tsb + ".all"
			output_file_matrix_simple = output_matrix_INDEL + file_prefix_simple + ".all"

	if initial_chrom != None:
		output_file_matrix += ".chr" + initial_chrom
		output_file_matrix_tsb += ".chr" + initial_chrom
		output_file_matrix_simple += ".chr" + initial_chrom


	with open (output_file_matrix, 'w') as out:
		# Prints all of the sample names into the first line of the file
		print ('MutationType\t', end='', flush=False, file=out)  
		samples.sort()
		for sample in samples:
			print (sample + '\t', end='', flush=False, file=out)
		print(file=out)
		
		# Prints the mutation count for each INDEL type across every sample
		for indel in indel_types:
			print (indel + '\t', end='', flush =False, file=out)
			for sample in samples:
				if sample not in indel_dict.keys():
					indel_dict[sample] = {}
				if indel in indel_dict[sample].keys():
					print (str(indel_dict[sample][indel]) + '\t', end='', file=out)
				else:
					print ('0\t', end='', file=out)
			print(file=out)


	with open (output_file_matrix_tsb, 'w') as out:
		# Prints all of the sample names into the first line of the file
		print ('MutationType\t', end='', flush=False, file=out)  
		samples.sort()
		for sample in samples:
			print (sample + '\t', end='', flush=False, file=out)
		print(file=out)

		types = sorted(indel_types_tsb, key=lambda val: (bias_sort[val[0]], val[2:]))
		# Prints the mutation count for each INDEL type across every sample
		for indel in types:
			print (indel + '\t', end='', flush =False, file=out)
			for sample in samples:
				if sample not in indel_tsb_dict.keys():
					indel_tsb_dict[sample] = {}
				if indel in indel_tsb_dict[sample].keys():
					print (str(indel_tsb_dict[sample][indel]) + '\t', end='', file=out)
				else:
					print ('0\t', end='', file=out)
			print(file=out)

	with open (output_file_matrix_simple, 'w') as out:
		# Prints all of the sample names into the first line of the file
		print ('MutationType\t', end='', flush=False, file=out)  
		samples.sort()
		for sample in samples:
			print (sample + '\t', end='', flush=False, file=out)
		print(file=out)

		# Prints the mutation count for each INDEL type across every sample
		for indel in indel_types_simple:
			print (indel + '\t', end='', flush =False, file=out)
			for sample in samples:
				if sample not in indel_simple_dict.keys():
					indel_simple_dict[sample] = {}
				if indel in indel_simple_dict[sample].keys():
					print (str(indel_simple_dict[sample][indel]) + '\t', end='', file=out)
				else:
					print ('0\t', end='', file=out)
			print(file=out)


	if plot:
		output_path = output_matrix + "plots/"
		if not os.path.exists(output_path):
			os.mkdir(output_path)
		try:
			sigPlt.plotID(output_file_matrix, output_path, project, '94', False)
		except:
			pass
		try:
			sigPlt.plotID(output_file_matrix_tsb, output_path, project, '96ID', False)
		except:
			print("no")
			pass

def matrix_generator_DINUC (output_matrix, samples, bias_sort, all_dinucs, all_mut_types, dinucs, project, exome, bed, chrom_start=None, plot=False):
	'''
	Writes the final mutational matrix for INDELS given a dictionary of samples, INDEL types, and counts

	Parameters:
			output_matrix  -> path where the final mutational matrix is stored
				  samples  -> a list of all sample names
				bias_sort  -> dictionary that provides the sorting order for the TSB matrices
			   all_dinucs  -> dictionary that contains all of the  mutation counts for each DINUC context
			all_mut_types  -> dictionary that contains all of the mutation types for each DINUC context
				   dinucs  -> dictionary with the counts for each DINUC type for each sample
				  project  -> unique name given to the set of samples (ex. 'BRCA') 
					exome  -> Boolean for whether the catalogue should be generated across the whole
							  genome or just the exome
					  bed  -> parameter used to filter the mutations on a user-provided BED file
			  chrom_start  -> current chromosome to generate the matrix for
					 plot  -> flag that generates the DINUC plots for the provided samples

	Returns:
		None

	Output:
		Write the final mutational matrix for DINUCs

	'''

	revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])
	revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1','U':'T','T':'U','B':'B','N':'N'}[B] for B in x][::-1])
	if not any(dinucs):
		return()

	output_matrix_DINUC = output_matrix + "DINUC/"
	if not os.path.exists(output_matrix_DINUC):
		os.mkdir(output_matrix_DINUC)

	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)
	for cont in all_dinucs.keys():
		mutation_types = all_mut_types[cont]
		dinucs = all_dinucs[cont]
		file_prefix = project + ".DBS" + cont
		if exome:
			output_file_matrix = output_matrix_DINUC + file_prefix + ".exome"
		else:
			if bed:
				output_file_matrix = output_matrix_DINUC + file_prefix + ".region"
			else:
				output_file_matrix = output_matrix_DINUC + file_prefix + ".all"
		if chrom_start != None:
			output_file_matrix += ".chr" + chrom_start

		with open (output_file_matrix, 'w') as out:
			# Prints all of the sample names into the first line of the file
			print ('MutationType\t', end='', flush=False, file=out)  
			for sample in samples:
				print (sample + '\t', end='', flush=False, file=out)
			print(file=out)
		
			if cont == '312' or cont == '4992':
				try:
					mutation_types = sorted(mutation_types, key=lambda val: (bias_sort[val[0]], val[2:]))
				except:
					print(mutation_types)
			# Prints the mutation count for each INDEL type across every sample
			for dinuc in mutation_types:
				if cont == '312' and dinuc[0] == 'T':
					count = 0
					if dinuc[2] == revcompl(dinuc[3]) and dinuc[5] == revcompl(dinuc[6]):
						try:
							if dinuc in dinucs[sample]:
								count += dinucs[sample][dinuc]
						except:
							dinucs[sample] = {dinuc:0}
						rev_dinuc = revbias(dinuc[0]) + dinuc[1:]
						if rev_dinuc in dinucs[sample]:
							count += dinucs[sample][rev_dinuc]
						if count >= 2:
							dinucs[sample][dinuc] = int(count/2)
							if count%2 != 0:							
								dinucs[sample][rev_dinuc] = int(count/2) + 1
							else:
								dinucs[sample][rev_dinuc] = int(count/2)
				
				if cont == '4992' and dinuc[0] == 'T':
					count = 0
					if dinuc[2] == revcompl(dinuc[10]) and dinuc[4] == revcompl(dinuc[5]) and dinuc[7] == revcompl(dinuc[8]):
						try:
							if dinuc in dinucs[sample]:
								count += dinucs[sample][dinuc]
						except:
							dinucs[sample] = {dinuc:0}
						rev_dinuc = revbias(dinuc[0]) + dinuc[1:]
						if rev_dinuc in dinucs[sample]:
							count += dinucs[sample][rev_dinuc]
						if count >= 2:
							dinucs[sample][dinuc] = int(count/2)
							if count%2 != 0:							
								dinucs[sample][rev_dinuc] = int(count/2) + 1
							else:
								dinucs[sample][rev_dinuc] = int(count/2)

				print (dinuc + '\t', end='', flush =False, file=out)
				for sample in samples:
					try:
						if dinuc in dinucs[sample]:
							print(str(dinucs[sample][dinuc]) + "\t",end='',file=out)
						else:
							print('0\t',end='',file=out)
					except:
						print('0\t',end='',file=out)
				print(file=out)
	

		if plot:
			output_path = output_matrix + "plots/"
			if not os.path.exists(output_path):
				os.mkdir(output_path)
			try:
				if cont == '78' or cont == '312':
					sigPlt.plotDBS(output_file_matrix, output_path, project, cont, False)
			except:
				pass
