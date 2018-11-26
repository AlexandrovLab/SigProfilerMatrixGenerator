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
import sys
import re
import argparse
import itertools
import pandas as pd
from itertools import chain
import time
import logging
import datetime
from scipy import stats
import statsmodels.stats.multitest as sm
import convert_input_to_simple_files as convertIn
import sigProfilerPlotting as sigPlt
import uuid
from collections import defaultdict
from collections import OrderedDict


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
		next(f)
		for lines in f:
			line = lines.strip().split()
			chrom = line[0]
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




def catalogue_generator_single (vcf_path, vcf_files, chrom_path, project, output_matrix, context, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, transcript_path, gs):
	'''
	Generates the mutational matrix for 96, 1536, 192, and 3072 context using a single vcf file with all samples of interest.

	Parameters:
					vcf_path  -> path to vcf file of interest
				   vcf_files  -> actual vcf file
				  chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
								file name: '1.txt', '2.txt', etc.
					 project  -> unique name given to the set of samples (ex. 'BRCA') 
			   output_matrix  -> path where the final mutational matrix is stored
					 context  -> desired context (ex. 96, 1536, 192, 3072)
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
			 transcript_path  -> path to the transcript files
						  gs  -> flag that generates a file for the strand bias on a gene basis.

	Returns:
		None

	Outputs:
		Outputs a mutational matrix for the desired context within /references/matrix/[project]/

	'''

	# Small functions to provide reverse complements of TSB and sequence info:
	revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])
	revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1','U':'T','T':'U','B':'B','N':'N'}[B] for B in x][::-1])
	
	# Provides the sorting order for the TSB matrices
	bias_sort = {'T':0,'U':1,'N':3,'B':2}
	tsb = ['T','U','N','B']
	bases = ['A','C','G','T']
	
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



	types = []
	samples = []
	mutation_dict = {}
	flag = True
	i = 0
	file = vcf_files[0]
	sample_start = None
	gene_counts = {}
	range_index = 0
	sample_count = {}
	if exome:
		exome_temp_file = "exome_temp.txt"
		exome_file = open(exome_temp_file, 'w')

	with open (vcf_path + file) as f:
			if gs:
				out = open(output_matrix + "gene_strand_bias_counts_SNV.txt", "w")
				out_hot = open(output_matrix + "gene_strand_bias_counts_hotspots_SNV.txt", "w")
				gene_ranges, gene_counts, gene_names, sample_mut_counts_per_gene, sample_mut_counts_per_mut_type = gene_range(transcript_path)
			
			# Skips any header lines
			for lines in f:
				if lines[0] == '#':
					next(f)

				# Saves all relevant data from the file
				else:
					line = lines.strip().split('\t')
					sample = line[1]
					chrom = line[5]
					if chrom in ncbi_chrom.keys():
						chrom = ncbi_chrom[chrom]
					if len(chrom) > 1:
						if chrom[0:3].upper() == 'CHR':
							chrom = chrom[-1]

					start = int(line[6])
					ref = line[8][0]
					mut = line[9][0]

					# Skips any mutations not present within the custom bed file
					if bed:
						if start not in bed_ranges[chrom]:
							print(chrom, start)
							continue
 
					if flag:
						sample_start = sample
						chrom_start = chrom

						# Opens each chromosome string path and bias string path
						try:
							with open(chrom_path + chrom_start + ".txt", "rb") as f:
								chrom_string = f.read()

						except:

							print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
							continue
						flag = False

					# Saves the sample name for later reference
					if sample not in samples:
						samples.append(sample)
						mutation_dict[sample] = {}
						sample_count[sample] = [0,0]

					if sample not in mutation_dict.keys():
						mutation_dict[sample] = {}

					# Opens the next chromosome data once a new chromosome is reached
					if chrom != chrom_start:
						range_index = 0
						logging.info("Chromosome " + chrom_start + " done")
						if chrom_based:
							matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start)
							mutation_dict = {}
							mutation_dict[sample] = {}
						chrom_start = chrom
						try:
							with open(chrom_path + chrom_start + ".txt", "rb") as f:
								chrom_string = f.read()
							i += 1
						except:
							print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
							continue

					# Pulls out the relevant sequence depending on the context
					sequence = ''
					for i in range (start-3, start+2, 1):
						try:
							sequence += tsb_ref[chrom_string[i]][1]
						except:
							print(i, chrom)
					bias = tsb_ref[chrom_string[start]][0]
					char = sequence[int(len(sequence)/2)]

					# Prints the sequence and position if the pulled sequence doesn't match
					# the variant from the file
					if char != ref and revcompl(char) != ref:
						print (sequence + "\t" + chrom + "\t" + char + "\t" + ref + "\t" + str(start))
				
					
					# Saves the sequence/mutation type if it matched the reference/reverse strand  
					else:
						if char == ref:
							strand = '1'
							bias_before = bias
							ref_before = ref
							if ref == 'A' or ref == 'G':
								strand = '-1'
								bias = revbias(bias)
								ref = revcompl(ref)
								mut = revcompl(mut)
								sequence = revcompl(sequence)



						else:
							strand = '-1'
							if ref == 'A' or ref == 'G':
								strand = '1'
								ref = revcompl(ref)
								mut = revcompl(mut) 
								sequence = revcompl(sequence)
							else:
								#pass
								bias = revbias(bias)

						# Performs the gene strand bias test if desired
						if gs:	
							for ranges in gene_ranges[chrom_start][range_index:]:
								dict_key = ref + ">" + mut
								if start < ranges[0]:
									break
								if ranges[0] <= start <= ranges[1]:
									gene_index = gene_ranges[chrom_start].index(ranges)
									gene = gene_names[chrom_start][gene_index]
									#if strand == ranges[2]:
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
									#elif int(strand) + int(ranges[2]) == 0:
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
						mut_key = sequence[0:int(len(sequence)/2)] + '[' + ref + '>' + mut + ']' + sequence[int(len(sequence)/2+1):]
						mut_key = bias + ':' + mut_key

						if mut_key not in types:
							types.append(mut_key)
						try:
							if mut_key not in mutation_dict[sample].keys():
								mutation_dict[sample][mut_key] = 1

							else:
								mutation_dict[sample][mut_key] += 1
						except:
							pass

						# If exome is specified, it will write the variant to a temporary exome file.
						if exome:
							exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + mut_key + "\n")


			# Once all variants are accounted for, complete the gene strand bias test/output to the final file. 
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
						if mut_count > 3 and mut_count/total_count > 0.5:
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

	logging.info("Chromosome " + chrom_start + " done")

	# Generate the matrix for the final chromosome if specified by user.
	if chrom_based:
		matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start)
		mutation_dict = {}

	# Organizes the variants based upon the exome
	if exome:
		exome_file.close()
		sort_file = exome_temp_file
		with open(sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()
		mutation_dict, samples = exome_check(genome, exome_temp_file)


	# Calls the function to generate the final matrix
	if functionFlag:
		chrom_start = None
		chrom_start=None
		matrices = matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start, functionFlag, plot)
		return(matrices)
	else:
		if not chrom_based:
			chrom_start=None
			matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start, functionFlag, plot)


def catalogue_generator_DINUC_single (vcf_path, vcf_files, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref):
	'''
	Generates the mutational matrix for the dinucleotide context.

	Parameters:
				   vcf_path  -> path to the vcf file of interest
				  vcf_files  -> vcf file of interest
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


	Returns:
		None

	Output:
		Saves the mutational matrix for the dinucleotide context

	'''

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
	# mutation_types = ['AA>CC','AA>CG','AA>CT','AA>GC','AA>GG','AA>GT','AA>TC','AA>TG','AA>TT',
	# 				  'AC>CA','AC>CG','AC>CT','AC>GA','AC>GG','AC>GT','AC>TA','AC>TG','AC>TT',
	# 				  'AG>CA','AG>CC','AG>CT','AG>GA','AG>GC','AG>GT','AG>TA','AG>TC','AG>TT',
	# 				  'AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA','CA>AC','CA>AG','CA>AT',
	# 				  'CA>GC','CA>GG','CA>GT','CA>TC','CA>TG','CA>TT','CC>AA','CC>AG','CC>AT',
	# 				  'CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT','CG>AA','CG>AC','CG>AT',
	# 				  'CG>GA','CG>GC','CG>TA','GA>AC','GA>AG','GA>AT','GA>CC','GA>CG','GA>CT',
	# 				  'GA>TC','GA>TG','GA>TT','GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA',
	# 				  'TA>AC','TA>AG','TA>AT','TA>CC','TA>CG','TA>GC']

	mutation_types = ['AC>CA','AC>CG','AC>CT','AC>GA','AC>GG','AC>GT','AC>TA','AC>TG','AC>TT',
					  'AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA','CC>AA','CC>AG','CC>AT',
					  'CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT','CG>AA','CG>AC','CG>AT',
					  'CG>GA','CG>GC','CG>TA','CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG',
					  'CT>TA','CT>TC','CT>TG','GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA',
					  'TA>AC','TA>AG','TA>AT','TA>CC','TA>CG','TA>GC','TC>AA','TC>AG','TC>AT',
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

	types = []
	samples = []
	mutation_dict = {}
	flag = True
	i = 0
	file = vcf_files[0]
	sample_start = None

	if exome:
		exome_temp_file = "exome_temp.txt"
		exome_file = open(exome_temp_file, 'w')

	with open (vcf_path + file) as data:
			initial_line = data.readline()
			initial_line_data = initial_line.strip().split('\t')
			previous_sample = initial_line_data[1]
			previous_chrom = initial_line_data[5] 
			previous_start = int(initial_line_data[6])
			previous_ref = initial_line_data[8]
			previous_mut = initial_line_data[9]
					
			
			for lines in data:
				# Skips any header lines
				if lines[0] == '#':
					next(data)

				# Saves the relevant data from each line
				else:
					line = lines.strip().split('\t')
					sample = line[1]
					chrom = line[5]
					if chrom in ncbi_chrom.keys():
						chrom = ncbi_chrom[chrom]
					if len(chrom) > 1:
						if chrom[0:3].upper() == 'CHR':
							chrom = chrom[-1]

					start = int(line[6])
					ref = line[8]
					mut = line[9]

					if bed:
						if start not in bed_ranges[chrom]:
							continue

					# Saves the sample name for later reference
					if sample not in samples:
						samples.append(sample)

					if flag:
						sample_start = sample
						chrom_start = chrom
						with open (chrom_path + chrom + ".txt", "rb") as f:
							chrom_string = f.read()
						flag = False

					# Breaks the loop once a new chromosome is reached
					if chrom != chrom_start:
						logging.info("Chromosome " + chrom_start + " done")
						if chrom_based:
							all_dinucs = {'78':dinucs, '312':dinucs_tsb, '1248':dinucs_context, '4992':dinucs_context_tsb}
							all_mut_types = {'78':mutation_types, '312':mutation_types_tsb, '1248':mutation_types_context, '4992':mutation_types_tsb_context}
							matrix_generator_DINUC (output_matrix, samples, bias_sort, all_dinucs, all_mut_types, dinucs, project, exome, bed, chrom_start, plot)
							dinucs = {}
						chrom_start = chrom
						with open (chrom_path + chrom + ".txt", "rb") as f:
							chrom_string = f.read()
						

					# Grabs the slice of interest
					else:
						if start == previous_start + 1:
							dinuc = previous_ref + ref + ">" + previous_mut + mut

							dinuc_seq = "".join([tsb_ref[chrom_string[start-3]][1],"[",dinuc,"]",tsb_ref[chrom_string[start]][1]])
							bias = tsb_ref[chrom_string[start-1]][0]
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

							if exome:
								exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc + "\n")

					previous_sample = sample
					previous_chrom = chrom 
					previous_start = start
					previous_ref = ref
					previous_mut = mut

	# Organizes the required dictionaries for the final matrix generation.
	all_dinucs = {'78':dinucs, '312':dinucs_tsb, '1248':dinucs_context, '4992':dinucs_context_tsb}
	all_mut_types = {'78':mutation_types, '312':mutation_types_tsb, '1248':mutation_types_context, '4992':mutation_types_tsb_context}
	logging.info("Chromosome " + chrom_start + " done")
	if chrom_based:
		matrix_generator_DINUC (output_matrix, samples, bias_sort, all_dinucs, all_mut_types, dinucs, project, exome, bed, chrom_start, plot)
		dinucs = {}

	if exome:
		exome_file.close()

		sort_file = exome_temp_file
		with open(sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		all_dinucs, samples = exome_check(genome, exome_temp_file)

	# Calls the function to generate the final mutational matrix
	if functionFlag:
		return(pd.DataFrame.from_dict(dinucs))
	else:
		if not chrom_based:
			chrom_start=None
			matrix_generator_DINUC (output_matrix, samples, bias_sort, all_dinucs, all_mut_types, dinucs, project, exome, bed, chrom_start, plot)
		

		
	

def catalogue_generator_INDEL_single (vcf_path, vcf_files, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, limited_indel, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, transcript_path, gs):
	'''
	Generates the mutational matrix for the INDEL context.

	Parameters:
				   vcf_path  -> path to the vcf file of interest
				  vcf_files  -> vcf file of interest
				 chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
								file name: '1.txt', '2.txt', etc.
					project  -> unique name given to the set of samples (ex. 'BRCA') 
			  output_matrix  -> path where the final mutational matrix is stored
				 ncbi_chrom  -> dictionary that allows for the converstion of ncbi chromosome names to standard format
								for the mm10 assembly.
						bed  -> parameter used to filter the mutations on a user-provided BED file
				 bed_ranges  -> dictionary that contains all of the ranges for each chromosome dictated by the user's input BED file
			transcript_path  -> path to the transcript files
						 gs  -> flag that generates a file for the strand bias on a gene basis.


	Returns:
		None

	Output:
		Saves the mutational matrix for the INDEL context

	'''

	# Instantiates the necessary variables/data structures
	indel_dict = {}
	indel_tsb_dict = {}
	samples = []
	indel_types_tsb = []
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

	# Instantiates the remaining varibales/data structures
	i = 0
	chrom_string = None
	count = 0
	non_matching = 0
	complex_muts = 0

	if exome:
		exome_temp_file = "exome_temp.txt"
		exome_file = open(exome_temp_file, 'w')

	with open (vcf_path + vcf_files[0]) as data:
		if gs:
			out = open(output_matrix + "gene_strand_bias_counts_indel.txt", "w")
			out_hot = open(output_matrix + "gene_strand_bias_counts_hotspots_indel.txt", "w")
			gene_ranges, gene_counts, gene_names, sample_mut_counts_per_gene, sample_mut_counts_per_mut_type = gene_range(transcript_path, True)

		first_flag = True

		# Saves the relevant data from each line
		for lines in data:
			line = lines.strip().split('\t')
			sample = line[1]
			chrom = line[5]
			if chrom in ncbi_chrom.keys():
				chrom = ncbi_chrom[chrom]
			if len(chrom) > 1:
				if chrom[0:3].upper() == 'CHR':
					chrom = chrom[-1]

			start = int(line[6])
			ref = line[8]
			mut = line[9]
			sub_type_length = 0

			mut_type = None

			# Skips the variant if it is not in the custom BED ranges
			if bed:
				if start not in bed_ranges[chrom]:
					continue

			if first_flag:
				initial_chrom = chrom
				try:
					with open (chrom_path + initial_chrom + '.txt', "rb") as f:
						chrom_string = f.read().strip()
					first_flag = False
				except:
					print(initial_chrom + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
					continue

			# Opens the next chromosome when a new chromosome is reached in the file
			if chrom != initial_chrom:
				if chrom_based:
					matrix_generator_INDEL(output_matrix, samples, indel_types, indel_types_tsb, indel_dict, indel_tsb_dict, project, exome, limited_indel, bed, initial_chrom)
					indel_dict = {}
					indel_tsb_dict = {}
				initial_chrom = chrom
				try:
					with open (chrom_path + initial_chrom + '.txt', "rb") as f:
						chrom_string = f.read().strip()
					i += 1
					logging.info("Chromosome "+ chrom + " complete")
				except:
					print(initial_chrom + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
					continue
			try:
				base = tsb_ref[chrom_string[start-1]][1]
			except:
				print(chrom, start)
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
					if sample not in indel_dict.keys():
						indel_dict[sample] = {}
						indel_tsb_dict[sample] = {}

					if not functionFlag:
						if 'complex' not in indel_dict[sample].keys():
							indel_dict[sample]['complex'] = 1
						else:
							indel_dict[sample]['complex'] += 1
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
				if limited_indel and indel_key_2 == 'Ins' and indel_key_3 == 'M':
						indel_key = str(indel_key_1) + ':' + indel_key_2 + ':' + 'R' + ':' + '0'

				else:        
					indel_key = str(indel_key_1) +':'+indel_key_2+':'+indel_key_3+':'+str(indel_key_4)

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

				if exome:
					exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + indel_key + "\n")
			else:
				if not functionFlag:
					if 'non_matching' not in indel_dict[sample].keys():
						indel_dict[sample]['non_matching'] = 1
					else:
						indel_dict[sample]['non_matching'] += 1


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
	logging.info("Non-matching mutations: " + str(non_matching))
	if chrom_based:
		matrix_generator_INDEL(output_matrix, samples, indel_types, indel_types_tsb, indel_dict, indel_tsb_dict, project, exome, limited_indel, bed, initial_chrom, plot)

	# Performs the final filter on the variants base upon the exome if desired
	if exome:
		exome_file.close()
		sort_file = exome_temp_file
		with open(sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'].index(x[1]), int(x[2]))):
			print('\t'.join(line), file=output)
		output.close()

		indel_dict, samples = exome_check(genome, exome_temp_file)

	# Calls the function to generate the final mutational matrix
	if functionFlag:
		return(pd.DataFrame.from_dict(indel_dict))
	else:
		if not chrom_based:
			initial_chrom=None
			matrix_generator_INDEL(output_matrix, samples, indel_types, indel_types_tsb, indel_dict, indel_tsb_dict, project, exome, limited_indel, bed, initial_chrom, plot)

def exome_check (genome, exome_temp_file):
	'''
	Filters the variants for those present within the exome. 

	Parameters:
				 genome  -> name of the genome of interest (ex: GRCh37)
		exome_temp_file  -> The temporary file that contains all of the variants used for filtering

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
	exome_file = ref_dir + "/references/chromosomes/exome/" + genome + "/" + genome + "_exome.interval_list"


	with open(exome_temp_file) as f, open(exome_file) as exome:
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

			# Saves a value for the x and y chromosomes as a numeric reference
			if chrom == 'X':
				chrom_value = -1
			elif chrom == 'Y':
				chrom_value = 0
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

	logging.info("Exome filtering is complete. Proceeding with the final catalogue generation...")
	print("Exome filtering is complete. Proceeding with the final catalogue generation...")
	return(mutation_dict, samples)




def matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed, chrom_start=None, functionFlag=False, plot=False):
	'''
	Writes the final mutational matrix given a dictionary of samples, mutation types, and counts

	Parameters:
					context  -> desired context (ex. 96, 1536, 192, 3072)
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



	Returns:
		None

	Output:
		Write the final mutational matrix for 96, 192, 1536, 3072 contexts

	'''

	# Prepares all of the required data structures and files
	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)

	contexts = ['96', '192', '1536', '6', '12']
	mut_types_all = {'96':[], '192':[], '1536':[], '6':[], '12':[], '6_pvalue':[], '7_pvalue':[]}
	mut_count_all = {'96':{}, '192':{}, '1536':{}, '6':{}, '12':{}, '6_pvalue':{}, '7_pvalue':{}, '6_pvalue_temp':{}, '7_pvalue_temp':{}}


	significant_tsb = open(output_matrix + "significantResults_strandBiasTest.txt", 'w')

	file_prefix = project + ".SBS" + context
	if exome:
		output_file_matrix = output_matrix + file_prefix + ".exome"
	else:
		if bed:
			output_file_matrix = output_matrix + file_prefix + ".region"
		else:
			output_file_matrix = output_matrix + file_prefix + ".all"

	if chrom_start != None:
		output_file_matrix += ".chr" + chrom_start

	with open (output_file_matrix, 'w') as out:

		# Prints all of the sample names into the first line of the file
		print ('MutationType\t', end='', flush=False, file=out)  
		for sample in samples:
			print (sample + '\t', end='', flush=False, file=out)
		print(file=out)

		# Sorts the mutation types to standardize the output file
		try:
			types = sorted(mut_types, key=lambda val: (bias_sort[val[0]], val[2:]))
		except:
			print(mut_types)

		# Prints the mutation count for each mutation type across every sample
		for mut_type in types:
			print (mut_type + '\t', end='', flush =False, file=out)
			if mut_type[3:-1] not in mut_types_all['96']:
				mut_types_all['96'].append(mut_type[3:-1])
			if mut_type[0:2] + mut_type[3:-1] not in mut_types_all['192']:
				mut_types_all['192'].append(mut_type[0:2] + mut_type[3:-1])
			if mut_type[2:] not in mut_types_all['1536']:
				mut_types_all['1536'].append(mut_type[2:])
			if mut_type[5:8] not in mut_types_all['6']:
				mut_types_all['6'].append(mut_type[5:8])
			if mut_type[0:2] + mut_type[5:8] not in mut_types_all['12']:
				mut_types_all['12'].append(mut_type[0:2] + mut_type[5:8])
			if mut_type[5:8] not in mut_types_all['6_pvalue']:
				mut_types_all['6_pvalue'].append(mut_type[5:8])
			if mut_type[5:8] not in mut_types_all['7_pvalue']:
				mut_types_all['7_pvalue'].append(mut_type[5:8])



			for sample in samples:
				if sample not in mutation_dict:
					mutation_dict[sample] = {}
				if sample not in mut_count_all['96']:
					mut_count_all['96'][sample] = {}
					mut_count_all['192'][sample] = {}
					mut_count_all['1536'][sample] = {}
					mut_count_all['6'][sample] = {}
					mut_count_all['12'][sample] = {}
					mut_count_all['6_pvalue_temp'][sample] = {}
					mut_count_all['7_pvalue_temp'][sample] = {}
					mut_count_all['6_pvalue'][sample] = {}
					mut_count_all['7_pvalue'][sample] = {}

				if mut_type[3:-1] not in mut_count_all['96'][sample]:
					mut_count_all['96'][sample][mut_type[3:-1]] = 0
				if (mut_type[0:2]+mut_type[3:-1]) not in mut_count_all['192'][sample]:
					mut_count_all['192'][sample][mut_type[0:2] + mut_type[3:-1]] = 0
				if mut_type[2:] not in mut_count_all['1536'][sample]:
					mut_count_all['1536'][sample][mut_type[2:]] = 0
				if mut_type[5:8] not in mut_count_all['6'][sample]:
					mut_count_all['6'][sample][mut_type[5:8]] = 0
				if (mut_type[0:2] + mut_type[5:8]) not in mut_count_all['12'][sample]:
					mut_count_all['12'][sample][mut_type[0:2] + mut_type[5:8]] = 0
					mut_count_all['6_pvalue_temp'][sample][mut_type[0:2] + mut_type[5:8]] = 0
					mut_count_all['7_pvalue_temp'][sample][mut_type[0:2] + mut_type[5:8]] = 0
					if mut_type[3:9] == 'A[T>C]':
						if mut_type[0:2] + 'T>CpN' not in mut_count_all['7_pvalue_temp'][sample]:
							mut_count_all['7_pvalue_temp'][sample][mut_type[0:2] + 'T>CpN'] = 0
					else:
						if (mut_type[0:2] + mut_type[5:8]) not in mut_count_all['7_pvalue_temp'][sample]:
							mut_count_all['7_pvalue_temp'][sample][mut_type[0:2] + mut_type[5:8]] = 0
						

				if mut_type in mutation_dict[sample]:
					mut_count_all['96'][sample][mut_type[3:-1]] += mutation_dict[sample][mut_type]
					mut_count_all['192'][sample][mut_type[0:2] + mut_type[3:-1]] += mutation_dict[sample][mut_type]
					mut_count_all['1536'][sample][mut_type[2:]] += mutation_dict[sample][mut_type]
					mut_count_all['6'][sample][mut_type[5:8]] += mutation_dict[sample][mut_type]
					mut_count_all['12'][sample][mut_type[0:2] + mut_type[5:8]] += mutation_dict[sample][mut_type]                                  
					mut_count_all['6_pvalue_temp'][sample][mut_type[0:2] + mut_type[5:8]] += mutation_dict[sample][mut_type]
					if mut_type[3:9] == 'A[T>C]':
						mut_count_all['7_pvalue_temp'][sample][mut_type[0:2] + 'T>CpN'] += mutation_dict[sample][mut_type]
					else:
						mut_count_all['7_pvalue_temp'][sample][mut_type[0:2] + mut_type[5:8]] += mutation_dict[sample][mut_type]

				if mut_type in mutation_dict[sample]:
					print (str(mutation_dict[sample][mut_type]) + '\t', end='', file=out)
				else:
					print ('0\t', end='', file=out)
			print(file=out)

		# Writes a strand bias file for all TSB mutaiton types and outputs any significant ones to a 
		# separate file
		with open (output_matrix + "strandBiasTest_3072.txt", 'w') as out2:
			print("Sample\tMutationType\tEnrichment[Trans/UnTrans]\tp.value\tFDR_q.value",file=out2)
			current_tsb = pd.DataFrame.from_dict(mutation_dict)
			for sample in samples:
				pvals = []
				enrichment = []
				for mut_type in types:
					if mut_type[0] == 'T':
						if mut_type not in mutation_dict[sample]:
							num1 = 0
						else:
							num1 = current_tsb.loc[mut_type][sample]

						if 'U:'+mut_type[2:] not in mutation_dict[sample]:
							num2 = 0
						else:
							num2 = current_tsb.loc['U:'+mut_type[2:]][sample]

						pvals.append(stats.binom_test([num1, num2]))

						if 'U:'+mut_type[2:] not in mutation_dict[sample] or current_tsb.loc['U:'+mut_type[2:]][sample] == 0:
							enrichment.append(0)
						else:
							enrichment.append(round(num1/num2, 4))
				qvals = sm.fdrcorrection(pvals)[1]
				p_index = 0
				for mut_type in types:
					if mut_type[0] == 'T':
						print(sample + "\t" + mut_type[2:] + "\t" + str(enrichment[p_index]) + "\t" + str(pvals[p_index]) + "\t" + str(qvals[p_index]), file=out2)
						if qvals[p_index] < 0.01:
							 print(sample + "\t" + mut_type[2:] + "\t" + str(enrichment[p_index]) + "\t" + str(pvals[p_index]) + "\t" + str(qvals[p_index]), file=significant_tsb)
					p_index += 1
  
		# If this code is run as an imported function, delete the physcial matrix.
		if functionFlag:
			os.system("rm " + output_file_matrix) 
			mut_count_all['3072'] = mutation_dict
			function_tsb_test = ['6_pvalue_temp', '7_pvalue_temp']
			for cont in function_tsb_test:
				cont_save = cont[:8]
				if cont == '6_pvalue_temp':
					types = mut_types_all['12']
				else:
					types = mut_types_all['12']
					types.append('T:T>CpN')
					types.append('U:T>CpN')
					print(types)
				types = list(set(types))
				current_tsb = pd.DataFrame.from_dict(mut_count_all[cont])
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

					#qvals = sm.fdrcorrection(pvals)[1]
			
			return(mut_count_all)

	strandBias_test = ['12','192', '3072']
	strandBias_test = set(strandBias_test)

	# Generates the matrices for the remaining matrices (1536, 192, 96, 12, 6) by 
	# summing the counts from the 3072 matrix
	for cont in contexts:
		types = mut_types_all[cont]
		types = list(set(types))
		mutation_dict = mut_count_all[cont]
		file_prefix = project + ".SBS" + cont
		if exome:
			output_file_matrix = output_matrix + file_prefix + ".exome"
		else:
			if bed:
				output_file_matrix = output_matrix + file_prefix + ".region"
			else:
				output_file_matrix = output_matrix + file_prefix + ".all"

		if chrom_start != None:
			output_file_matrix += ".chr" + chrom_start


		with open (output_file_matrix, 'w') as out:

			# Prints all of the sample names into the first line of the file
			print ('MutationType\t', end='', flush=False, file=out)  
			for sample in samples:
				print (sample + '\t', end='', flush=False, file=out)
			print(file=out)

			if cont == '192' or cont == '3072' or cont == '12':
				try:
					types = sorted(types, key=lambda val: (bias_sort[val[0]], val[2:]))
				except:
					print(mut_types)

			# Prints the mutation count for each mutation type across every sample
			for mut_type in types:


				print (mut_type + '\t', end='', flush =False, file=out)
				for sample in samples:
					if mut_type in mutation_dict[sample].keys():
						print (str(mutation_dict[sample][mut_type]) + '\t', end='', file=out)
					else:
						print ('0\t', end='', file=out)
				print(file=out)

			if cont in strandBias_test:
				with open (output_matrix + "strandBiasTest_" + cont + ".txt", 'w') as out2:
					print("Sample\tMutationType\tEnrichment[Trans/UnTrans]\tp.value\tFDR_q.value",file=out2)
					current_tsb = pd.DataFrame.from_dict(mut_count_all[cont])
					for sample in samples:
						pvals = []
						enrichment = []
						for mut_type in types:
							if mut_type[0] == 'T':
								if cont in strandBias_test:
									pvals.append(stats.binom_test([current_tsb.loc[mut_type][sample], current_tsb.loc['U:'+mut_type[2:]][sample]]))
									if current_tsb.loc['U:'+mut_type[2:]][sample] == 0:
										enrichment.append(0)
									else:
										enrichment.append(round(current_tsb.loc[mut_type][sample]/current_tsb.loc['U:'+mut_type[2:]][sample], 4))
						qvals = sm.fdrcorrection(pvals)[1]
						p_index = 0
						for mut_type in types:
							if mut_type[0] == 'T':
								print(sample + "\t" + mut_type[2:] + "\t" + str(enrichment[p_index]) + "\t" + str(pvals[p_index]) + "\t" + str(qvals[p_index]), file=out2)
								if qvals[p_index] < 0.01:
									 print(sample + "\t" + mut_type[2:] + "\t" + str(enrichment[p_index]) + "\t" + str(pvals[p_index]) + "\t" + str(qvals[p_index]), file=significant_tsb)
							p_index += 1
						


		# sorts the 96 and 1536 matrices by mutation type
		if cont == '96' or cont == '1536' or cont == '6':
			command1 = "head -1 " + output_file_matrix + " > " + output_matrix + "a.tmp;"
			command2 = "tail -n+2 " + output_file_matrix + " | sort -n >> " + output_matrix + "a.tmp;"

			
			os.system(command1)
			os.system(command2)
			os.system("cat " + output_matrix + "a.tmp > " + output_file_matrix)
			os.system("rm " + output_matrix + "a.tmp")

		if plot:
			output_path = ref_dir + "/plots/" + project + "/"
			if not os.path.exists(output_path):
				os.makedirs(output_path)
			if cont == '96':
				sigPlt.plot96(output_file_matrix, output_path, False, project, False)
			elif cont == '192':
				sigPlt.plot192(output_file_matrix, output_path, False, project, False)


	significant_tsb.close()

def matrix_generator_INDEL (output_matrix, samples, indel_types, indel_types_tsb, indel_dict, indel_tsb_dict, project, exome, limited_indel, bed, initial_chrom=None, plot=False):
	'''
	Writes the final mutational matrix for INDELS given a dictionary of samples, INDEL types, and counts

	Parameters:
			  output_matrix  -> path where the final mutational matrix is stored
					samples  -> a list of all sample names
				indel_types  -> list of the INDEL types 
			indel_types_tsb  -> list of the INDEL types for the TSB matrix
				 indel_dict  -> dictionary with the counts for each INDEL type for each sample
			 indel_tsb_dict  -> dictionary with the TSB counts for each INDEL type (only 1bp INDELs)
					project  -> unique name given to the set of samples (ex. 'BRCA') 
					  exome  -> Boolean for whether the catalogue should be generated across the whole
								genome or just the exome
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

	file_prefix = project + ".INDEL94"
	file_prefix_tsb = project + ".INDEL96"
	if exome:
		output_file_matrix = output_matrix + file_prefix + ".exome"
		output_file_matrix_tsb = output_matrix + file_prefix_tsb + ".exome"
	else:
		if bed:
			output_file_matrix = output_matrix + file_prefix + ".region"
			output_file_matrix_tsb = output_matrix + file_prefix_tsb + ".region"
		else:
			output_file_matrix = output_matrix + file_prefix + ".all"
			output_file_matrix_tsb = output_matrix + file_prefix_tsb + ".all"

	if initial_chrom != None:
		output_file_matrix += ".chr" + initial_chrom
		output_file_matrix_tsb += ".chr" + initial_chrom


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


	if plot:
		output_path = ref_dir + "/plots/" + project + "/"
		sigPlt.plotINDEL(output_file_matrix, output_path, False, project, False)

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
	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)
	for cont in all_dinucs.keys():
		mutation_types = all_mut_types[cont]
		dinucs = all_dinucs[cont]
		file_prefix = project + ".DBS" + cont
		if exome:
			output_file_matrix = output_matrix + file_prefix + ".exome"
		else:
			if bed:
				output_file_matrix = output_matrix + file_prefix + ".region"
			else:
				output_file_matrix = output_matrix + file_prefix + ".all"
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
						count += dinucs[sample][dinuc]
						rev_dinuc = revbias(dinuc[0]) + dinuc[1:]
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
						count += dinucs[sample][dinuc]
						rev_dinuc = revbias(dinuc[0]) + dinuc[1:]
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
			output_path = ref_dir + "/plots/" + project + "/"
			sigPlt.plotDINUC(output_file_matrix, output_path, False, project, False)

def main():
	start = time.time()

	ncbi_chrom = {'NC_000067.6':'1', 'NC_000068.7':'2', 'NC_000069.6':'3', 'NC_000070.6':'4', 
				  'NC_000071.6':'5', 'NC_000072.6':'6', 'NC_000073.6':'7', 'NC_000074.6':'8',
				  'NC_000075.6':'9', 'NC_000076.6':'10', 'NC_000077.6':'11', 'NC_000078.6':'12',
				  'NC_000079.6':'13', 'NC_000080.6':'14', 'NC_000081.6':'15', 'NC_000082.6':'16', 
				  'NC_000083.6':'17', 'NC_000084.6':'18', 'NC_000085.6':'19', 'NC_000086.7':'X', 
				  'NC_000087.7':'Y'}

	tsb_ref = {0:['N','A'], 1:['N','C'], 2:['N','G'], 3:['N','T'],
			   4:['T','A'], 5:['T','C'], 6:['T','G'], 7:['T','T'],
			   8:['U','A'], 9:['U','C'], 10:['U','G'], 11:['U','T'],
			   12:['B','A'], 13:['B','C'], 14:['B','G'], 15:['B','T'],
			   16:['N','N'], 17:['T','N'], 18:['U','N'], 19:['B','N']}

	contexts = ['3072', 'DINUC']

	exome = False
	indel = False
	limited_indel = False
	functionFlag = False
	bed = False
	bed_file = None
	bed_ranges = None
	chrom_based = False
	plot = False
	matrix_suffix = ''
	SNVs = False
	gs = False

	parser = argparse.ArgumentParser(description="Provide the necessary arguments to create the desired catalogue.")
	parser.add_argument("--project", "-p",help="Provide a unique name for your samples. (ex: BRCA)")
	parser.add_argument("--genome", "-g",help="Provide a reference genome. (ex: GRCh37, GRCh38, mm10)")
	parser.add_argument("-e", "--exome", help="Optional parameter instructs script to create the catalogues using only the exome regions. Whole genome context by default", action='store_true')
	parser.add_argument( "-snv", "--SNV", help="Optional parameter instructs script to create the catalogue for SNVs", action='store_true')
	parser.add_argument( "-i", "--indel", help="Optional parameter instructs script to create the catalogue for limited INDELs", action='store_true')
	parser.add_argument( "-ie", "--extended_indel", help="Optional parameter instructs script to create the catalogue for extended INDELs", action='store_true')
	parser.add_argument("-b", "--bed", nargs='?', help="Optional parameter instructs script to simulate on a given set of ranges (ex: exome). Whole genome context by default")
	parser.add_argument("-ch", "--chromosome", help="Optional parameter instructs script to generate the matrices per chromoosome", action='store_true')
	parser.add_argument("-pl", "--plot", help="Optional parameter instructs script to generate the plots for 96, 192, INDEL, and DINUC contexts", action='store_true')
	parser.add_argument("-gs", "--geneStrandBias", help="Optional parameter instructs script to generate the strand bias per gene", action='store_true')

	args=parser.parse_args()
	project = args.project
	genome = args.genome

	if args.exome:
		exome = True

	if args.extended_indel:
		indel = True

	if args.indel:
		indel = True
		limited_indel = True

	if args.bed:
		bed = True
		bed_file = args.bed

	if args.chromosome:
		chrom_based = True

	if args.plot:
		limited_indel = True
		indel = True
		plot = True

	if args.SNV:
		SNVs = True

	if args.geneStrandBias:
		gs = True


	# Organizes all of the reference directories for later reference:
	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)
	chrom_path =ref_dir + '/references/chromosomes/tsb/' + genome + "/"
	transcript_path = ref_dir + '/references/chromosomes/transcripts/' + genome + "/"
	
	time_stamp = datetime.date.today()
	output_logs = ref_dir + "/logs/"
	if not os.path.exists(output_logs):
		os.mkdir(output_logs)
	error_file = output_logs + 'SigProfilerMatrixGenerator_' + project + "_" + genome + "_" + str(time_stamp) + ".err"
	log_file = output_logs + 'SigProfilerMatrixGenerator_' + project + "_" + genome + "_" + str(time_stamp) + ".out"
	if os.path.exists(error_file):
		os.system("rm " + error_file)
	if os.path.exists(log_file):
		 os.system("rm " + log_file)
	sys.stderr = open(error_file, 'w')
	logging.basicConfig(filename=log_file, level=logging.INFO)


	# Organizes all of the input and output directories:
	output_matrix = ref_dir + "/references/matrix/" + project + "/"
	vcf_path = ref_dir + '/references/vcf_files/' + project + "/"
	bed_path = ref_dir + '/references/vcf_files/BED/' + project + "/"
	if not os.path.exists(output_matrix):
		os.system("mkdir " + output_matrix)

	# Gathers all of the vcf files:
	if SNVs:
		vcf_files_snv_temp = os.listdir(vcf_path + "SNV/")
	if indel:
		vcf_files_indel_temp = os.listdir(vcf_path + "INDEL/")

	vcf_files2 = [[],[]]

	if SNVs:
		for file in vcf_files_snv_temp:
			# Skips hidden files
			if file[0:3] == '.DS':
				pass
			else:
				vcf_files2[0].append(file)
	if indel:
		for file in vcf_files_indel_temp:
			# Skips hidden files
			if file[0:3] == '.DS':
				pass
			else:
				vcf_files2[1].append(file)   

	for i in range(0, len(vcf_files2), 1):
		if i == 1 and indel:
			contexts = ['INDEL']
		elif i ==1 and not indel:
			break
		elif i == 0 and not SNVs:
			continue
		vcf_path = ref_dir + '/references/vcf_files/' + project + "/"
		file_name = vcf_files2[i][0].split(".")
		file_extension = file_name[-1]

		unique_folder = project + str(uuid.uuid4())
		output_path = ref_dir + "/references/vcf_files/" + unique_folder + "/"
		if os.path.exists(output_path):
			os.system("rm -r " + output_path)

		os.makedirs(output_path)


		
		# Converts the input files into a temporary, single simple text file
		if file_extension == 'genome':
			if i == 1:
				convertIn.convertTxt(project, vcf_path + "INDEL/", genome, output_path, 'INDEL')
			else:
				convertIn.convertTxt(project, vcf_path + "SNV/", genome, output_path, 'SNV')
		else:
			if i == 1:
				if file_extension == 'txt':
					convertIn.convertTxt(project, vcf_path + "INDEL/",  genome, output_path, 'INDEL')
				elif file_extension == 'vcf':
					convertIn.convertVCF(project, vcf_path + "INDEL/", genome, output_path, 'INDEL')
				elif file_extension == 'maf':
					convertIn.convertMAF(project, vcf_path + "INDEL/", genome, output_path, 'INDEL')
				elif file_extension == 'tsv':
					convertIn.convertICGC(project, vcf_path + "INDEL/", genome, output_path, 'INDEL')
				else:
					print("File format not supported")
			else:
				if file_extension == 'txt':
					convertIn.convertTxt(project, vcf_path + "SNV/", genome, output_path, 'SNV')
				elif file_extension == 'vcf':
					convertIn.convertVCF(project, vcf_path + "SNV/", genome, output_path, 'SNV')
				elif file_extension == 'maf':
					convertIn.convertMAF(project, vcf_path + "SNV/", genome, output_path, 'SNV')
				elif file_extension == 'tsv':
					convertIn.convertICGC(project, vcf_path + "SNV/", genome, output_path, 'SNV')
				else:
					print("File format not supported")


				
		vcf_files = os.listdir(output_path)
		vcf_path = output_path

		# Include some kind of a flag for the INDEL option 
		sort_file = vcf_files[0]
		with open(vcf_path + sort_file) as f:
			lines = [tuple(t) for t in [line.strip().split() for line in f]]
			lines = set(lines)
		output = open(vcf_path + sort_file, 'w')

		for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT'].index(x[5]), x[1], int(x[6]))):
			print('\t'.join(line), file=output)
			

		output.close()

		print("Sorting complete...\nDetermining mutation type for each variant, one chromosome at a time. Starting catalogue generation...")

		# Calls the relevant functions for the current type of mutation/context
		if bed:
			bed_file_path = ref_dir + "/references/vcf_files/BED/" + project + "/" + bed_file
			bed_ranges = BED_filtering(bed_file_path)

		for context in contexts:
			if context != 'DINUC' and context != 'INDEL':
				catalogue_generator_single (vcf_path, vcf_files, chrom_path, project, output_matrix, context, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, transcript_path, gs)

			elif context == 'DINUC':
				catalogue_generator_DINUC_single (vcf_path, vcf_files, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref)

			elif context == 'INDEL':
				catalogue_generator_INDEL_single (vcf_path, vcf_files, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, limited_indel, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, transcript_path, gs)

			logging.info("Catalogue for " + context + " context is complete.")
			print("Catalogue for " + context + " context is complete.")
		os.system("rm -r " + vcf_path)

	end = time.time()
	logging.info("Job took " + str(end-start) + " seconds.")
	print("Job took ",str(end-start), " seconds.")



if __name__ == '__main__':
	main()