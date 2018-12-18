#!/usr/bin/env python3
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

from . import SigProfilerMatrixGenerator as matGen
import os
import re
import sys
import pandas as pd
import datetime
from SigProfilerMatrixGenerator.scripts import convert_input_to_simple_files as convertIn
import uuid
import shutil
import time

def SigProfilerMatrixGeneratorFunc (project, genome, vcfFiles, exome=False, bed_file=None, chrom_based=False, plot=False, gs=False):
	'''
	Allows for the import of the sigProfilerMatrixGenerator.py function. Returns a dictionary
	with each context serving as the first level of keys. 

	Parameters:
			   	  project  -> unique name given to the current samples
				   genome  -> reference genome 
				 vcfFiles  -> path where the input vcf files are located.
				 	exome  -> flag to use only the exome or not
			  	 bed_file  -> BED file that contains a list of ranges to be used in generating the matrices
			  chrom_based  -> flag to create the matrices on a per chromosome basis
			  		 plot  -> flag to generate the plots for each context
			  		   gs  -> flag that performs a gene strand bias test

	Returns:
			  	 matrices  -> dictionary (nested) of the matrices for each context

		example:
			matrices = {'96': {'PD1001a':{'A[A>C]A':23,
										 'A[A>G]A':10,...},
							  'PD1202a':{'A[A>C]A':23,
										 'A[A>G]A':10,...},...},
						'192':{'PD1001a':{'T:A[A>C]A':23,
										 'T:A[A>G]A':10,...},
							  'PD1202a':{'T:A[A>C]A':23,
										 'T:A[A>G]A':10,...},...},...}
	'''

	# Instantiates all of the required variables and references
	if gs:
		print("The Gene Strand Bias is not yet supported! Continuing with the matrix generation.")
		gs = False
	functionFlag = True
	bed = False
	bed_ranges = None
	limited_indel = True
	exome = exome
	plot = plot


	# Instantiates the final output matrix
	matrices = {'96':None, '1536':None, '384':None, '6144':None, 'DINUC':None, '6':None, '24':None, 'INDEL':None}

	# Provides a chromosome conversion from NCBI notation
	ncbi_chrom = {'NC_000067.6':'1', 'NC_000068.7':'2', 'NC_000069.6':'3', 'NC_000070.6':'4', 
				  'NC_000071.6':'5', 'NC_000072.6':'6', 'NC_000073.6':'7', 'NC_000074.6':'8',
				  'NC_000075.6':'9', 'NC_000076.6':'10', 'NC_000077.6':'11', 'NC_000078.6':'12',
				  'NC_000079.6':'13', 'NC_000080.6':'14', 'NC_000081.6':'15', 'NC_000082.6':'16', 
				  'NC_000083.6':'17', 'NC_000084.6':'18', 'NC_000085.6':'19', 'NC_000086.7':'X', 
				  'NC_000087.7':'Y'}

	# Provides the reference file conversion from binary to base information			  
	tsb_ref = {0:['N','A'], 1:['N','C'], 2:['N','G'], 3:['N','T'],
			   4:['T','A'], 5:['T','C'], 6:['T','G'], 7:['T','T'],
			   8:['U','A'], 9:['U','C'], 10:['U','G'], 11:['U','T'],
			   12:['B','A'], 13:['B','C'], 14:['B','G'], 15:['B','T'],
			   16:['N','N'], 17:['T','N'], 18:['U','N'], 19:['B','N']}

	# Instantiates the initial contexts to generate matrices for
	contexts = ['6144', 'DINUC']

	# Organizes all of the reference directories for later reference:
	current_dir = os.path.realpath(__file__)
	ref_dir = re.sub('\/scripts/SigProfilerMatrixGeneratorFunc.py$', '', current_dir)
	chrom_path =ref_dir + '/references/chromosomes/tsb/' + genome + "/"
	transcript_path = ref_dir + '/references/chromosomes/transcripts/' + genome + "/"

	# Terminates the code if the genome reference files have not been created/installed
	if not os.path.exists(chrom_path):
		print("The specified genome: " + genome + " has not been installed\nRun the following command to install the genome:\n\tpython3 sigProfilerMatrixGenerator/install.py -g " + genome)
		sys.exit()

	# Organizes all of the input and output directories:
	if vcfFiles[-1] != "/":
		vcfFiles += "/"
	vcf_path = vcfFiles + "input/"
	vcf_path_original = vcf_path
	if not os.path.exists(vcf_path) or len(os.listdir(vcf_path)) < 1:
		os.mkdir(vcf_path)
		input_files = os.listdir(vcfFiles)
		input_files.remove("input")
		input_files.remove(".DS_Store")
		for files in input_files:
			shutil.copy(vcfFiles + files, vcf_path + files)
	output_matrix = vcfFiles + "output/"
	if not os.path.exists(output_matrix):
		os.system("mkdir " + output_matrix)


	# Organizes the error and log files
	time_stamp = datetime.date.today()
	output_log_path = vcfFiles + "logs/"
	if not os.path.exists(output_log_path):
		os.system("mkdir " + output_log_path)
	error_file = output_log_path + 'SigProfilerMatrixGenerator_' + project + "_" + genome + str(time_stamp) + ".err"
	log_file = output_log_path + 'SigProfilerMatrixGenerator_' + project + "_" + genome + str(time_stamp) + ".out"

	if os.path.exists(error_file):
		os.system("rm " + error_file)
	if os.path.exists(log_file):
		 os.system("rm " + log_file)
	sys.stderr = open(error_file, 'w')
	#out = open(log_file, 'w')



	# Gathers all of the vcf files:
	vcf_files_temp = os.listdir(vcf_path)
	vcf_files = []
	first_extenstion = True
	for file in vcf_files_temp:
		# Skips hidden files
		if file[0:3] == '.DS' or file[0:2] == '__':
			pass
		else:
			vcf_files.append(file)

	# Creates a temporary folder for sorting and generating the matrices
	file_name = vcf_files[0].split(".")
	file_extension = file_name[-1]
	unique_folder = project + "_"+ str(uuid.uuid4())
	output_path = output_matrix + "temp/" + unique_folder + "/"
	if os.path.exists(output_path):
		os.system("rm -r " + output_path)
	os.makedirs(output_path)

	# Converts the input files to standard text in the temporary folder
	if file_extension == 'genome':
			convertIn.convertTxt(project, vcf_path, genome, output_path)
	else:
		if file_extension == 'txt':
			snv, indel = convertIn.convertTxt(project, vcf_path,  genome,  output_path)
		elif file_extension == 'vcf':
			snv, indel = convertIn.convertVCF(project, vcf_path,  genome, output_path)
		elif file_extension == 'maf':
			snv, indel = convertIn.convertMAF(project, vcf_path,  genome, output_path)
		elif file_extension == '.tsv':
			snv, indel = convertIn.convertICGC(project, vcf_path,  genome, output_path)
		else:
			print("File format not supported")


	# Instantiates variables for final output statistics
	analyzed_muts = [0, 0, 0]
	skipped_muts = 0
	sample_count_high = 0

	# Begins matrix generation for all possible contexts
	for i in range(0, 2, 1):
		if i == 0 and snv:
			output_path_snv = output_path + "SNV/"
			vcf_files = os.listdir(output_path_snv)
			vcf_path = output_path_snv
			print("Starting matrix generation for SNVs and DINUCs...", end='', flush=True)
			start = time.time()

		# Skips SNVs if none are present
		elif i == 0 and not snv:
			continue
		elif i == 1 and indel:
			contexts = ['INDEL']
			output_path_indel = output_path + "INDEL/"
			vcf_files = os.listdir(output_path_indel)
			vcf_path = output_path_indel
			print("Starting matrix generation for INDELs...", end='', flush=True)
			start = time.time()

		# Skips INDELs if none are present and deletes the temp folder
		elif i ==1 and not indel:
			os.system("rm -r " + output_matrix + "temp/")
			continue

		# Removes hidden files generated in macos
		if ".DS_Store" in vcf_files:
			vcf_files.remove(".DS_Store")

		# Sorts files based on chromosome, sample, and start position
		sort_file = vcf_files[0]
		with open(vcf_path + sort_file) as f:
			lines = [line.strip().split() for line in f]
		output = open(vcf_path + sort_file, 'w')
		try:
			for line in sorted(lines, key = lambda x: (['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'M'].index(x[5]), x[1], x[6])):
				print('\t'.join(line), file=output)
		except ValueError as e:
			err = str(e).split("'")
			err = err[1]
			print("There appears to be unrecognized chromosomes within the file. Please remove the lines that contain the chromosome: " + err)
			sys.exit()
		output.close()

		# Generates the bed regions if a bed file was provided
		if bed_file != None:
			bed = True
			bed_file_path = bed_file
			bed_ranges = matGen.BED_filtering(bed_file_path)
		else:
			bed_file_path = None

		# Creates the matrix for each context
		for context in contexts:
			if context != 'DINUC' and context != 'INDEL':
				matrix, skipped_mut, total, sample_count = matGen.catalogue_generator_single (vcf_path, vcf_path_original, vcf_files, bed_file_path, chrom_path, project, output_matrix, context, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, transcript_path, gs, log_file)
				matrices = matrix
				skipped_muts += skipped_mut
				analyzed_muts[0] = total
				if sample_count > sample_count_high:
					sample_count_high = sample_count

			elif context == 'DINUC':
				if analyzed_muts[0] > 1:
					matrix, skipped_mut, total, sample_count = matGen.catalogue_generator_DINUC_single (vcf_path, vcf_path_original, vcf_files, bed_file_path, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, log_file)
					matrices[context] = matrix
					analyzed_muts[1] = total
					if sample_count > sample_count_high:
						sample_count_high = sample_count
				else:
					pass

			elif context == 'INDEL':
				matrix, skipped_mut, total, sample_count = matGen.catalogue_generator_INDEL_single (vcf_path, vcf_path_original, vcf_files, bed_file_path, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, limited_indel, functionFlag, bed, bed_ranges, chrom_based, plot, tsb_ref, transcript_path, gs, log_file)
				matrices[context] = matrix
				skipped_muts += skipped_mut
				analyzed_muts[2] = total
				if sample_count > sample_count_high:
					sample_count_high = sample_count

			#logging.info("Catalogue for " + context + " context is complete.")
			with open(log_file, 'a') as f:
				print("Catalogue for " + context + " context is complete.", file=f)

		if i == 1:
			os.system("rm -r " + output_matrix + "temp/")
		end = time.time() - start
		print("Completed! Elapsed time: " + str(round(end, 2)) + " seconds.")

	# Prints a summary for the given run (total samples, skipped mutations, etc.)
	print("Matrices generated for " + str(sample_count_high) + " samples with " + str(skipped_muts) + " errors. Total of " + str(analyzed_muts[0]) + " SNVs, " + str(analyzed_muts[1]) + " DINUCs, and " + str(analyzed_muts[2]) + " INDELs were successfully analyzed.")
	final_matrices = {}
	for conts in matrices.keys():
		final_matrices[conts] = pd.DataFrame.from_dict(matrices[conts])
		for column in final_matrices[conts].columns:
			final_matrices[conts][column] = final_matrices[conts][column].fillna(0)
	# Returns a dictionary of panda data frames for each matrix context.
	return(final_matrices)

	