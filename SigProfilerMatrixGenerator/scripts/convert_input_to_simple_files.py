#!/usr/bin/env python3
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

import os
import sys
import re



def convertVCF (project, vcf_path, genome, output_path, ncbi_chrom, log_file):
	'''
	Converts input vcf files into a single simple text format.

	Parameters:
		 project  -> unique name given to the current samples
		vcf_path  -> path to the input vcf files
		  genome  -> reference genome
	 output_path  -> path to the temporary folder 

	Returns:
			 snv  -> Boolean that informs whether there are SNVs present in the 
					 input files
		   indel  -> Boolean that informs whether there are INDELs present
					 in the input files

	Ouput:
		Saves a single text file to the temporary file folder

	'''

	# Collect all input file names and instantiate flags
	files = os.listdir(vcf_path)
	first_indel = True
	first_SNV = True
	snv = False
	indel = False
	first_incorrect_file = True
	out = open(log_file, 'a')
	prev_line = None
	first_chrom = ''
	
	# Iterates through each file 
	for file in files:
		file_name = file.split(".")
		sample = file_name[0]
		if file == '.DS_Store':
			continue
		with open (vcf_path + file) as f:
			for lines in f:
				# Skips any header lines
				if lines[0] == "#":
					continue
				else:
					try:
						line = lines.strip().split('\t')
						chrom = line[0]
						if len(chrom) > 2:
							chrom = chrom[3:]
						if chrom in ncbi_chrom:
							chrom = ncbi_chrom[chrom]

						start = line[1]
						ref = line[3]
						mut = line[4]
						int(start)


					except:
						if first_incorrect_file:
							print("The given input files do not appear to be in the correct vcf format. Skipping this file: ", file)
							first_incorrect_file = False
						continue

					# Saves SNV mutations into an SNV simple text file
					if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
						snv = True

						if ref not in 'ACGT-':
							print("The ref base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
							out.flush()
							skipped_count += 1
							continue
						
						if mut not in 'ACGT-':
							print("The mutation base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
							out.flush()
							skipped_count += 1
							continue
						
						if ref == mut:
							print("The ref base appears to match the mutated base. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
							out.flush()
							skipped_count += 1
							continue

						if line == prev_line:
							print("There appears to be a duplicate single base substitution. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
							out.flush()
							skipped_count += 1
							continue

						if first_SNV:
							if not os.path.exists(output_path + "SNV/"):
								os.mkdir(output_path + 'SNV/')
							outputFile_X = open(output_path + "SNV/X_" + project + ".genome", "w")
							outputFile_Y = open(output_path + "SNV/Y_" + project + ".genome", "w")
							outputFile_1 = open(output_path + "SNV/1_" + project + ".genome", "w")
							outputFile_2 = open(output_path + "SNV/2_" + project + ".genome", "w")
							outputFile_3 = open(output_path + "SNV/3_" + project + ".genome", "w")
							outputFile_4 = open(output_path + "SNV/4_" + project + ".genome", "w")
							outputFile_5 = open(output_path + "SNV/5_" + project + ".genome", "w")
							outputFile_6 = open(output_path + "SNV/6_" + project + ".genome", "w")
							outputFile_7 = open(output_path + "SNV/7_" + project + ".genome", "w")
							outputFile_8 = open(output_path + "SNV/8_" + project + ".genome", "w")
							outputFile_9 = open(output_path + "SNV/9_" + project + ".genome", "w")
							outputFile_10 = open(output_path + "SNV/10_" + project + ".genome", "w")
							outputFile_11 = open(output_path + "SNV/11_" + project + ".genome", "w")
							outputFile_12 = open(output_path + "SNV/12_" + project + ".genome", "w")
							outputFile_13 = open(output_path + "SNV/13_" + project + ".genome", "w")
							outputFile_14 = open(output_path + "SNV/14_" + project + ".genome", "w")
							outputFile_15 = open(output_path + "SNV/15_" + project + ".genome", "w")
							outputFile_16 = open(output_path + "SNV/16_" + project + ".genome", "w")
							outputFile_17 = open(output_path + "SNV/17_" + project + ".genome", "w")
							outputFile_18 = open(output_path + "SNV/18_" + project + ".genome", "w")
							outputFile_19 = open(output_path + "SNV/19_" + project + ".genome", "w")
							outputFile_20 = open(output_path + "SNV/20_" + project + ".genome", "w")
							outputFile_21 = open(output_path + "SNV/21_" + project + ".genome", "w")
							outputFile_22 = open(output_path + "SNV/22_" + project + ".genome", "w")

							outFiles = {'X': outputFile_X, 'Y':outputFile_Y, '1':outputFile_1, '2':outputFile_2, '3':outputFile_3, '4':outputFile_4,
										'5':outputFile_5, '6':outputFile_6, '7':outputFile_7, '8':outputFile_8, '9':outputFile_9, '10':outputFile_10,
										'11':outputFile_11, '12':outputFile_12, '13':outputFile_13, '14':outputFile_14, '15':outputFile_15, '16':outputFile_16,
										'17':outputFile_17, '18':outputFile_18, '19':outputFile_19, '20':outputFile_20, '21':outputFile_21, '22':outputFile_22}
							#out_snv = open(outputFile, "w") 
							first_SNV = False

						#print("\t".join([project, sample, ".", genome, "SNV", chrom, start, start, ref, mut, "SOMATIC"]), file=out_snv)
						if chrom in outFiles:
							print("\t".join([sample, chrom, start, ref, mut]), file=outFiles[chrom])
						else:
							print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...", file=out)
							out.flush()



					# Saves INDEL mutations into an INDEL simple text file
					else:
						indel = True
						if first_indel:
							outputFile = output_path + "INDEL/" + project + ".genome"
							os.system("rm -f " + outputFile)
							if not os.path.exists(output_path + "INDEL/"):
								os.mkdir(output_path + "INDEL/")
							out_indel = open(outputFile, "w") 
							first_indel = False

						#print("\t".join([project, sample, ".", genome, "INDEL", chrom, start, start, ref, mut, "SOMATIC"]), file=out_indel)
						print("\t".join([project, sample, ".", genome, "I", chrom, start, start, ref, mut]), file=out_indel)

					prev_line = line


		first_incorrect_file = True

	# Closes the output files and returns the boolean flags
	if snv:
		for files in outFiles.values():
			files.close()
		#out_snv.close()
	if indel:
		out_indel.close()
	out.close()
	return(snv, indel)


def convertTxt (project, vcf_path, genome, output_path, ncbi_chrom, log_file):
	'''
	Converts input text files into a single simple text format.

	Parameters:
		 project  -> unique name given to the current samples
		vcf_path  -> path to the input text files
		  genome  -> reference genome
	 output_path  -> path to the temporary folder 

	Returns:
			 snv  -> Boolean that informs whether there are SNVs present in the 
					 input files
		   indel  -> Boolean that informs whether there are INDELs present
					 in the input files

	Ouput:
		Saves a single text file to the temporary file folder

	'''

	# Collect all input file names and instantiate flags
	out = open(log_file, 'a')
	files = os.listdir(vcf_path)
	first_indel = True
	first_SNV = True	
	snv = False
	indel = False
	first_incorrect_file = True
	prev_line = None

	# Iterates through each file 
	for file in files:
		if file == '.DS_Store':
			continue
		with open (vcf_path + file) as f:
			for lines in f:
				try:
					line = lines.strip().split('\t')
					sample = line[1]
					genome = line[3]
					chrom = line[5]
					if len(chrom) > 2:
						chrom = chrom[3:]
					if chrom in ncbi_chrom:
						chrom = ncbi_chrom[chrom]

					start = line[6]
					end = line[7]
					ref = line[8]
					mut = line[9]
					int(start)
					int(end)



				except:
					if first_incorrect_file:
						print("The given input files do not appear to be in the correct simple text format. Skipping this file: ", file)
						first_incorrect_file = False
					continue

				# Saves SNV mutations into an SNV simple text file
				if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
					snv = True
					if ref not in 'ACGT-':
						print("The ref base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue
					
					if mut not in 'ACGT-':
						print("The mutation base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue
					
					if ref == mut:
						print("The ref base appears to match the mutated base. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue

					if line == prev_line:
						print("There appears to be a duplicate single base substitution. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue

					if first_SNV:
						outputFile = output_path + "SNV/" + project + ".genome"
						os.system("rm -f " + outputFile)
						if not os.path.exists(output_path + "SNV/"):
							os.mkdir(output_path + 'SNV/')
						out_snv = open(outputFile, "w") 
						first_SNV = False

					print("\t".join([project, sample, ".", genome, "SNV", chrom, start, start, ref, mut, "SOMATIC"]), file=out_snv)
				
				# Saves INDEL mutations into an INDEL simple text file
				else:
					indel = True
					if first_indel:
						outputFile = output_path + "INDEL/" + project + ".genome"
						os.system("rm -f " + outputFile)
						if not os.path.exists(output_path + "INDEL/"):
							os.mkdir(output_path + "INDEL/")
						out_indel = open(outputFile, "w") 
						first_indel = False

					print("\t".join([project, sample, ".", genome, "INDEL", chrom, start, start, ref, mut, "SOMATIC"]), file=out_indel)
				
				prev_line = line
		first_incorrect_file = True

	# Closes the output files and returns the boolean flags
	if snv:
		out_snv.close()
	if indel:
		out_indel.close()
	out.close()
	return(snv, indel)

def convertMAF (project, vcf_path, genome, output_path, ncbi_chrom, log_file):
	'''
	Converts input MAF files into a single simple text format.

	Parameters:
		 project  -> unique name given to the current samples
		vcf_path  -> path to the input MAF files
		  genome  -> reference genome
	 output_path  -> path to the temporary folder 

	Returns:
			 snv  -> Boolean that informs whether there are SNVs present in the 
					 input files
		   indel  -> Boolean that informs whether there are INDELs present
					 in the input files

	Ouput:
		Saves a single text file to the temporary file folder

	'''

	# Collect all input file names and instantiate flags
	out = open(log_file, 'a')
	files = os.listdir(vcf_path)
	first_indel = True
	first_SNV = True	
	snv = False
	indel = False
	first_incorrect_file = True
	prev_line = None

	# Iterates through each file 
	for file in files:
		if file == '.DS_Store':
			continue
		name = file.split(".")
		sample = name[0]
		with open (vcf_path + file) as f:
			for lines in f:
				try:
					line = lines.strip().split('\t')
					genome = line[3]
					chrom = line[4]
					if len(chrom) > 2:
						chrom = chrom[3:]
					if chrom in ncbi_chrom:
						chrom = ncbi_chrom[chrom]

					start = line[5]
					end = line[6]
					ref = line[10]
					mut = line[47]
					int(start)
					int(end)


				except:
					if first_incorrect_file:
						print("The given input files do not appear to be in the correct MAF format. Skipping this file: ", file)
						first_incorrect_file = False
					continue

				# Saves SNV mutations into an SNV simple text file
				if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
					snv = True
					if ref not in 'ACGT-':
						print("The ref base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue
					
					if mut not in 'ACGT-':
						print("The mutation base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue
					
					if ref == mut:
						print("The ref base appears to match the mutated base. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue

					if line == prev_line:
						print("There appears to be a duplicate single base substitution. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue

					if first_SNV:
						outputFile = output_path + "SNV/" + project + ".genome"
						os.system("rm -f " + outputFile)
						if not os.path.exists(output_path + "SNV/"):
							os.mkdir(output_path + 'SNV/')
						out_snv = open(outputFile, "w") 
						first_SNV = False

					print("\t".join([project, sample, ".", genome, "SNV", chrom, start, start, ref, mut, "SOMATIC"]), file=out_snv)

				# Saves INDEL mutations into an INDEL simple text file				
				else:
					indel = True
					if first_indel:
						outputFile = output_path + "INDEL/" + project + ".genome"
						os.system("rm -f " + outputFile)
						if not os.path.exists(output_path + "INDEL/"):
							os.mkdir(output_path + "INDEL/")
						out_indel = open(outputFile, "w") 
						first_indel = False

					print("\t".join([project, sample, ".", genome, "INDEL", chrom, start, start, ref, mut, "SOMATIC"]), file=out_indel)
				
				prev_line = line
		first_incorrect_file = True

	# Closes the output files and returns the boolean flags
	if snv:
		out_snv.close()
	if indel:
		out_indel.close()
	out.close()
	return(snv, indel)


def convertICGC (project, vcf_path, genome, output_path, ncbi_chrom, log_file):
	'''
	Converts input ICGC files into a single simple text format.

	Parameters:
		 project  -> unique name given to the current samples
		vcf_path  -> path to the input ICGC files
		  genome  -> reference genome
	 output_path  -> path to the temporary folder 

	Returns:
			 snv  -> Boolean that informs whether there are SNVs present in the 
					 input files
		   indel  -> Boolean that informs whether there are INDELs present
					 in the input files

	Ouput:
		Saves a single text file to the temporary file folder

	'''

	# Collect all input file names and instantiate flags
	out = open(log_file, 'a')
	files = os.listdir(vcf_path)
	first_indel = True
	first_SNV = True	
	snv = False
	indel = False
	first_incorrect_file = True
	prev_line = None

	# Iterates through each file 
	for file in files:
		if file == '.DS_Store':
			continue
		with open (vcf_path + file) as f:
			for lines in f:
				try:
					line = lines.strip().split('\t')
					sample = line[1]
					icgc_sample_id = line[4]
					chrom = line[8]
					if len(chrom) > 2:
						chrom = chrom[3:]
					if chrom in ncbi_chrom:
						chrom = ncbi_chrom[chrom]

					start = line[9]
					end = line[10]
					genome = line[12]
					ref = line[15]
					mut = line[16]
					if ref == '-':
						mut = '-' + mut
					elif mut == '-':
						start -= 1
						ref = '-' + ref
					int(start)
					int(end)



				except:
					if first_incorrect_file:
						print("The given input files do not appear to be in the correct ICGC format.")
						first_incorrect_file = False
					continue

				# Saves SNV mutations into an SNV simple text file
				if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
					snv = True
					if ref not in 'ACGT-':
						print("The ref base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue
					
					if mut not in 'ACGT-':
						print("The mutation base is not recognized. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue
					
					if ref == mut:
						print("The ref base appears to match the mutated base. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue

					if line == prev_line:
						print("There appears to be a duplicate single base substitution. Skipping this mutation: " + chrom + " " + str(start) + " " + ref + " " + mut, file=out)
						out.flush()
						skipped_count += 1
						continue

					if first_SNV:
						outputFile = output_path + "SNV/" + project + ".genome"
						os.system("rm -f " + outputFile)
						if not os.path.exists(output_path + "SNV/"):
							os.mkdir(output_path + 'SNV/')
						out_snv = open(outputFile, "w") 
						first_SNV = False

					print("\t".join([project, sample, ".", genome, "SNV", chrom, start, start, ref, mut, "SOMATIC"]), file=out_snv)
				
				# Saves INDEL mutations into an INDEL simple text file				
				else:
					indel = True
					if first_indel:
						outputFile = output_path + "INDEL/" + project + ".genome"
						os.system("rm -f " + outputFile)
						if not os.path.exists(output_path + "INDEL/"):
							os.mkdir(output_path + "INDEL/")
						out_indel = open(outputFile, "w") 
						first_indel = False

					print("\t".join([project, sample, ".", genome, "INDEL", chrom, start, start, ref, mut, "SOMATIC"]), file=out_indel)
				prev_line = line

		first_incorrect_file = True

	# Closes the output files and returns the boolean flags
	if snv:
		out_snv.close()
	if indel:
		out_indel.close()
	out.close()
	return(snv, indel)



