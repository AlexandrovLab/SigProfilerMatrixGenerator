#!/usr/bin/env python3

#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

import os
import sys
import re
import subprocess
import argparse
import time
from scipy import spatial
import pandas as pd
import shutil
import logging
from SigProfilerMatrixGenerator.scripts import convert_input_to_simple_files as convertIn

def install_chromosomes (genomes, ref_dir, custom):
	if custom:
		for genome in genomes:
			os.system("gunzip references/chromosomes/fasta/" + genome + "/*.gz")
			chromosome_fasta_path = "references/chromosomes/fasta/" + genome + "/"
			os.system("python3 scripts/save_chrom_strings.py -g " + genome)
			print("Chromosome string files for " + genome + " have been created. Continuing with installation.")
			#os.system("rm -r " + chromosome_fasta_path)
	else:

		for genome in genomes:
			species = None
			chrom_number = None
			if genome == 'GRCh37' or genome == 'GRCh38': 
				species = "homo_sapiens"
				chrom_number = 24
			elif genome == 'mm10' or genome == 'mm9':
				species = "mus_musculus"
				chrom_number = 21
			else:
				print(genome + " is not supported. The following genomes are supported:\nGRCh37, GRCh38, mm10")
			
			chromosome_string_path = "references/chromosomes/chrom_string/" + genome + "/"
			chromosome_fasta_path = "references/chromosomes/fasta/" + genome + "/"

			wget_flag = True
			if os.path.exists(chromosome_string_path) == False or len(os.listdir(chromosome_string_path)) <= chrom_number:
				if os.path.exists(chromosome_fasta_path) == False or len(os.listdir(chromosome_fasta_path)) <= chrom_number:
					print("Chromosomes are not currently saved as individual text files for " + genome + ". Downloading the files now...")
					try:
						p = subprocess.Popen("wget", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					except:
						proceed = input("You may not have wget or homebrew installed. Download those dependencies now?[Y/N]").upper()
						if proceed == 'Y':
							try:
								os.system("brew install wget")
							except:
								os.system('/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"')
								os.system("brew install wget")
						else:
							print("Installation has stopped. Please download the chromosome files before proceeding with the installation.")
							wget_flag = False
							sys.exit()
					if wget_flag:
						try:
							if genome == 'GRCh37':
								os.system("wget -r -l1 -c -nc --no-parent -A '*.dna.chromosome.*' -nd -P " + chromosome_fasta_path + " ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/dna/ 2>> install.log ")
							elif genome == 'mm9':
								os.system("wget -r -l1 -c -nc --no-parent -A '*.dna.chromosome.*' -nd -P " + chromosome_fasta_path + " ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/ 2>> install.log")
							else:
								os.system("wget -r -l1 -c -nc --no-parent -A '*.dna.chromosome.*' -nd -P " + chromosome_fasta_path + " ftp://ftp.ensembl.org/pub/release-93/fasta/"+species+"/dna/ 2>> install.log")
							os.system("gunzip references/chromosomes/fasta/" + genome + "/*.gz")
						except:
							print("The ensembl ftp site is not currently responding.")
							sys.exit()

				print("Chromosome fasta files for " + genome + " have been installed. Creating the chromosome string files now...")
				os.system("python3 scripts/save_chrom_strings.py -g " + genome)
				print("Chromosome string files for " + genome + " have been created. Continuing with installation.")
				os.system("rm -r " + chromosome_fasta_path)

			else:
				print("Chromosome reference files exist for " + genome + ". Continuing with installation.")


def install_chromosomes_tsb (genomes, ref_dir, custom):
	for genome in genomes:
		chrom_string_path = "references/chromosomes/chrom_string/" + genome + "/"
		chrom_number = len(chrom_string_path)

		chromosome_TSB_path = "references/chromosomes/tsb/" + genome + "/"
		transcript_files = "references/chromosomes/transcripts/" + genome + "/"

		if os.path.exists(transcript_files) == False or len(os.listdir(transcript_files)) < 1:
			print("Please download the transcript files before proceeding. You can download the files from 'http://www.ensembl.org/biomart/martview'.")
			print("Follow the format presented in the README file:\n\n\tGene stable ID  Transcript stable ID    Chromosome/scaffold name    Strand  Transcript start (bp)   Transcript end (bp)\n\n\n")
			sys.exit()
		if os.path.exists(chromosome_TSB_path) == False or len(os.listdir(chromosome_TSB_path)) < chrom_number:
			print("The transcriptional reference data for " + genome + " has not been saved. Creating these files now")
			os.system("python3 scripts/save_tsb_192.py -g " + genome)

		print("The transcriptional reference data for " + genome + " has been saved.")

def install_chromosomes_tsb_BED (genomes, custom, ref_dir):
	for genome in genomes:
		if not os.path.exists(ref_dir + "chromosomes/tsb_BED/" + genome + "/") or len(os.listdir(ref_dir + "chromosomes/tsb_BED/" + genome + "/")) < 19:
			os.system("python3 scripts/save_chrom_tsb_separate.py -g " + genome)
			print("The TSB BED files for " + genome + " have been saved.")

def benchmark (ref_dir):
	if os.path.exists("scripts/Benchmark/BRCA_bench/"):
		shutil.move("scripts/Benchmark/BRCA_bench/", "references/vcf_files/")
	current_dir = os.path.realpath(__file__)
	ref_dir = re.sub('\/install.py$', '', current_dir)
	vcf_path = ref_dir + "/references/vcf_files/BRCA_bench/"
	output_path = ref_dir + "/references/matrix/"

	start_time = time.time()
	os.system("python3 " + ref_dir + "/scripts/sigProfilerMatrixGenerator.py -g GRCh37 -p BRCA_bench -vf " + vcf_path + " -op " + output_path + " -snv")
	end_time = time.time()

	original_matrix_96 = ref_dir + "/scripts/Benchmark/BRCA_bench_orig_96.txt"
	original_matrix_3072 = ref_dir + "/scripts/Benchmark/BRCA_bench_orig_3072.txt"
	new_matrix_96 = ref_dir + "/references/matrix/BRCA_bench/SBS/BRCA_bench.SBS96.all"
	new_matrix_3072 = ref_dir + "/references/matrix/BRCA_bench/SBS/BRCA_bench.SBS3072.all"

	genome = "GRCh37"

	############# Cosine Test ###################################################
	data_orig = pd.read_csv(original_matrix_96, sep='\t', header=0)
	data_new = pd.read_csv(new_matrix_96, sep='\t', header=0)
	count = 0
	range_count = min(len(data_orig.loc[0]), len(data_new.loc[0]))
	for i in range (1, range_count, 1):
	    orig_list = list(data_orig[data_orig.columns[i]])
	    new_list = list(data_new[data_new.columns[i]])
	    cosine_result = (1-spatial.distance.cosine(orig_list,new_list))
	    if cosine_result != 1:
	        count += 1
	if count != 0:
	    print("There seems to be some errors in the newly generated matrix. The installation may not have been successful.")


	data_orig = pd.read_csv(original_matrix_3072, sep='\t', header=0)
	data_new = pd.read_csv(new_matrix_3072, sep='\t', header=0)
	count = 0
	range_count = min(len(data_orig.loc[0]), len(data_new.loc[0]))
	for i in range (1, range_count, 1):
	    orig_list = data_orig[data_orig.columns[i]]
	    new_list = data_new[data_new.columns[i]]
	    cosine_result = (1-spatial.distance.cosine(orig_list,new_list))
	    if cosine_result <= 0.85:
	        count += 1
	if count != 0:
	    print("There seems to be some errors in the newly generated matrix. The installation may not have been successful.")

	end_time = time.time()
	print("Installation was succesful.\nSigProfilerMatrixGenerator took " + str(end_time-start_time) + " seconds to complete.")


def install (genome, custom=False):
	print("Beginning installation. This may take up to 20 minutes to complete.")
	first_path = os.getcwd()
	current_dir = os.path.realpath(__file__)
	ref_dir = re.sub('\/install.py$', '', current_dir)
	os.chdir(ref_dir)
	#genomes = ['mm9', 'mm10','GRCh37', 'GRCh38' ]
	#genomes = ['GRCh37']

	genomes = [genome]

	if os.path.exists("install.log"):
		os.system("rm install.log")

	ref_dir = "references/"
	chrom_string_dir = ref_dir + "chromosomes/chrom_string/"
	chrom_fasta_dir = ref_dir + "chromosomes/fasta/"
	chrom_tsb_dir = ref_dir + "chromosomes/tsb/"
	matrix_dir = ref_dir + "matrix/"
	vcf_dir = ref_dir + "vcf_files/"
	bed_dir = ref_dir + "vcf_files/BED/"
	log_dir = "logs/"
	new_dirs = [ref_dir, chrom_string_dir, chrom_fasta_dir, chrom_tsb_dir, matrix_dir, vcf_dir, bed_dir, log_dir]

	current_dir = os.getcwd()
	for dirs in new_dirs:
		if not os.path.exists(dirs):
			os.makedirs(dirs)

	install_chromosomes(genomes, ref_dir, custom)
	install_chromosomes_tsb (genomes, ref_dir, custom)

	if os.path.exists("BRCA_example/"):
		os.system("mv BRCA_example/ references/vcf_files/")
	if os.path.exists("example_test"):
		os.system("mv example_test/ references/vcf_files/")
	if os.path.exists("context_distributions/"):
		os.system("mv context_distributions/ references/chromosomes/")

	print("All reference files have been created.\nVerifying and benchmarking installation now...")
	if genome == "GRCh37":
		benchmark(ref_dir)
	# if os.path.exists(chrom_tsb_dir + "GRCh37/"):
	# 	benchmark(ref_dir)
	print ("To proceed with matrix_generation, please provide the path to your vcf files and an appropriate output path.")
	os.system("rm -r " + chrom_string_dir)
	print("Installation complete.")
	os.chdir(first_path)

def main ():

	first_path= os.getcwd()
	os.chdir(first_path + "/sigProfilerMatrixGenerator/")
	genomes = ['mm9', 'mm10','GRCh37', 'GRCh38' ]
	#genomes = ['GRCh37']
	custom = False
	parser = argparse.ArgumentParser(description="Provide the necessary arguments to install the reference files.")
	parser.add_argument("-g", "--genome", nargs='?', help="Optional parameter instructs script to install the custom genome.")
	parser.add_argument("-ct", "--custom", help="Optional parameter instructs script to create the reference files for a custom genome", action='store_true')
	args = parser.parse_args()

	if args.genome:
		genomes = [args.genome]
	if args.custom:
		custom = True

	if os.path.exists("install.log"):
		os.system("rm install.log")

	ref_dir = "references/"
	chrom_string_dir = ref_dir + "chromosomes/chrom_string/"
	chrom_fasta_dir = ref_dir + "chromosomes/fasta/"
	chrom_tsb_dir = ref_dir + "chromosomes/tsb/"
	matrix_dir = ref_dir + "matrix/"
	vcf_dir = ref_dir + "vcf_files/"
	bed_dir = ref_dir + "vcf_files/BED/"
	log_dir = "logs/"
	new_dirs = [ref_dir, chrom_string_dir, chrom_fasta_dir, chrom_tsb_dir, matrix_dir, vcf_dir, bed_dir, log_dir]

	current_dir = os.getcwd()
	for dirs in new_dirs:
		if not os.path.exists(dirs):
			os.makedirs(dirs)

	#initial_transcripts = os.listdir("transcripts_original/")
	#for file in initial_transcripts:
	#	name = file.split("_")
	#	if not os.path.exists("references/chromosomes/transcripts/"+name[0]+"/"):
	#		os.makedirs("references/chromosomes/transcripts/"+name[0]+"/")
	#	if name != ".DS":
	#		os.system("cp transcripts_original/"+file +" references/chromosomes/transcripts/" + name[0]+"/")

	#if os.path.exists("exome/"):
	#	os.system("mv exome/ references/chromosomes/")

	install_chromosomes(genomes, ref_dir, custom)
	install_chromosomes_tsb (genomes, ref_dir, custom)
	#install_chromosomes_tsb_BED (genomes, custom, ref_dir)
	if os.path.exists("BRCA_example/"):
		os.system("mv BRCA_example/ references/vcf_files/")
	if os.path.exists("example_test"):
		os.system("mv example_test/ references/vcf_files/")
	if os.path.exists("context_distributions/"):
		os.system("mv context_distributions/ references/chromosomes/")

	print("All reference files have been created.\nVerifying and benchmarking installation now...")
	if os.path.exists(chrom_tsb_dir + "GRCh37/"):
		benchmark(ref_dir)
	print ("Please place your vcf files for each sample into the 'references/vcf_files/[test]/[mutation_type]/' directory. Once you have done that, you can proceed with the matrix generation.")
	#os.system("rm -r " + chrom_string_dir)
	print("Installation complete.")
	os.chdir(first_path)

if __name__ == '__main__':
	main()