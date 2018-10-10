#!/usr/bin/env python3

#This file is part of Mutational Signatures Project.

#Mutational Signatures Project: need info on project

#Copyright (C) 2018 Erik Bergstrom

#

#Mutational Signatures is free software: need to include distribtution

#rights and so forth

#

#Mutational Signatures is distributed in the hope that it will be useful,

#but WITHOUT ANY WARRANTY; without even the implied warranty of

#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

#GNU General Public License for more details [change to whatever group we should include.
 

#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

import os
import sys
import re



def convertVCF (project, vcf_path, genome, mutType=None):
	if mutType == None:
		mutType = 'SNP'

	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)

	output_path = ref_dir + "/references/vcf_files/single/"
	outputFile = ref_dir + "/references/vcf_files/single/" + project + "_indels.genome"
	os.system("rm -f " + outputFile)

	if not os.path.exists(output_path):
		os.mkdir(output_path)

	out = open(outputFile, "w") 

	files = os.listdir(vcf_path)
	for file in files:
		with open (vcf_path + file) as f:
			for lines in f:
				line = lines.strip().split()
				chrom = line[0]
				start = line[1]
				sample = line[2]
				ref = line[3]
				mut = line[4]

				print("\t".join([project, sample, ".", genome, mutType, chrom, start, start, ref, mut, "SOMATIC"]), file=out)

	out.close()


def convertTxt (project, vcf_path, genome, mutType=None):
	if mutType == None:
		mutType = 'SNP'

	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)

	output_path = ref_dir + "/references/vcf_files/single/"
	outputFile = ref_dir + "/references/vcf_files/single/" + project + "_indels.genome"
	os.system("rm -f " + outputFile)

	if not os.path.exists(output_path):
		os.mkdir(output_path)

	out = open(outputFile, "w") 

	files = os.listdir(vcf_path)
	for file in files:
		with open (vcf_path + file) as f:
			for lines in f:
				line = lines.strip().split()
				sample = line[1]
				genome = line[3]
				#mutType = line[4]
				chrom = line[5]
				start = line[6]
				end = line[7]
				ref = line[8]
				mut = line[9]
				print("\t".join([project, sample, ".", genome, mutType, chrom, start, end, ref, mut, "SOMATIC"]), file=out)

	out.close()

def convertMAF (project, vcf_path, genome, mutType=None):
	if mutType == None:
		mutType = 'SNP'

	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)

	output_path = ref_dir + "/references/vcf_files/single/"
	outputFile = ref_dir + "/references/vcf_files/single/" + project + "_indels.genome"
	os.system("rm -f " + outputFile)

	if not os.path.exists(output_path):
		os.mkdir(output_path)

	out = open(outputFile, "w") 

	files = os.listdir(vcf_path)
	for file in files:
		name = file.split(".")
		sample = name[0]
		with open (vcf_path + file) as f:
			for lines in f:
				line = lines.strip().split()
				genome = line[3]
				chrom = line[4]
				start = line[5]
				end = line[6]
				ref = line[10]
				mut = line[47]
				print("\t".join([project, sample, ".", genome, mutType, chrom, start, end, ref, mut, "SOMATIC"]), file=out)

	out.close()


