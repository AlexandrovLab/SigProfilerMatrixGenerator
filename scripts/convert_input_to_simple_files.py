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



def convertVCF (project, vcf_path, genome, output_path, mutType=None):
	if mutType == None:
		mutType = 'SNP'

	outputFile = output_path + project + "_indels.genome"
	os.system("rm -f " + outputFile)

	if not os.path.exists(output_path):
		os.mkdir(output_path)

	out = open(outputFile, "w") 

	files = os.listdir(vcf_path)
	for file in files:
		file_name = file.split(".")
		sample = file_name[0]
		if file == '.DS_Store':
			continue
		with open (vcf_path + file) as f:
			for lines in f:
				if lines[0] == "#":
					continue
				else:
					line = lines.strip().split('\t')
					chrom = line[0]
					if len(chrom) > 2:
						chrom = chrom[3:]
					start = line[1]
					#sample = line[2]
					ref = line[3]
					mut = line[4]

					print("\t".join([project, sample, ".", genome, mutType, chrom, start, start, ref, mut, "SOMATIC"]), file=out)

	out.close()


def convertTxt (project, vcf_path, genome, output_path, mutType=None):
	if mutType == None:
		mutType = 'SNP'

	outputFile = output_path + project + "_indels.genome"
	os.system("rm -f " + outputFile)

	if not os.path.exists(output_path):
		os.mkdir(output_path)

	out = open(outputFile, "w") 

	files = os.listdir(vcf_path)
	for file in files:
		if file == '.DS_Store':
			continue
		with open (vcf_path + file) as f:
			for lines in f:
				line = lines.strip().split('\t')
				sample = line[1]
				genome = line[3]
				chrom = line[5]
				if len(chrom) > 2:
					chrom = chrom[3:]
				start = line[6]
				end = line[7]
				ref = line[8]
				mut = line[9]
				print("\t".join([project, sample, ".", genome, mutType, chrom, start, end, ref, mut, "SOMATIC"]), file=out)

	out.close()

def convertMAF (project, vcf_path, genome, output_path, mutType=None):
	if mutType == None:
		mutType = 'SNP'

	outputFile = output_path + project + "_indels.genome"
	os.system("rm -f " + outputFile)

	if not os.path.exists(output_path):
		os.mkdir(output_path)

	out = open(outputFile, "w") 

	files = os.listdir(vcf_path)
	for file in files:
		if file == '.DS_Store':
			continue
		name = file.split(".")
		sample = name[0]
		with open (vcf_path + file) as f:
			for lines in f:
				line = lines.strip().split('\t')
				genome = line[3]
				chrom = line[4]
				if len(chrom) > 2:
					chrom = chrom[3:]
				start = line[5]
				end = line[6]
				ref = line[10]
				mut = line[47]
				print("\t".join([project, sample, ".", genome, mutType, chrom, start, end, ref, mut, "SOMATIC"]), file=out)

	out.close()

def convertICGC (project, vcf_path, genome, output_path, mutType=None):
	if mutType == None:
		mutType = 'SNP'

	outputFile = output_path + project + "_indels.genome"
	os.system("rm -f " + outputFile)

	if not os.path.exists(output_path):
		os.mkdir(output_path)

	out = open(outputFile, "w") 

	files = os.listdir(vcf_path)
	for file in files:
		if file == '.DS_Store':
			continue
		with open (vcf_path + file) as f:
			for lines in f:
				line = lines.strip().split('\t')
				sample = line[1]
				icgc_sample_id = line[4]
				chrom = line[8]
				if len(chrom) > 2:
					chrom = chrom[3:]
				start = int(line[9])
				end = line[10]
				genome = line[12]
				ref = line[15]
				mut = line[16]
				if ref == '-':
					mut = '-' + mut
				elif mut == '-':
					start -= 1
					ref = '-' + ref
				print("\t".join([project, sample, ".", genome, mutType, chrom, str(start), end, ref, mut, "SOMATIC",icgc_sample_id]), file=out)

	out.close()


