#!/usr/bin/env python3

#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

from __future__ import print_function
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
import hashlib
from SigProfilerMatrixGenerator.scripts import convert_input_to_simple_files as convertIn
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import glob

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return (hash_md5.hexdigest())


def install_chromosomes (genomes, ref_dir, custom, rsync, bash):
	if custom:
		for genome in genomes:
			os.system("gzip -d references/chromosomes/fasta/" + genome + "/*.gz")
			chromosome_fasta_path = "references/chromosomes/fasta/" + genome + "/"
			os.system("python scripts/save_chrom_strings.py -g " + genome + " -c " + str(custom))
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
			elif genome == 'rn6':
				species = 'rattus_norvegicus'
				chrom_number = 22
			else:
				print(genome + " is not supported. The following genomes are supported:\nGRCh37, GRCh38, mm10")
				sys.exit()
			
			chromosome_string_path = "references/chromosomes/chrom_string/" + genome + "/"
			chromosome_fasta_path = "references/chromosomes/fasta/" + genome + "/"

			if os.path.exists(ref_dir + "chromosomes/tsb/" + genome) and len(os.listdir(ref_dir + "chromosomes/tsb/" + genome)) >= chrom_number:
				break
			wget_flag = True
			if os.path.exists(chromosome_string_path) == False or len(os.listdir(chromosome_string_path)) <= chrom_number:
				print("[DEBUG] Chromosome string files found at: " + ref_dir + chromosome_string_path)
				if os.path.exists(chromosome_fasta_path) == False or len(os.listdir(chromosome_fasta_path)) <= chrom_number:
					print("[DEBUG] Chromosome fasta files found at: " + ref_dir + chromosome_fasta_path)
					print("Chromosomes are not currently saved as individual text files for " + genome + ". Downloading the files now...")
					if not rsync:
					#os.system("rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/grch37/update/fasta/homo_sapiens/dna/ " + chromosome_fasta_path + " 2>&1>> install.log")
						# try:
						# 	p = subprocess.Popen("wget", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
						# except:
						# 	proceed = input("You may not have wget or homebrew installed. Download those dependencies now?[Y/N]").upper()
						# 	if proceed == 'Y':
						# 		try:
						# 			os.system("brew install wget")
						# 		except:
						# 			os.system('/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"')
						# 			os.system("brew install wget")
						# 	else:
						# 		print("Installation has stopped. Please download the chromosome files before proceeding with the installation.")
						# 		wget_flag = False
						# 		sys.exit()
						if wget_flag:
							try:
								if genome == 'GRCh37':
									if bash:
										os.system("bash -c '" + 'wget -r -l1 -c -nc --no-parent -A "*.dna.chromosome.*" -nd -P ' + chromosome_fasta_path + ' ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/ 2>> install.log' + "'")
									else:
										os.system('wget -r -l1 -c -nc --no-parent -A "*.dna.chromosome.*" -nd -P ' + chromosome_fasta_path + ' ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/ 2>> install.log')										
									#os.system("wget -r -l1 -c -nc --no-parent -A '*.dna.chromosome.*' -nd -P " + chromosome_fasta_path + " ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/dna/ 2>> install.log")
								elif genome == 'mm9':
									if bash:
										os.system("bash -c '" + 'wget -r -l1 -c -nc --no-parent -A "*.dna.chromosome.*" -nd -P ' + chromosome_fasta_path + ' ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/ 2>> install.log' + "'")
									else:
										os.system('wget -r -l1 -c -nc --no-parent -A "*.dna.chromosome.*" -nd -P ' + chromosome_fasta_path + ' ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/ 2>> install.log')

								elif genome == 'rn6':
									if bash:
										os.system("bash -c '" + 'wget -r -l1 -c -nc --no-parent -A "*.dna.chromosome.*" -nd -P ' + chromosome_fasta_path + ' ftp://ftp.ensembl.org/pub/release-96/fasta/rattus_norvegicus/dna/ 2>> install.log' + "'")									
									else:
										os.system('wget -r -l1 -c -nc --no-parent -A "*.dna.chromosome.*" -nd -P ' + chromosome_fasta_path + ' ftp://ftp.ensembl.org/pub/release-96/fasta/rattus_norvegicus/dna/ 2>> install.log')									
								else:
									if bash:
										os.system("bash -c '" + 'wget -r -l1 -c -nc --no-parent -A "*.dna.chromosome.*" -nd -P ' + chromosome_fasta_path + ' ftp://ftp.ensembl.org/pub/release-93/fasta/' +species+'/dna/ 2>> install.log' + "'")
									else:
										os.system('wget -r -l1 -c -nc --no-parent -A "*.dna.chromosome.*" -nd -P ' + chromosome_fasta_path + ' ftp://ftp.ensembl.org/pub/release-93/fasta/' +species+'/dna/ 2>> install.log')

								#os.system("gunzip references/chromosomes/fasta/" + genome + "/*.gz")
								os.system("gzip -d references/chromosomes/fasta/" + genome + "/*.gz")

							except:
								print("The ensembl ftp site is not currently responding.")
								sys.exit()
					else:
						try:
							if genome == 'GRCh37':
								if bash:
									os.system("bash -c '" + "rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/grch37/current/fasta/homo_sapiens/dna/ " + chromosome_fasta_path + " 2>&1>> install.log" + "'")
								else:
									os.system("rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/grch37/current/fasta/homo_sapiens/dna/ " + chromosome_fasta_path + " 2>&1>> install.log")
							elif genome == 'mm9':
								if bash:
									os.system("bash -c '" + "rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/release-67/fasta/mus_musculus/dna/ " + chromosome_fasta_path + " 2>&1>> install.log" + "'")
								else:
									os.system("rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/release-67/fasta/mus_musculus/dna/ " + chromosome_fasta_path + " 2>&1>> install.log")									
							elif genome == 'rn6':
								if bash:
									os.system("bash -c '" + "rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/release-96/fasta/rattus_norvegicus/dna/ " + chromosome_fasta_path + " 2>> install.log" + "'")									
								else:
									os.system("rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/release-96/fasta/rattus_norvegicus/dna/ " + chromosome_fasta_path + " 2>> install.log")									
							else:
								if bash:
									os.system("bash -c '" + "rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/release-93/fasta/"+species+"/dna/ " + chromosome_fasta_path + " 2>&1>> install.log" + "'")
								else:
									os.system("rsync -av -m --include='*/' --include='*.dna.chromosome.*' --exclude='*' rsync://ftp.ensembl.org/ensembl/pub/release-93/fasta/"+species+"/dna/ " + chromosome_fasta_path + " 2>&1>> install.log")

							#os.system("gunzip references/chromosomes/fasta/" + genome + "/*.gz")
							os.system("gzip -d references/chromosomes/fasta/" + genome + "/*.gz")

						except:
							print("The ensembl ftp site is not currently responding.")
							sys.exit()

				print("Chromosome fasta files for " + genome + " have been installed. Creating the chromosome string files now...")
				os.system("python scripts/save_chrom_strings.py -g " + genome)
				print("Chromosome string files for " + genome + " have been created. Continuing with installation.")
				# os.system("rm -r " + chromosome_fasta_path)
				# os.remove(chromosome_fasta_path)
				shutil.rmtree(chromosome_fasta_path)

			else:
				print("Chromosome reference files exist for " + genome + ". Continuing with installation.")


def install_chromosomes_tsb (genomes, ref_dir, custom):
	check_sum = {'GRCh37':
							{'1':'a7d51305e943cf06ff2029146bd91bca','2':'d24d0185af89356d44614ab0d6fd6a68','3':'ea5e033147dcaf77bfd4c70f50688d37',
							 '4':'00d7797c7184f1802367e33f6e2bc3da','5':'f74b1eeb329088242a9f22a16322b325','6':'b353cc4c4abc90340e7747509fe7b457',
							 '7':'bbadde91b3ef958c7d20e2b1088c8cd2','8':'0ff695692e3efebaf00c7905d0d536d7','9':'40b75a18acb66748a37888c53a76dcdb',
							 '10':'557881b744a932b4ceee8a06d6de85a4','11':'f8b8f118101d7cb04164b29a7acadca4','12':'52c18e9fefc3ed3e35c1d8771d1247de',
							 '13':'a241d1cdcadccfd94db792300ab000bf','14':'ed3907128336795669bc19d77c0aa409','15':'bfc66ad087c4e9025076d7571cffa30e',
							 '16':'bd251fddc42400bb54ef95d5e1002ece','17':'fcd36b1bf5c4bd74328dc4caaae244ae','18':'e015d4324c36374827582c5b1214a736',
							 '19':'5cfa7d47e2d73dbdbf8d68f97c8e8b23','20':'2fa0717bf4e8dddac64cd393f4134ff5','21':'ba5559776d4601b80ca42c82f02102a4',
							 '22':'ba762b6ae493df40d04d1ff63d9b2933','Y':'0303100be91874b966a998273cd7c8eb','X':'14e331d82736f6cfc177ff8c90f7bd78',
							 'MT':'dfd6db5743d399516d5c8dadee5bee78'},

				'GRCh38':
							{'1':'ebe083105e7703a49581a36d73732a96','2':'cd65e36dbdf12a8ac3d2c70ebac8cad4','3':'6c20a7008394f2fa9c304d231a1f391b',
							 '4':'5c7443e1678868adadeac0e57558f6e8','5':'45573232c8097c679503a6598f61e60b','6':'cfc137c7434d3a9a872332d405b5c553',
							 '7':'9d8210c22c1962db837e7b62a578975c','8':'665134fd44f21915cbeef955addf89ba','9':'758d0c0c71d8bafbe1ede86587191730',
							 '10':'397bb21acff1ca3052ac802f2aee06e0','11':'07707ff8a2a964656469a7be7bb3e576','12':'506d02539075e080ee12ebdf63908080',
							 '13':'03ed22f01ab43145733c0b6a647e0560','14':'8b93447086549e476c65699ed813a567','15':'cd0dfe9fa78cae2fc7becf8f8ec6c693',
							 '16':'e17bbb66eb4d6b62b7b0e2fbf062b6a6','17':'8fc95bb3101d024d890aa3543eb454c5','18':'a4870628045bb033a90e8c89f818e24d',
							 '19':'6a9d0c8298f0ba2fa13180e02b969f16','20':'aa75d35969cf3956bb4ace7bdc57b34e','21':'5d55f5ad6271d6a0d8806876924990f7',
							 '22':'efdb4e1d23ab7964302b828062a33447','Y':'3b38c639ad164d60f1a055b46fcd2748','X':'d5edbea3cf5d1716765dd4a7b41b7656',
							 'MT':'dfd6db5743d399516d5c8dadee5bee78'},

				'mm9':
							{'1':'c5afc4b3f7f2119696214511d7a04341','2':'a7b467475a1b032d2c893dac1c419a28','3':'f922bc529a17324f1cd858f9a8723d65',
							'4':'f3d6b74e3c04dbd229e2f1e363607506','5':'5fee4f1889c9fe20f7f8562c62bbeb0a','6':'481d47b87da45f3a20181c780fd796c2',
							'7':'454ef2bf49a5ba8cfea3d16dfcfc7f25','8':'2f4162d4c824db78a2a2a820cb4fec81','9':'0649e6aec61af1ab8ab4797ea8e54119',
							'10':'38296256bcfe886c8ae771418e4fd824','11':'b31cb0ce693e35eaa77031d44b12e474','12':'d2b3e4b015742b6aea30ceec5a972968',
							'13':'df77b6d0ed1b133224b128c189736372','14':'0ec3c0e6b3fa2cdb957541f19792e130','15':'44fcaf2ec9b82dae910f85ce41c3cfad',
							'16':'ad7a8dbdf46fa7077e0982a54eab70b7','17':'71aee1dee3cd2078e4619c485d88817e','18':'727ec4ed3128ecacd6cd2f7558083553',
							'19':'461a7119781ab7f4b654fdd9ef76e0ec','Y':'471ff3bbb4520c020cfaa7ca8371c543','X':'9ccadf96cd3aa0ed9d299894a3d7fde0',
							'MT':'a1d56043ed8308908965dd080a4d0c8d'},

				'mm10':
							{'1':'ef88c5ac276a32a2865c0408f92acd55','2':'ced7325ef9e2dfedea3fbe26428a6059','3':'9cd1794eeea27553077a018038303908',
							'4':'da616d7ed6c67f824487eb2ed09cd33b','5':'b327b82da6986bf947105d07c0ad6d2e','6':'fb9a8fa0b85561f8d4de633c22d5157a',
							'7':'12457fd80f6806779fc0d4cc8d36fbad','8':'5d98d86bd22bee1cb226406f49ee7caf','9':'b2f26613fcc622a4003e4c945ae55e25',
							'10':'e9f3589529e258ede66d2e77bb87d21d','11':'76bcd285c3c66471ad6fccfabe42294c','12':'ac34fc3616c9609d8e75a59069e9007a',
							'13':'f81b976e4e4617b25945d06f9aa30846','14':'95dc042eb2aa7d4cc0abe071d4d7966e','15':'fbf2477833aff73ae085537cd7ee0f85',
							'16':'77cbcd009ba50891571f785595717ec1','17':'cd9e4dfdd168ed3de05dac4d44c6e692', '18':'945e83694c7c8f69d6186e1a2abc9771',
							'19':'e57b25f8869de31a9dbce06510711db6','Y':'c2146ba4ab1ec262f5e38b2a1ebc5f5b','X':'9af543088be046fdc63976c2d41de94c',
							'MT':'a1d56043ed8308908965dd080a4d0c8d'},
				'rn6':
							{'1':'003723513cbdb3708fcc5d737c05199c','2':'53e52c5facc7f05462be533845f37425','3':'8d157a9b71fe9770cf783ea5459b19d7',
							'4':'a66dc1999bcc960ff11fe0b24c0d7b14','5':'601cf83411234adbdd9f911b89509564','6':'03b1f4af58fffdf213466ea85b570b3d',
							'7':'4ed05ddf9502ef79e121c02e391660e6','8':'3e2458daaf1b3e8ab4d0e0a9e60c067b','9':'8f83caeccec7ea6e35e404737138ee67',
							'10':'9c1af453a5facc9bfa821457bcfc4d30','11':'ef0480a905c55d76a3c58e295a85bc75','12':'643b6fe4a3a6363ffe64a6c316fa3e1a',
							'13':'102bb3fb420a4104c216bcdf99870374','14':'e26b8b63fba0ea7ced4f0330e93a8cdc','15':'da747616a1362d374d4786102fab6f9f',
							'16':'54e4f932eb0eda4cbf31156f96ef7235','17':'46c2facf5415e4eff8b0804161db722d', '18':'f1cb84f002967854b83bf266ec59a7a3',
							'19':'b85ca155fd1780fe5c327a4589c212a6','20':'899d3511352d78b9b9dc63f063d91b31','Y':'6a7a3539c329dc540dfa6db006003bb1',
							'X':'7a06bafab97c59a819f03633f0a6b7a2'},
				'c_elegans':
							{'I':'5a3ea8cf3dfbc641716b7bc805edcaae','II':'bf82edaa92809dd2fea2b791c38c9728','III':'d2df34b6743f41d3964549fc76c5f1a2',
							'IV':'23396bb57145d3acde2888947b5b8c3a','V':'09df3c53b12e5fd7d9035cc98ca221a3','X':'988046456f1409dfdb5e26444d84d238',
							'MtDNA':'48983f530959780de0125f74a87d4fc1'},
				'dog':
							{'1':'bef8283c1a36f9aef0e407de2ff6af00','2':'9cc961192bb5e58b3847060c3e9c1cfc','3':'d33263fa2de6666b41e140cb7a8da66c',
							 '4':'cd4ed39ebac1c04800ccf30466ec69f5','5':'c0f48a4a764e58388b48835aca2ec0a4','6':'4b472a2f8d0a53ac75cce04e7dc9279a',
							 '7':'12a61573a0da2c9306fff705bb1c39c1','8':'e22cf22a27560aa8523dc959ddcf6e25','9':'c079a73d719145cdd5c7c93969a1c392',
							 '10':'45805a518147f7846bd0457ca038c8df','11':'f38cda8508463a7607dff14a581ee7b0','12':'adb5de197f58bb827fa01fe924eb3a1d',
							 '13':'055a845ba97baad3b13d4d3359f88290','14':'27f0ba8e47996a058807a3827cf8e4a8','15':'2e9565c687a593eb0acbdd0962bb9255',
							 '16':'89b2225bb78d88b0fd1d38d9514ab0cb','17':'f0378253e2f083e42b665ea202fde3b0','18':'04d124e273f3b54a685ad6526223cd03',
							 '19':'67bae093919e6bb5ab6b9806c739d539','20':'5588387165a2e19c4533012cfb4998f3','21':'371cdf18a545728f7964b9db2fc72d5e',
							 '22':'fbf76865f88a018d93506e036f6a68bc','23':'085145e01d9fd9f0f999fb9e8e8d4400','24':'69b75a9962fb766b447e7d1252cb31ac',
							 '25':'12d5c6677b3e17170c317c1f5532d2a8','26':'13937d18e56b2b93d12fa5fcba48a138','27':'1d03d8ca5f201f4d156f5e1b38f7a67c',
							 '28':'c33395dec7fdc13e9d8f10afaa946f8c','29':'174f2db104ecaa5efef770f44241e3b0','30':'047d420ef9aecb933a7d83b6af820b23',
							 '31':'5be61f0c9944a5f2d7d1a5b2e75fb000','32':'212dcb867e95a642277a243fed8d8e41','33':'08a217b02cdd778cfdb0005dff4828b1',
							 '34':'4245d6fc370d9049ef4c25314fbef239','35':'1344aba8755b8a4e304629180fc0591a','36':'e4fff6ed84777905dc999ca6d6bc2557',
							 '37':'60d51ea6ae9e3f2fa316e3d03aff96b2','38':'4090ff76d94e6b38920916ae3ff2441c','X':'bce1372df64037d79b0995311d8ff971'},
				
				'ebv':
							{'gi_82503188_ref_NC_007605':'7f8894e52b0cac1f968c0402838bea07'}}


	for genome in genomes:
		chrom_number = None
		if genome == 'GRCh37' or genome == 'GRCh38': 
			chrom_number = 24
		elif genome == 'mm10' or genome == 'mm9':
			chrom_number = 21
		elif genome == 'rn6':
			chrom_number = 22

		if custom:
			chrom_number = len([x for x in os.listdir("references/chromosomes/chrom_string/" + genome + "/") if x != ".DS_Store"])


		chromosome_TSB_path = "references/chromosomes/tsb/" + genome + "/"
		transcript_files = "references/chromosomes/transcripts/" + genome + "/"
		print("[DEBUG] Chromosome tsb files found at: " + ref_dir +  chromosome_TSB_path)


		if os.path.exists(transcript_files) == False or len(os.listdir(transcript_files)) < 1:
			print("Please download the transcript files before proceeding. You can download the files from 'http://www.ensembl.org/biomart/martview'.")
			print("Follow the format presented in the README file:\n\n\tGene stable ID  Transcript stable ID    Chromosome/scaffold name    Strand  Transcript start (bp)   Transcript end (bp)\n\n\n")
			sys.exit()
		if os.path.exists(chromosome_TSB_path) == False or len(os.listdir(chromosome_TSB_path)) < chrom_number:
			print("The transcriptional reference data for " + genome + " has not been saved. Creating these files now")
			os.system("python scripts/save_tsb_192.py -g " + genome)

		if not custom:
			corrupt = False
			for files in os.listdir(chromosome_TSB_path):
				if "proportions" in files:
					continue
				if ".DS_Store" in files:
					continue
				chrom = files.split(".")
				chrom = chrom[0]
				check = md5(chromosome_TSB_path + files)
				if check_sum[genome][chrom] != check:
					corrupt = True
					os.remove(chromosome_TSB_path + files)
					print("[DEBUG] Chromosome " + chrom + " md5sum did not match => reference md5sum: " + str(check_sum[genome][chrom]) + "    new file md5sum: " + str(check))
			if corrupt:
				print("The transcriptional reference data appears to be corrupted. Please reinstall the " + genome + " genome.")
				sys.exit()
			
		print("The transcriptional reference data for " + genome + " has been saved.")

def install_chromosomes_tsb_BED (genomes, custom, ref_dir):
	for genome in genomes:
		if not os.path.exists(ref_dir + "chromosomes/tsb_BED/" + genome + "/") or len(os.listdir(ref_dir + "chromosomes/tsb_BED/" + genome + "/")) < 19:
			os.system("python scripts/save_chrom_tsb_separate.py -g " + genome)
			print("The TSB BED files for " + genome + " have been saved.")

def benchmark (genome, ref_dir):
	#current_dir = os.path.realpath(__file__)
	#ref_dir = re.sub('\/install.py$', '', current_dir)
	ref_dir = os.path.dirname(os.path.abspath(__file__))
	vcf_path = ref_dir + "/references/vcf_files/" + genome + "_bench/"

	start_time = time.time()
	matGen.SigProfilerMatrixGeneratorFunc(genome + "_bench", genome, vcf_path)
	end_time = time.time()

	original_matrix_96 = ref_dir + "/scripts/Benchmark/" + genome + "_bench_orig_96.txt"
	original_matrix_3072 = ref_dir + "/scripts/Benchmark/" + genome + "_bench_orig_3072.txt"
	new_matrix_96 = vcf_path + "output/SBS/" + genome + "_bench.SBS96.all"
	new_matrix_3072 = vcf_path + "output/SBS/" + genome + "_bench.SBS6144.all"

	#genome = "GRCh37"

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


def install (genome, custom=False, rsync=False, bash=True, ftp=True, fastaPath=None, transcriptPath=None, exomePath=None):
	if custom:
		ftp=False
	first_path= os.getcwd()
	ref_dir = os.path.dirname(os.path.abspath(__file__))
	os.chdir(ref_dir)

	if not custom:
		wget_install = shutil.which("wget") is not None
		rsync_install = shutil.which("rsync") is not None
		if not wget_install and not rsync_install:
			print("You do not have wget nor rsync installed. Please install wget or rsync and then reattempt to install your desired genome.")
			sys.exit()
		elif not wget_install and rsync_install:
			print("You do not have wget installed. Please install wget and then reattempt to install your desired genome.")
			sys.exit()
	tar_install = shutil.which("tar") is not None
	if not tar_install:
		print("You do not have tar installed. Please install tar and then reattempt to install your desired genome.")
		sys.exit()


	if os.path.exists("install.log"):
		os.remove("install.log")

	chrom_string_dir = ref_dir + "/references/chromosomes/chrom_string/"
	chrom_fasta_dir = ref_dir + "/references/chromosomes/fasta/"
	chrom_tsb_dir = ref_dir + "/references/chromosomes/tsb/"
	matrix_dir = ref_dir + "/references/matrix/"
	vcf_dir = ref_dir + "/references/vcf_files/"
	bed_dir = ref_dir + "/references/vcf_files/BED/"
	log_dir = "logs/"
	new_dirs = [ref_dir, chrom_string_dir, chrom_fasta_dir, chrom_tsb_dir, matrix_dir, vcf_dir, bed_dir, log_dir]


	for dirs in new_dirs:
		if not os.path.exists(dirs):
			os.makedirs(dirs)


	if custom:
		transcript_files = "references/chromosomes/transcripts/" + genome + "/"
		if os.path.exists(chrom_fasta_dir + genome + "/"):
			shutil.rmtree(chrom_fasta_dir + genome + "/")
		os.makedirs(chrom_fasta_dir + genome + "/")
		fastaFiles = glob.iglob(os.path.join(fastaPath, "*.gz"))
		for file in fastaFiles:
			if os.path.isfile(file):
				shutil.copy(file, chrom_fasta_dir + genome + "/")

		if os.path.exists(ref_dir + "/references/chromosomes/transcripts/" + genome + "/"):
			shutil.rmtree(ref_dir + "/references/chromosomes/transcripts/" + genome + "/")
		os.makedirs(ref_dir + "/references/chromosomes/transcripts/" +  genome + "/")
		shutil.copy(transcriptPath,ref_dir + "/references/chromosomes/transcripts/" +  genome + "/")

		if os.path.exists(ref_dir + "/references/chromosomes/exome/" + genome + "/"):
			shutil.rmtree(ref_dir + "/references/chromosomes/exome/" + genome + "/")
		os.makedirs(ref_dir + "/references/chromosomes/exome/" +  genome + "/")
		shutil.copy(exomePath,ref_dir + "/references/chromosomes/exome/" +  genome + "/")


	if ftp:
		check_sum = {'GRCh37':
								{'1':'a7d51305e943cf06ff2029146bd91bca','2':'d24d0185af89356d44614ab0d6fd6a68','3':'ea5e033147dcaf77bfd4c70f50688d37',
								 '4':'00d7797c7184f1802367e33f6e2bc3da','5':'f74b1eeb329088242a9f22a16322b325','6':'b353cc4c4abc90340e7747509fe7b457',
								 '7':'bbadde91b3ef958c7d20e2b1088c8cd2','8':'0ff695692e3efebaf00c7905d0d536d7','9':'40b75a18acb66748a37888c53a76dcdb',
								 '10':'557881b744a932b4ceee8a06d6de85a4','11':'f8b8f118101d7cb04164b29a7acadca4','12':'52c18e9fefc3ed3e35c1d8771d1247de',
								 '13':'a241d1cdcadccfd94db792300ab000bf','14':'ed3907128336795669bc19d77c0aa409','15':'bfc66ad087c4e9025076d7571cffa30e',
								 '16':'bd251fddc42400bb54ef95d5e1002ece','17':'fcd36b1bf5c4bd74328dc4caaae244ae','18':'e015d4324c36374827582c5b1214a736',
								 '19':'5cfa7d47e2d73dbdbf8d68f97c8e8b23','20':'2fa0717bf4e8dddac64cd393f4134ff5','21':'ba5559776d4601b80ca42c82f02102a4',
								 '22':'ba762b6ae493df40d04d1ff63d9b2933','Y':'0303100be91874b966a998273cd7c8eb','X':'14e331d82736f6cfc177ff8c90f7bd78',
								 'MT':'dfd6db5743d399516d5c8dadee5bee78'},
					'GRCh37_havana':
								{'1':'a954558ed4a92c33454f8539d5e47d43','2':'f2a1702c668d531656ffbeadb5bb1b66','3':'e6a4d10c433d8707a1a336e310a779bc',
								 '4':'709d6d8ab80b8f83f2f8892e23a9adba','5':'fa4fe82dce241cf4eeaf84ef9b44030f','6':'161d68063347fff0fe5d82a1c60cdef4',
								 '7':'5f866a10820d3d48d3b96fb8432ade8e','8':'585c7100d1de569df883910bbf01463c','9':'5667fc1a44a93b6f162bf60d4f6a21c1',
								 '10':'7684e3198ac4c484b56d6e91b8f2fc15','11':'e3b5501a11be7d5fcf9219933467c86a','12':'511a8b89d19e47df4633da5f277373b5',
								 '13':'7b621f0c18db0592c0f82db9b760ceef','14':'16b383fa1dbf30feb8a5d1d89a9ac7af','15':'f3295833a4dfd7a628ca9befe1972739',
								 '16':'9985d8cf725dcd29b3ded0a4eeb426e2','17':'9118e1200b7cd95a4d5c95ff0e039c7c','18':'116ae0ba70ea43b10983ddd07387d6dd',
								 '19':'930c568fe249b78a602a3400ce094b83','20':'fa24f27c160dffd5d8b7c6270dff4a0d','21':'4f6f61d43c698e2c66dc62f21f7ad6cd',
								 '22':'77661703e7b65ad2b6adcb75e7e3e53b','Y':'b86042fd443490fb0061478037392fc0','X':'02b7328d7d74704d571fd38149bbf814',
								 },

					'GRCh38':
								{'1':'ebe083105e7703a49581a36d73732a96','2':'cd65e36dbdf12a8ac3d2c70ebac8cad4','3':'6c20a7008394f2fa9c304d231a1f391b',
								 '4':'5c7443e1678868adadeac0e57558f6e8','5':'45573232c8097c679503a6598f61e60b','6':'cfc137c7434d3a9a872332d405b5c553',
								 '7':'9d8210c22c1962db837e7b62a578975c','8':'665134fd44f21915cbeef955addf89ba','9':'758d0c0c71d8bafbe1ede86587191730',
								 '10':'397bb21acff1ca3052ac802f2aee06e0','11':'07707ff8a2a964656469a7be7bb3e576','12':'506d02539075e080ee12ebdf63908080',
								 '13':'03ed22f01ab43145733c0b6a647e0560','14':'8b93447086549e476c65699ed813a567','15':'cd0dfe9fa78cae2fc7becf8f8ec6c693',
								 '16':'e17bbb66eb4d6b62b7b0e2fbf062b6a6','17':'8fc95bb3101d024d890aa3543eb454c5','18':'a4870628045bb033a90e8c89f818e24d',
								 '19':'6a9d0c8298f0ba2fa13180e02b969f16','20':'aa75d35969cf3956bb4ace7bdc57b34e','21':'5d55f5ad6271d6a0d8806876924990f7',
								 '22':'efdb4e1d23ab7964302b828062a33447','Y':'3b38c639ad164d60f1a055b46fcd2748','X':'d5edbea3cf5d1716765dd4a7b41b7656',
								 'MT':'dfd6db5743d399516d5c8dadee5bee78'},

					'GRCh38_havana':
								{'1':'c4ef4ee14a4f0f7b319e9ed01f2a9742','2':'189cb33da673afcf161b7724d78d314b','3':'19fb72b2f4bebca66969d9eaaab0b06b',
								 '4':'04069627d35462add3095cd14299c181','5':'e3fbab13e69a5ba48c9bf706215aef19','6':'7de52b0136f65c25b2f3652e8888449a',
								 '7':'0071cc7e4024b0d81941df353ce99fdc','8':'c19f322ee6a786e94e665b1ab3f459ea','9':'00920f091df48d3dba8b270432bb38a8',
								 '10':'dfcdfcc0a2524888e91f9a40aa2289e5','11':'01118568e33f876f31f4d0585085b844','12':'11b6ffc45b5d198344e8d0f2e51fa08a',
								 '13':'e1241284e80f9a3afe204131c9eb6192','14':'b8d68bd919dc87483679314d0dae9514','15':'ca2c6d19dadc4d905013b37dfaa75fa6',
								 '16':'ccf7fe495c7ff7c11c430a935050bb4d','17':'ec09c6c67cc0496117da52bc584cab8f','18':'e8edcd9275c659101745570b0ddc436d',
								 '19':'293befa3316b6bb8b5a2a2c75e5e8030','20':'4065b4403e1b9659f9f049c4c886ed7f','21':'3f661bcfd29bf7153b4bfb95e8d112d8',
								 '22':'ea31a6d2cdac7d8aeb2e3be0688b4159','Y':'a68939c10c3ebc16340785bf99366328','X':'57623734b88f441f5499955d2c83a6f9',
								 },

					'mm9':
								{'1':'c5afc4b3f7f2119696214511d7a04341','2':'a7b467475a1b032d2c893dac1c419a28','3':'f922bc529a17324f1cd858f9a8723d65',
								'4':'f3d6b74e3c04dbd229e2f1e363607506','5':'5fee4f1889c9fe20f7f8562c62bbeb0a','6':'481d47b87da45f3a20181c780fd796c2',
								'7':'454ef2bf49a5ba8cfea3d16dfcfc7f25','8':'2f4162d4c824db78a2a2a820cb4fec81','9':'0649e6aec61af1ab8ab4797ea8e54119',
								'10':'38296256bcfe886c8ae771418e4fd824','11':'b31cb0ce693e35eaa77031d44b12e474','12':'d2b3e4b015742b6aea30ceec5a972968',
								'13':'df77b6d0ed1b133224b128c189736372','14':'0ec3c0e6b3fa2cdb957541f19792e130','15':'44fcaf2ec9b82dae910f85ce41c3cfad',
								'16':'ad7a8dbdf46fa7077e0982a54eab70b7','17':'71aee1dee3cd2078e4619c485d88817e','18':'727ec4ed3128ecacd6cd2f7558083553',
								'19':'461a7119781ab7f4b654fdd9ef76e0ec','Y':'471ff3bbb4520c020cfaa7ca8371c543','X':'9ccadf96cd3aa0ed9d299894a3d7fde0',
								'MT':'a1d56043ed8308908965dd080a4d0c8d'},

					'mm10':
								{'1':'ef88c5ac276a32a2865c0408f92acd55','2':'ced7325ef9e2dfedea3fbe26428a6059','3':'9cd1794eeea27553077a018038303908',
								'4':'da616d7ed6c67f824487eb2ed09cd33b','5':'b327b82da6986bf947105d07c0ad6d2e','6':'fb9a8fa0b85561f8d4de633c22d5157a',
								'7':'12457fd80f6806779fc0d4cc8d36fbad','8':'5d98d86bd22bee1cb226406f49ee7caf','9':'b2f26613fcc622a4003e4c945ae55e25',
								'10':'e9f3589529e258ede66d2e77bb87d21d','11':'76bcd285c3c66471ad6fccfabe42294c','12':'ac34fc3616c9609d8e75a59069e9007a',
								'13':'f81b976e4e4617b25945d06f9aa30846','14':'95dc042eb2aa7d4cc0abe071d4d7966e','15':'fbf2477833aff73ae085537cd7ee0f85',
								'16':'77cbcd009ba50891571f785595717ec1','17':'cd9e4dfdd168ed3de05dac4d44c6e692', '18':'945e83694c7c8f69d6186e1a2abc9771',
								'19':'e57b25f8869de31a9dbce06510711db6','Y':'c2146ba4ab1ec262f5e38b2a1ebc5f5b','X':'9af543088be046fdc63976c2d41de94c',
								'MT':'a1d56043ed8308908965dd080a4d0c8d'},


					'mm10_havana':
								{'1':'8351fa274b2e49b159fdfd1a17b8e522','2':'afa7459e3fa1a40f96279d318e3f6fc4','3':'106064694826ed8b9cfd29c5cfe39f1d',
								'4':'807aafbf65a28995b2ff64c040f54ff8','5':'afc3c4640788db5b8866d7644be4c49d','6':'2eff5f397f61735e90dc4089705751bc',
								'7':'f62fe7d720c219bacbe91fb0624fdb33','8':'1e523de22029fc4feea1e3222638d3ec','9':'3b75a1f2360ef49cbdc3da16a1e40d23',
								'10':'2d3cbd3aa1b5f6c78a11bc799787d226','11':'9b5e62ac65c2e05e7c09d7b3ac9aa507','12':'07e2bfe001933efb9052458f203fd8d4',
								'13':'4df084dfa4d2ea80519da4fa825d7c6e','14':'c4c01a321c2e8b67ed619dc91588d1cc','15':'6ccda29f891ac9d17567c69a1b47a40a',
								'16':'7df7554d6b3d7ca8331fe163e4930a64','17':'60ecd93b145f36226c246356701422f2', '18':'b52767b6bbc3a293993652e78db7ddac',
								'19':'5a246e5963dc0a300765b5cdd87e07d9','Y':'da73d749cd65735db5fd10556ae2376a','X':'ab9c1ad518e71b15ff543cda415a463f',
								},

					'rn6':
								{'1':'003723513cbdb3708fcc5d737c05199c','2':'53e52c5facc7f05462be533845f37425','3':'8d157a9b71fe9770cf783ea5459b19d7',
								'4':'a66dc1999bcc960ff11fe0b24c0d7b14','5':'601cf83411234adbdd9f911b89509564','6':'03b1f4af58fffdf213466ea85b570b3d',
								'7':'4ed05ddf9502ef79e121c02e391660e6','8':'3e2458daaf1b3e8ab4d0e0a9e60c067b','9':'8f83caeccec7ea6e35e404737138ee67',
								'10':'9c1af453a5facc9bfa821457bcfc4d30','11':'ef0480a905c55d76a3c58e295a85bc75','12':'643b6fe4a3a6363ffe64a6c316fa3e1a',
								'13':'102bb3fb420a4104c216bcdf99870374','14':'e26b8b63fba0ea7ced4f0330e93a8cdc','15':'da747616a1362d374d4786102fab6f9f',
								'16':'54e4f932eb0eda4cbf31156f96ef7235','17':'46c2facf5415e4eff8b0804161db722d', '18':'f1cb84f002967854b83bf266ec59a7a3',
								'19':'b85ca155fd1780fe5c327a4589c212a6','20':'899d3511352d78b9b9dc63f063d91b31','Y':'6a7a3539c329dc540dfa6db006003bb1',
								'X':'7a06bafab97c59a819f03633f0a6b7a2'},

					'c_elegans':
								{'I':'5a3ea8cf3dfbc641716b7bc805edcaae','II':'bf82edaa92809dd2fea2b791c38c9728','III':'d2df34b6743f41d3964549fc76c5f1a2',
								'IV':'23396bb57145d3acde2888947b5b8c3a','V':'09df3c53b12e5fd7d9035cc98ca221a3','X':'988046456f1409dfdb5e26444d84d238',
								'MtDNA':'48983f530959780de0125f74a87d4fc1'},

					'dog':
								{'1':'bef8283c1a36f9aef0e407de2ff6af00','2':'9cc961192bb5e58b3847060c3e9c1cfc','3':'d33263fa2de6666b41e140cb7a8da66c',
								 '4':'cd4ed39ebac1c04800ccf30466ec69f5','5':'c0f48a4a764e58388b48835aca2ec0a4','6':'4b472a2f8d0a53ac75cce04e7dc9279a',
								 '7':'12a61573a0da2c9306fff705bb1c39c1','8':'e22cf22a27560aa8523dc959ddcf6e25','9':'c079a73d719145cdd5c7c93969a1c392',
								 '10':'45805a518147f7846bd0457ca038c8df','11':'f38cda8508463a7607dff14a581ee7b0','12':'adb5de197f58bb827fa01fe924eb3a1d',
								 '13':'055a845ba97baad3b13d4d3359f88290','14':'27f0ba8e47996a058807a3827cf8e4a8','15':'2e9565c687a593eb0acbdd0962bb9255',
								 '16':'89b2225bb78d88b0fd1d38d9514ab0cb','17':'f0378253e2f083e42b665ea202fde3b0','18':'04d124e273f3b54a685ad6526223cd03',
								 '19':'67bae093919e6bb5ab6b9806c739d539','20':'5588387165a2e19c4533012cfb4998f3','21':'371cdf18a545728f7964b9db2fc72d5e',
								 '22':'fbf76865f88a018d93506e036f6a68bc','23':'085145e01d9fd9f0f999fb9e8e8d4400','24':'69b75a9962fb766b447e7d1252cb31ac',
								 '25':'12d5c6677b3e17170c317c1f5532d2a8','26':'13937d18e56b2b93d12fa5fcba48a138','27':'1d03d8ca5f201f4d156f5e1b38f7a67c',
								 '28':'c33395dec7fdc13e9d8f10afaa946f8c','29':'174f2db104ecaa5efef770f44241e3b0','30':'047d420ef9aecb933a7d83b6af820b23',
								 '31':'5be61f0c9944a5f2d7d1a5b2e75fb000','32':'212dcb867e95a642277a243fed8d8e41','33':'08a217b02cdd778cfdb0005dff4828b1',
								 '34':'4245d6fc370d9049ef4c25314fbef239','35':'1344aba8755b8a4e304629180fc0591a','36':'e4fff6ed84777905dc999ca6d6bc2557',
								 '37':'60d51ea6ae9e3f2fa316e3d03aff96b2','38':'4090ff76d94e6b38920916ae3ff2441c','X':'bce1372df64037d79b0995311d8ff971'},
					'ebv':
								{'gi_82503188_ref_NC_007605':'7f8894e52b0cac1f968c0402838bea07'},

					'yeast':
							{'I':'c7cd4d6148ff54d1d003d7fe3c39bfdf', 'II':'3303648ecf106df7c6c1449cb55980e8', 'III':'12f784db9c1cd2bc632fbefdb5329504',
							 'IV':'e5bbd701bae1ab3965da90e1131f5b5c', 'V':'b4447ed088b10a404925c153ca3055d4', 'VI':'3285962c283977064340b62a7f843601', 
							 'VII':'72e1438a61fde52bd66cef12c53ae4db', 'VIII':'0e00648307fd4d64f50d0feb7e428e13', 'IX':'43466267caa5ae7aec23df022b583564',
							 'X':'dac706139aeb5cd3d75f24d0e405131e', 'XI':'5d4b415e1e313918bba6cd061d416d58', 'XII':'e26e59a4b4944ab96b819a92b60c9868',
							 'XIII':'a12878f65ed7464e2824cb0c291d359e', 'XIV':'807431901ba993e6ba87861401da1cb2', 'XV':'44fd7223d9d30bdcf94cf627e5a172ca',
							 'XVI':'670d7aad863043ea7a125f004e281483'}}

		chromosome_fasta_path = ref_dir + "/references/chromosomes/tsb/"
		print("Beginning installation. This may take up to 40 minutes to complete.")
		if not rsync:
			try:
				if bash:
					try:
						os.system("bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chromosome_fasta_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerMatrixGenerator/' + genome + '.tar.gz 2>> install.log' + "'")
					except:
						print("The UCSD ftp site is not responding...pulling from sanger ftp now.")
					try:
						os.system("bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chromosome_fasta_path + ' ftp://ngs.sanger.ac.uk/scratch/project/mutographs/SigProf/' + genome + '.tar.gz 2>> install.log' + "'")
					except:
						print("The Sanger ftp site is not responding. Please check your internet connection/try again later.")
				else:
					os.system('wget -r -l1 -c -nc --no-parent -nd -P ' + chromosome_fasta_path + ' ftp://ngs.sanger.ac.uk/scratch/project/mutographs/SigProf/' + genome + '.tar.gz 2>> install.log')										
				os.system("tar -xzf " + ref_dir + "/references/chromosomes/tsb/" + genome + ".tar.gz -C " + ref_dir + "/references/chromosomes/tsb/")
				os.remove(ref_dir + "/references/chromosomes/tsb/" + genome + ".tar.gz")
			except:
				print("The ensembl ftp site is not currently responding.")
				sys.exit()
		else:
			print("Direct download for RSYNC is not yet supported")
			sys.exit()
			# try:
			# 	if bash:
			# 		os.system("bash -c '" + "rsync -av -m --include='*/' rsync://ftp.ngs.sanger.ac.uk/scratch/project/mutographs/SigProf/" + genome + ".tar.gz " + chromosome_fasta_path + " 2>&1>> install.log" + "'")
			# 	else:
			# 		os.system("rsync -av -m rsync://ftp://ngs.sanger.ac.uk/scratch/project/mutographs/SigProf/" + genome + ".tar.gz " + chromosome_fasta_path + " 2>&1>> install.log")
			# 	os.system("tar -xzf " + ref_dir + "/references/chromosomes/tsb/" + genome + ".tar.gz -C " + ref_dir + "/references/chromosomes/tsb/")
			# 	os.remove(ref_dir + "/references/chromosomes/tsb/" + genome + ".tar.gz")
			# except:
			# 	print("The ensembl ftp site is not currently responding.")
			# 	sys.exit()

		chromosome_TSB_path = chromosome_fasta_path + genome + "/"
		corrupt = False
		for files in os.listdir(chromosome_TSB_path):
			if "proportions" in files:
				continue
			if ".DS_Store" in files:
				continue
			chrom = files.split(".")
			chrom = chrom[0]
			check = md5(chromosome_TSB_path + files)
			if check_sum[genome][chrom] != check:
				corrupt = True
				os.remove(chromosome_TSB_path + files)
				print("[DEBUG] Chromosome " + chrom + " md5sum did not match => reference md5sum: " + str(check_sum[genome][chrom]) + "    new file md5sum: " + str(check))
		if corrupt:
			print("The transcriptional reference data appears to be corrupted. Please reinstall the " + genome + " genome.")
			sys.exit()
		print("The transcriptional reference data for " + genome + " has been saved.")


	else:
		print("Beginning installation. This may take up to 20 minutes to complete.")
		first_path = os.getcwd()
		ref_dir = os.path.dirname(os.path.abspath(__file__))
		os.chdir(ref_dir)

		print("[DEBUG] Path to SigProfilerMatrixGenerator used for the install: ", ref_dir)

		genomes = [genome]

		if os.path.exists("install.log"):
			os.remove("install.log")

		install_chromosomes(genomes, ref_dir, custom, rsync, bash)
		install_chromosomes_tsb (genomes, ref_dir, custom)

	if os.path.exists("BRCA_example/"):
		shutil.copy("BRCA_example/", "references/vcf_files/")
	if os.path.exists("example_test"):
		shutil.copy("example_test/", "references/vcf_files/")
	if os.path.exists("context_distributions/"):
		shutil.copy("context_distributions/", "references/chromosomes/")

	print("All reference files have been created.")
	if 'havana' in genome:
		genome = genome.split("_")[0]
	if genome != "rn6" and genome != 'dog' and genome != 'c_elegans' and not custom:
		print("Verifying and benchmarking installation now...")
		benchmark(genome, ref_dir)

	print ("To proceed with matrix_generation, please provide the path to your vcf files and an appropriate output path.")
	shutil.rmtree(chrom_string_dir)
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


	install_chromosomes(genomes, ref_dir, custom)
	install_chromosomes_tsb (genomes, ref_dir, custom)
	#install_chromosomes_tsb_BED (genomes, custom, ref_dir)
	if os.path.exists("BRCA_example/"):
		os.system("mv BRCA_example/ references/vcf_files/")
	if os.path.exists("example_test"):
		os.system("mv example_test/ references/vcf_files/")
	if os.path.exists("context_distributions/"):
		os.system("mv context_distributions/ references/chromosomes/")

	
	if os.path.exists(chrom_tsb_dir + "GRCh37/"):
		print("All reference files have been created.\nVerifying and benchmarking installation now...")
		benchmark(ref_dir)
	else:
		print("All reference files have been created.")
	print ("Please place your vcf files for each sample into the 'references/vcf_files/[test]/[mutation_type]/' directory. Once you have done that, you can proceed with the matrix generation.")
	#os.system("rm -r " + chrom_string_dir)
	print("Installation complete.")
	os.chdir(first_path)

if __name__ == '__main__':
	main()