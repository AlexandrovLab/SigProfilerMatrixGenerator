#!/usr/bin/env python3
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

from __future__ import print_function
import os
import sys
import time
import re
import argparse

start_time = time.time()

def save_tsb (chromosome_string_path, transcript_path, output_path):
	'''
	Creates binary files that contain the transcriptional information at a given 
	base. The transcriptional information entails Transcribed, Untranscribed, 
	Bi-directionally transcribed, and Non-transcribed. These files are required to
	simulate transcriptional strand bias.

	Input:
		chromosome_string_path -> path to chromosomes saved in individual text files
			   transcript_path -> path to transcript data. The following should be used as the format
								  (additional information may be saved, however it will not
								  be used when creating the files):

		Gene stable ID  Transcript stable ID    Chromosome/scaffold name    Strand  Transcript start (bp)   Transcript end (bp)

				   output_path -> path where the binary transcript files are saved.

	Output:
		Binary files that contain the transciptional information at each base for
		each chromosome. The files are saved as "[chrom]_192.txt"

	'''

	# Retrieves the transcript file 
	transcript_files = os.listdir(transcript_path)

	if len(transcript_files) < 3:
		# Instantiates required flag variables
		initial_flag = True
		initial_chrom = None
		for file in transcript_files:
			file_name = file.split("_")
			name = file_name[0]

			# Skips hidden files
			if name == '.DS':
				pass
			else:

				# Sorts the transcript file by chromosome and transcript range
				#os.system("sort -t $'\t' -k 3,3n -k 3,3 -k 5,5 -k 6,6 " + transcript_path+file + " -o " + transcript_path+file)
				with open(transcript_path + file) as f:
					lines = [line.strip().split() for line in f]
				output = open(transcript_path + file, 'w')
				for line in sorted(lines, key = lambda x: (['I','II','III','IV','V','chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI','X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23','24',
															'25','26','27','28','29','30','31','32','33','34','35','36','37','38','39', 'MT', 'M', 'MtDNA', 'chrM'].index(x[2]), int(x[4]), int(x[5]))):
					print('\t'.join(line), file=output)

				with open (transcript_path + file) as f:
					#next(f)

					# Saves all transcriptional information for a chromosome into a unique file
					for lines in f:
						line = lines.strip()
						line_split = line.split()
						chrom = line_split[2]

						if initial_flag:
							initial_chrom = chrom
							initial_flag = False
							out = open (transcript_path+chrom + "_transcripts.txt", 'w')

						if chrom == initial_chrom:
							print(line, file=out)

						else:
							out.close()
							initial_chrom = chrom
							out = open(transcript_path+chrom+"_transcripts.txt",'w')
					out.close()

			# Removes the initial transcript file
			# os.system("rm " + transcript_path + file)
			os.remove(transcript_path + file)

	# Retrieves the chromosome based list of transcript files 
	transcript_files = os.listdir(transcript_path)
	if '.DS_Store' in transcript_files:
		transcript_files.remove('.DS_Store')
	for files in transcript_files:
		with open(transcript_path + files) as f:
			lines = [line.strip().split() for line in f]

		output = open(transcript_path + files, 'w')

		for line in sorted(lines, key = lambda x: (int(x[4]))):
			print('\t'.join(line), file=output)

		output.close()

	print("Creating the transcriptional reference files now. This may take awhile...")

	first = True
	chrom_total = len(transcript_files)
	z = 1
	for files in transcript_files:
		file_name = files.split("_")
		chrom = file_name[0]
		if "chr" in chrom:
			chrom = chrom[3:]
		
		if chrom == '.DS':
			pass
		else:
			# try:
			if True:
				if os.path.exists(output_path + chrom + ".txt"):
					continue

				else:
					# Create a binary file for the output.
					outFile = open (output_path + chrom + ".txt", "wb")


					# Instantiates all of the required parameters.
					transcripts = None
					bi_start = 0
					bi_end= None
					I_chrom = None
					I_T_start = None
					I_T_end = None
					I_strand = None
					location = 0
					chrom_length = None

					# Save the length for the current chromosome.
					with open (chromosome_string_path + chrom + ".txt") as f:
						chrom_string = f.readline().strip()
						chrom_length = len(chrom_string)

					# Establishes a unique binary number for each transcript type.
					Non_transcribed_A = '00000000' # 0
					Non_transcribed_C = '00000001' # 1
					Non_transcribed_G = '00000010' # 2
					Non_transcribed_T = '00000011' # 3
					Transcribed_A = '00000100'     # 4
					Transcribed_C = '00000101'     # 5
					Transcribed_G = '00000110'     # 6
					Transcribed_T = '00000111'     # 7
					Untranscribed_A = '00001000'   # 8
					Untranscribed_C = '00001001'   # 9
					Untranscribed_G = '00001010'   # 10
					Untranscribed_T = '00001011'   # 11
					BiDirectional_A = '00001100'   # 12
					BiDirectional_C = '00001101'   # 13
					BiDirectional_G = '00001110'   # 14
					BiDirectional_T = '00001111'   # 15
					Non_transcribed_N = '00010000' # 16
					Transcribed_N = '00010001'     # 17
					Untranscribed_N = '00010010'   # 18 
					BiDirectional_N = '00010011'   # 19

					# Organizes the binary/integer formats for the data structures.
					nN_A = int(Non_transcribed_A, 2)
					nN_C = int(Non_transcribed_C, 2)
					nN_G = int(Non_transcribed_G, 2)
					nN_T = int(Non_transcribed_T, 2)
					nT_A = int(Transcribed_A, 2)
					nT_C = int(Transcribed_C, 2)
					nT_G = int(Transcribed_G, 2)
					nT_T = int(Transcribed_T, 2)
					nU_A = int(Untranscribed_A, 2)
					nU_C = int(Untranscribed_C, 2)
					nU_G = int(Untranscribed_G, 2)
					nU_T = int(Untranscribed_T, 2)
					nB_A = int(BiDirectional_A, 2)
					nB_C = int(BiDirectional_C, 2)
					nB_G = int(BiDirectional_G, 2)
					nB_T = int(BiDirectional_T, 2)
					
					nN_N = int(Non_transcribed_N, 2)
					nT_N = int(Transcribed_N, 2)
					nU_N = int(Untranscribed_N, 2)
					nB_N = int(BiDirectional_N, 2)

					byte_ref = {'N_A':bytes([nN_A]), 'N_C':bytes([nN_C]), 'N_G':bytes([nN_G]), 'N_T':bytes([nN_T]),
								'T_A':bytes([nT_A]), 'T_C':bytes([nT_C]), 'T_G':bytes([nT_G]), 'T_T':bytes([nT_T]),
								'U_A':bytes([nU_A]), 'U_C':bytes([nU_C]), 'U_G':bytes([nU_G]), 'U_T':bytes([nU_T]),
								'B_A':bytes([nB_A]), 'B_C':bytes([nB_C]), 'B_G':bytes([nB_G]), 'B_T':bytes([nB_T]),
								'N_N':bytes([nN_N]), 'T_N':bytes([nT_N]), 'U_N':bytes([nU_N]), 'B_N':bytes([nB_N])}

					stringTest = ''


					l = 0
					line_count = 1
					pointer = 0
					with open(transcript_path + files, 'r') as f:
						all_lines2 = f.readlines()
						all_lines = all_lines2[:]
			#            all_lines = all_lines2[1:]

						for lines in all_lines:
							
							# Handles the first line separately. 
							if pointer == 0:
								first_line = lines.split()
								location = int(first_line[4]) -1
								I_chrom = first_line[2][3:]
								I_T_start = int(first_line[4]) -1
								I_T_end = int(first_line[5]) - 1
								I_strand = first_line[3]
								
								# Saves Non-transcribed data up until
								# the first transcript.
								for i in range(0, location, 1):
									nuc = chrom_string[i].upper()
									outFile.write(byte_ref['N_' + nuc])
									l += 1
								pointer = 1

							# Handles all other lines after the first.
							else:
								first_line = lines.split()
								I_chrom = first_line[2][3:]
								I_T_start = int(first_line[4]) - 1
								I_T_end = int(first_line[5]) - 1
								I_strand = first_line[3]

								# Saves Non-transcribed data up until the 
								# next transcript.
								if I_T_start > location:
									for i in range(location, I_T_start, 1):
										nuc = chrom_string[i].upper()
										outFile.write(byte_ref['N_' + nuc])
										l += 1
									location = I_T_start
								
							# Reads the subsequent line to look for bi-directional
							# transcription sites.
							for line2 in all_lines[line_count:]:
								next_line = line2.split()
								c_chrom = next_line[2][3:]
								c_start = int(next_line[4]) - 1
								c_end= int(next_line[5]) - 1
								c_strand = next_line[3]
								
								# Breaks to the next line if the two transcripts 
								# don't overlap.
								if c_start > I_T_end:
									break

								# Checks if the transcripts are on opposite strands
								# if they overlap.    
								else:
									if c_strand != I_strand:
										bi_start = c_start
										if c_end > I_T_end:
											bi_end = I_T_end
										else:
											bi_end = c_end

									# Saves Un/Transcribed data up until the start of
									# the bi-directional transcription.
									if bi_start != 0 and bi_start > location:
										for i in range (location, bi_start, 1):
											if I_strand == '-1':
												nuc = chrom_string[i].upper()
												outFile.write(byte_ref['T_' + nuc])
												l += 1
											else:
												nuc = chrom_string[i].upper()
												outFile.write(byte_ref['U_' + nuc])
												l += 1
										location = bi_start

										# Saves Bi-directional data for the length
										# of the opposing transcripts overlap
										for i in range (location, bi_end, 1):
											nuc = chrom_string[i].upper()
											outFile.write(byte_ref['B_' + nuc])
											l += 1
										location = bi_end
									bi_start = 0
									bi_end = 0

							# Saves Un/Transcribed data up to the end of the current
							# transcript if data has not already been saved for these bases.
							if I_T_end > location:
								for i in range(location, I_T_end, 1):
									if I_strand == '-1':
										nuc = chrom_string[i].upper()
										outFile.write(byte_ref['T_' + nuc])
										l += 1
									else:
										nuc = chrom_string[i].upper()
										outFile.write(byte_ref['U_' + nuc])
										l +=1
								location = I_T_end
							line_count += 1

					# Save Non-transcribed data for the remainder of the chromosome after
					# the final transcript is analyzed         
					if location < chrom_length:
						for i in range(location, chrom_length,1):
							nuc = chrom_string[i].upper()
							outFile.write(byte_ref['N_' + nuc]) 
							l += 1        
						
					outFile.close()
					if first:
						print("The transcript reference file has been created for Chromosome: " +chrom + " (" + str(z) + "/" + str(chrom_total) + ")")
						first = False
					else:
						print("                                                               " +chrom + " (" + str(z) + "/" + str(chrom_total) + ")")
					z += 1
			#except:
			 #   print(chrom + " has been skipped. This chromosome is not supported.")

	end_time = time.time()
	print("Transcript files created.\n Job took: ", end_time-start_time, " seconds")

def main ():
	parser = argparse.ArgumentParser(description="Provide the necessary arguments to create the transcriptional strand bias reference strings.")
	parser.add_argument("--genome", "-g",help="Provide a reference genome. (ex: GRCh37, GRCh38, mm10)")
	args=parser.parse_args()
	genome = args.genome

	#script_dir = os.getcwd()
	#ref_dir = re.sub('\/scripts$', '', script_dir)
	ref_dir, tail = os.path.split(os.path.dirname(os.path.abspath(__file__)))

	chromosome_string_path = ref_dir + "/references/chromosomes/chrom_string/" + genome + "/"
	transcript_path = ref_dir + "/references/chromosomes/transcripts/" + genome + "/"
	output_path = ref_dir + "/references/chromosomes/tsb/" + genome + "/"
	if os.path.exists(output_path) == False:
		os.makedirs(output_path)

	save_tsb(chromosome_string_path, transcript_path, output_path)

if __name__ == '__main__':
	main()


