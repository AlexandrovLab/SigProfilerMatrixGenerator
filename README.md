# SigProfilerMatrixGenerator
SigProfilerMatrixGenerator creates mutational matrices for all types of somatic mutations. It allows downsizing the generated mutations only to parts for the genome (e.g., exome or a custom BED file). The tool seamlessly integrates with other SigProfiler tools. 

**INTRODUCTION**

The purpose of this document is to provide a guide for using the SigProfilerMatrixGenerator framework to generate mutational matrices for a set of samples with associated mutational catalogues. 

**PREREQUISITES**

The framework is written in PYTHON, however, it also requires the following software with the given versions (or newer):

  * PYTHON          version 3.4 or newer
  * PANDAS [MODULE]        any version
  * WGET                   version 1.9
  * SigProfilerPlotting    newest version (https://github.com/AlexandrovLab/SigProfilerPlotting)

By default the installation process will save the FASTA files for all chromosomes for the default genome 
assemblies (GRCh37, GRCH38, mm10, mm9). As a result, ~3 Gb of storage must be available for the downloads for each genome.

**QUICK START GUIDE**

This section will guide you through the minimum steps required to create mutational matrices:
1. Install the python package using pip:
```
                          pip install SigProfilerMatrixGenerator
```
2. Install your desired reference genome from the command line/terminal as follows (available reference genomes are: GRCh37, GRCh38, mm9, and mm10):
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37')
```
    This will install the human 37 assembly as a reference genome. You may install as many genomes as you wish.

3. Place your vcf files in your desired output folder. It is recommended that you name this folder based on your project's name
4. From within a python session, you can now generate the matrices as follows:
```
$ python3
>>from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
>>matrices = matGen.SigProfilerMatrixGeneratorFunc("test", "GRCh37", "/Users/ebergstr/Desktop/test",plot=True)
```
  The layout of the parameters are as follows:
  
      SigProfilerMatrixGeneratorFunc(project, reference_genome, path_to_input_files, plot)
      
  where project, reference_genome, and path_to_input_files must be strings (surrounded by quotation marks, ex: "test") and plot is a boolean argument (True or False)


**INPUT FILE FORMAT**

This tool currently supports maf, vcf, simple text file, and ICGC formats. The user must provide variant data adhering to one of these four formats. If the user’s files are in vcf format, each sample must be saved as a separate files. 


**Output File Structure**

The output structure is divided into three folders: input, output, and logs. The input folder contains copies of the user-provided input files. The outputfolder contains
a DINUC, SBS, INDEL, and TSB folder (there will also be a plots folder if this parameter is chosen). The matrices will be saved into the appropriate folders. The logs 
folder contains the error and log files for the submitted job.


**SUPPORTED GENOMES**

This tool currently supports the following genomes:

GRCh38.p12 [GRCh38] (Genome Reference Consortium Human Reference 37), INSDC
Assembly GCA_000001405.27, Dec 2013. Released July 2014. Last updated January 2018. This genome was downloaded from ENSEMBL database version 93.38.

GRCh37.p13 [GRCh37] (Genome Reference Consortium Human Reference 37), INSDC
Assembly GCA_000001405.14, Feb 2009. Released April 2011. Last updated September 2013. This genome was downloaded from ENSEMBL database version 93.37. 

GRCm38.p6 [mm10] (Genome Reference Consortium Mouse Reference 38), INDSDC
Assembly GCA_000001635.8, Jan 2012. Released July 2012. Last updated March 2018. This genome was downloaded from ENSEMBL database version 93.38. 

GRCm37 [mm9] (Release 67, NCBIM37), INDSDC Assembly GCA_000001635.18.
Released Jan 2011. Last updated March 2012. This genome was downloaded from ENSEMBL database version release 67.


**LOG FILES**

All errors and progress checkpoints are saved into *sigProfilerMatrixGenerator_[project]_[genome].err* and *sigProfilerMatrixGenerator_[project]_[genome].out*, respectively. 
For all errors, please email the error and progress log files to the primary contact under CONTACT INFORMATION.

**COPYRIGHT**

This software and its documentation are copyright 2018 as a part of the sigProfiler project. The sigProfilerMatrixGenerator framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Erik Bergstrom at ebergstr@eng.ucsd.edu