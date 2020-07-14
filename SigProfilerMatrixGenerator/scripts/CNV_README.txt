**COPY NUMBER MATRIX GENERATION**

In order to generate a copy number matrix, provide the an abosolue path to a multi-sample segmentation file obtained from one of the following copy number calling tools (if you have individual sample files, please combine them into one file with the first column corresponding to the sample name):

1. ASCAT
2. ASCAT_NGS
3. SEQUENZA
4. ABSOLUTE
5. FACETS

In addition, provide the output directory for the resulting matrix containing the 48 channels found in CNV_features.txt as rows and samples as columns.

An example comman to generate the CNV matrix is as follows:

$ python3
>>from SigProfilerMatrixGenerator.scripts import CNVMatrixGenerator as scna
>>scna.generateCNVMatrix(file_type, input_file, output_path)

