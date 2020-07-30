import pandas as pd
import numpy as np
import os
import shutil

def generateCNVMatrix(file_type, input_file, project, output_path):

    super_class = ['het', 'LOH', "homdel"]
    # het_sub_class = ['amp+', 'amp', 'gain', 'neut']
    # loh_subclass = ['amp+', 'amp', 'gain', 'neut', "del"]
    hom_del_class = ['0-100kb', '100kb-1Mb', '>1Mb']
    x_labels = ['>40Mb', '10Mb-40Mb', '1Mb-10Mb', '100kb-1Mb', '0-100kb']

    df = pd.read_csv(input_file, sep='\t')

    #make sample by feature matrix for nmf with rows as features and samples as columns
    features = []
    with open('SigProfilerMatrixGenerator/references/CNV/CNV_features.tsv') as f:
        for line in f:
            features.append(str(line.strip()))
    assert(len(features) == 48)
    columns = list(df[df.columns[0]].unique())
    arr = np.zeros((48, len(columns)), dtype='int')
    nmf_matrix = pd.DataFrame(arr, index=features, columns=columns)


    # 2 - total copy number {del=0-1; neut=2; gain=3-4; amp=5-8; amp+=9+}
    CN_classes = ["1","2","3-4","5-8","9+"] # different total CN states
    CN_class = []
    if file_type == 'ASCAT_NGS':
        for tcn in df['Tumour TCN']:
            if tcn == 2:
                CN_class.append("2")
            elif tcn == 0 or tcn == 1:
                CN_class.append("1")
            elif tcn == 3 or tcn == 4:
                CN_class.append("3-4")
            elif tcn >= 5 and tcn <= 8:
                CN_class.append("5-8")
            else:
                CN_class.append("9+")

    elif file_type == 'SEQUENZA':
        for tcn in df['CNt']:
            if tcn == 2:
                CN_class.append("2")
            elif tcn == 0 or tcn == 1:
                CN_class.append("1")
            elif tcn == 3 or tcn == 4:
                CN_class.append("3-4")
            elif tcn >= 5 and tcn <= 8:
                CN_class.append("5-8")
            else:
                CN_class.append("9+")
    elif file_type == "ASCAT":
        for acn, bcn in zip(df['nMajor'], df['nMinor']):
            tcn = acn + bcn
            if tcn == 2:
                CN_class.append("2")
            elif tcn == 0 or tcn == 1:
                CN_class.append("1")
            elif tcn == 3 or tcn == 4:
                CN_class.append("3-4")
            elif tcn >= 5 and tcn <= 8:
                CN_class.append("5-8")
            else:
                CN_class.append("9+")
    elif file_type == 'ABSOLUTE':
        for acn, bcn in zip(df['Modal_HSCN_1'], df['Modal_HSCN_2']):
            tcn = acn + bcn
            if tcn == 2:
                CN_class.append("2")
            elif tcn == 0 or tcn == 1:
                CN_class.append("1")
            elif tcn == 3 or tcn == 4:
                CN_class.append("3-4")
            elif tcn >= 5 and tcn <= 8:
                CN_class.append("5-8")
            else:
                CN_class.append("9+")

    else:
        pass

    df['CN_class'] = CN_class

    # 1 - LOH status {hom del; heterozygous; LOH}.
    LOH_status = []

    if file_type == 'ASCAT':
       for acn, bcn in zip(df['nMajor'], df['nMinor']):
            t = acn + bcn
            if t == 0:
                LOH_status.append("homdel")
            elif acn == 0 or bcn == 0:
                LOH_status.append("LOH")
            else:
                LOH_status.append("het")

    elif file_type == 'SEQUENZA':
        for t, a, b in zip(df['CNt'], df['A'], df['B']):
            if t == 0:
                LOH_status.append("homdel")
            elif a == 0 or b == 0:
                LOH_status.append("LOH")
            else:
                LOH_status.append("het")
    elif file_type == 'ASCAT_NGS':
        normal_ACN = np.asarray(df['Normal TCN']) - np.asarray(df['Normal BCN'])
        tumour_ACN = np.asarray(df['Tumour TCN']) - np.asarray(df['Tumour BCN'])
        df["Normal ACN"] = list(normal_ACN)
        df["Tumour ACN"] = list(tumour_ACN)

        A_CN = np.asarray(df['Tumour ACN']) #copy number of A allele in tumor
        B_CN = np.asarray(df['Tumour BCN']) #copy number of B allele in tumor
        loh = np.minimum(A_CN, B_CN) #minimum copy number when considering both A and B alleles
        df['loh'] = list(loh)
        for t, a in zip(df['Tumour TCN'], df['loh']):
            if t == 0:
                LOH_status.append("homdel")
            elif a == 0:
                LOH_status.append("LOH")
            else:
                LOH_status.append("het")

    elif file_type == 'ABSOLUTE':
        for acn, bcn in zip(df['Modal_HSCN_1'], df['Modal_HSCN_2']):
            t = acn + bcn
            if t == 0:
                LOH_status.append("homdel")
            elif acn == 0 or bcn == 0:
                LOH_status.append("LOH")
            else:
                LOH_status.append("het")
    else:
        print("Please provide a proper file type")


    df['LOH'] = LOH_status

    lengths = []

    #get chromosomal sizes
    if file_type == 'ASCAT_NGS':
        for start, end in zip(df['Start Position'], df['End Position']):
            lengths.append((end - start)/1000000) #megabases
    elif file_type == 'SEQUENZA':
        for start, end in zip(df['start.pos'], df['end.pos']):
            lengths.append((end - start)/1000000)
    elif file_type == 'ABSOLUTE': #Start End
        for start, end in zip(df['Start'], df['End']):
            lengths.append((end - start)/1000000)
    elif file_type == 'ASCAT':
        for start, end in zip(df['startpos'], df['endpos']):
            lengths.append((end - start)/1000000)
    else:
        pass


    df['length'] = lengths

    sizes = []
    size_bins = [] #features of matrix(matches Chris's classification)
    hom_del_class = ['0-100kb', '100kb-1Mb', '>1Mb']

    for l, s in zip(lengths, df['LOH']): #keep in mind the lengths are in megabases
        if s == 'homdel':
            if l > -0.01 and l <= 0.1:
                size = "0-100kb"
                size_bin = "(-0.01,0.1]"
            elif l > 0.1 and l <= 1:
                size = "100kb-1Mb"
                size_bin = "(0.1,1]"
            else:
                size = '>1Mb'
                size_bin = "(1,Inf]"
        else:
            if l > -0.01 and l <= 0.1:
                size = "0-100kb"
                size_bin = "(-0.01,0.1]"
            elif l > 0.1 and l <= 1:
                size = "100kb-1Mb"
                size_bin = "(0.1,1]"
            elif l > 1 and l <= 10:
                size = "1Mb-10Mb"
                size_bin = "(1,10]"
            elif l > 10 and l <= 40:
                size = "10Mb-40Mb"
                size_bin = "(10,40]"
            else:
                size = ">40Mb"
                size_bin = "(40,Inf]"
        sizes.append(size)
        size_bins.append(size_bin)
    df['size_classification'] = sizes


    for sample, tcn, loh, size in zip(df[df.columns[0]], df['CN_class'], df['LOH'], df['size_classification']):
        if loh == "homdel":
            channel = "0" + ":" + loh + ":" + size
        else:
            channel = tcn + ":" + loh + ":" + size
        nmf_matrix.at[channel, sample] += 1
 
    nmf_matrix.index.name = 'classification'
    output_path = output_path + project + "/"
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)
    nmf_matrix.to_csv(output_path + file_type + '.CNV.matrix.tsv', sep='\t')
    #return nmf_matrix
