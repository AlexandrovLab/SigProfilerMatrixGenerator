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

    features = ['0:homdel:0-100kb', '0:homdel:100kb-1Mb', '0:homdel:>1Mb', '1:LOH:0-100kb', 
            '1:LOH:100kb-1Mb', '1:LOH:1Mb-10Mb', '1:LOH:10Mb-40Mb', '1:LOH:>40Mb', 
            '2:LOH:0-100kb', '2:LOH:100kb-1Mb', '2:LOH:1Mb-10Mb', '2:LOH:10Mb-40Mb', '2:LOH:>40Mb', 
            '3-4:LOH:0-100kb', '3-4:LOH:100kb-1Mb', '3-4:LOH:1Mb-10Mb', '3-4:LOH:10Mb-40Mb', '3-4:LOH:>40Mb', 
            '5-8:LOH:0-100kb', '5-8:LOH:100kb-1Mb', '5-8:LOH:1Mb-10Mb', '5-8:LOH:10Mb-40Mb', '5-8:LOH:>40Mb', 
            '9+:LOH:0-100kb', '9+:LOH:100kb-1Mb', '9+:LOH:1Mb-10Mb', '9+:LOH:10Mb-40Mb', '9+:LOH:>40Mb', 
            '2:het:0-100kb', '2:het:100kb-1Mb', '2:het:1Mb-10Mb', '2:het:10Mb-40Mb', '2:het:>40Mb', 
            '3-4:het:0-100kb', '3-4:het:100kb-1Mb', '3-4:het:1Mb-10Mb', '3-4:het:10Mb-40Mb', '3-4:het:>40Mb', 
            '5-8:het:0-100kb', '5-8:het:100kb-1Mb', '5-8:het:1Mb-10Mb', '5-8:het:10Mb-40Mb', '5-8:het:>40Mb', 
            '9+:het:0-100kb', '9+:het:100kb-1Mb', '9+:het:1Mb-10Mb', '9+:het:10Mb-40Mb', '9+:het:>40Mb']
    
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
    elif file_type == 'PCAWG':
        for tcn in df["copy_number"]:
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
    elif file_type == 'FACETS':
        for tcn in df["tcn.em"]:
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
    elif file_type == 'PCAWG':
        for tcn, m in zip(df['copy_number'], df['mutation_type']):
            tcn = int(tcn)
            if m == 'copy neutral LOH' or m == 'amp LOH' or m == "hemizygous del LOH" or (m == "loss" and tcn == 1):
                LOH_status.append("LOH")
            elif m == "copy neutral" or m == "gain"                                                      :
                LOH_status.append("het")
            elif m == "loss" and tcn == 0:
                LOH_status.append("homdel")
            else:
                print(tcn, m)
                #raise ValueError('Unable to determine zygosity')
    elif file_type == 'FACETS':
        for t, a in zip(df['tcn.em'], df['lcn.em']):
            if (t == 0 and a == 0) or (t == 0 and a == "NA") or (t == 0 and pd.isnull(a)):
                LOH_status.append("homdel")
            elif t >= 1 and (a == 0 or a == "NA" or pd.isnull(a)):
                LOH_status.append("LOH")
            else:
                LOH_status.append("het")
                if t == 1:
                    print(t, a) 

    else:
        print("Please provide a proper file type")


# chrom   seg     num.mark        nhet    cnlr.median     mafR    segclust        cnlr.median.clust       mafR.clust      start   end     cf.em   tcn.em  lcn.em  cnlr.median-dipLogR
# 1       1       148     2       -0.524597703662466      1.13764233698258        21      -0.561362815875377      NA      10700   588600  0.72333842611374        1       0       -0.643682513766257
# 1       2       72      0       -2.07028566583277       0       5       -2.07028566583277       NA      589500  687700  0.72333842611374        0       0       -2.18937047593656
# 1       3       6485    10      -0.534644008651447      0.412452808996456       21      -0.561362815875377      NA      690300  7305200 0.72333842611374        1       0       -0.653728818755238


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
    elif file_type == 'PCAWG':
        for start, end in zip(df['chromosome_start'], df['chromosome_end']):
            lengths.append((end - start)/1000000)
    elif file_type == 'FACETS':
        for start, end in zip(df['start'], df['end']):
            lengths.append((end - start)/1000000)
    else:
        pass


    df['length'] = lengths

    sizes = []
    size_bins = []
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
    nmf_matrix.reindex([features])
    nmf_matrix.to_csv(output_path + file_type + '.CNV.matrix.tsv', sep='\t')