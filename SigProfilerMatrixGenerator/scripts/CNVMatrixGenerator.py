import pandas as pd
import numpy as np
import os
from pathlib import Path
import shutil
from typing import Union


def _bucketize_total_copy_number(copy_number: Union[float,int] ) -> str:
    """Transform copy number into bucket label."""
    copy_number = round(copy_number)

    if copy_number == 2:
        return "2"
    elif copy_number in (0, 1):
        return "1"
    elif copy_number in (3, 4):
        return "3-4"
    elif 5 <= copy_number <= 8:
        return "5-8"
    return "9+"


def _bucketize_segment_length(length: float, homozygous_deletion: bool) -> str:
    """Transform segment length (in megabases) into bucket label.

    Args:
        length: Copy number segment length in megabases (mb).
        homozygous_deletion: Is this copy number variant a homozygous deletion?
    """
    if -0.01 < length <= 0.1:
        return "0-100kb"
    elif 0.1 < length <= 1:
        return "100kb-1Mb"
    # Bucketization for homozygous deletion is coarser.
    if homozygous_deletion:
        return '>1Mb'

    if 1 < length <= 10:
        return "1Mb-10Mb"
    elif 10 < length <= 40:
        return "10Mb-40Mb"
    return ">40Mb"


def _classify_heterozygosity_state(minor_allele: Union[int, float], major_allele: Union[int, float]) -> str:
    """Transform minor and major allele copy number into status.

    As described on p. 9 in Steele et al., Nature ('22).

    Returns:
        homdel: Homozygous deletion.
        het: Heterozygous segment.
        LOH: Loss of heterozygocity.
    """
    # Round floating point copy numbers to nearest integer.
    minor_allele, major_allele = round(minor_allele), round(major_allele)

    total_copy_number = minor_allele + major_allele
    if total_copy_number == 0:
        return "homdel"
    elif minor_allele == 0 or major_allele == 0:
        return "LOH"
    else:
        return "het"


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

    sample_names = df[df.columns[0]]
    # For filetype PURPLE, only one sample per file.
    if file_type == 'PURPLE':
        # Use filename as sample name, stripping trailing suffix.
        name = Path(input_file).stem
        sample_names = [name] * df.shape[0]  # One for each record.
        columns = [name]
    # Other formats accomodate multiple samples per file.
    else:
        columns = list(sample_names.unique())
    arr = np.zeros((48, len(columns)), dtype='int')
    nmf_matrix = pd.DataFrame(arr, index=features, columns=columns)


    # 2 - total copy number {del=0-1; neut=2; gain=3-4; amp=5-8; amp+=9+}
    CN_classes = ["1","2","3-4","5-8","9+"] # different total CN states
    CN_class = []
    if file_type == 'ASCAT_NGS':
        for tcn in df['Tumour TCN']:
            CN_class.append(_bucketize_total_copy_number(tcn))

    elif file_type == 'SEQUENZA':
        for tcn in df['CNt']:
            CN_class.append(_bucketize_total_copy_number(tcn))

    elif file_type == "ASCAT":
        for acn, bcn in zip(df['nMajor'], df['nMinor']):
            tcn = acn + bcn
            CN_class.append(_bucketize_total_copy_number(tcn))

    elif file_type == 'ABSOLUTE':
        for acn, bcn in zip(df['Modal_HSCN_1'], df['Modal_HSCN_2']):
            tcn = acn + bcn
            CN_class.append(_bucketize_total_copy_number(tcn))

    elif file_type == 'PCAWG':
        for tcn in df["copy_number"]:
            CN_class.append(_bucketize_total_copy_number(tcn))

    elif file_type == 'FACETS':
        for tcn in df["tcn.em"]:
            CN_class.append(_bucketize_total_copy_number(tcn))

    elif file_type == 'PURPLE':
        CN_class = df['copyNumber'].apply(_bucketize_total_copy_number)

    else:
        pass


    df['CN_class'] = CN_class

    # 1 - LOH status {hom del; heterozygous; LOH}.
    LOH_status = []

    if file_type == 'ASCAT':
       for acn, bcn in zip(df['nMajor'], df['nMinor']):
            LOH_status.append(_classify_heterozygosity_state(bcn, acn))

    elif file_type == 'SEQUENZA':
        for t, a, b in zip(df['CNt'], df['A'], df['B']):
            assert t == a + b 
            LOH_status.append(_classify_heterozygosity_state(b, a))
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
            LOH_status.append(_classify_heterozygosity_state(t - a, a))

    elif file_type == 'ABSOLUTE':
        for acn, bcn in zip(df['Modal_HSCN_1'], df['Modal_HSCN_2']):
            LOH_status.append(_classify_heterozygosity_state(bcn, acn))
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

    elif file_type == 'PURPLE':
        for minor, major in zip(df['minorAlleleCopyNumber'], df['majorAlleleCopyNumber']):
            LOH_status.append(_classify_heterozygosity_state(minor, major))

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
    elif file_type in ('FACETS', 'PURPLE'):
        for start, end in zip(df['start'], df['end']):
            lengths.append((end - start)/1000000)
    else:
        pass


    df['length'] = lengths

    sizes = []
    hom_del_class = ['0-100kb', '100kb-1Mb', '>1Mb']

    for l, s in zip(lengths, df['LOH']): #keep in mind the lengths are in megabases
        is_homozygous_deletion = s == 'homdel'
        size = _bucketize_segment_length(length=l, homozygous_deletion=is_homozygous_deletion)
        sizes.append(size)
    df['size_classification'] = sizes

    for sample, tcn, loh, size in zip(sample_names, df['CN_class'], df['LOH'], df['size_classification']):
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