import pandas as pd
import numpy as np

def generateCNVMatrix(file_type, input_matrix):
    
    super_class = ['het', 'LOH', "Hom del"]
    het_sub_class = ['amp+', 'amp', 'dup', 'neut']
    loh_subclass = ['amp+', 'amp', 'dup', 'neut', "del"]
    hom_del_class = ['0-100kb', '100kb-1Mb', '>1Mb']
    x_labels = ['>40Mb', '10Mb-40Mb', '1Mb-10Mb', '100kb-1Mb', '0-100kb']

    if file_type == 'ASCAT':
        df = pd.read_csv(input_matrix, sep='\t')
        
        #make sample by feature matrix for nmf with rows as features and samples as columns
        features = []
        with open('CNV_features1.tsv') as f:
            for line in f:
                features.append(line.strip())           
        # features = pd.read_csv('CNV_features1.tsv', sep='\t', header=None) #this file determines the ordering the x-axis for plots
        # features = features.iloc[:, 0]       
        columns = df['Sample'].unique()
        nmf_matrix = pd.DataFrame(index=features, columns=columns)
        
        # 2 - total copy number {del=0-1; neut=2; gain=3-4; amp=5-8; amp+=9+}.
        CN_classes = ["del","neut","dup","amp","amp+"] # different total CN states
        CN_class = []
        for tcn in df['Tumour TCN']:
            if tcn == 2:
                CN_class.append("neut")
            elif tcn == 0 or tcn == 1:
                CN_class.append("del")
            elif tcn == 3 or tcn == 4:
                CN_class.append("gain")
            elif tcn >= 5 and tcn <= 8:
                CN_class.append("amp")
            else:
                CN_class.append("amp+")
        df['CN_class'] = CN_class
        
        # 1 - LOH status {hom del; heterozygous; LOH}.
        LOH_status = []
        normal_ACN = np.asarray(df['Normal TCN']) - np.asarray(df['Normal BCN'])
        tumour_ACN = np.asarray(df['Tumour TCN']) - np.asarray(df['Tumour BCN'])
        df["Normal ACN"] = list(normal_ACN)
        df["Tumour ACN"] = list(tumour_ACN)

        a = np.asarray(df['Tumour ACN'])
        b = np.asarray(df['Tumour BCN'])
        loh = np.minimum(a, b)
        df['loh'] = list(loh)
                
        for t, a in zip(df['Tumour TCN'], df['loh']):
            if t == 0:
                LOH_status.append("homdel") 
            elif a == 0:
                LOH_status.append("LOH")
            else:
                LOH_status.append("het")
                              
        df['LOH'] = LOH_status
        
        lengths = []
        
        for start, end in zip(df['Start Position'], df['End Position']):
            lengths.append((end - start)/1000000) #megabases
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
        df['size_bin'] = size_bins
        
        counts = {} #dictionary that maps (sample, feature) to frequency     
        for a, c1 in enumerate(super_class):
            df1 = df[df['LOH'] == c1]
            if c1 == 'het':
                for b, c2 in enumerate(het_sub_class): #amp+, amp, etc.
                    df2 = df1[df1['CN_class'] == c2]                   
                    for c, x in enumerate(x_labels):
                        df3 = df2[df2['size_classification'] == x]
                        for s in df3['Sample'].unique(): 
                            sample_df = df3[df3['Sample'] == s]
                            for a, b, c in zip(sample_df['CN_class'], sample_df['LOH'], sample_df['size_classification']):
                                f = a+":"+b+":"+c
                                key = (s, f)
                                value = sample_df.shape[0]
                                if f not in set(features):
                                    print (f)
                                else:
                                    counts[key] = value
                                break
                                
            elif c1 == 'LOH':
                for b, c2 in enumerate(loh_subclass): #amp+, amp, etc.
                    df2 = df1[df1['CN_class'] == c2]
                    for c, x in enumerate(x_labels):
                        df3 = df2[df2['size_classification'] == x]
                        for s in df3['Sample'].unique(): 
                            sample_df = df3[df3['Sample'] == s]
                            for a, b, c in zip(sample_df['CN_class'], sample_df['LOH'], sample_df['size_classification']):
                                f = a+":"+b+":"+c
                                key = (s, f)
                                value = sample_df.shape[0]
                                if f not in set(features):
                                    print (f)
                                else:
                                    counts[key] = value
                                                  
            else: #Hom del
                for b, c2 in enumerate(hom_del_class):
                    df3 = df1[df1['size_classification'] == c2]
                    for s in df3['Sample'].unique(): 
                        sample_df = df3[df3['Sample'] == s]
                        for a, b, c in zip(sample_df['CN_class'], sample_df['LOH'], sample_df['size_classification']):
                            f = "del:homdel:" + c
                            key = (s, f)
                            value = sample_df.shape[0]
                        if f not in set(features):
                            print (f)
                        else:
#                             print (key)
                            counts[key] = value
                        
      
        #use counts dictionary(which maps (sample, CNV feature) to frequency observed) to populate matrix
        for i, row in enumerate(nmf_matrix.index):
            for j, sample in enumerate(nmf_matrix.columns):
                if (sample, row) in counts:
                    nmf_matrix.iat[i, j] = counts[(sample, row)]
                else:
                    nmf_matrix.iat[i, j] = 0
                    
                    
        nmf_matrix.index.name = 'classification'     
        nmf_matrix.to_csv('CNV.matrix.tsv', sep='\t')
        

if __name__ == "__main__":
    generateCNVMatrix("ASCAT", "ASCAT-seg.tsv")