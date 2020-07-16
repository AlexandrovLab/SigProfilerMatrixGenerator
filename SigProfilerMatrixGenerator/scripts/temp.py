#make sample by feature matrix for nmf with rows as features and samples as columns
features = []
with open('./SigProfilerMatrixGenerator/scripts/CNV_features.tsv') as f:
    for line in f:
        features.append(line.strip())
assert(len(features) == 48)
print(features)
