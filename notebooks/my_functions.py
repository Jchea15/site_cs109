def import_gene_info():
    
    import pandas as pd

    # Read in gene expression data
    GENE_EXPRESSION = pd.read_csv('ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv', low_memory=False, index_col=0, header=2)
    GENE_EXPRESSION = GENE_EXPRESSION[6:]
    GENE_EXPRESSION = GENE_EXPRESSION.rename(columns={GENE_EXPRESSION.columns[0]: 'ProbeSet', GENE_EXPRESSION.columns[1]: 'LocusLink'})
    GENE_EXPRESSION = GENE_EXPRESSION.transpose()

    # Save gene expression data and info on genes separately
    GENE_EXPRESSION_DATA = GENE_EXPRESSION[2:]
    GENE_EXPRESSION_GENE_INFO = GENE_EXPRESSION[:2]
    GENE_EXPRESSION_GENE_INFO.columns = GENE_EXPRESSION_GENE_INFO.columns.values
    GENE_INFO_T = GENE_EXPRESSION_GENE_INFO.transpose()
    
    return GENE_INFO_T


def get_gene_ID(Name, GENE_INFO_T):
    return GENE_INFO_T.loc[GENE_INFO_T['LocusLink'] == Name].index[0]

def get_gene_name(ID, GENE_INFO_T):
    return GENE_INFO_T.loc[GENE_INFO_T.index == ID].LocusLink.values[0]

def for_loop_status(length, index = 0):
    import sys
    sys.stdout.write('\r%f%%' % ((index/length)*100))
    sys.stdout.flush()
    index += 1
    return index

def for_loop_status2(length, length2, index = 0, index2 = 0):
    import sys
    sys.stdout.write('\r%f%% of %f%%' % ((index/length)*100, (index2/length2)*100))
    sys.stdout.flush()
    index += 1
    index2 += 1
    return index, index2

def scale_predictor(df, predictor):
    mean = df[predictor].mean() #get mean
    std = df[predictor].std() #get training standard deviation
    df[predictor] = (df[predictor] - mean)/std #scale df set
    return df

def make_dict(train, test, dep):
    data = {'xtrain' : train.drop(dep, axis = 1).values,
               'ytrain' : train[dep].values,
               'xtest' : test.drop(dep, axis = 1).values,
               'ytest' : test[dep].values}
    return data


