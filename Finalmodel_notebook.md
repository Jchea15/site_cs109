Cleaning: Jupyter notebook
===================

[Download this notebook here.](https://raw.githubusercontent.com/Pagel56/site_cs109/master/notebooks/Final%20Model.ipynb)	

```python
cd ~/Desktop/proj_cs109/Data
```

    /Users/Nick/Desktop/proj_cs109/Data
    


```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn.metrics import accuracy_score
from numpy.random import randint
import sklearn.metrics as metrics
from sklearn.preprocessing import PolynomialFeatures
from sklearn.decomposition import PCA
from collections import OrderedDict
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV
from my_functions import for_loop_status
from my_functions import for_loop_status2
import matplotlib
from matplotlib import cm
from sklearn.linear_model import LogisticRegressionCV
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import export_graphviz
from sklearn.model_selection import cross_val_score
from sklearn.tree import export_graphviz
from sklearn.model_selection import KFold
import matplotlib.patches as mpatches
from IPython.display import Image
from IPython.display import display
import sys as sys
from sklearn.ensemble import GradientBoostingClassifier
%matplotlib inline
```


```python
# import data
gene_corr_base = pd.read_csv('gene_corr_base.csv', index_col = None, header=None)
gene_corr_cut = pd.read_csv('gene_corr_cut.csv', index_col = None, header=None)
gene_corr_p = pd.read_csv('gene_corr_p.csv', index_col = None, header=None)
gene_names = pd.read_csv('gene_names.csv', index_col = None, header=None)
df_full = pd.read_csv('post_imputation.csv', index_col = 0)

# just gene expression data
X_full = df_full.iloc[:,14:49400]

# make correlations dataframe
corrs_df = pd.DataFrame()
corrs_df['Gene_Name'] = gene_names[0]
corrs_df['Base'] = gene_corr_base[0]/15
corrs_df['Cut'] = gene_corr_cut[0]/15
corrs_df['p'] = gene_corr_p[0]/15

# choose most significant genes
Used = corrs_df[corrs_df.p < 0.01]
used_genes = list(Used.Gene_Name.values)

# interaction and polynomial terms
predictors_expo_list = []
for_index = 0
index = 0
for predictor in used_genes:
    for_index = for_loop_status(len(used_genes), for_index)
    predictors_temp_list = []
    for expo in [2,3,4]:
        #if expo = 2 and the predictor is atemp
        #then you'll get "atemp2" which is atemp^2
        predictors_temp_list.append('{0}{1}'.format(predictor, expo))
        df_full['{0}{1}'.format(predictor, expo)] = df_full['{0}'.format(predictor)]**int('{0}'.format(expo))
    predictors_expo_list.append(predictors_temp_list)
    for predictor2 in used_genes[int(index + 1):]:
        if predictor != predictor2:
            df_full['{0}_x_{1}'.format(predictor, predictor2)] = df_full['{0}'.format(predictor)] * df_full['{0}'.format(predictor2)]
            #for predictor3 in predictors[int(index + 2):]:
             #   if predictor2 != predictor3:
              #      df_full['{0}_x_{1}_x_{2}'.format(predictor, predictor2, predictor3)] = df_full['{0}'.format(predictor)] * df_full['{0}'.format(predictor2)] * df_full['{0}'.format(predictor3)]
    index += 1

# output processed data
df_full.to_csv('poly_interactions.csv')
df_full.dtypes.to_csv('poly_interactions_columntypes.txt')

# find index of first polynomial term
firstpoly = df_full.columns.get_loc('DX_Final_Rate') + 1

# remove rows where diagnosis changes
df_full_removed = df_full[df_full['DX_Final_Progression'] != 0]
df_full_removed.to_csv('df_full_removed.csv')
df_full = df_full[df_full['DX_Final_Progression'] == 0]

#original genes and new poly/interaction terms
temp_genes = used_genes.copy()
full_genes = temp_genes + list(df_full.columns[firstpoly:].values)
X_full = df_full[full_genes]
y_full = df_full['DX']
X_full_removed = df_full_removed[full_genes]
y_full_removed = df_full_removed['DX']
X_full_removed.to_csv('X_full_removed.csv')
y_full_removed.to_csv('y_full_removed.csv')

# split into train and test
np.random.seed(9001)
msk = np.random.rand(len(X_full)) < 0.5
X_train = X_full[msk]
X_test = X_full[~msk]
y_train = y_full[msk]
y_test = y_full[~msk]

# making classifiers to test our model against
# random classifier
yrand_train = randint(3, size = len(X_train))
yrand_test = randint(3, size = len(X_test))
# all zeros classifier
zeros = {}
zeros['train_pred'] = [0] * len(X_train)
zeros['test_pred'] = [0] * len(X_test)
zeros['model'] = np.zeros(len(X_test))
# all ones classifier
ones = {}
ones['train_pred'] = [1] * len(X_train)
ones['test_pred'] = [1] * len(X_test)
ones['model'] = np.ones(len(X_test))

# all twos classifier
twos = {}
twos['train_pred'] = [2] * len(X_train)
twos['test_pred'] = [2] * len(X_test)
twos['model'] = np.ones(len(X_test))*2

pre_X_train = X_train
pre_X_test = X_test
pre_X_train['class'] = y_train
pre_X_test['class'] = y_test

# PCA with all components
pca_full_fit = PCA()
pca_full = {}
pca_full['Xtrain'] = pca_full_fit.fit_transform(pre_X_train)

# find number of components that explain 90% of predictor variance
n_components = (np.argwhere((np.cumsum(pca_full_fit.explained_variance_ratio_)) > 0.9)[0] + 1)[0]

# PCA with 90% of variance explained
pca90 = {}
pca90_fit = PCA(n_components, random_state = 9001)
pca90_fit.fit(pre_X_train)
X_train = pca90_fit.transform(pre_X_train)
X_test = pca90_fit.transform(pre_X_test)

# run cv on multiple parameters in gradient boosting
param_dict = OrderedDict(
    max_depth = range(1,20),
    n_estimators = range(1,20),
    learning_rate = np.arange(0.05,1,0.05)
)
est = GradientBoostingClassifier(random_state = 9001)
gb_cv = GridSearchCV(est, param_grid = param_dict, cv=3, n_jobs=-1)
gb_cv.fit(X_train, y_train)
opt_depth = gb_cv.best_estimator_.max_depth
opt_n_est = gb_cv.best_estimator_.n_estimators
opt_lr = gb_cv.best_estimator_.learning_rate

# output accuracy for our new method, and comparison to benchmark models
print('Random class. accuracy, train: ', accuracy_score(y_train, yrand_train))
print('Random class. accuracy, test: ', accuracy_score(y_test, yrand_test))
print('All zeros class. accuracy, train: ', accuracy_score(y_train, zeros['train_pred']))
print('All zeros class. accuracy, test: ', accuracy_score(y_test, zeros['test_pred']))
print('All ones class. accuracy, train: ', accuracy_score(y_train, ones['train_pred']))
print('All ones class. accuracy, test: ', accuracy_score(y_test, ones['test_pred']))
print('All twos class. accuracy, train: ', accuracy_score(y_train, twos['train_pred']))
print('All twos class. accuracy, test: ', accuracy_score(y_test, twos['test_pred']))
print('Gradient Boost with PCA and CV class. accuracy, train: ', metrics.accuracy_score(y_train, gb_cv.best_estimator_.fit(X_train, y_train).predict(X_train)))
print('Gradient Boost with PCA and CV class. accuracy, test: ', metrics.accuracy_score(y_test, gb_cv.best_estimator_.fit(X_train, y_train).predict(X_test)))
```

    99.481865%

    /Users/Nick/anaconda/lib/python3.6/site-packages/ipykernel_launcher.py:96: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
    /Users/Nick/anaconda/lib/python3.6/site-packages/ipykernel_launcher.py:97: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
    

    Random class. accuracy, train:  0.305882352941
    Random class. accuracy, test:  0.373684210526
    All zeros class. accuracy, train:  0.411764705882
    All zeros class. accuracy, test:  0.415789473684
    All ones class. accuracy, train:  0.4
    All ones class. accuracy, test:  0.421052631579
    All twos class. accuracy, train:  0.188235294118
    All twos class. accuracy, test:  0.163157894737
    Gradient Boost with PCA and CV class. accuracy, train:  0.888235294118
    Gradient Boost with PCA and CV class. accuracy, test:  0.421052631579
    


```python
np.linspace(0.53333333, 0.66666667, 7)
np.random.seed(9001)
df = pd.read_csv('poly_interactions.csv', index_col = 0, dtype = {'DX': np.int32})
cols = df.dtypes
cols.to_csv('columntypes.txt')
df['DX'] = df['DX'] + 1
firstgene = df.columns.get_loc('DX') + 1
middle_start = df.columns.get_loc('SubjectID') - 1
end_start = df.columns.get_loc('DX_Final_Rate') + 1
df_start = df.iloc[:,firstgene:(middle_start+1)]
df_end = df.iloc[:,end_start:]
df = pd.concat([df_start, df_end, df.DX], axis = 1)
msk = np.random.rand(len(df)) < 0.5
data_train = df[msk]
data_test = df[~msk]
data = {}
models = {}
data['alz'] = {'xtrain' : data_train.drop('DX', axis = 1).values,
               'ytrain' : data_train['DX'].values,
               'xtest' : data_test.drop('DX', axis = 1).values,
               'ytest' : data_test['DX'].values}
```


```python
# define classifiction accuracy formula
def class_acc(y, y_pred):
    mis = 0
    length = len(y_pred)
    for prediction in range(0,len(y)):
        if y_pred[prediction] == 0:
            length -= 1
        elif y[prediction] != y_pred[prediction]:
            mis += 1
    score = (1 - mis/len(y_pred)) #calculate 1 - misclassification rate
    return score
```


```python
"""
Function
--------
calc_cost

Inputs
------
y = the true class
y_pred = the predicted class (or abstain)
ab_cost = the cost of abstaining
mis_cost = the cost of a misclassification
ab_class = what class is abstain?

Returns
-------
the cost of this model on the test data per patient
note that a misclassification = 5000
and an abstain (3) = 1000

"""
def calc_cost(y, y_pred, ab_cost = 6850, mis_cost = 12000, ab_class = 0):
    mis = 0
    ab = 0
    for prediction in range(0,len(y)):
        if y[prediction] != y_pred[prediction] and y_pred[prediction] != ab_class:
            mis += 1
        if y_pred[prediction] == ab_class:
            ab += 1
    return ((ab * ab_cost + mis * mis_cost)/len(y))
```


```python
"""
Function
--------
abstain

Inputs
------
preds = a list of the the predictions
probs = a list of the probability for each prediciton
mci_thresh = the threshold for mci
dementia_thresh = the threshold for dementia
cn_thresh = the threshold for cn

Returns
-------
ab_pred = the predictions with abstains

"""
def abstain(preds, probs, mci_thresh = 0, dementia_thresh = 0, cn_thresh = 0):
    cog_norm = 1
    mci = 2
    dementia = 3
    ab = 0
    
    ab_pred = []
    
    for i in range(0,len(preds)):
        if preds[i] == cog_norm and probs[i] < cn_thresh or preds[i] == mci and probs[i] < mci_thresh or preds[i] == dementia and probs[i] < dementia_thresh:
            ab_pred.append(ab) #abstain if the probability isn't high enough
        else:
            ab_pred.append(preds[i]) #accept the original probability
    
    return ab_pred
```


```python
"""
Function
--------
ab_model

Inputs
------
name = the name to use for the model
data_ = a dictionary of the data to use to train and test
    which should contain 'xtest' 'xtrain' 'ytest' and 'ytrain'
    if using with a validation set, make sure the data
    passed in is still labelled as 'xtest' and 'ytest'
    (i.e. not 'xvalid' and 'yvalid')
mci_thresh = the threshold for abstaining from mci
    a lower threshold = less abstaining
dementia_thresh = threshold for abstaining from dementia
cn_thresh = threshold for abstaining from cog_normal

Returns
-------
a dictionary containing this model's
(1) mci, dementia, and cog_normal thresholds
(2) the cost of this model
(3) the acc of this model (for those on which it does not abstain)


"""
def ab_model_original(name, data_, mci_thresh = 0, dementia_thresh = 0, cn_thresh = 0): 
    bm = gb_cv.best_estimator_
    models[name] = bm
    models[name].fit(data[name]['xtrain'], data[name]['ytrain'])

    y_pred = models[name].predict(data[name]['xtest'])
    y_probs = models[name].predict_proba(data[name]['xtest'])
    y_prob = []
    
    for i in range(0,len(y_pred)):
        y_prob.append(y_probs[i][y_pred[i] - 1])
    
    y_ab = abstain(y_pred, y_prob, mci_thresh = mci_thresh, dementia_thresh = dementia_thresh, cn_thresh = cn_thresh)
    ab_percent = (np.count_nonzero(y_ab))/len(y_ab)
    
    cost = calc_cost(data[name]['ytest'], y_ab)
    acc = class_acc(data[name]['ytest'], y_ab)
    
    ab_models[name] = {'model' : models[name],
                       'mci_thresh' : mci_thresh,
                       'dementia_thresh' : dementia_thresh,
                       'cn_thresh' : cn_thresh,
                       'raw predictions' : y_pred,
                       'predictions' : y_ab,
                       'probabilities' : y_prob,
                       'cost' : cost,
                       'acc' : acc,
                       'percent' : ab_percent,
                      }
    
    return ab_models[name]
```


```python
"""
Function
--------
ab_model

Inputs
------
name = the name to use for the model
data_ = a dictionary of the data to use to train and test
    which should contain 'xtest' 'xtrain' 'ytest' and 'ytrain'
    if using with a validation set, make sure the data
    passed in is still labelled as 'xtest' and 'ytest'
    (i.e. not 'xvalid' and 'yvalid')
mci_thresh = the threshold for abstaining from mci
    a lower threshold = less abstaining
dementia_thresh = threshold for abstaining from dementia
cn_thresh = threshold for abstaining from cog_normal

Returns
-------
a dictionary containing this model's
(1) mci, dementia, and cog_normal thresholds
(2) the cost of this model
(3) the acc of this model (for those on which it does not abstain)


"""
def ab_model(name, y_pred, y_prob, mci_thresh, dementia_thresh, cn_thresh): 
    bm = gb_cv.best_estimator_
    models[name] = bm
    
    y_ab = abstain(y_pred, y_prob, mci_thresh = mci_thresh, dementia_thresh = dementia_thresh, cn_thresh = cn_thresh)
    ab_percent = (np.count_nonzero(y_ab))/len(y_ab)
    
    cost = calc_cost(data[name]['ytest'], y_ab)
    acc = class_acc(data[name]['ytest'], y_ab)
    
    ab_model_dict = {'model' : models[name],
                       'mci_thresh' : mci_thresh,
                       'dementia_thresh' : dementia_thresh,
                       'cn_thresh' : cn_thresh,
                       'raw predictions' : y_pred,
                       'predictions' : y_ab,
                       'probabilities' : y_prob,
                       'cost' : cost,
                       'acc' : acc,
                       'percent' : ab_percent,
                      }
    
    return ab_model_dict
```


```python
#now we're going to actually fine-tune the threshold values
#using reasonable estimates from the above examination
#of the data
#note that this can be used to find threshold values for
#literally any combination of abstain costs and misclassification costs
#since you would simply need to adjust the cost in the function
#calc_cost above as that is used for the accuracy

#define lists to use for the different thresholds checks to get a model that Werks
iteration_max = 5
models = {}

#get best model
bm = gb_cv.best_estimator_

index = 0
n_models = 5
fold = 0
ab_models = []

#first make all of the models
#splits the training and valid
#note that it does the same every time
for train, valid in KFold(n_models, shuffle = True, random_state = 9001).split(data_train):
    #prepare the train and valid data
    fold_name = 'fold {0}'.format(fold) #name of the model
    fold += 1
    data[fold_name] = {'xtrain' : data_train.drop('DX', axis = 1).iloc[train].values,
                                  'xtest' : data_train.drop('DX', axis = 1).iloc[valid].values,
                                  'ytrain' : data_train['DX'].iloc[train].values,
                                  'ytest' : data_train['DX'].iloc[valid].values,
                      }
    #make and fit the model
    models[fold_name] = bm 
    models[fold_name].fit(data[fold_name]['xtrain'], data[fold_name]['ytrain'])
    #predict and get probabilities
    y_pred = models[fold_name].predict(data[fold_name]['xtest'])
    y_probs = models[fold_name].predict_proba(data[fold_name]['xtest'])
    
    #get the probabilities corresponding to the class of each prediction
    #i.e. if the first prediction is "1", then get the probability of prediction of class "1",
    #discarding "0" and "2" for the first prediction
    #since the original probabilities are stored in a matrix
    y_prob = []
    for i in range(0,len(y_pred)):
        y_prob.append(y_probs[i][y_pred[i] - 1])
    
    temp_model = {'model' : models[fold_name],
                  'name' : fold_name,
                  'predictions' : y_pred,
                  'probabilities' : y_prob,
                 }
    
    ab_models.append(temp_model)
    index = for_loop_status(n_models, index)
```

    80.000000%


```python
#now we will use the above fitted models for the threshold!
start_cn = 0
stop_cn = 1
start_mci = 0
stop_mci = 1
start_dementia = 0.6
stop_dementia = 1
n_thresh = 11
cn_thresholds = list(np.linspace(start_cn, stop_cn, n_thresh))
mci_thresholds = list(np.linspace(start_mci, stop_mci, n_thresh))
dementia_thresholds = list(np.linspace(start_dementia, stop_dementia, n_thresh))
index = 0
print('done!\t\tmci\tdem\tcn\tcost\tacc\tpercent')
length = length = len(cn_thresholds) * len(mci_thresholds) * len(dementia_thresholds) * (iteration_max)
for iteration in range(1,iteration_max + 1):
    
    thresh_models = {} #prepare dictionary to hold the separate models

    #iterates through every single possibility
    for cn_thresh in cn_thresholds:
        for mci_thresh in mci_thresholds:
            for dementia_thresh in dementia_thresholds:
                temp_costs = []
                temp_accs = []
                temp_percents = []
                thresh_name = 'threshold set {}'.format(index) 
                for ab_cv_model in ab_models:
                    temp_thresh_model = (ab_model(ab_cv_model['name'], ab_cv_model['predictions'], ab_cv_model['probabilities'], mci_thresh = mci_thresh, dementia_thresh = dementia_thresh, cn_thresh = cn_thresh)) #get the model
                    temp_costs.append(temp_thresh_model['cost']) #get the cost
                    temp_accs.append(temp_thresh_model['acc']) #get the acc
                    temp_percents.append(temp_thresh_model['percent']) #get abstainp ercent
                avg_cost = sum(temp_costs)/len(temp_costs) #average the costs
                avg_acc = sum(temp_accs)/len(temp_accs) #average the accs
                avg_percent = sum(temp_percents)/len(temp_percents) #get percent abstained
                
                #make an entry
                thresh_models[thresh_name] = {'mci_thresh' : mci_thresh,
                                              'dementia_thresh' : dementia_thresh,
                                              'cn_thresh' : cn_thresh,
                                              'cost' : avg_cost,
                                              'acc' : avg_acc,
                                              'percent' : avg_percent,
                                             }
                index = for_loop_status(length, index)
                
    #find the model with the minimum cost
    #and store the name of that model in min_name
    min_cost = 100000000 #arbitrarily high value
    min_name = ''

    for key, model in thresh_models.items():
        if min_cost > model['cost']: #if it's lower than the current minimum
            min_cost = model['cost']
            min_name = key

    print('\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.2f}\t{4:.4f}\t{5:.4f}'.format(thresh_models[min_name]['mci_thresh'],
                                                                                  thresh_models[min_name]['dementia_thresh'],
                                                                                 thresh_models[min_name]['cn_thresh'],
                                                                                 thresh_models[min_name]['cost'],
                                                                                 thresh_models[min_name]['acc'],
                                                                                 thresh_models[min_name]['percent']))
    #create thresholds based on cost-minimizing model
    thresh_min_dementia = thresh_models[min_name]['dementia_thresh']
    thresh_min_mci = thresh_models[min_name]['mci_thresh']
    thresh_min_cn = thresh_models[min_name]['cn_thresh']
    cn_thresh_index = cn_thresholds.index(thresh_min_cn)
    mci_thresh_index = mci_thresholds.index(thresh_min_mci)
    dementia_thresh_index = dementia_thresholds.index(thresh_min_dementia)
    
    #prevent thresholds from going above 1 or below 0
    if cn_thresholds.index(thresh_min_cn) == 0:
        cn_thresh_index = 1
    elif cn_thresholds.index(thresh_min_cn) == (len(cn_thresholds) - 1):
        cn_thresh_index = len(cn_thresholds) - 2
    if mci_thresholds.index(thresh_min_mci) == 0:
        mci_thresh_index = 1
    elif mci_thresholds.index(thresh_min_mci) == (len(mci_thresholds) - 1):
        mci_thresh_index = len(mci_thresholds) - 2
    if dementia_thresholds.index(thresh_min_dementia) == 0:
        dementia_thresh_index = 1
    elif dementia_thresholds.index(thresh_min_dementia) == (len(dementia_thresholds) - 1):
        dementia_thresh_index = len(dementia_thresholds) - 2
        
    #change thresholds to be one decimal point more
    start_cn = cn_thresholds[cn_thresh_index - 1]
    stop_cn = cn_thresholds[cn_thresh_index + 1]
    start_mci = mci_thresholds[mci_thresh_index - 1]
    stop_mci = mci_thresholds[mci_thresh_index + 1]
    start_dementia = dementia_thresholds[dementia_thresh_index - 1]
    stop_dementia = dementia_thresholds[dementia_thresh_index + 1]
    
    #create new thresholds    
    cn_thresholds = list(np.linspace(start_cn, stop_cn, n_thresh))
    mci_thresholds = list(np.linspace(start_mci, stop_mci, n_thresh))
    dementia_thresholds = list(np.linspace(start_dementia, stop_dementia, n_thresh))
```

    done!		mci	dem	cn	cost	acc	percent
    19.984974%	0.7000	1.0000	1.0000	6173.51	0.8502	0.3612
    39.984974%	0.7400	1.0000	0.9800	6133.70	0.8642	0.3424
    59.984974%	0.7320	0.9968	0.9800	6133.70	0.8642	0.3424
    79.984974%	0.7312	0.9958	0.9784	6133.70	0.8642	0.3424
    99.984974%	0.7312	0.9958	0.9781	6133.70	0.8642	0.3424
    


```python
#make the final model using the threshold values above
#and find its cost
name = 'final model'
models[name] = bm 
data[name] = data['alz']
models[name].fit(data['alz']['xtrain'], data['alz']['ytrain'])
y_pred = models[name].predict(data['alz']['xtest'])
y_probs = models[name].predict_proba(data['alz']['xtest'])
    
#get the probabilities corresponding to the class of each prediction
#i.e. if the first prediction is "1", then get the probability of prediction of class "1",
#discarding "0" and "2" for the first prediction
#since the original probabilities are stored in a matrix
y_prob = []
for i in range(0,len(y_pred)):
    y_prob.append(y_probs[i][y_pred[i] - 1])
final_model = ab_model(name, y_pred, y_prob, thresh_models[min_name]['mci_thresh'], thresh_models[min_name]['dementia_thresh'], thresh_models[min_name]['cn_thresh'])
print('The cost per patient for the model set with thresholds is: ${0:.2f} and classification accuracy is {1:.2f} and % diagnosed is {2:.2f}.'.format(final_model['cost'], final_model['acc'], final_model['percent']))
```

    The cost per patient for the model set with thresholds is: $6776.32 and classification accuracy is 0.89 and % diagnosed is 0.21.
    
