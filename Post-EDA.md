

```python
cd ~/Desktop/proj_cs109/Data
```

    /Users/Nick/Desktop/proj_cs109/Data
    


```python
%matplotlib inline
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LassoCV
from sklearn.linear_model import RidgeCV
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from numpy.random import randint
from sklearn.metrics import mean_squared_error
import seaborn as sns
import sklearn.metrics as metrics
from sklearn.metrics import r2_score
from sklearn.metrics import auc
from sklearn.preprocessing import PolynomialFeatures
```


```python
# define classifiction accuracy formula
def class_ac(y, y_pred):
    mis = 0
    for prediction in range(0,len(y)):
        if y[prediction] != y_pred[prediction]:
            mis += 1
    score = (1 - mis/len(y_pred)) #calculate 1 - misclassification rate
    return score
```


```python
#load dataset
df_full = pd.read_csv('FINAL_GENE_EXPRESSION.csv', index_col = 0)

#drop columns with high missingness
df_full = df_full.drop(['EcogPtMem', 'EcogPtLang', 'EcogPtVisspat', 'EcogPtPlan', 'EcogPtOrgan', 'EcogPtDivatt', 'EcogPtTotal', 'RAVLT_forgetting'], axis = 1)
```


```python
#locate rows with missing values (found earlier)
missing_rows = df_full.apply(lambda x: 49431-x.count(), axis=0)
```


```python
#drop rows with DX missing
df_full = df_full.dropna(subset = ['DX'])  
```


```python
#get names of missing columns in each row
lst = []
for index, row in df_full.iterrows():
    mask = row.isnull()
    lst += [row[mask].index.tolist()]
missing_list = []
for listy in lst:
    if len(listy) != 0:
        missing_list.append(listy)
```


```python
#calculate number of missing values in each row
missing_rows = df_full.isnull().sum(axis=1)
missing_r = pd.DataFrame(missing_rows, columns = ['Count'])
missing_r = missing_r[missing_r['Count'] > 0]
missing_r['Missing_vals'] = missing_list
missing_r = missing_r[missing_r['Count'] > 1]
```


```python
#set alpha values
alphas = np.zeros(11)
for i in np.arange(-5,6,1):
    alphas[i] = 10.0** i
```


```python
#imputation through a model
#adapted from lab7
df_miss = pd.read_csv('FINAL_GENE_EXPRESSION.csv', index_col = 0)
df_miss = df_miss[['CDRSB', 'ADAS11', 'ADAS13',
       'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting',
       'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang', 'EcogSPVisspat', 'EcogSPPlan',
       'EcogSPOrgan', 'EcogSPDivatt', 'EcogSPTotal', 'DX']]
df_miss_cut = df_miss.dropna(how='all')
df_filled = df_miss.copy()
#using the intact data
#build a model to use to impute
#onto the dataset that has missing values
#for each given predictors
df_miss_cut = df_miss_cut.drop('DX', axis =1 )
null_predictors = df_miss_cut.columns[df_miss_cut.isnull().any()] #get null predictors
print(null_predictors)
df_miss_cut.shape
```

    Index(['ADAS11', 'ADAS13', 'RAVLT_immediate', 'RAVLT_learning',
           'RAVLT_perc_forgetting', 'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang',
           'EcogSPVisspat', 'EcogSPPlan', 'EcogSPOrgan', 'EcogSPDivatt',
           'EcogSPTotal'],
          dtype='object')
    




    (441, 16)




```python
#imputation using LASSO
#define model dict
model = {}

#impute
for predictor in null_predictors:
    print('Now working on {}.'.format(predictor))
    #first the missing values
    #from the dataset
    null_miss = df_miss_cut[pd.DataFrame(data = df_miss_cut[predictor], index = df_miss_cut.index).isnull().any(axis=1)] #get all rows with missing values in that predict
    #and the missing index from earlier
    xmiss = null_miss.drop(null_predictors, axis=1)
    miss_index = df_miss_cut[predictor][df_miss_cut[predictor].isnull()].index #index of missing rows
    #now use the intact database
    #to create the model
    xnomiss = df_miss_cut.drop(null_predictors, axis=1)
    ynomiss = pd.DataFrame(data = df_miss_cut[predictor], index = df_miss_cut[predictor].index) #note that ynomiss is looking at that particular predictor, not "is_cancer
    xnomiss = xnomiss.drop(miss_index)
    ynomiss = ynomiss.drop(miss_index)
    #create the model
    model['Lasso'] = LassoCV(alphas = alphas) #linear regression
    model['Lasso'].fit(xnomiss, ynomiss) #fit
    ynomiss_pred = model['Lasso'].predict(xnomiss) #predict on not-missing data
    print('The r2 score is {}.'.format(r2_score(ynomiss.values, ynomiss_pred)))
    ymiss_pred = model['Lasso'].predict(xmiss) #predict on missing data
    #now include noise in the model
    #using mean squared error
    #calculated by the not-missing data
    #so we know how generally accurate
    #the model is and can thus
    #add appropriate noise
    ymiss_pred_noisy = ymiss_pred + np.random.normal(loc = 0, scale = np.sqrt(mean_squared_error(ynomiss.values, ynomiss_pred)), size = ymiss_pred.shape)
    #put together the imputed values
    miss_series = pd.Series(data = ymiss_pred.flatten(), index = miss_index)
    df_filled[predictor] = df_filled[predictor].fillna(miss_series)
null_predictors_check = df_filled.columns[df_filled.isnull().any()] #get null predictors
print(null_predictors_check) #check to make sure this is empty!
```

    Now working on ADAS11.
    The r2 score is 0.6832152845413335.
    Now working on ADAS13.
    The r2 score is 0.6822339396037445.
    Now working on RAVLT_immediate.
    The r2 score is 0.39138497985134724.
    Now working on RAVLT_learning.
    The r2 score is 0.2021633388767654.
    Now working on RAVLT_perc_forgetting.
    The r2 score is 0.2684739057011736.
    Now working on FAQ.
    The r2 score is 0.7531252677282848.
    Now working on MOCA.
    

    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    

    The r2 score is 0.5850571375652185.
    Now working on EcogSPMem.
    The r2 score is 0.5300661421776194.
    Now working on EcogSPLang.
    The r2 score is 0.48388589143871863.
    Now working on EcogSPVisspat.
    The r2 score is 0.5554280667685116.
    Now working on EcogSPPlan.
    The r2 score is 0.5508304299920903.
    Now working on EcogSPOrgan.
    The r2 score is 0.5145290923953789.
    Now working on EcogSPDivatt.
    The r2 score is 0.4639033601271356.
    Now working on EcogSPTotal.
    The r2 score is 0.6189780074542853.
    Index([], dtype='object')
    

    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/Nick/anaconda/lib/python3.6/site-packages/sklearn/linear_model/coordinate_descent.py:1082: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    


```python
#now do the split and do the logistic regression as usual
np.random.seed(9001)
df_filled = df_filled.dropna(how='all')
msk = np.random.rand(len(df_filled)) < 0.75
train_df = df_filled[msk]
test_df = df_filled[~msk]
#prepare train and tests
xtrain = train_df.drop('DX', axis=1).values
ytrain = train_df['DX'].values
xtest = test_df.drop('DX', axis=1).values
ytest = test_df['DX'].values
```


```python
#logistic regression fine-tuned via crossvalidation
#using the L2 penalty
model['log_lasso'] = LogisticRegressionCV(penalty='l2')
model['log_lasso'].fit(xtrain,ytrain) 
print('The classification accuracy of the fitted logistic regression is {0}.'.format(class_ac(ytest, model['log_lasso'].predict(xtest))))
yzeros = np.zeros(ytest.shape[0]) #all zeroes model
print('The classification accuracy of the all-zeros classifier is {0}.'.format(class_ac(ytest, yzeros)))
```

    The classification accuracy of the fitted logistic regression is 0.9173553719008265.
    The classification accuracy of the all-zeros classifier is 0.42148760330578516.
    


```python
#imputation through Ridge
#adapted from lab7
df_miss = pd.read_csv('FINAL_GENE_EXPRESSION.csv', index_col = 0)
df_miss = df_miss[['CDRSB', 'ADAS11', 'ADAS13',
       'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting',
       'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang', 'EcogSPVisspat', 'EcogSPPlan',
       'EcogSPOrgan', 'EcogSPDivatt', 'EcogSPTotal', 'DX']]
df_miss_cut = df_miss.dropna(how='all')
df_filled = df_miss.copy()
#using the intact data
#build a model to use to impute
#onto the dataset that has missing values
#for each given predictors
df_miss_cut = df_miss_cut.drop('DX', axis =1 )
null_predictors = df_miss_cut.columns[df_miss_cut.isnull().any()] #get null predictors
print(null_predictors)
df_miss_cut.shape
#impute
for predictor in null_predictors:
    print('Now working on {}.'.format(predictor))
    #first the missing values
    #from the dataset
    null_miss = df_miss_cut[pd.DataFrame(data = df_miss_cut[predictor], index = df_miss_cut.index).isnull().any(axis=1)] #get all rows with missing values in that predict
    #and the missing index from earlier
    xmiss = null_miss.drop(null_predictors, axis=1)
    miss_index = df_miss_cut[predictor][df_miss_cut[predictor].isnull()].index #index of missing rows
    #now use the intact database
    #to create the model
    xnomiss = df_miss_cut.drop(null_predictors, axis=1)
    ynomiss = pd.DataFrame(data = df_miss_cut[predictor], index = df_miss_cut[predictor].index) #note that ynomiss is looking at that particular predictor, not "is_cancer
    xnomiss = xnomiss.drop(miss_index)
    ynomiss = ynomiss.drop(miss_index)
    #create the model
    model['Ridge'] = RidgeCV(alphas = alphas) #linear regression
    model['Ridge'].fit(xnomiss, ynomiss) #fit
    ynomiss_pred = model['Ridge'].predict(xnomiss) #predict on not-missing data
    print('The r2 score is {}.'.format(r2_score(ynomiss.values, ynomiss_pred)))
    ymiss_pred = model['Ridge'].predict(xmiss) #predict on missing data
    #now include noise in the model
    #using mean squared error
    #calculated by the not-missing data
    #so we know how generally accurate
    #the model is and can thus
    #add appropriate noise
    ymiss_pred_noisy = ymiss_pred + np.random.normal(loc = 0, scale = np.sqrt(mean_squared_error(ynomiss.values, ynomiss_pred)), size = ymiss_pred.shape)
    #put together the imputed values
    miss_series = pd.Series(data = ymiss_pred.flatten(), index = miss_index)
    df_filled[predictor] = df_filled[predictor].fillna(miss_series)
null_predictors_check = df_filled.columns[df_filled.isnull().any()] #get null predictors
print(null_predictors_check) #check to make sure this is empty!
#now do the split and do the logistic regression as usual
np.random.seed(9001)
df_filled = df_filled.dropna(how='all')
msk = np.random.rand(len(df_filled)) < 0.75
train_df = df_filled[msk]
test_df = df_filled[~msk]
#prepare train and tests
xtrain = train_df.drop('DX', axis=1).values
ytrain = train_df['DX'].values
xtest = test_df.drop('DX', axis=1).values
ytest = test_df['DX'].values
#logistic regression fine-tuned via crossvalidation
#using the L2 penalty
model['log_ridge'] = LogisticRegressionCV(penalty='l2')
model['log_ridge'].fit(xtrain,ytrain) 
print('The classification accuracy of the fitted logistic regression is {0}.'.format(class_ac(ytest, model['log_ridge'].predict(xtest))))
yzeros = np.zeros(ytest.shape[0]) #all zeroes model
print('The classification accuracy of the all-zeros classifier is {0}.'.format(class_ac(ytest, yzeros)))
```

    Index(['ADAS11', 'ADAS13', 'RAVLT_immediate', 'RAVLT_learning',
           'RAVLT_perc_forgetting', 'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang',
           'EcogSPVisspat', 'EcogSPPlan', 'EcogSPOrgan', 'EcogSPDivatt',
           'EcogSPTotal'],
          dtype='object')
    Now working on ADAS11.
    The r2 score is 0.6875508934320097.
    Now working on ADAS13.
    The r2 score is 0.6842067261601028.
    Now working on RAVLT_immediate.
    The r2 score is 0.39267190524296536.
    Now working on RAVLT_learning.
    The r2 score is 0.20223000434719196.
    Now working on RAVLT_perc_forgetting.
    The r2 score is 0.26844855500150255.
    Now working on FAQ.
    The r2 score is 0.753149669397961.
    Now working on MOCA.
    The r2 score is 0.5849235641842113.
    Now working on EcogSPMem.
    The r2 score is 0.5323540409795313.
    Now working on EcogSPLang.
    The r2 score is 0.4874731056225644.
    Now working on EcogSPVisspat.
    The r2 score is 0.555398751695281.
    Now working on EcogSPPlan.
    The r2 score is 0.5508378516572242.
    Now working on EcogSPOrgan.
    The r2 score is 0.5173184554735832.
    Now working on EcogSPDivatt.
    The r2 score is 0.46628578294422185.
    Now working on EcogSPTotal.
    The r2 score is 0.6225011510433546.
    Index([], dtype='object')
    The classification accuracy of the fitted logistic regression is 0.9090909090909091.
    The classification accuracy of the all-zeros classifier is 0.42148760330578516.
    


```python
#imputation through linear regression
#adapted from lab7
df_miss = pd.read_csv('FINAL_GENE_EXPRESSION.csv', index_col = 0)
df_miss = df_miss[['CDRSB', 'ADAS11', 'ADAS13',
       'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting',
       'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang', 'EcogSPVisspat', 'EcogSPPlan',
       'EcogSPOrgan', 'EcogSPDivatt', 'EcogSPTotal', 'DX']]
df_miss_cut = df_miss.dropna(how='all')
df_filled = df_miss.copy()
#using the intact data
#build a model to use to impute
#onto the dataset that has missing values
#for each given predictors
df_miss_cut = df_miss_cut.drop('DX', axis =1 )
null_predictors = df_miss_cut.columns[df_miss_cut.isnull().any()] #get null predictors
print(null_predictors)
df_miss_cut.shape
#impute
for predictor in null_predictors:
    print('Now working on {}.'.format(predictor))
    #first the missing values
    #from the dataset
    null_miss = df_miss_cut[pd.DataFrame(data = df_miss_cut[predictor], index = df_miss_cut.index).isnull().any(axis=1)] #get all rows with missing values in that predict
    #and the missing index from earlier
    xmiss = null_miss.drop(null_predictors, axis=1)
    miss_index = df_miss_cut[predictor][df_miss_cut[predictor].isnull()].index #index of missing rows
    #now use the intact database
    #to create the model
    xnomiss = df_miss_cut.drop(null_predictors, axis=1)
    ynomiss = pd.DataFrame(data = df_miss_cut[predictor], index = df_miss_cut[predictor].index) #note that ynomiss is looking at that particular predictor, not "is_cancer
    xnomiss = xnomiss.drop(miss_index)
    ynomiss = ynomiss.drop(miss_index)
    #create the model
    model['Linear'] = LinearRegression() #linear regression
    model['Linear'].fit(xnomiss, ynomiss) #fit
    ynomiss_pred = model['Linear'].predict(xnomiss) #predict on not-missing data
    print('The r2 score is {}.'.format(r2_score(ynomiss.values, ynomiss_pred)))
    ymiss_pred = model['Linear'].predict(xmiss) #predict on missing data
    #now include noise in the model
    #using mean squared error
    #calculated by the not-missing data
    #so we know how generally accurate
    #the model is and can thus
    #add appropriate noise
    ymiss_pred_noisy = ymiss_pred + np.random.normal(loc = 0, scale = np.sqrt(mean_squared_error(ynomiss.values, ynomiss_pred)), size = ymiss_pred.shape)
    #put together the imputed values
    miss_series = pd.Series(data = ymiss_pred.flatten(), index = miss_index)
    df_filled[predictor] = df_filled[predictor].fillna(miss_series)
null_predictors_check = df_filled.columns[df_filled.isnull().any()] #get null predictors
print(null_predictors_check) #check to make sure this is empty!
#now do the split and do the logistic regression as usual
np.random.seed(9001)
df_filled = df_filled.dropna(how='all')
msk = np.random.rand(len(df_filled)) < 0.75
train_df = df_filled[msk]
test_df = df_filled[~msk]
#prepare train and tests
xtrain = train_df.drop('DX', axis=1).values
ytrain = train_df['DX'].values
xtest = test_df.drop('DX', axis=1).values
ytest = test_df['DX'].values
#logistic regression fine-tuned via crossvalidation
#using the L2 penalty
model['log_lin'] = LogisticRegressionCV(penalty='l2')
model['log_lin'].fit(xtrain,ytrain) 
print('The classification accuracy of the fitted logistic regression is {0}.'.format(class_ac(ytest, model['log_lin'].predict(xtest))))
yzeros = np.zeros(ytest.shape[0]) #all zeroes model
print('The classification accuracy of the all-zeros classifier is {0}.'.format(class_ac(ytest, yzeros)))
```

    Index(['ADAS11', 'ADAS13', 'RAVLT_immediate', 'RAVLT_learning',
           'RAVLT_perc_forgetting', 'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang',
           'EcogSPVisspat', 'EcogSPPlan', 'EcogSPOrgan', 'EcogSPDivatt',
           'EcogSPTotal'],
          dtype='object')
    Now working on ADAS11.
    The r2 score is 0.6878303984629959.
    Now working on ADAS13.
    The r2 score is 0.684485626711936.
    Now working on RAVLT_immediate.
    The r2 score is 0.39289383486922924.
    Now working on RAVLT_learning.
    The r2 score is 0.20249048997663088.
    Now working on RAVLT_perc_forgetting.
    The r2 score is 0.2686892288562157.
    Now working on FAQ.
    The r2 score is 0.753179103890844.
    Now working on MOCA.
    The r2 score is 0.5852076887985204.
    Now working on EcogSPMem.
    The r2 score is 0.532377400974708.
    Now working on EcogSPLang.
    The r2 score is 0.4874871929067146.
    Now working on EcogSPVisspat.
    The r2 score is 0.5554284491898449.
    Now working on EcogSPPlan.
    The r2 score is 0.5508647799918345.
    Now working on EcogSPOrgan.
    The r2 score is 0.5173512568531835.
    Now working on EcogSPDivatt.
    The r2 score is 0.46630500577535494.
    Now working on EcogSPTotal.
    The r2 score is 0.6225275914495195.
    Index([], dtype='object')
    The classification accuracy of the fitted logistic regression is 0.9090909090909091.
    The classification accuracy of the all-zeros classifier is 0.42148760330578516.
    


```python
#imputation through linear regression with polynomial features

#include polynomial terms
quad_terms = PolynomialFeatures(degree = 2) #include quadratic terms

#adapted from lab7
df_miss = pd.read_csv('FINAL_GENE_EXPRESSION.csv', index_col = 0)
df_miss = df_miss[['CDRSB', 'ADAS11', 'ADAS13',
       'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting',
       'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang', 'EcogSPVisspat', 'EcogSPPlan',
       'EcogSPOrgan', 'EcogSPDivatt', 'EcogSPTotal', 'DX']]
df_miss_cut = df_miss.dropna(how='all')
df_filled = df_miss.copy()
#using the intact data
#build a model to use to impute
#onto the dataset that has missing values
#for each given predictors
df_miss_cut = df_miss_cut.drop('DX', axis =1 )
null_predictors = df_miss_cut.columns[df_miss_cut.isnull().any()] #get null predictors
print(null_predictors)
df_miss_cut.shape
#impute
for predictor in null_predictors:
    print('Now working on {}.'.format(predictor))
    #first the missing values
    #from the dataset
    null_miss = df_miss_cut[pd.DataFrame(data = df_miss_cut[predictor], index = df_miss_cut.index).isnull().any(axis=1)] #get all rows with missing values in that predict
    #and the missing index from earlier
    xmiss = null_miss.drop(null_predictors, axis=1)
    xmiss = quad_terms.fit_transform(xmiss)
    miss_index = df_miss_cut[predictor][df_miss_cut[predictor].isnull()].index #index of missing rows
    #now use the intact database
    #to create the model
    xnomiss = df_miss_cut.drop(null_predictors, axis=1)
    ynomiss = pd.DataFrame(data = df_miss_cut[predictor], index = df_miss_cut[predictor].index) #note that ynomiss is looking at that particular predictor, not "is_cancer
    xnomiss = xnomiss.drop(miss_index)
    ynomiss = ynomiss.drop(miss_index)
    xnomiss = quad_terms.fit_transform(xnomiss)
    #create the model
    model['Poly'] = LinearRegression() #linear regression
    model['Poly'].fit(xnomiss, ynomiss) #fit
    ynomiss_pred = model['Poly'].predict(xnomiss) #predict on not-missing data
    print('The r2 score is {}.'.format(r2_score(ynomiss.values, ynomiss_pred)))
    ymiss_pred = model['Poly'].predict(xmiss) #predict on missing data
    #now include noise in the model
    #using mean squared error
    #calculated by the not-missing data
    #so we know how generally accurate
    #the model is and can thus
    #add appropriate noise
    ymiss_pred_noisy = ymiss_pred + np.random.normal(loc = 0, scale = np.sqrt(mean_squared_error(ynomiss.values, ynomiss_pred)), size = ymiss_pred.shape)
    #put together the imputed values
    miss_series = pd.Series(data = ymiss_pred.flatten(), index = miss_index)
    df_filled[predictor] = df_filled[predictor].fillna(miss_series)
null_predictors_check = df_filled.columns[df_filled.isnull().any()] #get null predictors
print(null_predictors_check) #check to make sure this is empty!
#now do the split and do the logistic regression as usual
np.random.seed(9001)
df_filled = df_filled.dropna(how='all')
msk = np.random.rand(len(df_filled)) < 0.75
train_df = df_filled[msk]
test_df = df_filled[~msk]
#prepare train and tests
xtrain = train_df.drop('DX', axis=1).values
ytrain = train_df['DX'].values
xtest = test_df.drop('DX', axis=1).values
ytest = test_df['DX'].values
#logistic regression fine-tuned via crossvalidation
#using the L2 penalty
model['log_poly'] = LogisticRegressionCV(penalty='l2')
model['log_poly'].fit(xtrain,ytrain) 
print('The classification accuracy of the fitted logistic regression is {0}.'.format(class_ac(ytest, model['log_poly'].predict(xtest))))
yzeros3 = np.zeros(ytest.shape[0]) #all zeroes model
print('The classification accuracy of the all-zeros classifier is {0}.'.format(class_ac(ytest, yzeros)))
```

    Index(['ADAS11', 'ADAS13', 'RAVLT_immediate', 'RAVLT_learning',
           'RAVLT_perc_forgetting', 'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang',
           'EcogSPVisspat', 'EcogSPPlan', 'EcogSPOrgan', 'EcogSPDivatt',
           'EcogSPTotal'],
          dtype='object')
    Now working on ADAS11.
    The r2 score is 0.7154929492672664.
    Now working on ADAS13.
    The r2 score is 0.7171846220305957.
    Now working on RAVLT_immediate.
    The r2 score is 0.44372227239655027.
    Now working on RAVLT_learning.
    The r2 score is 0.23365963126397093.
    Now working on RAVLT_perc_forgetting.
    The r2 score is 0.3204344253774243.
    Now working on FAQ.
    The r2 score is 0.7621861374078318.
    Now working on MOCA.
    The r2 score is 0.6008333153764035.
    Now working on EcogSPMem.
    The r2 score is 0.6466117970433768.
    Now working on EcogSPLang.
    The r2 score is 0.5293601837756068.
    Now working on EcogSPVisspat.
    The r2 score is 0.5710668891738866.
    Now working on EcogSPPlan.
    The r2 score is 0.5801879166443826.
    Now working on EcogSPOrgan.
    The r2 score is 0.5526735381250583.
    Now working on EcogSPDivatt.
    The r2 score is 0.5231586599060172.
    Now working on EcogSPTotal.
    The r2 score is 0.683518526659628.
    Index([], dtype='object')
    The classification accuracy of the fitted logistic regression is 0.9338842975206612.
    The classification accuracy of the all-zeros classifier is 0.42148760330578516.
    


```python
#drop poorly imputed (r^2 < .5) columns
df_filled = df_filled.drop(['RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting'], axis = 1)
```


```python
#now do the split and do the logistic regression for final imputed model
np.random.seed(9001)
df_filled = df_filled.dropna(how='all')
msk = np.random.rand(len(df_filled)) < 0.75
train_df = df_filled[msk]
test_df = df_filled[~msk]
#prepare train and tests
xtrain = train_df.drop('DX', axis=1).values
ytrain = train_df['DX'].values
xtest = test_df.drop('DX', axis=1).values
ytest = test_df['DX'].values
model['final'] = LogisticRegressionCV(penalty='l2')
model['final'].fit(xtrain,ytrain) 
print('The classification accuracy of the fitted logistic regression is {0}.'.format(class_ac(ytest, model['final'].predict(xtest))))
yzeros = np.zeros(ytest.shape[0]) #all zeroes model
print('The classification accuracy of the all-zeros classifier is {0}.'.format(class_ac(ytest, yzeros)))
```

    The classification accuracy of the fitted logistic regression is 0.9338842975206612.
    The classification accuracy of the all-zeros classifier is 0.42148760330578516.
    


```python
#drop imputed columns
df_full = df_full.drop(['RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting', 'CDRSB', 'ADAS11', 'ADAS13',
       'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting',
       'FAQ', 'MOCA', 'EcogSPMem', 'EcogSPLang', 'EcogSPVisspat', 'EcogSPPlan',
       'EcogSPOrgan', 'EcogSPDivatt', 'EcogSPTotal', 'DX'], axis = 1)
```


```python
#add imputed columns
df_final = df_filled.join(df_full)
```


```python
#export csv
df_final.to_csv('Post-Imputation.csv')
```


```python
#check for any null
df_final.isnull().values.any()
```




    True


