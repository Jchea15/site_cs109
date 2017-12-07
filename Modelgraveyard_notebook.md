---
title: Discarded models: Jupyter notebook
notebook: Model%20Graveyard.ipynb
---

[Download this notebook here.](https://raw.githubusercontent.com/Pagel56/site_cs109/master/notebooks/Model%20Graveyard.ipynb)

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
import sklearn.metrics as metrics
from sklearn.model_selection import cross_val_score
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import KFold
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn import discriminant_analysis

from my_functions import make_dict
```


```python
cd ~/Documents/GitHub/Data_Proj
```

    C:\Users\Jackie\Documents\GitHub\Data_Proj
    


```python
gene_corr_base = pd.read_csv('gene_corr_base.csv', index_col = None, header=None)
gene_corr_cut = pd.read_csv('gene_corr_cut.csv', index_col = None, header=None)
gene_corr_p = pd.read_csv('gene_corr_p.csv', index_col = None, header=None)
gene_names = pd.read_csv('gene_names.csv', index_col = None, header=None)
df_full = pd.read_csv('post_imputation.csv', index_col = 0)
```


```python
df_full.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CDRSB</th>
      <th>ADAS11</th>
      <th>ADAS13</th>
      <th>MMSE</th>
      <th>FAQ</th>
      <th>MOCA</th>
      <th>EcogSPMem</th>
      <th>EcogSPLang</th>
      <th>EcogSPVisspat</th>
      <th>EcogSPPlan</th>
      <th>...</th>
      <th>M</th>
      <th>update_stamp</th>
      <th>First_DX</th>
      <th>First_date</th>
      <th>Final_DX</th>
      <th>Final_date</th>
      <th>First_Delta_Time</th>
      <th>Final_Delta_Time</th>
      <th>DX_Final_Progression</th>
      <th>DX_Final_Rate</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>7</th>
      <td>0.0</td>
      <td>3.0</td>
      <td>4.0</td>
      <td>30.0</td>
      <td>1.0</td>
      <td>28.0</td>
      <td>1.125</td>
      <td>1.00000</td>
      <td>1.00000</td>
      <td>1.00</td>
      <td>...</td>
      <td>60</td>
      <td>2017-10-06 23:19:46.0</td>
      <td>0</td>
      <td>2011-06-16</td>
      <td>0</td>
      <td>2015-06-02</td>
      <td>0.0</td>
      <td>1447.0</td>
      <td>0.0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>14</th>
      <td>1.5</td>
      <td>16.0</td>
      <td>24.0</td>
      <td>28.0</td>
      <td>5.0</td>
      <td>20.0</td>
      <td>3.250</td>
      <td>3.66667</td>
      <td>1.28571</td>
      <td>2.25</td>
      <td>...</td>
      <td>0</td>
      <td>2017-10-06 23:19:54.0</td>
      <td>1</td>
      <td>2011-08-25</td>
      <td>2</td>
      <td>2013-10-29</td>
      <td>0.0</td>
      <td>796.0</td>
      <td>1.0</td>
      <td>0.001256</td>
    </tr>
    <tr>
      <th>22</th>
      <td>0.0</td>
      <td>7.0</td>
      <td>8.0</td>
      <td>28.0</td>
      <td>0.0</td>
      <td>21.0</td>
      <td>1.000</td>
      <td>1.11111</td>
      <td>1.00000</td>
      <td>1.00</td>
      <td>...</td>
      <td>0</td>
      <td>2017-10-06 23:19:55.0</td>
      <td>0</td>
      <td>2011-09-13</td>
      <td>0</td>
      <td>2015-09-24</td>
      <td>0.0</td>
      <td>1472.0</td>
      <td>0.0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>27</th>
      <td>1.0</td>
      <td>7.0</td>
      <td>9.0</td>
      <td>30.0</td>
      <td>0.0</td>
      <td>22.0</td>
      <td>1.000</td>
      <td>1.00000</td>
      <td>1.00000</td>
      <td>1.00</td>
      <td>...</td>
      <td>0</td>
      <td>2017-10-06 23:19:55.0</td>
      <td>1</td>
      <td>2011-09-27</td>
      <td>1</td>
      <td>2012-10-02</td>
      <td>0.0</td>
      <td>371.0</td>
      <td>0.0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>35</th>
      <td>0.0</td>
      <td>6.0</td>
      <td>11.0</td>
      <td>28.0</td>
      <td>0.0</td>
      <td>25.0</td>
      <td>1.250</td>
      <td>1.00000</td>
      <td>1.00000</td>
      <td>1.00</td>
      <td>...</td>
      <td>0</td>
      <td>2017-10-06 23:19:55.0</td>
      <td>0</td>
      <td>2011-10-04</td>
      <td>0</td>
      <td>2015-10-20</td>
      <td>0.0</td>
      <td>1477.0</td>
      <td>0.0</td>
      <td>0.000000</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 49435 columns</p>
</div>




```python
# just gene expression data
X_full = df_full.iloc[:,14:49400]
X_full.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>11715100_at</th>
      <th>11715101_s_at</th>
      <th>11715102_x_at</th>
      <th>11715103_x_at</th>
      <th>11715104_s_at</th>
      <th>11715105_at</th>
      <th>11715106_x_at</th>
      <th>11715107_s_at</th>
      <th>11715108_x_at</th>
      <th>11715109_at</th>
      <th>...</th>
      <th>AFFX-r2-TagH_at</th>
      <th>AFFX-r2-TagIN-3_at</th>
      <th>AFFX-r2-TagIN-5_at</th>
      <th>AFFX-r2-TagIN-M_at</th>
      <th>AFFX-r2-TagJ-3_at</th>
      <th>AFFX-r2-TagJ-5_at</th>
      <th>AFFX-r2-TagO-3_at</th>
      <th>AFFX-r2-TagO-5_at</th>
      <th>AFFX-r2-TagQ-3_at</th>
      <th>AFFX-r2-TagQ-5_at</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>7</th>
      <td>-0.004618</td>
      <td>0.696194</td>
      <td>-1.317020</td>
      <td>0.315948</td>
      <td>-0.034504</td>
      <td>-0.378302</td>
      <td>0.612318</td>
      <td>-0.203451</td>
      <td>-0.206527</td>
      <td>-0.301305</td>
      <td>...</td>
      <td>-0.334344</td>
      <td>-0.096450</td>
      <td>-0.054353</td>
      <td>0.639972</td>
      <td>-1.911019</td>
      <td>-0.309482</td>
      <td>-1.040793</td>
      <td>-0.536355</td>
      <td>-1.040062</td>
      <td>-0.749751</td>
    </tr>
    <tr>
      <th>14</th>
      <td>-1.200963</td>
      <td>-0.476742</td>
      <td>0.331783</td>
      <td>1.228517</td>
      <td>-0.655580</td>
      <td>-0.343305</td>
      <td>0.378144</td>
      <td>1.305417</td>
      <td>-1.739029</td>
      <td>1.347586</td>
      <td>...</td>
      <td>-0.154246</td>
      <td>0.594417</td>
      <td>0.998769</td>
      <td>0.032531</td>
      <td>-0.757878</td>
      <td>-1.838261</td>
      <td>-0.209825</td>
      <td>-0.108993</td>
      <td>0.539111</td>
      <td>0.087313</td>
    </tr>
    <tr>
      <th>22</th>
      <td>-0.419879</td>
      <td>-0.319452</td>
      <td>-0.338550</td>
      <td>0.713567</td>
      <td>-0.544255</td>
      <td>0.125656</td>
      <td>-0.095186</td>
      <td>0.399328</td>
      <td>-1.401597</td>
      <td>1.846088</td>
      <td>...</td>
      <td>-0.381326</td>
      <td>-0.711159</td>
      <td>-0.454836</td>
      <td>0.657837</td>
      <td>-0.397119</td>
      <td>1.068100</td>
      <td>-0.025166</td>
      <td>-0.240972</td>
      <td>-0.780795</td>
      <td>0.287931</td>
    </tr>
    <tr>
      <th>27</th>
      <td>-0.563242</td>
      <td>0.210841</td>
      <td>-0.192590</td>
      <td>0.782010</td>
      <td>-1.018851</td>
      <td>-0.511291</td>
      <td>-1.609842</td>
      <td>0.215039</td>
      <td>1.818063</td>
      <td>-0.509470</td>
      <td>...</td>
      <td>-0.608406</td>
      <td>-1.314987</td>
      <td>0.509290</td>
      <td>0.693569</td>
      <td>0.105368</td>
      <td>0.396109</td>
      <td>-0.080563</td>
      <td>-1.749311</td>
      <td>0.931940</td>
      <td>-1.510718</td>
    </tr>
    <tr>
      <th>35</th>
      <td>-1.922725</td>
      <td>-1.231736</td>
      <td>-1.333238</td>
      <td>-1.043128</td>
      <td>-0.456367</td>
      <td>0.650612</td>
      <td>-0.289500</td>
      <td>-2.088576</td>
      <td>0.074666</td>
      <td>-1.265440</td>
      <td>...</td>
      <td>-0.905960</td>
      <td>-0.787317</td>
      <td>-1.389296</td>
      <td>-0.414116</td>
      <td>-0.197413</td>
      <td>0.152512</td>
      <td>-0.966929</td>
      <td>0.645177</td>
      <td>0.067716</td>
      <td>-0.888108</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 49386 columns</p>
</div>




```python
corrs_df = pd.DataFrame()
corrs_df['Gene_Name'] = gene_names[0]
corrs_df['Base'] = gene_corr_base[0]/15
corrs_df['Cut'] = gene_corr_cut[0]/15
corrs_df['p'] = gene_corr_p[0]/15
corrs_df.sort_values('p').head(20)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Gene_Name</th>
      <th>Base</th>
      <th>Cut</th>
      <th>p</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>42374</th>
      <td>11757474_x_at</td>
      <td>0.234427</td>
      <td>0.234427</td>
      <td>0.000019</td>
    </tr>
    <tr>
      <th>31151</th>
      <td>11746251_x_at</td>
      <td>0.212479</td>
      <td>0.212479</td>
      <td>0.000037</td>
    </tr>
    <tr>
      <th>6909</th>
      <td>11722009_a_at</td>
      <td>0.198839</td>
      <td>0.198839</td>
      <td>0.000216</td>
    </tr>
    <tr>
      <th>47835</th>
      <td>11762935_x_at</td>
      <td>0.187445</td>
      <td>0.187445</td>
      <td>0.000357</td>
    </tr>
    <tr>
      <th>15665</th>
      <td>11730765_at</td>
      <td>0.174991</td>
      <td>0.174991</td>
      <td>0.000755</td>
    </tr>
    <tr>
      <th>10789</th>
      <td>11725889_at</td>
      <td>0.166186</td>
      <td>0.166186</td>
      <td>0.000824</td>
    </tr>
    <tr>
      <th>7203</th>
      <td>11722303_x_at</td>
      <td>0.169622</td>
      <td>0.169622</td>
      <td>0.000936</td>
    </tr>
    <tr>
      <th>4518</th>
      <td>11719618_a_at</td>
      <td>0.161189</td>
      <td>0.161189</td>
      <td>0.001258</td>
    </tr>
    <tr>
      <th>26090</th>
      <td>11741190_a_at</td>
      <td>0.160000</td>
      <td>0.160000</td>
      <td>0.001670</td>
    </tr>
    <tr>
      <th>6810</th>
      <td>11721910_at</td>
      <td>0.160999</td>
      <td>0.160999</td>
      <td>0.001846</td>
    </tr>
    <tr>
      <th>32506</th>
      <td>11747606_a_at</td>
      <td>0.151940</td>
      <td>0.151940</td>
      <td>0.001913</td>
    </tr>
    <tr>
      <th>46475</th>
      <td>11761575_at</td>
      <td>0.155691</td>
      <td>0.155691</td>
      <td>0.002033</td>
    </tr>
    <tr>
      <th>9183</th>
      <td>11724283_a_at</td>
      <td>0.162228</td>
      <td>0.162228</td>
      <td>0.002094</td>
    </tr>
    <tr>
      <th>10788</th>
      <td>11725888_at</td>
      <td>0.154584</td>
      <td>0.154584</td>
      <td>0.002121</td>
    </tr>
    <tr>
      <th>18119</th>
      <td>11733219_x_at</td>
      <td>0.153854</td>
      <td>0.153854</td>
      <td>0.002261</td>
    </tr>
    <tr>
      <th>32487</th>
      <td>11747587_a_at</td>
      <td>0.153094</td>
      <td>0.153094</td>
      <td>0.002452</td>
    </tr>
    <tr>
      <th>26628</th>
      <td>11741728_x_at</td>
      <td>0.152116</td>
      <td>0.152116</td>
      <td>0.002461</td>
    </tr>
    <tr>
      <th>12501</th>
      <td>11727601_a_at</td>
      <td>0.155935</td>
      <td>0.155935</td>
      <td>0.002570</td>
    </tr>
    <tr>
      <th>43104</th>
      <td>11758204_s_at</td>
      <td>0.156271</td>
      <td>0.156271</td>
      <td>0.002714</td>
    </tr>
    <tr>
      <th>22057</th>
      <td>11737157_x_at</td>
      <td>0.150977</td>
      <td>0.150977</td>
      <td>0.002723</td>
    </tr>
  </tbody>
</table>
</div>




```python
Used = corrs_df[corrs_df.p < 0.01]
used_genes = list(Used.Gene_Name.values)
```


```python
def for_loop_status(length, index = 0):
    sys.stdout.write('\r%f%%' % ((index/length)*100))
    sys.stdout.flush()
    index += 1
    return(index)

def for_loop_status2(length, length2, index = 0, index2 = 0):
    sys.stdout.write('\r%f%% of %f%%' % ((index/length)*100, (index2/length2)*100))
    sys.stdout.flush()
    index += 1
    index2 += 1
    return(index, index2)
```


```python
# interaction and polynomial terms
predictors_expo_list = []
for_index = 0
index = 0
for predictor in used_genes:
    for_index = for_loop_status(len(used_genes), for_index)
    predictors_temp_list = []
    for expo in [2,3,4]:
        #if expo = 2 and the predictor is atemp, then you'll get "atemp2" which is atemp^2
        predictors_temp_list.append('{0}{1}'.format(predictor, expo))
        df_full['{0}{1}'.format(predictor, expo)] = df_full['{0}'.format(predictor)]**int('{0}'.format(expo))
    predictors_expo_list.append(predictors_temp_list)
    for predictor2 in used_genes[int(index + 1):]:
        if predictor != predictor2:
            df_full['{0}_x_{1}'.format(predictor, predictor2)] = df_full['{0}'.format(predictor)] * df_full['{0}'.format(predictor2)]
    index += 1
```

    99.481865%


```python
# find index of first polynomial term
firstpoly = df_full.columns.get_loc('DX_Final_Rate') + 1
df_full.columns[firstpoly:]
```




    Index(['11715176_at2', '11715176_at3', '11715176_at4',
           '11715176_at_x_11715574_a_at', '11715176_at_x_11715665_a_at',
           '11715176_at_x_11716095_s_at', '11715176_at_x_11716149_a_at',
           '11715176_at_x_11716279_x_at', '11715176_at_x_11716702_a_at',
           '11715176_at_x_11717003_s_at',
           ...
           '11762935_x_at4', '11762935_x_at_x_11762998_x_at',
           '11762935_x_at_x_200064_PM_at', '11762998_x_at2', '11762998_x_at3',
           '11762998_x_at4', '11762998_x_at_x_200064_PM_at', '200064_PM_at2',
           '200064_PM_at3', '200064_PM_at4'],
          dtype='object', length=19107)




```python
df_full_removed = df_full[df_full['DX_Final_Progression'] != 0]
df_full = df_full[df_full['DX_Final_Progression'] == 0]
```


```python
temp_genes = used_genes.copy()
#significant genes and new poly/interaction terms
full_genes = temp_genes + list(df_full.columns[firstpoly:].values)
X_full = df_full[full_genes]
y_full = df_full['DX']
```


```python
# split into train and test
np.random.seed(9001)
msk = np.random.rand(len(X_full)) < 0.5

X_train = X_full[msk]
X_test = X_full[~msk]
y_train = y_full[msk]
y_test = y_full[~msk]
```

# PCA


```python
pre_X_train = X_full[msk]
pre_X_test = X_full[~msk]
pre_y_train = y_full[msk]
pre_y_test = y_full[~msk]

pre_X_train['class'] = pre_y_train
pre_X_test['class'] = pre_y_test
```

    C:\Users\Jackie\Anaconda3\lib\site-packages\ipykernel_launcher.py:6: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      
    C:\Users\Jackie\Anaconda3\lib\site-packages\ipykernel_launcher.py:7: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      import sys
    


```python
# PCA with all components
pca_full_fit = PCA()
pca_full = {}
pca_full['Xtrain'] = pca_full_fit.fit_transform(pre_X_train)

# find number of components that explain 90% of predictor variance
n_components = (np.argwhere((np.cumsum(pca_full_fit.explained_variance_ratio_)) > 0.9)[0] + 1)[0]

# PCA with n_components
pca90 = {}
pca90_fit = PCA(n_components)
pca90_fit.fit(pre_X_train)
X_train = pca90_fit.transform(pre_X_train)
X_test = pca90_fit.transform(pre_X_test)
y_train = pre_y_train
y_test = pre_y_test
```


```python
# random classifier
rand = {}
rand['train_pred'] = randint(3, size = len(X_train))
rand['test_pred'] = randint(3, size = len(X_test))
```

# AdaBoost with PCA and CV


```python
from collections import OrderedDict
from sklearn.model_selection import GridSearchCV

param_dict = OrderedDict(
    n_estimators = range(1,20),
    learning_rate = np.arange(0.05,1,0.05)
)

est = AdaBoostClassifier(random_state = 9001)
gb_cv = GridSearchCV(est, param_grid = param_dict, cv=3, n_jobs=-1)
gb_cv.fit(X_train, y_train)
gb_cv.best_estimator_
```




    AdaBoostClassifier(algorithm='SAMME.R', base_estimator=None,
              learning_rate=0.60000000000000009, n_estimators=19,
              random_state=9001)




```python
print('AdaBoost with PCA & CV, train: ', accuracy_score(y_train, gb_cv.best_estimator_.fit(X_train, y_train).predict(X_train)))
print('AdaBoost with PCA & CV, test: ', accuracy_score(y_test, gb_cv.best_estimator_.fit(X_train, y_train).predict(X_test)))
print('Random class. accuracy, train: ', accuracy_score(y_train, rand['train_pred']))
print('Random class. accuracy, test: ', accuracy_score(y_test, rand['test_pred']))
```

    AdaBoost with PCA & CV, train:  0.664705882353
    AdaBoost with PCA & CV, test:  0.421052631579
    Random class. accuracy, train:  0.264705882353
    Random class. accuracy, test:  0.310526315789
    

# Logistic Regression, Discriminant Analysis, kNN


```python
# fit logistic regression model
lin_class = {}
lin_class['logit_ovr'] = LogisticRegressionCV(Cs=7, penalty='l2', random_state = 9001, multi_class='ovr')
lin_class['logit_ovr'].fit(X_train, y_train)

# classify
lin_class['logit_ovr_train'] = lin_class['logit_ovr'].predict(X_train)
lin_class['logit_ovr_test'] = lin_class['logit_ovr'].predict(X_test)

# fit logistic regression model
lin_class['logit_multinomial'] = LogisticRegressionCV(Cs=7, penalty='l2', random_state = 9001, multi_class='multinomial')
lin_class['logit_multinomial'].fit(X_train, y_train)

# classify
lin_class['logit_multinomial_train'] = lin_class['logit_multinomial'].predict(X_train)
lin_class['logit_multinomial_test'] = lin_class['logit_multinomial'].predict(X_test)

# LDA
lda = {}
lda['model'] = discriminant_analysis.LinearDiscriminantAnalysis()
lda['model'].fit(X_train, y_train)
lda['train_pred'] = lda['model'].predict(X_train)
lda['test_pred'] = lda['model'].predict(X_test)

# QDA
qda = {}
qda['model'] = discriminant_analysis.QuadraticDiscriminantAnalysis()
qda['model'].fit(X_train, y_train)
qda['train_pred'] = qda['model'].predict(X_train)
qda['test_pred'] = qda['model'].predict(X_test)

#Fit a knn model
knn = {}

# knn cv code from https://kevinzakka.github.io/2016/07/13/k-nearest-neighbor/#parameter-tuning-with-cross-validation
knn['cv_scores'] = []
knn['neighbours'] = list(range(1,50))
for k in knn['neighbours']:
    current_model = KNN(n_neighbors=k)
    scores = cross_val_score(current_model, X_train, y_train, cv=10, scoring='accuracy')
    knn['cv_scores'].append(scores.mean())
knn['MSE'] = [1 - x for x in knn['cv_scores']]
# determining best k
knn['optimal_k'] = knn['neighbours'][knn['MSE'].index(min(knn['MSE']))]

knn['model'] = KNN(n_neighbors=knn['optimal_k'])
knn['model'].fit(X_train, y_train)
knn['train_pred'] = knn['model'].predict(X_train)
knn['test_pred'] = knn['model'].predict(X_test)
```

    C:\Users\Jackie\Anaconda3\lib\site-packages\sklearn\discriminant_analysis.py:695: UserWarning: Variables are collinear
      warnings.warn("Variables are collinear")
    


```python
# print all classification accuracies
print('Logistic OvR class. accuracy, train: ', accuracy_score(y_train, lin_class['logit_ovr_train']))
print('Logistic OvR class. accuracy, test: ', accuracy_score(y_test, lin_class['logit_ovr_test']))
print('Logistic multinomial class. accuracy, train: ', accuracy_score(y_train, lin_class['logit_multinomial_train']))
print('Logistic multinomial class. accuracy, test: ', accuracy_score(y_test, lin_class['logit_multinomial_test']))
print('LDA class. accuracy, train: ', accuracy_score(y_train, lda['train_pred']))
print('LDA class. accuracy, test: ', accuracy_score(y_test, lda['test_pred']))
print('QDA class. accuracy, train: ', accuracy_score(y_train, qda['train_pred']))
print('QDA class. accuracy, test: ', accuracy_score(y_test, qda['test_pred']))
print('KNN class. accuracy, train: ', accuracy_score(y_train, knn['train_pred']))
print('KNN class. accuracy, test: ', accuracy_score(y_test, knn['test_pred']))
print('Random class. accuracy, train: ', accuracy_score(y_train, rand['train_pred']))
print('Random class. accuracy, test: ', accuracy_score(y_test, rand['test_pred']))
```

    Logistic OvR class. accuracy, train:  0.635294117647
    Logistic OvR class. accuracy, test:  0.468421052632
    Logistic multinomial class. accuracy, train:  0.770588235294
    Logistic multinomial class. accuracy, test:  0.347368421053
    LDA class. accuracy, train:  0.611764705882
    LDA class. accuracy, test:  0.473684210526
    QDA class. accuracy, train:  0.888235294118
    QDA class. accuracy, test:  0.426315789474
    KNN class. accuracy, train:  0.723529411765
    KNN class. accuracy, test:  0.421052631579
    Random class. accuracy, train:  0.264705882353
    Random class. accuracy, test:  0.310526315789
    

## Decision Tree CV for max depth


```python
df_train = pd.read_csv('Xtrain.csv')
df_test = pd.read_csv('Xtest.csv')
```


```python
orig_names = df_train.columns.values
for value in orig_names:
    if value[0] == ' ':
        df_train = df_train.rename(columns = {value: value[1:]})
        df_test = df_test.rename(columns = {value: value[1:]})
```


```python
# extract matrices from dataframe
dtrain_X_full = df_train.drop('class', axis = 1)
#dtrain_X = df_train.iloc[:, :-8]
dtest_X_full = df_test.drop('class', axis = 1)
#dtest_X = df_test.iloc[:, :-8]

dt_class = {}
dt_class['Xtrain'] = dtrain_X_full.values
dt_class['Xtest'] = dtest_X_full.values
dt_class['ytrain'] = df_train['class'].values
dt_class['ytest'] = df_test['class'].values
```


```python
# plot training and test accuracies as function of tree depth
dt = {}
depths = list(range(2,11))
train_scores = []
test_scores = []
for i in depths:
    dt[str(i)] = DecisionTreeClassifier(max_depth = i)
    dt[str(i)].fit(dt_class['Xtrain'], dt_class['ytrain'])
    train_scores.append(dt[str(i)].score(dt_class['Xtrain'], dt_class['ytrain']))
    test_scores.append(dt[str(i)].score(dt_class['Xtest'], dt_class['ytest']))
plt.plot(depths, train_scores, "o-", label = "Training")
plt.plot(depths, test_scores, "o-", label = "Test")
plt.xlabel("Tree depth")
plt.ylabel("Accuracy")
plt.legend(loc='best')

# use 5-fold cross validation to find optimal tree depth
cv_scores = []
for i in depths:
    cv_scores.append(np.mean(cross_val_score(dt[str(i)], dt_class["Xtrain"], dt_class["ytrain"], cv=5)))
    cv_max = max(cv_scores)
    opt_depth = depths[cv_scores.index(cv_max)]
print("Optimal tree depth is %d with %f accuracy" %(opt_depth, cv_max))
max_depth = opt_depth

# refit tree w/optimal tree depth
cv_dt = DecisionTreeClassifier(max_depth = opt_depth)
cv_dt.fit(dt_class["Xtrain"], dt_class["ytrain"])

print('Tree class. accuracy, train: ', accuracy_score(dt_class["ytrain"], cv_dt.predict(dt_class['Xtrain'])))
print('Tree class. accuracy, test: ', accuracy_score(dt_class["ytest"], cv_dt.predict(dt_class['Xtest'])))
```

    Optimal tree depth is 9 with 0.347324 accuracy
    Tree class. accuracy, train:  1.0
    Tree class. accuracy, test:  0.363157894737
    


![png](Model%20Graveyard_files/Model%20Graveyard_27_1.png)
