Longitudinal attempt: Jupyter notebook
===================

[Download this notebook here.](https://raw.githubusercontent.com/Pagel56/site_cs109/master/notebooks/Longitudinal%20Attempt.ipynb)

```python
cd ~/Desktop/CS109/Data
```

    /Users/Leah/Desktop/CS109/Data
    


```python
%matplotlib inline
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
import statsmodels.api as sm
from statsmodels.api import OLS
from sklearn.metrics import r2_score

from my_functions import for_loop_status
from my_functions import for_loop_status2
```

    /anaconda/lib/python3.6/site-packages/statsmodels/compat/pandas.py:56: FutureWarning: The pandas.core.datetools module is deprecated and will be removed in a future version. Please use the pandas.tseries module instead.
      from pandas.core import datetools
    


```python
# import data
gene_corr_base = pd.read_csv('gene_corr_base.csv', index_col = None, header=None)
gene_corr_cut = pd.read_csv('gene_corr_cut.csv', index_col = None, header=None)
gene_corr_p = pd.read_csv('gene_corr_p.csv', index_col = None, header=None)
gene_names = pd.read_csv('gene_names.csv', index_col = None, header=None)
df_full = pd.read_csv('post_imputation.csv', index_col = 0)
# drop rows where response variable is NaN
df_full = df_full[np.isfinite(df_full['DX_Final_Rate'])]

# make correlations dataframe
corrs_df = pd.DataFrame()
corrs_df['Gene_Name'] = gene_names[0]
corrs_df['Base'] = gene_corr_base[0]/14
corrs_df['Cut'] = gene_corr_cut[0]/14
corrs_df['p'] = gene_corr_p[0]/14

# choose most significant genes
Used = corrs_df[corrs_df.p < 0.01]
used_genes = list(Used.Gene_Name.values)

# just gene expression data
X_full = df_full[used_genes]
y_full = df_full['DX_Final_Rate']

# get polynomial terms up to degree 3, since that is computationally feasible
gen_poly_terms = PolynomialFeatures(degree=3, interaction_only=False)
X_full_with_poly = gen_poly_terms.fit_transform(X_full)

# split into train and test
np.random.seed(9001)
msk = np.random.rand(len(X_full)) < 0.5
X_train = X_full[msk]
X_test = X_full[~msk]
y_train = y_full[msk]
y_test = y_full[~msk]

# PCA with all components
pca_full_fit = PCA()
pca_full = {}
pca_full['X_train'] = pca_full_fit.fit_transform(X_train)

# find number of components that explain 90% of predictor variance
n_components = (np.argwhere((np.cumsum(pca_full_fit.explained_variance_ratio_)) > 0.9)[0] + 1)[0]

# PCA with 90% of variance explained
pca90 = {}
pca90_fit = PCA(n_components)
pca90_fit.fit(X_train)
X_train = pca90_fit.transform(X_train)
X_test = pca90_fit.transform(X_test)

# fit linear regression model
X_train = sm.add_constant(X_train)
X_test = sm.add_constant(X_test)
model = sm.OLS(y_train, X_train).fit()
y_train_pred = model.predict(X_train)
y_test_pred = model.predict(X_test)
y_train_r2 = r2_score(y_train, y_train_pred)
y_test_r2 = r2_score(y_test, y_test_pred)
y_train_r2, y_test_r2
```




    (0.29874232067403128, -0.019842587361204478)


