---
title: The Data
nav_include: 2
---

blah blah blah.

----------


The ADNI database
-------------

explanation of what the FUCK it is


Gene expression and ADNIMERGE
-------------
what data we using and why, why ADNI2 and not ADNIGO/1, why not WSG


Cognitive testing, and 50,000 genes, and Alzheimer's diagnosis, oh my!
-------------
"how preliminary EDA defined the question"



![dx to final](images/EDA_24_1.png)

A diagnosis of 0 is cognitively normal, 1 is mild cognitive impairment, and 2 is dementia. It appears that most of the patients who were diagnosed as cognitively normal during the visit where the gene expression data was measured stayed cognitively normal until their final visit, while a significant number of patients diagnosed with mild cognitive impairment were later diagnosed with dementia.



![corr_matrix](images/EDA_34_2.png)

By looking at the last row of the correlation matrix, we can see how each cognitive test correlates to the final diagnosis. RAVLT_forgetting and self-reported Ecog are least correlated (study-partner-reported Ecog was more correlated), so we won’t use them as response variables.

![correlations](images/EDA_11_0.png)

The correlation of genes with the final diagnosis is approximately normally distributed, centered at 0 with the most extreme values having magnitudes of approximately 0.15 or higher. By looking at the genes with the most extreme correlations to the final diagnosis, we tried to find genes for which the expression levels could be good predictors, and ideally suggest pathways involved in Alzheimer’s.