---
title: The Cleaning
nav_include: 4
---

blah blah blah.

----------


Preliminary data cleaning
-------------

just got ADNI2, removed longitudinal data, got final_dx, etc...basically the first email we sent

In order build our data set, we extracted gene expression data for all ADNI2 patients, as well as their corresponding cognitive tests from ADNIMERGE. ADNI2 was used as it was the largest sample size for any single protocol in the ADNI database. The data was merged using the patient IDs, a value unique to each patient. 

To control for differing numbers of visits, only the cognitive tests from the visit where gene expression data was taken were used. This allowed us to combine each patient's tests and gene expression into a single row respresenting one visit, vastly improving consistency in the collection protocols for the data.

From there, we recorded the dates of the patient's first and last visits, the days in between the visits, the final diagnosis (of Alzheimer's dimentia), and the diagnosis on the day of the gene expression collection (or most recent diagnosis, since a null value was present for patients with no updated diagnosis). This allowed us to have a complete picture of the patients vists to the doctor, and the progression of their disease, while minimizing the amount of logitundinal data required.

Finally, standard data processing was performed. Types were corrected, features were standardized and diagnoses were converted to numerical representations.

Featuring: too many features!
-------------
blah blah brief intro plz

| Feature  | what the heck it mean | purpose |
| ------------- | ------------- | ------------- |
| FDG  | Average 18F-fluorodeoxyglucose position emission tomography (PET)  | imaging for diagnosis |
| PIB | Pittsburg compound B (PET radio traer for beta amyloid plaque) standard uptake value ratio | imaging for diagnosis |
| AV45 | 18F-AV-45 (florbetapir) (PET radio tracer for beta amyloid plaques) standard uptake value ratio | imaging for diagnosis |
| CDRSB | Clinical Dementia Rating Scale: Sum of Boxes, a measurement of dementia | cognitive tests for dementia |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |


Missingness and mindfulness
-------------

Given the small sample size, it is critical that missing cognitive tests be imputed. However, for 8 patients, all cognitive tests were missing. This would make imputing the value of cognitive tests completely random, and thus useless. All of these patients were dropped. Following this, patients with one or two cognitive tests
