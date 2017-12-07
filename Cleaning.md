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
| FDG, FDG_bl | average 18F-fluorodeoxyglucose position emission tomography (PET); 'bl' stands for baseline, meaning the value from the very first visit | imaging for beta amyloid plaques |
| PIB, PIB_bl | Pittsburg compound B (PET radio tracer) standard uptake value ratio | imaging for beta amyloid plaques |
| AV45, AV45_bl | 18F-AV-45 (florbetapir) (PET radio tracer) standard uptake value ratio | imaging for beta amyloid plaques |
| CDRSB | Clinical Dementia Rating Scale: Sum of Boxes, a measurement of dementia | cognitive tests for dementia |
| ADAS11, ADAS13, ADAS11_bl, ADAS13_bl | Alzheimer's Disease Assessment Scale, 11 or 13 items questionnaire; 'bl' stands for baseline | cognitive tests for dementia |
| MMSE, MMSE_bl | mini mental state examination, tests memory, attention, and language | cognitive tests for dementia |
| RAVLT_immediate, RAVLT_learning, RAVLT_forgetting, RALVT_perc_forgetting, all aforementioned \_bl | Rey Auditory Verbal Learning Test, examining short-term and longer-term verbal memory | cognitive tests for dementia |
| FAQ, FAQ_bl | function activities questionnaires, which tests daily activities such as remembering appointments | cognitive tests for dementia |
| MOCA, MOCA_bl | Montreal Cognitive Assessment | cognitive tests for dementia |
| EcogPtMem, EcogPtLang, EcogPtVisspat, EcogPtPlan, EcogPtOrgan, EcogPtDivatt, EcogPtTotal, all aforementioned \_bl | Every Cognition tests as reported by the patient, which respectively test abilities in memory, language, visual-spatial, planning, organising, dividing attention | cognitive tests for dementia |
| EcogSPMem, EcogSPLang, EcogSPVisspat, EcogSPPlan, EcogSPOrgan, EcogSPDivatt, EcogSPTotal, all aforementioned \_bl | Every Cognition, as reported by a study partner rather than the patient themself | cognitive tests for dementia |
| FLDSTRENG | the field strength of the MRI used, either 1.5T or 3T | MRI-related |
| FSVERSION | another measure of field strength; FreeSurfer version | MRI-related |
| Ventricles, Hippocampus, WholeBrain, Entorhinal, Fusiform, MidTemp, ICV, all aforementioned \_bl | a metric of volume and cortical thickness for these regions of the brain | MRI-related |
| DX | diagnosis (0 = cognitively normal, 1 = mild cognitive impaired, 2 = dementia) | the official diagnosis |


Missingness and mindfulness
-------------

Given the small sample size, it is critical that missing cognitive tests be imputed. However, for 8 patients, all cognitive tests were missing. This would make imputing the value of cognitive tests completely random, and thus useless. All of these patients were dropped. Following this, patients with one or two cognitive tests
