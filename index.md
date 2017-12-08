---
title: Alzheimer's through Gene Expression
---

Group number: 69

Group members: Jackie Chea, Nike Izmaylov, Nick Pagel, Leah Rosen

## Project Statement

>*"Did we just cure Alzheimer's?" -- our group, right before we realized we had not.*

For some patients, gene expression tests, particularly for a small number of targeted genes, may be relatively inexpensive and less invasive compared to traditional methods such as radiolabeled PET scans. This project seeks to use targetted gene expression data to help reduce the costs of diagnosing Alzheimer's dementia and aid in earlier detection. Given that Alzheimer's disease is the third leading cause of death in the elderly, early detection is critical to limiting its rate of progression. 

We  first used the ADNIMERGE database, as well as gene expression data from ADNI2 in order to isolate genes with a strong correlation to Alzheimer's dementia. Using these genes, we constructed a model that aims to diagnose Alzheimer's dementia based soley upon gene expression data, a far cheaper method of diagnosis. Furhtermore, the model was altered to be able to abstain from a diagnosis in certain cases in order to minimize the overall cost of diagnosis. This allows the model to not only minimize costs in testing, but also in treatment for those individuals who are wrongly diagnosed. 

In the end, the model managed to diagnose 46% of patients accurately; such a custom gene expression test for slightly more than a hundred genes would cost about a hundred dollars. Counting the possible problems of abstaining and misdiagnosis, the model, once thresholded, has an average cost of $6801.47 per patient. This is below the normal average cost per patient of $10.0042.11 for the next-best model. With further research, it is not unrealistic to eventually expect a model based upon gene expression data to be much more accurate, significantly reducing the difficulty in predicting a disease that currently only be diagnosed post-humously.

Our original project goal was "Using the ADNI database, we will analyze the genetic data of patients with and without Alzheimer’s, seeking genes/biomarkers that are correlated with AD, aside from TOMM40 and APOE, in the hopes of finding genetic markers and tests for early detection of Alzheimer’s, and potentially even risk prediction prior to onset." This is pretty much what our goal still is. However, as we worked with the data we realised that the gene sequencing data was a lot more unwieldy than the gene expression data, while both allowed for the same level of analysis. Therefore we decided to go for the dataset where we would be able to do more analysis in the given scope. We decided to not only look at diagnosis as correlated to gene expression data, but knew from the start that we wanted to look at several ADNI tests. As we became more familiar with the data, we decided that we would use all of the reliable cognitive tests, as well as APOE4 allele to decide relevant genes. We weren't sure how to integrate multiple variables into a single response variable, and in the end decided to only use diagnosis as the final response variable, but select which of the >45000 genes to use by making sure that they were correlated with every reliable cognitive test, as well as diagnosis, and APOE4 allele. We reasoned that this would minimise random correlation and make as sure as possible that the selected genes were actually causally connected with cognitive impairment.
