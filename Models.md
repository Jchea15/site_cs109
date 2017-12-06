---
title: The Models
nav_include: 5
---

The ultimate goal of this project was to create a model that is more predictive of Alzheimer's dementia than current methods. Using the most correlated genes, a model was constructed that outperforms both random chance and a model including all genes and models.

----------


Baseline model for comparison
-------------

The baseline model is a simple logistic regression using all of the availible data. This model managed to achieve a classifciation accuracy of ______. However, it is computationally complex and, as it fits a model with tens of thousands of predictors to a dataset of only 285 samples, prone to overfitting. For this reason, it was necessary to construct a simpler model, one that can predcit Alzheimer's dementia more accurately.


Beyond the baseline
-------------
other models that we tried


A star performance
-------------
summarize performance for each one


Abstinence is the best policy
-------------

While the model performed well, it can be improved. The model predicts based upon the probabilities that a given patient is of a certain diagnosis. This means that a certain patient can be misdiagnosed if they are even just above the probability threshold for a given diagnosis. As the model has a decent variance, it is likely that a small number of patients with diagnosis proabilities' near the threshold will be misdiagnosed. When it comes to medical misdaignoses, this can lead to significant uncessary expenses, and should be minimized. 

In order to do just this, we implemented a version of the model that is able to abstain in order to minimize the cost of a diagnosis. A misdiagnosis of positive results in follow-up tests and visits to specialists, carrying a pricetag of 12000, while false negative would result in later testing, a cost of $6850 (cost of traditional tests minus gene expression data). The model tests models with varying thresholds and then calculates the costs, outputting the optimum thresholds to minimize the costs. 

This model had a **insert bullshit I don't wanna write here**
