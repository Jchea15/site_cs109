---
title: The Models
nav_include: 5
---

The ultimate goal of this project was to create a model that is more predictive of Alzheimer's dementia than current methods. Using the most correlated genes, a model was constructed that outperforms both random chance and a model including all genes and models.

----------


Baseline model for comparison
-------------

The baseline model is a simple logistic regression using the 10 highest and 10 lowest correlated genes with fianl diagnosis on the training dataset. This model managed to achieve a classifciation accuracy of 0.360 on the testing dataset, and 0.648 on the training dataset. However, this model is very basic, as it is a linear model, and uses only 20 genes, correlated to a single response variable. This means that it may not include all genes that could help us predict diagnosis, and also it is likely that some of these genes are randomly correlated.


Random and no-work models
-------------
We made 4 models to benchmark how our models are doing: predicting 0, 1, 2 randomly, all 0s, all 1s, and all 2s. The random model varied, but would do around 0.33 in accuracy, usually fluctuating within a 0.04 margin. For the other models, we got the following accuracies:

| Model | Dataset | Accuracy |
| ------------- | ------------- | ------------- |
| Random class. | train | 0.305882352941 |
| Random class. | test | 0.373684210526 |
| All zeros class. | train | 0.411764705882 |
| All zeros class. | test | 0.415789473684 |
| All ones class. | train | 0.4 |
| All ones class. | test | 0.421052631579 |
| All twos class. | train | 0.188235294118 |
| All twos class. | test | 0.163157894737 |


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
