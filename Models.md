---
title: The Models
nav_include: 5
---

The ultimate goal of this project was to create a model that is more predictive of Alzheimer's dementia than current methods. Using the most correlated genes, a model was constructed that outperforms both random chance and a model including all genes and models.

----------


Back to the basics: a baseline model for comparison
-------------

The baseline model is a simple logistic regression using the 10 highest and 10 lowest correlated genes with final diagnosis on the training dataset. This model managed to achieve a classification accuracy of 0.360 on the testing dataset, and 0.648 on the training dataset. However, this model is very basic, as it is a linear model, and uses only 20 genes, correlated to a single response variable. This means that it may not include all genes that could help us predict diagnosis, and also it is likely that some of these genes are randomly correlated.

To view our process, please see [our baseline model immediately after preliminary EDA](EDA_notebook.md).

On your benchmarks!
-------------
We made four models to benchmark how our models are doing: (1) predicting 0, 1, 2 randomly; (2) all 0s; (3) all 1s; and (4) all 2s. We observed the following accuracies, using a random seed of 9001 for the random class:

| Model     | Train Accuracy | Test Accuracy |
| --------- | -------------- | ------------- |
| Random    | 0.306          | 0.374         |
| All zeros | 0.412          | 0.416         |
| All ones  | 0.400          | 0.421         |
| All twos  | 0.188          | 0.163         |


Ready, set, go: beyond the baseline
-------------

| Model                                    | Train Accuracy | Test Accuracy |
| ---------------------------------------- | -------------- | ------------- |
| Decision Tree with CV for max depth      | 1.000          | 0.363         |
| Logistic One-Vs-Rest                     | 0.635          | 0.468         |
| Logistic Multinomial                     | 0.771          | 0.347         |
| LDA                                      | 0.612          | 0.474         |
| QDA                                      | 0.888          | 0.426         |
| kNN                                      | 0.724          | 0.421         |
| AdaBoost with PCA and CV for number of estimators | 0.665          | 0.421         |


And the award goes to...
------------------
Our best model was a gradient boost algorithm with PCA and CV for number of estimators, with a test accuracy of 0.489.


Abstinence is the best policy
-------------

While the model performed well, it can be improved. The model predicts based upon the probabilities that a given patient is of a certain diagnosis. This means that a certain patient can be misdiagnosed if they are even just above the probability threshold for a given diagnosis. As the model has a decent variance, it is likely that a small number of patients with diagnosis proabilities' near the threshold will be misdiagnosed. When it comes to medical misdaignoses, this can lead to significant uncessary expenses, and should be minimized. 

In order to do just this, we implemented a version of the model that is able to abstain in order to minimize the cost of a diagnosis. A misdiagnosis of positive results in follow-up tests and visits to specialists, carrying a pricetag of \$12000, while false negative would result in later testing, a cost of \$6850 (cost of traditional tests minus gene expression data). The model tests models with varying thresholds and then calculates the costs, outputting the optimum thresholds to minimize the costs. 

This model had a **insert bullshit I don't wanna write here**

To view our process, please see [our documented creation and tuning of the final model](Finalmodel_notebook.md). You may also be interested in [our alternative models that were discarded](Modelgraveyard_notebook.md).
