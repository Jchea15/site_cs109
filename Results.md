---
title: The Results
nav_include: 6
---

After experimenting with many different models and combinations of models and dimensionality reduction techniques, the model with the best performance was a gradient boosted Decision Tree using PCA.

----------


A gene-ius look at most correlated genes
-------------
Leah’s spreadsheet
t-test/p-value notebook
alpha=0.01 (116 genes)
Table to categorize all genes
More detailed description of top 10 genes (like Milestone 2)

| Gene Name | Average P value            | Connection to Alzheimer's | 
|-----------|----------------------------|---------------------------| 
| CLIC1     | 1.9E-05, 3.67E-05, 2.1E-04 | "Amyloid β-peptide (Aβ) accumulation in plaques is a hallmark of familial Alzheimer disease. The truncated Aβ25-35 species was shown previously to increase the expression of CLIC1 chloride conductance in cortical microglia and to provoke microglial neurotoxicity." <sup>[1](#myfootnote1)</sup> (Skaper, S D, et al. 2013) | 
| HLA-DQB1  | 3.60E-04                   | MHCII proteins, many MHCII proteins associated with Alzheimer's: "Marked increases in MHC class II-expressing microglia have been shown in many neuropathologic disorders, including Alzheimer's disease (AD)." (Perlmutter et all 1992)                                                          | 
| CD177     | 7.55E-04                   | "In contrast, the CD177+Êpopulation was significantly increased in mAD (= mild stage Alzheimer's disease) patients (p < 0.05)." (Page et al. 2015)                                                                                                                                                | 
| MCEMP1    | 8.2E-04, 2.0E-03           | Associated with strokes, (Wood, 2016)                                                                                                                                                                                                                                                             | 
| LRRFIP1   | 9.36E-04                   | The DNA sensorLRRFIP1Êmediates the production of type I IFN via a _-catenin-dependent pathway (Yang et al, 2010)                                                                                                                                                                                  | 
| PEX5      | 1.26E-03                   | The molecular biological analysis showed that the changes of these TF activities and their target genes in the interactions of signaling proteins in cell cycle, chronic inflammation and immune response play important roles in the deterioration of AD. (Kong et al., 2017)                    | 
| P2RY10    | 1.67E-03                   | Expressed in lymphoid cells, immunological, involved in Calcium influx ativation of diacylglyceride-dependent protein kinases. AlzheimerÕs disease IS a dementia characterised by aberrant calcium signalling (Adrian et al., 2000; Ghosh et al., 2015)                                           | 
| CRAMP1L   | 1.85E-03                   | No information except RNA& protein expressed highly in male tissues, and protein expressed highly in brain according to human protein atlas                                                                                                                                                       | 
| NCAPD2    | 1.91E-03                   | Regulates chromatin during the cell cycle. Gene polymorphisms are associated with Alzheimer's disease (Zhang et al. 2014)                                                                                                                                                                         | 
| NTHL1     | 2.03E-03                   | DNA damage repair, specifically a pathway that is known to be involved in AD (Ray et al., 2007)                                                                                                                                                                                                   | 
| ANXA3     | 2.09E-03                   | Leads to white matter inflammation, as is the case in Alzheimer's (Raj et al., 2017) Also associated with Huntington's disease a different neurological disease.                                                                                                                                  | 
| TFAP4     | 2.26E-03                   | Transcribes AP4 Complex subunits, which mediate sorting of APP (Alzheimer's disease amyloid precursor protein)                                                                                                                                                                                    | 
| VWA9      | 2.45E-03                   | Probably involved in transcription of snRNAs U1, and U2, and "aggregates of U1 snRNA and U1 small nuclear ribonucleoproteins represent a new pathological hallmark of AD" (Hales et al., 2014)                                                                                                    | 
| NFATC3    | 2.46E-03                   | Inhibition of NFAT pathway found to alleviate amyloid beta neurotoxicity in a mouse model of Alzheimer's disease                                                                                                                                                                                  | 
| WFDC1     | 2.57E-03                   | Previously found negatively correlated with AD (Miller et al., 2013)                                                                                                                                                                                                                              | 
| PPP2R5A   | 2.71E-03                   | Can modulate PP2A  catalytic activity, and "alterations in PP2A regulators and PP2A catalytic activity, subunit expression, methylation and/or phosphorylation, have been reported in AD-affected brain regions" (Sontag et al., 2014)                                                            | 
| TMEM241   | 2.72E-03                   | a ubiquitous sugar transporter protein, gene variant suggested to contribute to increased triglyceride levels in Mexicans. Observed relation between high triglycerides, diabetes, and vascular dementia (Raffaitin et al., 2009; Rodr’guez et al., 2016)                                         | 



Diagnosing and detecting Alzheimer's
-------------
Summarize model results
To view our process, please see [our final model building documentation](Finalmodel_notebook.md).


Strengths and shortcomings
-------------
take a wild guess bitch
For more information, please see [our attempt at working with longitudinal data](Long_notebook.md).


Where do we go from here?
-------------
future directions
Down. Straight down to hell.

Given the sample size was so small, this model can easily be improved with a larger sample size. Future research should seek to replicate these methods with the ADNI1/GO/3 databases. Moreover, it would be useful to attempt to create a new database with more standardized collection protocols to allow for easier analysis in future research. Finally, the failure of our attempt at a longitundinal model to predict the progression of Alzheimer's dementia hihglights a sorely lacking aspect of the databases: longitudinal data. Progression of the disease could very well be modeled by gene expression as well, but it is impossible to tell given the limited availible data.


Footnones
-------------
<a name="myfootnote1">1</a>: Skaper, S D, et al. “Intracellular Ion Channel CLIC1: Involvement in Microglia-Mediated β-Amyloid Peptide(1-42) Neurotoxicity.” Neurochemical Research., U.S. National Library of Medicine, Sept. 2013, www.ncbi.nlm.nih.gov/pubmed/23743620.
