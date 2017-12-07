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
| HLA-DQB1  | 3.60E-04                   | This is a type of MHCII protein, MHCII proteins have previously been associated with Alzheimer's: "Marked increases in MHC class II-expressing microglia have been shown in many neuropathologic disorders, including Alzheimer's disease (AD)." <sup>[2](#myfootnote2)</sup> (Perlmutter et all 1992) | 
| CD177     | 7.55E-04                   | Study found that CD177+ neutrophil population significantly increased in mAD (= mild stage Alzheimer's disease) patients (p < 0.05)." <sup>[3](#myfootnote3)</sup> (Page et al. 2015) | 
| MCEMP1    | 8.2E-04, 2.0E-03           | Associated with strokes so this may be a correlation unrelated to Alzheimer's, <sup>[4](#myfootnote4)</sup> (Wood, 2016) | 
| LRRFIP1   | 9.36E-04                   | The DNA sensor LRRFIP1 mediates the production of type I IFN via a β-catenin-dependent pathway, which mediates neuro-inflammatory events in models of Alzheimer's disease. <sup>[5](#myfootnote5)</sup> (Yang et al, 2010; Taylor et al, 2014) | 
| PEX5      | 1.26E-03                   | "The molecular biological analysis showed that the changes of these transcription factor activities and their target genes in the interactions of signaling proteins in cell cycle, chronic inflammation and immune response play important roles in the deterioration of Alzheimer's Disease." (on group of genes including PEX5) <sup>[6](#myfootnote6)</sup> (Kong et al., 2017)                    | 
| P2RY10    | 1.67E-03                   | Expressed in lymphoid cells, immunological, involved in Calcium influx ativation of diacylglyceride-dependent protein kinases. "Alzheimer's disease is a dementia characterised by aberrant calcium signalling." <sup>[7](#myfootnote7)</sup> (Adrian et al., 2000; Ghosh et al., 2015)                                           | 
| CRAMP1L   | 1.85E-03                   | No information except RNA & protein expressed highly in male tissues, and protein expressed highly in brain according to human protein atlas. Can't infer connection to AD. | 
| NCAPD2    | 1.91E-03                   | Regulates chromatin during the cell cycle. Gene polymorphisms in NCAPD2 are associated with Alzheimer's disease <sup>[8](#myfootnote8)</sup> (Zhang et al. 2014) | 
| NTHL1     | 2.03E-03                   | DNA damage repair, specifically involved in a DNA repair pathway that is known to be involved in AD <sup>[9](#myfootnote9)</sup> (Ray et al., 2007) | 
| ANXA3     | 2.09E-03                   | Leads to white matter inflammation, as is the case in Alzheimer's <sup>[10](#myfootnote10)</sup> (Raj et al., 2017) Also associated with Huntington's disease, a different neurological disease. | 
| TFAP4     | 2.26E-03                   | Transcribes AP4 Complex subunits, which mediate sorting of APP (Alzheimer's disease amyloid precursor protein) <sup>[11](#myfootnote11)</sup> (Burgos et al., 2010 | 
| VWA9      | 2.45E-03                   | Probably involved in transcription of snRNAs U1, and U2, and "aggregates of U1 snRNA and U1 small nuclear ribonucleoproteins represent a new pathological hallmark of AD" <sup>[12](#myfootnote12)</sup> (Hales et al., 2014) | 
| NFATC3    | 2.46E-03                   | "Inhibition of NFAT pathway found to alleviate amyloid beta neurotoxicity in a mouse model of Alzheimer's disease." <sup>[13](#myfootnote13)</sup> (Hudry et al., 2012) | 
| WFDC1     | 2.57E-03                   | Previously found negatively correlated with AD <sup>[14](#myfootnote14)</sup> (Miller et al., 2013) | 
| PPP2R5A   | 2.71E-03                   | Can modulate PP2A  catalytic activity, and "alterations in PP2A regulators and PP2A catalytic activity, subunit expression, methylation and/or phosphorylation, have been reported in AD-affected brain regions" <sup>[15](#myfootnote15)</sup> (Sontag et al., 2014) | 
| TMEM241   | 2.72E-03                   | a ubiquitous sugar transporter protein, gene variant suggested to contribute to increased triglyceride levels in Mexicans. Observed relation between high triglycerides, diabetes, and vascular dementia <sup>[16](#myfootnote16)</sup> (Raffaitin et al., 2009; Rodr’guez et al., 2016) | 



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

<a name="myfootnote2">2</a>: Perlmutter, L S, et al. “MHC Class II-Positive Microglia in Human Brain: Association with Alzheimer Lesions.”Journal of Neuroscience Research., U.S. National Library of Medicine, Dec. 1992, www.ncbi.nlm.nih.gov/pubmed/1484388.

<a name="myfootnote3">3</a>: Page, Aurélie LE, et al. “Immune Signatures of Alzheimer's Disease: Profiles of Neutrophils. (HUM1P.301).” The Journal of Immunology, American Association of Immunologists, 1 May 2015, www.jimmunol.org/content/194/1_Supplement/52.26.

<a name="myfootnote4">4</a>: Wood, Heather. “Stroke: MCEMP1 — a New Prognostic and Diagnostic Biomarker for Stroke?” Nature News, Nature Publishing Group, 19 Feb. 2016, www.nature.com/articles/nrneurol.2016.17.

<a name="myfootnote5">5</a>: Yang, P, et al. “The Cytosolic Nucleic Acid Sensor LRRFIP1 Mediates the Production of Type I Interferon via a Beta-Catenin-Dependent Pathway.” Nature Immunology., U.S. National Library of Medicine, June 2010, www.ncbi.nlm.nih.gov/pubmed/20453844. Taylor, J M, et al. “Type-1 Interferon Signaling Mediates Neuro-Inflammatory Events in Models of Alzheimer's Disease.” Neurobiology of Aging., U.S. National Library of Medicine, May 2014, www.ncbi.nlm.nih.gov/pubmed/24262201.

<a name="myfootnote6">6</a>: Kong, Wei, et al. “Differences of Immune Disorders between Alzheimer’s Disease and Breast Cancer Based on Transcriptional Regulation.” PLOS ONE, Public Library of Science, journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0180337.

<a name="myfootnote7">7</a>: Adrian, K., Bernhard, M. K., Breitinger, H.-G., Ogilvie, A. Expression of purinergic receptors (ionotropic P2X1-7 and metabotropic P2Y1-11) during myeloid differentiation of HL60 cells. Biochim. Biophys. Res. Acta 1492: 127-138, 2000. Ghosh, Anshua, and Karl Peter Giese. “Calcium/Calmodulin-Dependent Kinase II and Alzheimer's Disease.” Molecular Brain, BioMed Central, 24 Nov. 2015, molecularbrain.biomedcentral.com/articles/10.1186/s13041-015-0166-2.

<a name="myfootnote8">8</a>: Zhang, P, et al. “Non-SMC Condensin I Complex, Subunit D2 Gene Polymorphisms Are Associated with Parkinson's Disease: a Han Chinese Study.” Genome., U.S. National Library of Medicine, May 2014, www.ncbi.nlm.nih.gov/pubmed/25166511.

<a name="myfootnote9">9</a>: Ray, Monika and Zhang, Weixiong, "DNA repair in incipient Alzheimer's disease" Report Number: WUCSE-2007-31 (2007). All Computer Science and Engineering Research.

<a name="myfootnote10">10</a>: Raj, Divya, et al. “Increased White Matter Inflammation in Aging- and Alzheimer’s Disease Brain.” Frontiers in Molecular Neuroscience, Frontiers Media S.A., 2017, www.ncbi.nlm.nih.gov/pmc/articles/PMC5492660/.

<a name="myfootnote11">11</a>: Burgos, P V, et al. “Sorting of the Alzheimer's Disease Amyloid Precursor Protein Mediated by the AP-4 Complex.” Developmental Cell., U.S. National Library of Medicine, 16 Mar. 2010, www.ncbi.nlm.nih.gov/pubmed/20230749?dopt=Abstract.

<a name="myfootnote12">12</a>: Hales, Chadwick M., et al. “Aggregates of Small Nuclear Ribonucleic Acids (SnRNAs) in Alzheimer’s Disease.’” Brain Pathology (Zurich, Switzerland), U.S. National Library of Medicine, July 2014, www.ncbi.nlm.nih.gov/pmc/articles/PMC4096308/.

<a name="myfootnote13">13</a>: Hudry, Eloise, et al. “INHIBITION OF THE NFAT PATHWAY ALLEVIATES AMYLOID BETA NEUROTOXICITY IN A MOUSE MODEL OF ALZHEIMER’S DISEASE.” The Journal of Neuroscience, U.S. National Library of Medicine, 29 Feb. 2012, www.ncbi.nlm.nih.gov/pmc/articles/PMC3296329/.

<a name="myfootnote14">14</a>: Miller, Jeremy A, et al. “Genes and Pathways Underlying Regional and Cell Type Changes in Alzheimer's Disease.” Genome Medicine, BioMed Central, 25 May 2013, genomemedicine.biomedcentral.com/articles/10.1186/gm452.

<a name="myfootnote15">15</a>: Sontag, Jean-Marie, and Estelle Sontag. “Protein Phosphatase 2A Dysfunction in Alzheimer’s Disease.” Frontiers in Molecular Neuroscience, Frontiers Media S.A., 2014, www.ncbi.nlm.nih.gov/pmc/articles/PMC3949405/.

<a name="myfootnote16">16</a>: Raffaitin, Christelle, et al. “Metabolic Syndrome and Risk for Incident Alzheimer's Disease or Vascular Dementia: The Three-City Study.” Diabetes Care, American Diabetes Association, Jan. 2009, www.ncbi.nlm.nih.gov/pmc/articles/PMC2606808/. Rodríguez, A, et al. “Molecular Characterization of the Lipid Genome-Wide Association Study Signal on Chromosome 18q11.2 Implicates HNF4A-Mediated Regulation of the TMEM241 Gene.” Arteriosclerosis, Thrombosis, and Vascular Biology., U.S. National Library of Medicine, July 2016, www.ncbi.nlm.nih.gov/pubmed/27199446.
