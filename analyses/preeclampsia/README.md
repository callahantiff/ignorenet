# Preeclampsia Ignorome

## Background  
Preeclampsia is a leading cause of maternal and fetal morbidity and mortality. Currently, the only definitive treatment of preeclampsia is delivery of the placenta, which is central to the pathogenesis of the disease. Transcriptional profiling of human placenta from pregnancies complicated by preeclampsia has been extensively performed to identify differentially expressed genes (DEGs). DEGs are identified using unbiased assays, however the decisions to investigate DEGs experimentally are biased by many factors, causing many DEGs to remain uninvestigated. A set of DEGs which are associated with a disease experimentally, but which have no known association to the disease in the literature are known as the ignorome. Preeclampsia has an extensive body of scientific literature, a large pool of DEG data, and only one definitive treatment. Tools facilitating knowledge-based analyses, which are capable of combining disparate data from many sources in order to suggest underlying mechanisms of action, may be a valuable resource to support discovery and improve our understanding of this disease.

#### Objective  
Demonstrate how a large-scale heterogeneous biomedical knowledge graph (KG) could be used to identify novel preeclampsia mechanisms from previously analyzed transcriptomic experiments.  

___

## Analysis Workflow  
Our analysis had two analysis steps:  

### Step 1 - Indentify the Preeclampsia Ignorome
#### Identify the preeclampsia molecular signature     
A meta-analysis was performed to identify relevant transcriptomic data on the Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)). Using the keyword `“preeclampsia”`, publicly available human experiments deposited in GEO were examined. The initial set of identified studies were further reviewed for the following criteria to ensure:  
1. Processed samples were from a human placenta biopsy (i.e., chorionic villi, decidua basalis, and placenta);  
2. Samples were processed using Agilent, Affymetrix, Applied Biosystems, Illumina, or NimbleGen;  
3. Studies provided normalized data and/or DEG lists.

Each study’s normalized data were processed using standard R pipelines using the ignorenet library (https://github.com/callahantiff/ignorenet). To make our analysis transparent and reproducbile, we created a R Markdown notebook that includes the code that was run and describes each step in detail. The R Notebook can be accessed here: https://github.com/callahantiff/ignorenet/blob/master/notebooks/Biohackathon_GEO_Full_pipeline.nb.html.  

The final gene list was assembled by selecting significant  DEGs (p<0.05) in at least 50% of the studies and included (listed by GEO study identifier):  
- [`GSE4707`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4707)  
- [`GSE10588`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10588)  
- [`GSE25906`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25906)    
- [`GSE14722`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14722)    
- [`GSE24129`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24129)    
- [`GSE30186`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30186)    
- [`GSE43942`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43942)    
- [`GSE74341`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74341)    
- [`GSE35574`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35574)    
- [`GSE44711`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44711)    
- [`GSE60438`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60438)    
- [`GSE73374`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73374)    


#### Idenitfify the genes with a known association with preeclampsia.  
To identify known preeclampsia genes two strategies were employed:  
**Literature-Driven.** Identify relevant genes via keyword search against PubTator, DisGeNET, and Malacards. For this step, all queried results were manually verified for accuracy (i.e., verified that hits obtained were actually to preeclampsia and the associated keywords and were not errors or mismatches to closely associated synonyms or acronyms) and all valid associations were used to create a final unique list of genes.  
*Keywords:*  `“Preeclampsia”`, `“HELLP Syndrome”`, `“Severe Preeclampsia”`, and `“Placenta Disease”`

**Gene-Driven.** Implemented during the 10th annual BioHackathon, this strategy aimed to identify relevant articles by querying 18 keywords (similar to the `Literature-Driven Approach`) in addition to the  the preeclampsia molecular signature DEGs against PubAnnotation. Similar to the Literature-Driven Approach, all results were manually verified for accuracy and all associations were used to create a final unique list of genes.  

*Keywords:* `"HELLP Syndrome"`, `"gestational hypertension"`, `"hypertensive disorders of pregnancy"`, `"pre-eclampsia"`, `"pre-eclamptic patient"`, `"pre-eclamptic toxaemia"`, `"pre-eclamptic toxemia"`, `"preeclampsia"`, `"preeclamptic toxaemia"`, `"preeclamptic toxemia"`, `"pregnancy hypertension"`, `"pregnancy-related hypertensive disorder"`, `"toxaemia of pregnancy"`, `"toxaemic mother"`, `"toxaemic pregnancy"`, `"toxemia of pregnancy"`, `"toxemic mother"`, `"toxemic pregnancy”`

The final list of known preeclampsia genes was derived as the unique set of genes produced from both strategies.

#### Identify the Preeclampsia Ignorome  
The preeclampsia ignorome was generated from the set difference of the gene lists generated from the preeclampsia molecular signature and the list of genes known to be associated with preeclampsia in the literature. The complete list of these genes can be found here (Entrez identifiers): [Supplemental Table 3).  

<br>  

### Step 2 - Annotate the Preeclampsia Ignorome  
#### Knowledge Graph Node Embeddings
Knowledge Graph Data and Files: https://github.com/callahantiff/PheKnowLator/wiki/v1.0.0

#### Generate Node Embeddings
**Data:**  
- Node Embeddings are publicly available on Zenodo: https://zenodo.org/record/3830982/files/embeddings.zip


#### Knowledge-based Mechanistic Enrichment   
Enrichment Analysis Code: https://gist.github.com/callahantiff/86d174f27838b5a6d243859fdd3b8e1b


<br>  

___

## Resource List  
- [Transcriptomic Metanalysis Notebook (Jupyter Notebook)](https://github.com/callahantiff/ignorenet/blob/master/notebooks/Biohackathon_GEO_Full_pipeline.nb.html)   
- [Knowledge Graph Node Embeddings](https://zenodo.org/record/3830982/files/embeddings.zip)  
- [Enrichment Analysis Script (GitHub Gist)](https://gist.github.com/callahantiff/86d174f27838b5a6d243859fdd3b8e1b)  

**PheKnowLator Knowledge Graph (`v1.0`):** [https://github.com/callahantiff/PheKnowLator/wiki/v1.0.0](https://github.com/callahantiff/PheKnowLator/wiki/v1.0.0)  



<br>  

___

## Additional Information  
For more information on the analysis that was performed including a description of the results, please see the published manuscript on arXiv:  
> Callahan TJ, Stefanski AL, Kim JD, Baumgartner Jr WA, Wyrwa JM, Hunter LE. Knowledge-Driven Mechanistic Enrichment of the Preeclampsia Ignorome. 2022. [Preprint]. [arXiv:2207.14294](https://arxiv.org/abs/2207.14294)  
