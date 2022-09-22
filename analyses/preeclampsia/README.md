# Preeclampsia Ignorome

## Background  
Preeclampsia is a leading cause of maternal and fetal morbidity and mortality. Currently, the only definitive treatment of preeclampsia is delivery of the placenta, which is central to the pathogenesis of the disease. Transcriptional profiling of human placenta from pregnancies complicated by preeclampsia has been extensively performed to identify differentially expressed genes (DEGs). DEGs are identified using unbiased assays, however the decisions to investigate DEGs experimentally are biased by many factors, causing many DEGs to remain uninvestigated. A set of DEGs which are associated with a disease experimentally, but which have no known association to the disease in the literature are known as the ignorome. Preeclampsia has an extensive body of scientific literature, a large pool of DEG data, and only one definitive treatment. Tools facilitating knowledge-based analyses, which are capable of combining disparate data from many sources in order to suggest underlying mechanisms of action, may be a valuable resource to support discovery and improve our understanding of this disease.

#### Objective  
Demonstrate how a large-scale heterogeneous biomedical knowledge graph (KG) could be used to identify novel preeclampsia mechanisms from previously analyzed transcriptomic experiments.  

___

## Analysis Workflow  
All of the data files references in the [analysis_script.py](https://github.com/callahantiff/ignorenet/blob/master/analyses/preeclampsia/analysis_script.py) can be downloaded from Zenodo (https://zenodo.org/record/7105487/files/analysis_data.zip). 

### Step 1 - Identify the Preeclampsia Ignorome
#### Identify the preeclampsia molecular signature     
A meta-analysis was performed to identify relevant transcriptomic data on the Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)). Using the keyword `“preeclampsia”`, publicly available human experiments deposited in GEO were examined. The initial set of identified studies were further reviewed for the following criteria to ensure:  
1. Processed samples were from a human placenta biopsy (i.e., chorionic villi, decidua basalis, and placenta);  
2. Samples were processed using Agilent, Affymetrix, Applied Biosystems, Illumina, or NimbleGen;  
3. Studies provided normalized data and/or DEG lists.

Each study’s normalized data were processed using standard R pipelines using the ignorenet library (https://github.com/callahantiff/ignorenet). To make our analysis transparent and reproducible, we created a R Markdown notebook that includes the code that was run and describes each step in detail. The R Notebook can be accessed here: https://github.com/callahantiff/ignorenet/blob/master/notebooks/Biohackathon_GEO_Full_pipeline.nb.html.  

The final gene list was assembled by selecting significant DEGs (p<0.05) in at least 50% of the studies and included (listed by GEO study identifier):  
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

The processed data for each study can be found here: https://github.com/callahantiff/ignorenet/tree/master/R/DE_Results/Final_studies

#### Identify the genes with a known association with preeclampsia.  
To identify known preeclampsia genes two strategies were employed:  
**Literature-Driven.** Identify relevant genes via keyword search against PubTator, DisGeNET, and Malacards. For this step, all queried results were manually verified for accuracy (i.e., verified that hits obtained were actually to preeclampsia and the associated keywords and were not errors or mismatches to closely associated synonyms or acronyms) and all valid associations were used to create a final unique list of genes.  
*Keywords:*  `“Preeclampsia”`, `“HELLP Syndrome”`, `“Severe Preeclampsia”`, and `“Placenta Disease”`

**Gene-Driven.** Implemented during the 10th annual BioHackathon, this strategy aimed to identify relevant articles by querying 18 keywords (similar to the `Literature-Driven Approach`) in addition to the  preeclampsia molecular signature DEGs against PubAnnotation. Similar to the Literature-Driven Approach, all results were manually verified for accuracy and all associations were used to create a final unique list of genes.  

*Keywords:* `"HELLP Syndrome"`, `"gestational hypertension"`, `"hypertensive disorders of pregnancy"`, `"pre-eclampsia"`, `"pre-eclamptic patient"`, `"pre-eclamptic toxaemia"`, `"pre-eclamptic toxemia"`, `"preeclampsia"`, `"preeclamptic toxaemia"`, `"preeclamptic toxemia"`, `"pregnancy hypertension"`, `"pregnancy-related hypertensive disorder"`, `"toxaemia of pregnancy"`, `"toxaemic mother"`, `"toxaemic pregnancy"`, `"toxemia of pregnancy"`, `"toxemic mother"`, `"toxemic pregnancy”`

The final list of known preeclampsia genes was derived as the unique set of genes produced from both strategies.

#### Identify the Preeclampsia Ignorome  
The preeclampsia ignorome was generated from the set difference of the gene lists generated from the preeclampsia molecular signature and the list of genes known to be associated with preeclampsia in the literature. The complete list of these genes can be found here (Entrez identifiers): https://github.com/callahantiff/ignorenet/blob/master/analyses/preeclampsia/preeclampsia_gene_lists.txt  

<br>  

### Step 2 - Annotate the Preeclampsia Ignorome  
#### Knowledge Graph
A v1.0 PheKnowLator KG was used for this project. This KG was built using publicly available Linked Open Data and Open Biological and Biomedical Ontology Foundry ontologies. The core set of ontologies included phenotypes (Human Phenotype Ontology [HP]), diseases (Human Disease Ontology [DOID]), and biological processes, molecular functions, and cellular components (Gene Ontology [GO]). Genes, pathways, and chemicals were added to the core set of ontologies to form the foundation of the KG which was extended by adding relations between phenotypes, diseases, and GO biological processes, molecular functions, and cellular components. As shown in the figure below, the [Basic Formal Ontology](http://basic-formal-ontology.org/) and [Relation Ontology](https://github.com/oborel/obo-relations/) ontologies were then used to create edges between the node types. The downloaded resource information for generating this information can be accessed [here](https://www.dropbox.com/s/lyegt57269mvey3/resource_info.txt?dl=0).

<img src="https://user-images.githubusercontent.com/8030363/63867692-88cb7280-c972-11e9-9963-a437c4b63d9f.png" width="800" height="500">

Details on the knowledge graph can be found from our preprint (https://zenodo.org/record/5893789) as well as on the PheKnowLator GitHub page: https://github.com/callahantiff/PheKnowLator/wiki/v1-Build-Details. The data used are available from Zenodo (https://zenodo.org/record/7030201).  

#### Knowledge Graph Node Embeddings
A [modified](https://github.com/xgfs/deepwalk-c) version of the [DeepWalk algorithm](https://github.com/phanein/deepwalk) was implemented to generate molecular mechanism embeddings from the biomedical knowledge graph. A t-SNE plot of the dimensionality reduced mechanism embeddings is shown in [Figure 2](https://github.com/callahantiff/PheKnowLator/wiki/v1.0.0/figure-2-t-sne-plot-of-molecular-mechanisms) below. For this release, the hyperparameters were set to 512 dimensions, 100 walks, walk length of 20, and a window of 10. Hyperparameter settings for t-SNE plots were determined through experimentation and for node embeddings, we used the recommended settings provided by the algorithm developers.  

Node Embeddings are publicly available on Zenodo: https://zenodo.org/record/7030189


#### Knowledge-based Mechanistic Enrichment   
Using the PheKnowLator node embeddings, the 100 nearest disease, drug, gene, GO concepts, pathway, and phenotype (i.e., domains) annotations as measured by pairwise cosine similarity of the node embeddings (normalized dot product: `K(X, Y)=<X, Y>/ (||X||*||Y||`) were obtained for each ignorome gene. These annotations were reviewed by a PhD molecular biologist specializing in reproductive science (ALS; 08-09/2021).

To determine if these associations occurred by chance, we:  
1. Examined the overlap between the top-100 closest associations to each ignorome gene in the expert-verified list and the associations generated when enriching the preeclampsia ignorome using ToppGene.  
2. Computed how often the reviewed associations occurred by chance in 1,000 ignorome-sized random samples drawn from all non-ignorome genes represented in the KG. For each sample, the top-100 closest annotations to each gene, by domain were obtained and the number of annotations that overlapped with the expert-verified list was recorded. P-values were obtained for each domain by dividing the number of overlapping annotations out of the 1,000 samples, where a p-value of 0.05 indicates there is a 50 in 1,000 chance of observing a sample annotation that overlaps with the expert-verified preeclampsia ignorome annotations.


<br>  

___

## Resource List  
- GEO Metanalysis Notebook (Jupyter Notebook): https://github.com/callahantiff/ignorenet/blob/master/notebooks/Biohackathon_GEO_Full_pipeline.nb.html   
- Analysis Script: https://github.com/callahantiff/ignorenet/blob/master/analyses/preeclampsia/analysis_script.py  
- Zipped Analysis Data: https://zenodo.org/record/7105487/files/analysis_data.zip    

**PheKnowLator Knowledge Graph (`v1.0`):**  
- GitHub: https://github.com/callahantiff/PheKnowLator/wiki/v1-Build-Details  
- Knowledge Graph: https://zenodo.org/record/7030201
- Node Embeddings: https://zenodo.org/record/7030189     

**Ignorome Gene Lists:**  
- Preeclampsia Gene Lists: https://github.com/callahantiff/ignorenet/blob/master/analyses/preeclampsia/preeclampsia_gene_lists.txt  



<br>  

___

## Additional Information  
For more information on the analysis that was performed including a description of the results, please see the published manuscript on arXiv:  
> Callahan TJ, Stefanski AL, Kim JD, Baumgartner Jr WA, Wyrwa JM, Hunter LE. Knowledge-Driven Mechanistic Enrichment of the Preeclampsia Ignorome. 2022. [Preprint]. [arXiv:2207.14294](https://arxiv.org/abs/2207.14294)  

and conference poster on Zenodo:  
> Callahan TJ, Stefanski AL, Kim JD, Baumgartner Jr WA, Wyrwa JM, Hunter LE. Knowledge-Driven Mechanistic Enrichment of the Preeclampsia Ignorome. 2022. [Poster]. [Zenodo:7105487](https://zenodo.org/record/7105487/files/2023_PSB_Poster_Callahan.pdf)  
