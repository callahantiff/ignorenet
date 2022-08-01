### Background
- Preeclampsia (PE) is life-threatening, acute-onset hypertension and proteinuria at > 20 weeks gestation.<sup>1</sup> PE accounts for 40% of fetal mortality.<sup>2</sup> The only known cure is placenta delivery.  
- Currently, more than one-third of all protein-coding genes have no known function or published literature (i.e., the ignorome).<sup>3,4</sup> Many disease-associated genes may provide important insight when examined within the context of of other diseases/phenotypes.  
- An extensive body of scientific literature and data exist for PE. The PheKnowLator Ecosystem<sup>5</sup> helps users construct large-scale knowledge graphs (KGs) from a wide variety of biomedical data.

**Objective:** Can PheKnowLator be used to identify novel and actionable molecular mechanisms from the PE ignorome?

<br>  

___

### Methods  
The experimental design is highlighted in `Figure 1`.

**Identification of the PE Molecular Signature**
- A meta-analysis of domain expert-selected [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) studies. 

**Identification of Known PE-Associated Genes**  
- **Literature-Driven.** Mine [PubTator](https://www.ncbi.nlm.nih.gov/research/pubtator/), [DisGeNET](https://www.disgenet.org/), and [Malacards](https://www.malacards.org/).  
- **Gene-Driven.** Differentially expressed genes (DEGs) queried against PubAnntotation.

The PE ignorome was identified as genes from the PE molecular signature with no known PE-association in the literature. 

**PheKnowLator KG Enrichment**
- Generate KG node embeddings using [Walking RDF/OWL](https://github.com/bio-ontology-research-group/walking-rdf-and-owl).  
- The 100 nearest KG concepts (gold circles, `Figure 1`) to each ignorome gene were identified, reviewed by domain experts, and compared to gene set enrichment results produced by [ToppGene](https://toppgene.cchmc.org/help/publications.jsp).

<img src="https://user-images.githubusercontent.com/8030363/177888222-0c5c5113-08d8-4603-86f9-49c8a29c61e2.png" width="800" height="550"/>

**Figure 1.** Overview of Results for Finding the Preeclampsia Ignorome. The figure provides an overview of the procedures utilized in order to obtain the preeclampsia ignorome. Acronyms - PE: Preeclampsia.

<br>

___

### Results  
- The PE ignorome contains 445 genes (Venn diagram, `Figure 1`). ToppGene Enrichment revealed that 90% of the PE ignorome genes were associated with a disease other than PE (`Figure 2`), most often neoplasms (48.7%).  
- PheKnowLator-derived enrichment of the 100 KG concepts nearest to each PE ignorome gene resulted in 2,227 unique annotations.  
- Expert reviewed reduced the 2,227 PheKnowLator annotations to 53 deemed worthy of experimental follow-up. None of the identified diseases, biological processes, cellular components, molecular functions, pathways, or phenotype associations overlapped with the ToppGene.  
- Mechanistic explanations were derived for the 53 expert-selected annotations. An example of a novel disease association and mechanistic explanation is shown below:  
<img width="500" src="https://user-images.githubusercontent.com/8030363/177889687-c32c8e88-e12d-4453-abe3-ebc1214b2fbd.png">

<img src="https://user-images.githubusercontent.com/8030363/177890243-50a40fe7-93a9-49f2-b7c2-a1ef57eaa32e.png" width="1200" height="350"/>

**Figure 2.** Preeclampsia Ignorome Gene Annotations in Other Diseases. (A) illustrates the literature coverage of the 445 preeclampsia ignorome genes to other diseases. The x-axis represents the number of disease-annotated articles for each gene as of December 2017. The left y-axis shows the number of genes as bars, where the red bar contains the number of genes with no literature annotations to any disease. The right y-axis shows the number of diseases annotated to each preeclampsia gene and the number of annotations to diseases other than preeclampsia that were found for each ignorome gene in the literature. (B) Plots the counts of literature annotations to high-level disease categories.

<br>

___

### Conclusions
- Expert-led multiplatform microarray meta-analysis and literature mining identified the PE ignorome (n=445). The majority of the ignorome genes were associated with a disease other than PE.  
- The KG-based enrichment strategy produced 53 highly relevant novel PE associations, thus potentially identifying additional targets for prevention/intervention.  
- The PheKnowLator Ecosystem can aid researchers and bench scientists in relevant and biologically-actionable discovery and provide new opportunities to leverage existing resources.

**Limitations:** Limited to transcriptionally-regulated molecules and should be considered with respect to the current lack of agreed upon standards for microarray meta-analysis.

<br><br>

____
____

### Data, Code, Gists, and Notebooks
- [Transcriptomic Metanalysis Notebook (Jupyter Notebook)](https://github.com/callahantiff/ignorenet/blob/master/notebooks/Biohackathon_GEO_Full_pipeline.nb.html)   
- [Knowledge Graph Node Embeddings](https://zenodo.org/record/3830982/files/embeddings.zip)  
- [Enrichment Analysis Script (GitHub Gist)](https://gist.github.com/callahantiff/86d174f27838b5a6d243859fdd3b8e1b)  

**PheKnowLator Knowledge Graph (`v1.0`):** [https://github.com/callahantiff/PheKnowLator/wiki/v1.0.0](https://github.com/callahantiff/PheKnowLator/wiki/v1.0.0)  


<br>

___

#### References
1. Chaiworapongsa T, Romero R, Whitten AE, et al. The use of angiogenic biomarkers in maternal blood to identify which SGA fetuses will require a preterm delivery and mothers who will develop pre-eclampsia. J Matern Fetal Neonatal Med. 2016;29(8):1214-28. [`PMID: 26303962`](https://pubmed.ncbi.nlm.nih.gov/26303962/)  
2. Anderson UD, Olsson MG, Kristensen KH, et al. Biochemical markers to predict preeclampsia. Placenta. 2012;33:S42-7. [`PMID: 22197626`](https://pubmed.ncbi.nlm.nih.gov/22197626/)   
3. Pandey AK, Lu L, Wang X, et al. Functionally enigmatic genes: a case study of the brain ignorome. PLoS One. 2014;9(2):e88889. [`PMID: 24523945`](https://pubmed.ncbi.nlm.nih.gov/24523945/)   
4. Riba M, Garcia Manteiga JM, Bošnjak B, et al. Revealing the acute asthma ignorome: characterization and validation of uninvestigated gene networks. Sci Rep. 2016;6:24647. [`PMID: 27097888`](https://pubmed.ncbi.nlm.nih.gov/27097888/) 
5. Callahan TJ, Tripodi IJ, Wyrwa JM, et al. Phenotype Knowledge Translator: A FAIR Ecosystem for Representing Large-Scale Biomedical Knowledge. 2022. [`Zenodo:5716383`](https://doi.org/10.5281/zenodo.5716383) 


____

📢 Preprint can be found here 👉 [https://arxiv.org/abs/2207.14294](https://arxiv.org/abs/2207.14294)

```
A preprint of this article submitted for consideration in Pacific Symposium on Biocomputing © 2022 copyright World Scientific Publishing Company https://psb.stanford.edu/

```
