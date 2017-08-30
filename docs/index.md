# ignorenet

### Project Description
Integrating biological knowledge with experimental data is vital for understanding the mechanisms underlying complex diseases and clinical phenotypes. Unfortunately, most genes differentially expressed in association with disease have no currently described disease-associated function in the literature (aka the ignorome) [1-2]. The current project aims to integrate existing linked open biomedical resources and software/tools with publicly available high-throughput data repositories to identify and annotate the function of uncharacterized disease-associated genes for a specific disease. This repository documents our efforts to reveal and annotate disease-specific ignoromes.
<br>
### Approach
The general approach we plan to use (described using diabetes as an example). Please note that this workflow is subject to change and currently under development.

#### Build an initial network representing our current knowledge of diabetes 
Steps to build the network include:
   - Identifying DOID and HPO concepts annotated to diabetes
   - Identify genes associated with diabetes (Orphanet, OMIM, GAD, and UMLS)
   - Using the list of disease-associated genes, find annotations to other biomedical ontologies and knowledge sources

#### Add evidence and metadata to the diabetes knowledge network
PubAnnotation concept annotations will be mapped to concepts in the network
  - Concept co-occurrence across all annotated articles related to diabetes
  - Metadata will be added to the network using Colil and other sources other sources [8]
  - *Additional measures are under development*

#### Identify the ignorome
We will use publicly available gene expression data to identify the ignorome (procedure similar to Riba et al. [2])
  - Identify, merge, and analyze relevant GEO datasets
  - Produce a list of genes that are differentially expressed (DEGs) in diabetes
  - Map differentially expressed genes to the genes in the diabetes knowledge network created
  - The genes that we are not able to be mapped are diabetes-specific ignorome candidates
  - Initial validation of gene list using genes identified in Riba et al. [2]

#### Annotate the ignorome gene candidates
Annotate the genes using the same resources used when building the initial diabetes knowledge network
  - Determine if any of the ignorome concept annotations can be used to connect the ignorome genes to the diabetes knowledge network

#### Details regarding network inference under development

<br><br>

### Resources: Data, Software, and Tools
#### SPARQL Endpoints and RDF Code/Repositories
  * [PubAnnotation](http://sparql.pubannotation.org/) [3]
  * [DisGeNET](http://rdf.disgenet.org/sparql/) [4]
  * [Orphanet](http://www.orpha.net/sparql) [5]
  * The Knowledge Base of Biomedicine [6]
  * PMC article metadata [7]

#### Software and Tools
  * [Colil](http://colil.dbcls.jp/browse/papers/) [7]
  * [PubCases](https://pubcases.dbcls.jp/)

#### Public Data and Repositories
  * [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)

#### Dependencies
  * SPARQL Endpoint APIs
      - Java
      - Python
  * Network Generation, Inference, and Visualization
      - Python - Networkx
      - Matlab
      - Cytoscape
  * Bioinformatic Analyses
      - R + Bioconductor

<br><br>

### Collaborators
  * University of Colorado Denver Anschutz Medical Campus (AMC)
  * Database Center for Life Science (DBCLS)

<br><br>

#### References
  1. Pandey AK et al. PLoS One. 2014:PMC3921226 
  2. Riba et al. Sci Rep. 2016:PMC4838989 
  3. [PubAnnotation](http://pubannotation.org/)
  4. [DisGeNET](http://www.disgenet.org/web/DisGeNET/menu)
  5. [Orphanet](http://www.orpha.net/consor/cgi-bin/index.php)
  6. Livingston et al. BMC Bioinformatics. 2015:PMC4448321
  7. Fujimara et al. J Biomed Semantics. 2015:PMC4617487
  8. https://github.com/callahantiff/pmc-metadata
