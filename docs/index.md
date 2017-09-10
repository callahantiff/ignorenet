# ignorenet

### Project Description
Integrating biological knowledge with experimental data is vital for understanding the mechanisms underlying complex diseases and clinical phenotypes. Most genes differentially expressed in association with disease have no currently described disease-associated function in the literature ("ignorome") [1-2]. The current project aims to leverage existing linked open biomedical resources and software/tools and publicly available high-throughput data repositories to identify and annotate the function of uncharacterized disease-associated genes for a specific disease. This repository documents our efforts to reveal and annotate disease-specific ignoromes.

#### BioHackathon 2017 Use Case
Preeclampsia affects 3-8% of all pregnancies accounting for 18% of maternal deaths and 40% of fetal mortality [4]. The only existing cure for preeclampsia is placental delivery and preterm birth [3]. Preeclampsia is good candidate for the ignorome because it is a substantial public health burden, has an largely unknown etiology, clear clinical characteristics, is well represented by publicly availible data, and has insufficient existing mouse models.
<br>

### Approach
<img src="https://github.com/callahantiff/ignorenet/blob/master/images/approach.png" width="800">

<br>

### Resources: Data, Software, and Tools
#### SPARQL Endpoints and RDF Code/Repositories
  * [PubAnnotation](http://sparql.pubannotation.org/) [5]
  * [DisGeNET](http://rdf.disgenet.org/sparql/) [6]
  * [Orphanet](http://www.orpha.net/sparql) [7]
  * The Knowledge Base of Biomedicine [8]
  * Human Variant Database [9]

#### Software and Tools
  * [Colil](http://colil.dbcls.jp/browse/papers/)
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

### Collaborators
  * University of Colorado Denver Anschutz Medical Campus
  * Database Center for Life Science
  * University of Maryland, Baltimore County

#### References
  1. Pandey AK et al. PLoS One. 2014:PMC3921226 
  2. Riba et al. Sci Rep. 2016:PMC4838989 
  3. Anderson et al. Placenta. 2012:PMID22197626 
  4. [Preeclampsia Foundation](https://www.preeclampsia.org/)
  5. [PubAnnotation](http://pubannotation.org/)
  6. [DisGeNET](http://www.disgenet.org/web/DisGeNET/menu)
  7. [Orphanet](http://www.orpha.net/consor/cgi-bin/index.php)
  8. Livingston et al. BMC Bioinformatics. 2015:PMC4448321
  9. Peterson et al. J Mol Biol. 2013:PMC3807015
