# Meeting Notes - Preparing for Hackathon

## 05/22/2017

1. What Tiffany needs to have complete before arriving at the hackathon:
    * For at least one disease, the expression data needs to be analyzed and ready for integration via list of differential expressed genes
    * Develop, build, and query KaBOB to build an initial knowledge networks for at least 1 common and 1 rare disease
    * Work on pipeline for mapping PubAnnotation concepts to the knowledge network. This is meant to provide us with a method to add evidence to the knowledge networks. 
    * Some initial ideas for how to add this evidence include:
      - Co-occurrence (low precision) - will have to decide if we do this on the sentence level or if we do it on the article level.
      - Parsing-Based Approach (medium precision) - Find sentences that include 2 interesting entities - if there is a syntactic connection then we can assume that this represents some relationship between the 2 entities
      - Explicit Relations (high precision) - when we have relation annotations these will be used.
      - I am considering leveraging the predicates from SemRep. 
      - Also discussed how it would be useful to be able to map relations to the Relation Ontology (this may not be within the scope of the hackathon)

2. PubAnnotation SPARQL endpoint update - been delayed because a more permanent solution is being implemented. Really looking forward to using this resource! :D

3. Alignment of the ignorome with PubCases work (https://pubcases.dbcls.jp/search). There are many different ways to think about an ignorome:
   * The way that I proposed - trying to identify genes that are ignored in the literature
   * Toyofumi’s method
     - There are many databases - cases that are published as case reports, only part of them are curated in the databases.
     - Want to better identify phenotypes underlying rare diseases to help physicians improve their diagnostic abilities (is this correct?) 
     - Many of these services use data from OMIM or Orphanet, but the data in OMIM is curated (human) - there is more information that can be harvested through literature mining. Looking for connections between phenotypes to diseases: gender; age; trying to find a good combination of ontologies to represent the cases (e.g., Orphanet, OMIM, HPO, DOID; ICD-10. 
     - Target user of the system: medical doctors

4. Tiffany will do the following:
   * Share the list of ontologies/sources proposed for use with the Ignorome with Jin-Dong
     - DOID
     - HPO
     - GO
     - Orphanet
     - Protein Ontology
     - Reactome
     - Uberon/pato
     - ChEBI

5. Develop and refine the knowledge representation of an initial knowledge network with examples/clear details regarding how to map it to PubAnnotation
