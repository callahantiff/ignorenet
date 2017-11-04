#########################################################################################
### GEO_Processing Pipeline
### version 2.1.0
### date: 10.24.17
#########################################################################################


## LIBRARIES
library(dplyr)
library(rentrez)
library(stats)
library(tidyr)

# load scripts
source("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/GEO_Processing.R")
source("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/PubTator.R")

############################################
## ANALYZE GSE STUDIES
# process all GSE ids in a list
gse_ids <- c("GSE4707", "GSE10588", "GSE25906", "GSE14722", "GSE24129", "GSE30186", "GSE47187", "GSE74341", "GSE35574","GSE44711", "GSE60436", "GSE73374")

full_results <- GEO_Pipeline(gse_ids, 'FULL')
comp_results <- GEO_Pipeline(gse_ids, 'COMP')

############################################
############################################
## PROCESS RESULTS
## full dataset
dir_full <- "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_studies"

# run post-processing and format data
full_dataset <- Post_Process(dir_full)
full_dataset <- Data_Format(full_dataset)
full_dataset <- Alias_Annotation(full_dataset)

# save the data objects
saveRDS(full_dataset, "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_studies/studies_results.rds")

## comparisons dataset
dir_comp <- "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_fullcomp"

# run post-processing and format data
comp_dataset <- Post_Process(dir_comp)
comp_dataset <- Data_Format(comp_dataset)
comp_dataset <- Alias_Annotation(comp_dataset)

# Save the data objects
saveRDS(comp_dataset, "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_fullcomp/pairwise_comparison_results.rds")


# # Load the data objects
# pairwise <- readRDS("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_fullcomp/pairwise_comparison_results.rds")
# studies <- readRDS("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_studies/studies_results.rds")


############################################
############################################
## GET GENE COUNTS
# get counts by study
grp_symb = as.data.frame(table(grp = full_dataset$group, symb = full_dataset$symbol))
tot_count <- stats::aggregate(grp_symb$Freq, by=list(symbol=grp_symb$symb), FUN=sum)
table(tot_count$x)

# get counts by group - second method
study.proportions <- group_by(full_dataset, symbol, group) %>% 
  summarize(n = length(group)) %>% # count per gender
  ungroup %>% group_by(symbol) %>% 
  mutate(count = sum(n)) 
study.proportions$n <- NULL

comp.proportions <- group_by(comp_dataset, symbol, group) %>% 
  summarize(n = length(group)) %>% # count per gender
  ungroup %>% group_by(symbol) %>% 
  mutate(count = sum(n)) 
comp.proportions$n <- NULL

# get count of unique genes across all comparisons
gene_count <- unique(comp_dataset$symbol)


############################################
############################################
## PUBTATOR ANNOTATION
# retrieve pubmed IDs matching specific search criteria
r_search <- rentrez::entrez_search(db="pubmed", term="preeclampsia AND Homo Sapiens[ORGN]", retmax = 90000)

###### ADRIANNE
# r1_search <- rentrez::entrez_search(db="pubmed", term="(asthma OR
#                                     chronic obstructive pulmonary disease OR
#                                     idiopatic pulmonary fibrosis OR
#                                     emphysema)
#                                     AND (Homo Sapiens[ORGN] OR Mus Musculus[ORGN])
#                                     AND 2007:2017[PDAT]", retmax = 90000)

dir = 'DE_Results/Annotations/PubTator_Results/Adrianne/'

# save the results
saveRDS(r_search, "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/PubMed_Results/pubmed_search_res.rds")
#adrianne
# saveRDS(r1_search, paste(dir, "pubmed_search_res.rds", sep = ""))

# load data
r_search <- readRDS("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/PubMed_Results/pubmed_search_res.rds")


# get the abstracts
pubmed_res <- getPubmedAbstracts(r1_search$ids)
pubmed_res$Abstract[pubmed_res$Abstract == ""] <- NA
pubmed_res_filt <- pubmed_res[complete.cases(pubmed_res$Abstract), ]

# save the results
saveRDS(pubmed_res_filt, "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/PubMed_Results/pubmed_abstracts.rds")
#adrianne
# saveRDS(pubmed_res_filt, paste(dir, "pubmed_abstracts.rds", sep = ""))

# load data
pubmed_res_filt <- readRDS("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/PubMed_Results/pubmed_abstracts.rds")
# adrianne
# pubmed_res_filt <- readRDS(paste(dir, "pubmed_abstracts.rds", sep = ""))


# return gene annotations from PubAnnotation
# chunk processesing option - NOT CURRENTLY USING
pubmed_res_filt$chunkid <- as.numeric(ggplot2::cut_width(1:nrow(pubmed_res_filt),5,boundary=0))
result = ddply(pubmed_res_filt, .(chunkid), getPubtatorMatches(pubmed_res_filt, info = c("Genes", "Diseases", "Mutations", "Chemicals", "Species")))

# not using chunk processing - USING
# dir <- "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Annotations/PubTator_Results/Adrianne/"
pubtator <- getPubtatorMatches(pubmed_res_filt, 
                               info = c("Genes", "Diseases", "Mutations", "Chemicals", "Species"),
                               dir)

# read in directory of output files
# set directory to location of files
dir_pt <- "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/PubMed_Results/Annotations/PubTator_Results/text_files"
dir_pt2 <- "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/PubMed_Results/Annotations/PubTator_Results/text_files/Old"

# adrianne
# dir_pt <- "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Annotations/PubTator_Results/Adrianne/"

# run post-processing and format data
pubtator <- Post_Process(dir_pt)
pubtator2 <- Post_Process(dir_pt2)

# remove rows missing a gene and disease annotation
pubtator_noNA <- pubtator %>% drop_na(gene, disease)
pubtator_noNA2 <- pubtator2 %>% drop_na(gene, disease)

# remove duplicate rows
pubtator_noNA <- pubtator_noNA[!duplicated(pubtator_noNA[,2:6]), ]
pubtator_noNA2 <- pubtator_noNA2[!duplicated(pubtator_noNA2[,2:6]), ]


write.table(pubtator_noNA2, 
            paste("~/Desktop/", 'Preeclampsia_results.txt', sep = ''), 
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(pubtator_noNA, 
            paste("~/Desktop/", 'results2.txt', sep = ''), 
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)

# ADRIANNE
# for (i in 1:nrow(dataset)) {
#   dis = dataset[i, 6]
#   dis_rate = dataset[i, 3]
#   
#   print(i)
#   
#   if(length(grep("idiopathic", dis))>0 && length(grep("idiopathic", dis_rate))>0){
#     dataset$disease_type[i] <- "IDIOPATHIC PULMONARY FIBROSIS"
#   }
# 
#   if(length(grep("asthma", dis)) >0 && length(grep("asthma", dis_rate))>0){
#     dataset$disease_type[i] <- "ASTHMA"
#   }
# 
#   if(length(grep("emphysema", dis))>0 && length(grep("emphysema", dis_rate))>0){
#     dataset$disease_type[i] <- "EMPHYSEMA"
#   }
#   if(length(grep("chronic obstructive pulmonary disease|COPD", dis))>0  && length(grep("chronic obstructive pulmonary disease", dis_rate))>0){
#     dataset$disease_type[i] <- "COPD"
#   }
# }

# create window to chunk data by
chunks <- seq(100, nrow(pubtator_noNA), 100)
res <- results_process(pubtator_noNA, chunks)

# Save the results
saveRDS(pubtator_noNA, "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/PubMed_Results/pubtator_co-occur_processed.rds")




