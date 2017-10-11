# import the library
# load needed libraries
library(rentrez); library(RISmed); library(RCurl); library(XML); library(plyr); library(data.table); library(dplyr); library(org.Hs.eg.db); library(plyr); library(quanteda); library(openNLP); library(NLP)


# retrieve pubmed abstracts and titles for pmids in data set
getPubmedAbstracts <- function(pmids) {
  # Load required package
  require(RISmed, quietly = TRUE)
  
  # get information from Pubmed articles included in the list of
  message("Retrieving Pubmed data...")
  records = RISmed::EUtilsGet(unique(pmids))
  
  # Create a data frame with PMID, Title and Abstract
  df <- data.frame('PMID' = as.numeric(RISmed::PMID(records)),
                   'Title' = RISmed::ArticleTitle(records),
                   'Abstract'= RISmed::AbstractText(records),
                   stringsAsFactors = FALSE)
  
  # Order data frame by PMIDs
  df <- df[order(df$PMID),]
  
  # Return data frame with abstracts extracted from Pubmed
  df
}


# search results and return genes list
getPubtatorMatches <- function(abstracts,
                               info = c("Genes", "Diseases", "Mutations", "Chemicals", "Species")) {
  
  #-------------------------------------------------
  # Function to perform regular expression to a sentence
  findInSentences <- function(sentence, symbols) {
    tryCatch({
      # match patterns and returns sentences that have at least two matches within the same sentence
      ifelse(!is.null(symbols) & length(which(gregexpr(paste(tolower(symbols), collapse = "|"), tolower(sentence))[[1]] > 0)) > 1, sentence[[1]], "")
    },error = function(e) {
      # if the regular expression has an error, return the following error message
      "*** REGEX ERROR ***"
    })
  }
  
  # Load required package
  require(pubmed.mineR, quietly = TRUE) && require(plyr, quietly = TRUE)
  
  # Loop through each abstract
  df <- plyr::ddply(abstracts, .(PMID), function(abstract) {
    
    # extract PubTator information from selected abstract
    # pubtator_output <- pubmed.mineR::pubtator_function(abstract["PMID"])
    pubtator_output = NULL
    
    while(!inherits(pubtator_output,'list')) {
      pubtator_output <- tryCatch(
        pubmed.mineR::pubtator_function(abstract["PMID"]),
        error=function(e) {Sys.sleep(10); message('retrying');
          return(e)})
    }
    
    # Check what PubTator information to use for list of symbols.
    # The user can define any or all types: Genes, Diseases, Mutations, Chemicals and Species
    symbols <- c()
    for(pubType in c("Genes", "Diseases", "Mutations", "Chemicals", "Species")) {
      if(tolower(pubType) %in% tolower(info) & !is.null(pubtator_output[pubType][[1]])){
        symbols <- c(symbols, pubtator_output[pubType][[1]])
      }
    }
    
    # Separate the abstract in different sentences
    sentences <- pubmed.mineR::SentenceToken(abstract["Abstract"])
    
    # Loop throught the sentences and find matches of symbols
    sent <- plyr::laply(sentences, findInSentences, symbols = symbols)
    
    # Remove sentence with no matches
    sent <- sent[sent!=""]
    
    # find matches in the title of the abstract
    title <- findInSentences(abstract["Title"], symbols)
    
    # Create a data frame with all symbols found, corresponding title and sentences
    data.frame(symbols=paste(symbols, collapse = "|"), title=title, sentences=paste(sent, collapse = "|"))
  }, .progress = "text")
  
  # Order data frame by PMIDs
  df <- df[order(df$PMID),]
  
  # 0 = no match, no symbols found by Pubtator, or REGEX ERROR in the sentences
  # 1 = Match at least 2 symbols in one sentence/title
  df$match <- ifelse(df$symbols == "" | df$sentences == "" | grepl("REGEX ERROR", df$sentences), 0, 1)
  
  # Return data frame with Pubtator matches
  df
}

# retrieve pubmed IDs matching specific search criteria
r_search <- rentrez::entrez_search(db="pubmed", term="preeclampsia AND Homo Sapiens[ORGN]", retmax = 90000)


# get the abstracts
pubmed_res <- getPubmedAbstracts(r_search$ids)
pubmed_res$Abstract[pubmed_res$Abstract == ""] <- NA
pubmed_res_filt <- pubmed_res[complete.cases(pubmed_res$Abstract), ]

# return gene annotations from PubAnnotation 
pubtator <- getPubtatorMatches(pubmed_res_filt, info = c("Genes"))

# Save the results
saveRDS(gene_data_pairwise, "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/PubMed_Results/pairwise_comparison_results.rds")

# Load the results
pairwise <- readRDS("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_fullcomp/pairwise_comparison_results.rds")



