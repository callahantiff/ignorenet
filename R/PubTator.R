#########################################################################################
### PubTator - code to retrieve abstracts from PubMed and run PubTator (adapted from Andrey Soares)
### version 1.2.0
### date: 10.11.17
#########################################################################################


# load needed libraries
library(rentrez)
library(RISmed)
library(RCurl)
library(XML)
library(plyr)
library(data.table)
library(dplyr)
library(plyr)
library(quanteda)
library(openNLP)
library(NLP)
library(annotate)
library(org.Hs.eg.db)

##########################################
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

##########################################
# search results and return genes list
getPubtatorMatches <- function(data,
                               info = c("Genes", "Diseases", "Mutations", "Chemicals", "Species"),
                               dir) {

  # Load required package
  require(pubmed.mineR, quietly = TRUE) && require(plyr, quietly = TRUE)
  
  # Loop through each abstract
  df <- data.frame(PMID = rep(NA, nrow(pubmed_res_filt)), 
                   gene = rep(NA, nrow(pubmed_res_filt)),
                   disease = rep(NA, nrow(pubmed_res_filt)),
                   mutation = rep(NA, nrow(pubmed_res_filt)),
                   species = rep(NA, nrow(pubmed_res_filt)), 
                   sentence = rep(NA, nrow(pubmed_res_filt))) 
  
  # create window to chunk data by
  pause <- seq(200, nrow(data), 200)
  
  # loop over data
  for (i in 1:nrow(data)) {
    print(i)
    
    # for chunk processing
    if(length(grep(paste("^",i,"$", sep=''), pause)) > 0) {
      # write dataset to file
      write.table(df, 
                  paste(dir, 'Pubtator_', i, '_results.txt', sep = ''), 
                  sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
      
      # pause system for 1 minute
      Sys.sleep(120)
    }
    
    # save row to variable 
    abstract <- data[i,]
    
    # extract PubTator information from selected abstract
    pubtator_output <- pubmed.mineR::pubtator_function(abstract["PMID"])

  # Create a data frame with all symbols found, corresponding title and sentences
  if (length(pubtator_output) != 1) {
    # Separate the abstract in different sentences
    sentences <- pubmed.mineR::SentenceToken(abstract["Abstract"])
  
    # Remove sentence with no matches
    sentences <- sentences[sentences!=""]
    keep_sent <- sentences

    # # loop over sentences and find those where a gene and disease co-occur
    # keep_sent <- c()
    # for (i in 1:length(sentences)) {
    #   sent <- sentences[[i]]
    # 
    #   if (any(lapply(pubtator_output[[1]], function(x) length(agrep(x, sent))) > 0)
    #       && any(lapply(pubtator_output[[2]], function(x) length(agrep(x, sent))) > 0)) {
    #     keep_sent <- c(keep_sent, sent)
    #   }
    # }
    
    # save output to data frame
    df[i, 1] <- abstract$PMID
    df[i, 2] <- ifelse(!is.null(pubtator_output[[1]]), paste(pubtator_output[[1]], collapse = " | "), NA)
    df[i, 3] <- ifelse(!is.null(pubtator_output[[2]]), paste(pubtator_output[[2]], collapse = " | "), NA)
    df[i, 4] <- ifelse(!is.null(pubtator_output[[3]]), paste(pubtator_output[[3]], collapse = " | "), NA)
    df[i, 5] <- ifelse(!is.null(pubtator_output[[5]]), paste(pubtator_output[[5]], collapse = " | "), NA)
    df[i, 6] <- ifelse(paste(keep_sent, collapse = "|") != "", paste(keep_sent, collapse = "|"), NA)
    
  }
  else {
    df[i, 1] <- abstract$PMID
    df[i, 2] <- NA
    df[i, 3] <- NA
    df[i, 4] <- NA
    df[i, 5] <- NA
    df[i, 6] <- NA
  }
  
  }
  
  # Order data frame by PMIDs
  df <- df[order(df$PMID),]
  
  # write out the last chunk of processed data
  write.table(df, 
              paste(dir, 'Pubtator_', nrow(data), '_results.txt', sep = ''), 
              sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
  }


##########################################
# Process gene list from PubTator results
results_process <- function(data, chunks) {
  # parse results to create a single gene list
  for (i in 1:nrow(data)) {
    keep_genes <- c()
    print(i)
    genes <- strsplit(pubtator_noNA[i,2], "[|]")
    
    # for chunk processing
    if(length(grep(paste("^",i,"$", sep=''), chunks)) > 0) {
      
      # pause system for 1 minute
      Sys.sleep(120)
    }
    
    for (j in 1:length(genes[[1]])) {
      # cat("j", j)
      gene2 <- sub('\\)|\\(', "", trimws(genes[[1]][j]))
      # print(gene2)
      
      if (length(grep(gene2, " ")) > 0 | nchar(gene2) > 7 | grepl("^[[:lower:]]+$", gene2)) {
        
        # convert protein name to gene identifier
        prot = paste(gene2, "[PROT]", sep = "")
        protein_number = rentrez::entrez_search(db="protein", term=prot)
        
        if(length(protein_number$ids) > 0) {
          protein_links <- rentrez::entrez_link(dbfrom='protein', id=as.numeric(protein_number$ids[1]), db='all')
          
          if(is.null(protein_links$links$protein_genome) && is.null(protein_links$links$protein_gene)) {
            protein_links <- gene2}
          else {
            protein_links <- ifelse(!is.null(protein_links$links$protein_gene),
                                    protein_links$links$protein_gene, protein_links$links$protein_genome)}
          
          # get annotations
          if(is.na(suppressWarnings(as.numeric(protein_links)))) {
            # check to convert roman numerals in gene names
            rom <- suppressWarnings(gsub(as.roman(strsplit(protein_links, " ")[[1]][length(strsplit(protein_links, " ")[[1]])]), 
                                         as.numeric(as.roman(strsplit(protein_links, " ")[[1]][length(strsplit(protein_links, " ")[[1]])])),
                                         protein_links))
            
            protein_links <- ifelse(is.na(rom), protein_links, rom)
            
            # annotate genes using either gene name or entrez id
            protein_links <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys=as.character(sub('-', ' ', protein_links)), columns=c('ENTREZID'), keytype='GENENAME'), error = function(e) NULL)
            # remove NA
            protein_links <- ifelse(is.null(protein_links), "", as.numeric(protein_links$ENTREZID))
            
            keep_genes <- c(keep_genes, protein_links)}
          else{keep_genes <- c(keep_genes, as.numeric(protein_links))}
        }
      }
      
      else {
        # convert protein name to gene identifier
        gene = paste(gene2, "[GENE]", sep = "")
        gene_number = rentrez::entrez_search(db="gene", term=gene)
        
        if(length(gene_number$ids) > 0) {
          gene_number <-  gene_number$ids[[1]]
          keep_genes <- c(keep_genes, as.numeric(gene_number))}
      }
    }
    
    # get gene symbols
    gene_annots <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys=as.character(unique(keep_genes)), columns=c('SYMBOL'),
                                                  keytype='ENTREZID'), error = function(e) NULL)
    
    # append solution to database as collapsed list
    data$entrez[i] <- paste(unique(keep_genes), collapse = "|")
    data$symbol[i] <- paste(unique(gene_annots$SYMBOL), collapse = "|")
  }
}
