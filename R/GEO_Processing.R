#########################################################################################
### GEO_Processing - code that contians functions needed to support ignorenets
### version 2.2.0
### date: 10.17.17
#########################################################################################


# LIBRARIES
require(pacman)
library(AnnotationDbi)
library(annotate)
library(Biobase)
library(data.table)
library(EMA)
library(GEOmetadb)
library(GEOquery)
library(ggplot2)
library(gplots)
library(knitr)
library(latex2exp)
library(limma)
library(mygene)
library(RColorBrewer)
library(R.utils)
library(GEOquery)
library(annotate)
library(hgug4112a.db) #GSE4707
library(org.Hs.eg.db) #GSE10588
library(illuminaHumanv2.db) #GSE25906, GSE35574
library(hgu133a.db) #GSE14722
library(hugene10sttranscriptcluster.db) #GSE24129
library(illuminaHumanv4.db) #GSE30186, GSE44711
library(hugene20sttranscriptcluster.db) #GSE73374


# # connect to database
# con <- DBI::dbConnect(SQLite(),'GEO_Database/GEOmetadb.sqlite')
# geo_tables <- dbListTables(con)
# dbListFields(con,'gpl')
# dbGetQuery(con,'select bioc_package from gpl where gpl = "GPL570"')

############################################
############################################
## QUERY GEO
# Download recent GEO data
getLatestGEO <- function(){
  ## Download and unzip database
  require(GEOmetadb, quietly = TRUE)
  
  # temporarily set working directory to sub-folder
  current_wd <- getwd()
  setwd(paste(current_wd, "/ignorenet/GEO_Database", sep=""))
  
  # download and compress the data using SQLite - checks for file before automatically downloading
  if(!file.exists('GEOmetadb.sqlite')) GEOmetadb::getSQLiteFile()
  
  # reset working directory
  setwd(current_wd)
}

############################################
# Get list of all organisms in GEO
GEOChoices <- function() {
  
  # load libraries
  pacman::p_load(GEOmetadb, DBI)
  
  # open database connection
  con <- DBI::dbConnect(SQLite(),'/ignorenet/GEO_Database/GEOmetadb.sqlite')
  
  # set query
  choices <- DBI::dbGetQuery(con,'SELECT DISTINCT organism FROM gpl')
  
  # close database
  DBI::dbDisconnect(con)
  return(choices)
}

############################################
# Query GEO using user input
GEOQuery <- function(text1, text2) {
  
  # load libraries
  pacman::p_load(GEOmetadb, DBI)
  
  # download GEO
  # getLatestGEO()
  
  # verify input has valid data
  disease <- ifelse(text1 == "", "*", text1)
  org <- ifelse(text2 == "", "*", text2)
  
  # open database connection
  con <- DBI::dbConnect(SQLite(),'/ignorenet/GEO_Database/GEOmetadb.sqlite')
  
  # set query
  sql <- paste("SELECT DISTINCT
               gse.gse AS GSE,
               gse.title AS Title,
               gpl.organism AS Organism,
               gse.type AS Type,
               gsm.gpl AS GPL,
               gpl.manufacturer AS Platform,
               gse.summary AS Summary,
               gse.status AS Status,
               gse.pubmed_id AS PMID
               FROM gsm
               JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm
               JOIN gse ON gse_gsm.gse=gse.gse
               JOIN gse_gpl ON gse_gpl.gse=gse.gse
               JOIN gpl ON gse_gpl.gpl=gpl.gpl
               WHERE",
               "LOWER(gse.type) like '%expression profiling by array%' AND",
               "LOWER(gse.title) LIKE", paste('"%', disease, '%"', sep=''), "AND",
               "LOWER(gpl.organism) LIKE", paste('"%', org, '%"', sep=''), sep=" ")
  
  # query GEO
  results <- DBI::dbGetQuery(con, sql)
  
  # close database
  DBI::dbDisconnect(con)
  
  return(results)
  
}

############################################
############################################
## GET GEO DATA
GEOData <- function(ID) {
  # load libraries
  pacman::p_load(Biobase, GEOquery, limma)
  
  # gse file
  # dirs <- grep("_output",list.dirs("./ignorenet/",recursive=FALSE),value=TRUE)
  setwd("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/ignorenet/GEO_output")
  # file_list <- list.files()
  # gse_file <- list.files(file_list, ID, recursive=TRUE, full.names=TRUE)
  # file <- paste(ID, "_series_matrix.txt.gz", sep = "")
  gse_file <- list.files(pattern=paste(ID, "_series_matrix.txt.gz", sep = ""), recursive=TRUE)
  
  # check length of files
  gse_file <- ifelse(length(gse_file) ==1, gse_file, gse_file[1])
  
  # pull GSE data
  gset = NULL
  if(!file.exists(gse_file)) {
    print(paste("Downloading", ID))
    
    while(!inherits(gset,'list')) {
      gset=tryCatch(
        suppressWarnings(GEOquery::getGEO(ID, destdir = './ignorenet/GEO_output/',
                                          GSEMatrix=TRUE)),
        error=function(e) {Sys.sleep(10); message('retrying');
          return(e)})
    }
    
    # keep only the first platform annotations
    if (length(gset) > 1) idx <- grep(gset[[1]]@annotation, attr(gset, 'names')) else idx <- 1
    
    # subset data
    gset <- gset[[1]]} 
  
  # if the data is already downloaded then load it
  else {
    gset <- suppressWarnings(GEOquery::getGEO(filename = gse_file))
  }
  
  # update column names
  Biobase::fvarLabels(gset) <- make.names(Biobase::fvarLabels(gset))
  
  setwd("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R")
  
  # remove file holders
  rm(gse_file)
  
  return(gset)
}


############################################
## REMOVE CONTROL PROBES
Probe_collapse <- function(gset){
  # convert gset to expression data
  exp <- Biobase::exprs(gset)
  
  # get annotation data
  annot_data <- Biobase::fData(gset)
  
  if(ID == 'GSE10588') {
    annot_data <- annot_data[annot_data$Status == 'current',]}
  
  if(ID == 'GSE25906' | ID == 'GSE14722' | ID == 'GSE30186' | ID == 'GSE35574' | ID == 'GSE44711' | ID == 'GSE60438' | ID == 'GSE73374') {
    annot_data <- annot_data}
  
  if(ID == 'GSE74341') {
    annot_data <- annot_data[as.character(annot_data$SPOT_ID) == as.character(annot_data$NAME),]}
  
  if(ID == 'GSE24129') {
    annot_data <- annot_data[annot_data$category == 'main',]
  }
  
  if(ID == "GSE4707") {
    # remove probes that are of control type
    annot_data <- annot_data[annot_data$CONTROL_TYPE == 'FALSE',]
  }
  
  # subest expression data to only include non-control type probes
  exp <- exp[rownames(exp) %in% annot_data$ID, ]
  
  # renames exp rownames
  row.names(exp) <- rownames(exp)
  
  return(exp)
}


############################################
## LOG-TRANSFORM AND FILTER DATA
GEO_preprocess <- function(exp) {
  # Log2 transformation - determine if data has been log2 transformed
  ex <- exp
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  
  if(LogC) { 
    ex[ex <= 1.0] <- 1.0
    exprs_data <- log2(ex) 
  }
  else {
    exprs_data <- ex
  }
  
  # Filtering - summarize average gene expression values
  avg <- rowMeans(exprs_data)
  sum <- summary(avg)
  
  # filter out genes that "aren't expressed" - will appear red in the plot
  cutoff <- sum[2][[1]] # remove genes below first quantile
  exp_filt <- expFilter(exprs_data, 
                        threshold = round(sum[2][[1]], 4), graph = FALSE)
  
  return(exp_filt)
}


############################################
## CREATE GROUPS FOR DE ANALYSIS
Targets_full <- function(ID, gset){
  
  # get phenotype data about samples
  target_data <- pData(gset)
  
  if(ID == 'GSE10588'| ID == 'GSE30186') {
    groups <- as.character(lapply(as.character(tolower(target_data$title)), function(x) ifelse(length(grep('normal|control', x)) == 1, 'Control', 'Case')))}
  
  if(ID == 'GSE14722') {
    groups <- as.character(lapply(as.character(tolower(target_data$title)), function(x) ifelse(length(grep('preterm', x)) == 1, 'Control', 'Case')))}
  
  if(ID == 'GSE73374') {
    groups <- as.character(lapply(as.character(tolower(target_data$title)), function(x) ifelse(length(grep('norm', x)) == 1, 'Control', 'Case')))}
  
  if(ID == 'GSE44711') {
    groups <- as.character(tolower(target_data$title))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('eo', x)) == 1, 'Case', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('control', x)) == 1, 'Control', x))) }
  
  if(ID == 'GSE74341') {
    groups <- as.character(tolower(target_data$title))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('eo|lo', x)) == 1, 'Case', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('preterm|at', x)) == 1, 'Control', x))) }
  
  if(ID =='GSE4707') {
    groups <- as.character(tolower(target_data$title))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('eo|lo', x)) == 1, 'Case', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('_n', x)) == 1, 'Control', x)))}
  
  if(ID =='GSE24129') {
    groups <- as.character(tolower(target_data$title))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('norm', x)) == 1, 'Control', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('fetal|pre', x)) == 1, 'Case', x)))}
  
  if(ID == 'GSE25906' | ID == 'GSE35574') {
    groups <- as.character(lapply(as.character(tolower(target_data$characteristics_ch1.2)), function(x) ifelse(length(grep('control', x)) == 1, 'Control', 'Case'))) }
  
  if(ID == 'GSE60438') {
    groups <- as.character(lapply(as.character(tolower(target_data$characteristics_ch1)), function(x) ifelse(length(grep('norm', x)) == 1, 'Control', 'Case'))) }
  
  if(ID == "GSE35574") {
    groups <- as.character(tolower(target_data$characteristics_ch1.2))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('control', x)) == 1, 'Control', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('iugr|pe', x)) == 1, 'Case', x)))}
  
  return(groups)
}


Targets_comparisons <- function(ID, gset){
  
  # get phenotype data about samples
  target_data <- pData(gset)
  
  if(ID == 'GSE10588'| ID == 'GSE30186' | ID == 'GSE44711') {
    groups <- as.character(lapply(as.character(tolower(target_data$title)), function(x) ifelse(length(grep('normal|control', x)) == 1, 'Control', 'Case')))}
  
  if(ID == 'GSE14722') {
    groups <- as.character(lapply(as.character(tolower(target_data$title)), function(x) ifelse(length(grep('preterm', x)) == 1, 'Control', 'Case')))}
  
  if(ID == 'GSE73374') {
    groups <- as.character(lapply(as.character(tolower(target_data$title)), function(x) ifelse(length(grep('norm', x)) == 1, 'Control', 'Case')))}
  
  if(ID == 'GSE74341') {
    groups <- as.character(tolower(target_data$title))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('eo', x)) == 1, 'EO', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('lo', x)) == 1, 'LO', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('at', x)) == 1, 'AT', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('preterm', x)) == 1, 'PT', x)))}
  
  if(ID =='GSE4707') {
    groups <- as.character(tolower(target_data$title))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('eo', x)) == 1, 'EO', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('lo', x)) == 1, 'LO', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('_n', x)) == 1, 'Control', x)))}
  
  if(ID =='GSE24129') {
    groups <- as.character(tolower(target_data$title))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('norm', x)) == 1, 'Control', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('fetal', x)) == 1, 'FGR', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('pre', x)) == 1, 'Case', x)))}
  
  if(ID == 'GSE25906') {
    groups <- as.character(lapply(as.character(tolower(target_data$characteristics_ch1.2)), function(x) ifelse(length(grep('control', x)) == 1, 'Control', 'Case'))) }
  
  if(ID == 'GSE60438') {
    groups <- as.character(lapply(as.character(tolower(target_data$characteristics_ch1)), function(x) ifelse(length(grep('norm', x)) == 1, 'Control', 'Case'))) }
  
  if(ID == "GSE35574") {
    groups <- as.character(tolower(target_data$characteristics_ch1.2))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('control', x)) == 1, 'Control', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('iugr', x)) == 1, 'IUGR', x)))
    groups <- as.character(lapply(groups, function(x) ifelse(length(grep('pe', x)) == 1, 'Case', x)))}
  
  return(groups)
}


############################################
## REFORMAT OUTPUT FILES
Output_Format <- function(ID, comp_type) {
  
  if (comp_type != 'full') {
    # read in data
    data <- read.table(paste('DE_Results/', ID ,'_DEgenes_pairwise.txt', sep = ''), header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
    
    if(ID == 'GSE10588' | ID == 'GSE25906' | ID == 'GSE14722' | ID == 'GSE30186' | ID == 'GSE44711' | ID == 'GSE60438' | ID == 'GSE73374'){
      # reformat data
      data1 <- data[, c(4,9:11)]
      data1$group <- rep(paste(ID, '_PEvC', sep = ''), nrow(data1))
      
      # rearrange columns
      data1 <- data1[, c(2,3,1,4,5)]
      
      # rename columns
      names(data1) <- c("probe", "id", "pvalue", "symbol", "group")
      
      # remove rows with p-value >= 0.05
      data1 <- data1[data1$pvalue <= 0.05, ]
      
      # write both datasets to file
      write.table(data1, 
                  paste('DE_Results/Final_fullcomp/', ID, '_PEvsC','.txt', sep = ''), 
                  sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
    }
    
    else {
      # reformat data
      data1 <- data[, c(6,14:16)]
      grp_name1 <- strsplit(names(data1[1]), '[.]')[[1]][3]
      data1$group <- rep(paste(ID, '_', grp_name1, sep = ''), nrow(data1))
      
      # rearrange columns
      data1 <- data1[, c(2,3,1,4,5)]
      
      # rename columns
      names(data1) <- c("probe", "id", "pvalue", "symbol", "group")
      
      # remove rows with p-value >= 0.05
      data1 <- data1[data1$pvalue <= 0.05, ]
      
      # data set 2
      data2 <- data[, c(7,14:16)]
      grp_name2 <- strsplit(names(data2[1]), '[.]')[[1]][3]
      data2$group <- rep(paste(ID, '_', grp_name2, sep = ''), nrow(data2))
      
      # rearrange columns
      data2 <- data2[, c(2,3,1,4,5)]
      
      # rename columns
      names(data2) <- c("probe", "id", "pvalue", "symbol", "group")
      
      # remove rows with p-value >= 0.05
      data2 <- data2[data2$pvalue <= 0.05, ]
      
      # write both datasets to file
      write.table(data1, 
                  paste('DE_Results/Final_fullcomp/', ID, '_', grp_name1,'.txt', sep = ''), 
                  sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
      
      write.table(data2, 
                  paste('DE_Results/Final_fullcomp/', ID, '_', grp_name2,'.txt', sep = ''), 
                  sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
    }
  }
  
  if (comp_type == 'full') {
    
    # read in data
    data <- read.table(paste('DE_Results/', ID ,'_DEgenes_studies.txt', sep = ''), header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
    
    # reformat dataå
    data1 <- data[, c(4,9:11)]
    data1$group <- rep(paste(ID, '_PEvC', sep = ''), nrow(data1))
    
    # rearrange columns
    data1 <- data1[, c(2,3,1,4,5)]
    
    # rename columns
    names(data1) <- c("probe", "id", "pvalue", "symbol", "group")
    
    # remove rows with p-value >= 0.05
    data1 <- data1[data1$pvalue <= 0.05, ]
    
    # write both datasets to file
    write.table(data1, 
                paste('DE_Results/Final_studies/', ID, '_PEvsC','.txt', sep = ''), 
                sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
}


############################################
## PERFORM DIFFERENTIAL EXPRESSION ANALYSIS
GEO_DEA_full <- function(targets, gset, gset_pre, ID) {
  
  # create design matrix
  design <- stats::model.matrix(~ targets + 0, data.frame(gset_pre))
  colnames(design) <- gsub('targets', '', colnames(design))
  
  # make contrasts
  if(length(unique(targets)) == 2) {
    contr.matrix <- limma::makeContrasts(
      CasevsControl = Case-Control,
      levels = colnames(design))
    contr.matrix }
  
  # run limma
  fit <- limma::lmFit(data.frame(gset_pre), design)
  fit <- limma::contrasts.fit(fit, contrasts=contr.matrix)
  efit <- limma::eBayes(fit)
  
  # print number of up and down regulated genes between the samples
  top_dec_test <- decideTests(efit)
  print(summary(top_dec_test))
  
  # get gene annotations
  annot <- Biobase::fData(gset)

  # get symbols
  if(ID == 'GSE4707') {
    annot_db = 'hgug4112a'
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$SPOT_ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE35574') {
    annot_db = 'illuminaHumanv2'
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE44711') {
    annot_db = 'illuminaHumanv4'
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE24129') {
    annot_db = 'hugene10sttranscriptcluster'
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE73374') {
    annot_db = 'hugene20sttranscriptcluster'
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE30186') {
    annot_db = 'illuminaHumanv4'
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE14722') {
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID),
                                       as.character(lapply(as.character(annot$Gene.Symbol), function(x) strsplit(as.character(x), " /// ")[[1]][1]))), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    row.names(genes) <- as.character(genes$V1)
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE74341') {
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), as.character(annot$NAME), as.character(annot$GENE_SYMBOL)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    row.names(genes) <- as.character(genes$V1)
  }
  
  if(ID == 'GSE60438') {
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), as.character(annot$Symbol)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE25906') {
    annot_db = 'illuminaHumanv2'
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE10588'){
    genes_annot <- annot[, c(17,1,6)]
    genes <- genes_annot[genes_annot$ID %in% row.names(gset_pre), ]
  }
  
  # subest annotation data to only include the same rows as exp data
  names(genes) <- c("probe", "id", "symbol")
  
  # fill missing probe identifiers with "NA"
  genes$probe[genes$probe == ""] <- NA
  
  # add gene information to efit
  efit$genes <- genes
  
  # process Results
  limma::write.fit(efit, top_dec_test, adjust='BH', file=paste('DE_Results/', ID ,'_DEgenes_studies.txt', sep = ''), sep = '\t')
  
}

GEO_DEA_comparisons <- function(targets, gset, gset_pre, ID) {
  
  # create design matrix
  design <- stats::model.matrix(~ targets + 0, data.frame(gset_pre))
  colnames(design) <- gsub('targets', '', colnames(design))
  
  # make contrasts
  if(ID == 'GSE4707') {
    contr.matrix <- limma::makeContrasts(
      EOvsC = EO-Control,
      LOvsC = LO-Control,
      levels = colnames(design))
    contr.matrix
  }
  
  if(ID == 'GSE74341') {
    contr.matrix <- limma::makeContrasts(
      EOvsPT = EO-PT,
      LOvsAT = LO-AT,
      levels = colnames(design))
    contr.matrix}
  
  if(ID == 'GSE35574') {
    contr.matrix <- limma::makeContrasts(
      IUGRvsC = IUGR-Control,
      CasevsControl = Case-Control,
      levels = colnames(design))
    contr.matrix }
  
  if(ID == 'GSE24129') {
    contr.matrix <- limma::makeContrasts(
      FGRvsC = FGR-Control,
      CasevsControl = Case-Control,
      levels = colnames(design))
    contr.matrix}
  
  if(length(unique(targets)) == 2) {
    contr.matrix <- limma::makeContrasts(
      CasevsControl = Case-Control,
      levels = colnames(design))
    contr.matrix}
  
  # run limma
  fit <- limma::lmFit(data.frame(gset_pre), design)
  fit <- limma::contrasts.fit(fit, contrasts=contr.matrix)
  efit <- limma::eBayes(fit)
  
  # print number of up and down regulated genes between the samples
  top_dec_test <- decideTests(efit)
  print(summary(top_dec_test))
  
  # get gene annotation
  annot <- Biobase::fData(gset)
  
  # get symbols
  if(ID == 'GSE4707') {
    annot_db = 'hgug4112a'
    genes_annot <- as.data.frame(cbind(annot$ID, annotate::getSYMBOL(as.character(annot$SPOT_ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE24129') {
    annot_db = 'hugene10sttranscriptcluster'
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE35574') {
    annot_db = 'illuminaHumanv2'
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE60438') {
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), as.character(annot$Symbol)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE44711') {
    annot_db = 'illuminaHumanv4'
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE74341') {
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), as.character(annot$NAME), as.character(annot$GENE_SYMBOL)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    row.names(genes) <- as.character(genes$V1)
  }
  
  if(ID == 'GSE25906') {
    annot_db = 'illuminaHumanv2'
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE14722') {
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID),
                                       as.character(lapply(as.character(annot$Gene.Symbol), function(x) strsplit(as.character(x), " /// ")[[1]][1]))), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    row.names(genes) <- as.character(genes$V1)
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE10588'){
    genes_annot <- annot[, c(17,1,6)]
    genes <- genes_annot[genes_annot$ID %in% row.names(gset_pre), ]
  }
  
  if(ID == 'GSE30186') {
    annot_db = 'illuminaHumanv4'
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  if(ID == 'GSE73374') {
    annot_db = 'hugene20sttranscriptcluster'
    annot <- annot[!is.na(annot$ID),]
    genes_annot <- as.data.frame(cbind(as.character(annot$ID), annotate::getSYMBOL(as.character(annot$ID), annot_db)), stringsAsFactors = FALSE)
    
    # reduce size of file
    genes <- genes_annot[genes_annot$V1 %in% row.names(gset_pre), ]
    genes <- cbind(row.names(genes), genes)
  }
  
  # subest annotation data to only include the same rows as exp data
  names(genes) <- c("probe", "id", "symbol")
  
  # fill missing probe identifiers with "NA"
  genes$probe[genes$probe == ""] <- NA
  
  # add gene information to efit
  efit$genes <- genes
  
  # process Results
  limma::write.fit(efit, top_dec_test, adjust='BH', file=paste('DE_Results/', ID ,'_DEgenes_pairwise.txt', sep = ''), sep = '\t')
  
}


############################################
## RUN FULL PIPELINE
GEO_Pipeline <- function(gse_ids, comp_type) {
  
  for (i in 1:length(gse_ids)) {
    # set ID
    ID <- gse_ids[i]
    print(ID)
    
    # pull data from GEO
    gset <- GEOData(ID)
    
    # remove control probes and collapse probes annotated to same gene
    exp <- Probe_collapse(gset)
    dim(exp)[1]
    
    # verify log-transformation & filter genes with exp below 1 quantile
    gset_pre <- GEO_preprocess(exp)
    dim(gset_pre)[1]
    
    # run at study-level
    if(comp_type == "full") {
      # get targets
      targets <- Targets_full(ID, gset)
      print(table(targets))
      
      # differential expression analysis
      DE <- GEO_DEA_full(targets, gset, gset_pre, ID)
      
      # format ourput
      Output_Format(ID, comp_type)
    }
    
    # run all possible comparisons
    if(comp_type != "full") {
      # get targets
      targets <- Targets_comparisons(ID, gset)
      print(table(targets))
      
      # differential expression analysis
      DE <- GEO_DEA_comparisons(targets, gset, gset_pre, ID)
      
      # format ourput
      Output_Format(ID, "comp")
    }
  }
}


############################################
############################################
## PROCESS RESULTS
# LOAD RESULTS INTO SINGLE DATAFRAME
Post_Process <- function(dir_loc) {
  
  # set directory to location containing files
  setwd(dir_loc)
  
  # get files in directory
  file_list <- list.files()
  
  for (file in file_list){
    
    if(length(grep('.txt', file)) != 0) {
      print(file)
      
      # if the merged dataset doesn't exist, create it
      if (!exists("dataset")){
        dataset <- read.table(file, header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")}
      
      # if the merged dataset does exist, append to it
      if (exists("dataset")){
        temp_dataset <- read.table(file, header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
        dataset <- rbind(dataset, temp_dataset)
        rm(temp_dataset)}}
  }
  # get entrez ids
  # remove genes with NA
  dataset <- dataset[complete.cases(dataset[3:5]),]
  dataset$id <- rownames(dataset)
  
  # remove duplicate genes
  dataset_updated = data.frame()
  for (study in unique(dataset$group)){
    # subset data by study
    study_sub <- dataset[which(dataset$group==study), ]
    
    # remove duplicates
    res <- study_sub[!duplicated(study_sub$symbol),]
    
    # append results to dataframe
    dataset_updated <- rbind(dataset_updated, res)
  }

  # re-set directory to original project location
  setwd("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R")
  return(dataset_updated)
  rm(dataset)
}


############################################
## CLEAN UP FORMATTING
Data_Format <- function(dataset) {
  
  for (i in 1:length(dataset$symbol)) {
    print(i)
    gene <- dataset$symbol[i]
    dataset$symbol[i] <- gene
    
    if(gene == "" | gene == "N/A" | is.na(gene) | gene == "NA"){
      dataset$symbol[i] <- NA}
    
    if(length(grep("\\d+\\-[A-Z]+", gene)) > 0){
      num <- gsub("\\-[A-Z]+", "", gene)
      mon <- gsub("\\d+\\-", "", gene)
      
      if(mon == "MAR" | mon == "Mar") {
        dataset$symbol[i] <- paste("MARCH", num, sep="")}
      if(mon == "SEP" | mon == "Sep") {
        dataset$symbol[i] <- paste("SEPT", num, sep="")}
      if(mon == "NOV" | mon == "DEC" | mon == "APR" | mon == "Nov" | mon == "Dec" | mon == "Apr") {
        dataset$symbol[i] <- paste(toupper(mon), num, sep="")}}
    
    if(length(grep('/', gene)) > 0 && length(grep('\\-', gene)) == 0 && gene != "N/A"){
      dataset$symbol[i] <- gsub('"', "", strsplit(gene, " // ")[[1]][1])}
    
    if(length(grep(' // ', gene)) > 0) {
      dataset$symbol[i] <- gsub('"', "", strsplit(gene, " // ")[[1]][1])}
    
    if(length(grep(' /// |\\,', gene)) > 0) {
      dataset$symbol[i] <- gsub('"', "", strsplit(gene, " /// ")[[1]][1])}
    
    if(length(grep('\\,', gene)) > 0) {
      dataset$symbol[i] <- gsub('"', "", strsplit(gene, ",")[[1]][2])}
  }
  dataset$entrez <- NULL
  
  return(dataset)
}



