#########################################################################################
### GEO_Preprocessing - code that contians functions needed to support ignorenets
### version 1.0.0
### date: 07.11.17
#########################################################################################


## LIBRARIES
require(pacman)
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

############################################
############################################
## QUERY GEO
# Download recent GEO data
getLatestGEO <- function(){
  ## Download and unzip database
  require(GEOmetadb, quietly = TRUE)
  
  # temporarily set working directory to sub-folder
  current_wd = getwd()
  setwd(paste(current_wd, "/ignorenet/GEO_Database", sep=""))
  
  # download and compress the data using SQLite - checks for file before automatically downloading
  if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
  
  # reset working directory
  setwd(current_wd)
}

############################################
# Get list of all organisms in GEO
GEOChoices <- function() {
  
  # load libraries
  pacman::p_load(GEOmetadb, DBI)
  
  # open database connection
  con <- DBI::dbConnect(SQLite(),'GEO_Database/GEOmetadb.sqlite')
  
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
  con <- DBI::dbConnect(SQLite(),'GEO_Database/GEOmetadb.sqlite')
  
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
## ANALYZE GSE STUDIES
# http://wiki.bits.vib.be/index.php/Analyse_GEO2R_data_with_R_and_Bioconductor

## USE CASE
gse_ids = c("GSE73374",
            "GSE24129")

# Get GEO data
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
  gse_file = ifelse(length(gse_file) ==1, gse_file, gse_file[1])
  
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
    gset <- gset[[1]]
    
  } # if the data is already downloaded then load it
  else {
    gset = suppressWarnings(GEOquery::getGEO(filename = gse_file)) }
  
  # update column names
  Biobase::fvarLabels(gset) <- make.names(Biobase::fvarLabels(gset))
  
  setwd("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R")
  
  # remove file holders
  rm(gse_file)
  
  return(gset)
}

############################################
# Annotations - also needs to be fixed
Gene_Annotation <- function(ID, gset){
  
  # get gene annotation information
  fset <- Biobase::fData(gset)
  
  # identify probes to remove (control probes)
  control = grep("control", names(fset), value = TRUE, ignore.case = TRUE)
  if(length(control) != 0) {
    fset <- fset[which(fset[, control]=='FALSE'),]
  }
  
  if(ID == 'GSE10588') {
    fset = fset[,c(1, 7, 6)]
    names(fset) = c('ID', 'Spot', 'Symbol')
  }
  
  if(ID == 'GSE24129') {
    # load libraries
    pacman::p_load(annotate, illuminaHumanv4.db)
    
    # extract accession ids
    fset$match <- mapply(unlist, strsplit(gsub("(N[A-Z]{1}_\\d+)|.", "\\1 ", fset$gene_assignment), "\\s+"))
    # convert to single list for annotation
    genes <- unlist(fset$match, recursive=FALSE)
    genes <- genes[!genes == '']
    
    # annotate accession ids
    annots<- AnnotationDbi::select(illuminaHumanv4.db, 
                                   keys=as.character(genes),
                                   columns=c('SYMBOL'), 
                                   keytype='ACCNUM')
    
    for (i in 1:length(fset$match)) {
      print(i)
      dat_list <- fset$match[i]
      annot <- c()
      
      for (j in 1:length(dat_list[[1]])) {
        if (any(grepl(dat_list[[1]][j], annots$ACCNUM))) {
          match = annots[which(annots$ACCNUM == dat_list[[1]][j]),]
          match = unique(match$SYMBOL)
          match <- match[! is.na(match)]
          annot <- c(annot, match)
          annot <- paste(unlist(annot), collapse=' // ')
        } else {
          annot = ""
        }
      }
      
      # append to data frame
      fset$symbol[i] <- annot
    }
  }
  
  if(ID == 'GSE73374') {
    # load libraries
    pacman::p_load(annotate, hugene20sttranscriptcluster.db)
    
    # extract accession ids
    fset$match <- mapply(unlist, strsplit(gsub("(N[A-Z]{1}_\\d+)|.", "\\1 ", fset$GB_ACC), "\\s+"))
    # convert to single list for annotation
    genes <- fset$GB_ACC
    genes <- genes[!genes == '']
    
    # annotate accession ids
    annots<- AnnotationDbi::select(hugene20sttranscriptcluster.db, 
                                   keys=as.character(genes),
                                   columns=c('SYMBOL'), 
                                   keytype='ACCNUM')
    
    for (i in 1:length(fset$GB_ACC)) {
      print(i)
      dat_list <- fset$GB_ACC[i]
      if (any(grepl(dat_list[[1]][j], annots$ACCNUM))) {
        match = annots[which(annots$ACCNUM == dat_list[[1]]),]
        match = unique(match$SYMBOL)
        match <- match[! is.na(match)]
        annot <- c(annot, match)
        annot <- paste(unlist(annot), collapse=' // ')
      } else {
        annot = ""
      }
      # append to data frame
      fset$symbol[i] <- annot
    }
  }

  # get column for gene symbol
  symbol = match(grep("symbol", names(fset), value = TRUE, ignore.case = TRUE),
                 names(fset))
  probe = match(grep("spot", names(fset), value = TRUE, ignore.case = TRUE),
                names(fset))
  probe = ifelse(sum(is.na(fset[,probe])) > 1, NA, probe)
  
  # spot ID
  if(!is.na(probe) && symbol > 1) {
    GSE_genes = fset[,c(1, probe,symbol)]
  }
  else {
    GSE_genes = fset[,c(1,1, symbol)]
  }
  
  # name data
  names(GSE_genes) = c('ID', 'Probe', 'Symbol')
  return(GSE_genes)
}

############################################
# Annotate probes to GSE-GPL Bioconductor Annotation libraries
Annotate <- function(genes) {
  
  # load libraries
  pacman::p_load(annotate, org.Hs.eg.db)
  
  # annotate genes
  genes_annot <- AnnotationDbi::select(org.Hs.eg.db, 
                                       keys=as.character(genes$Symbol),
                                       columns=c('ENTREZID',
                                                 'GENENAME',
                                                 'ALIAS'), 
                                       keytype='SYMBOL')
  
  
  # aggregate 1:many mappings into single list
  dt <- as.data.table(genes_annot)
  gene_data = dt[, lapply(.SD, function(x) toString(unique(x))), by = SYMBOL]
  
  # merge with probes
  genes_merge = merge(genes, gene_data, by.x = "Symbol", by.y = "SYMBOL", all = TRUE)
  genes_merge <- genes_merge[order(genes_merge$Probe),]
  genes_merge[] <- lapply(genes_merge, as.character)
  
  return(genes_merge)
}

############################################
# collapse genes with multiple probes
Probe_collapse <- function(genes, gset){
  # convert gset to expression data
  exp = exprs(gset)
  
  # subest expression data to only include non-control type probes
  exp <- exp[rownames(exp) %in% genes$ID, ]
  
  # transform matrix into data frame
  exp_m = cbind(genes$Probe, data.frame(exp), stringsAsFactors = FALSE)
  
  # collapse genes + label rownames as genes
  gene.max <- aggregate(. ~ exp_m[,1], data = exp_m, max)
  rownames(gene.max) <- gene.max[,1]
  
  # keep only samples
  keep <- grep('GSM', colnames(gene.max))
  gene.max = gene.max[,c(keep)]
  
}

############################################
# Pre-Processing
GEO_preprocess <- function(gset) {
  ## Log2 Transformation
  # determine if data has been log2 transformed
  ex <- gset
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  
  if (LogC) { 
    ex[ex <= 1.0] <- 1.0
    exprs_data <- log2(ex) 
  }
  else{
    exprs_data <- ex
  }
  
  ## Filtering
  # summarize average gene expression values
  avg <- rowMeans(exprs_data)
  sum = summary(avg)
  
  # filter out genes that "aren't expressed" - will appear red in the plot
  cutoff <- sum[2][[1]] # remove genes below first quantile
  exp_filt = expFilter(exprs_data, 
                       threshold = round(sum[2][[1]], 4), graph = FALSE)
  
  return(exp_filt)
}

############################################
# Access group assignments - harded coded for now, will need to be updated
Targets <- function(ID, gset){
  
  if(ID == 'GSE10588'|ID == 'GSE14722' | ID == 'GSE30186'| ID == "GSE44711") {
    groups = as.character(lapply(as.character(tolower(gset@phenoData@data$title)), function(x) ifelse(length(grep('pree|control', x)) == 1, 'Control', 'Case')))
  }
  
  if(ID == 'GSE73374') {
    groups = as.character(lapply(as.character(tolower(gset@phenoData@data$title)), function(x) ifelse(length(grep('norm', x)) == 1, 'Control', 'Case')))
  }
  
  # if(ID == 'GSE74341') {
  #   groups = as.character(tolower(gset@phenoData@data$title))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('eo', x)) == 1, 'EO', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('lo', x)) == 1, 'LO', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('at', x)) == 1, 'AT', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('preterm', x)) == 1, 'PT', x)))
  # }
  
  if(ID == 'GSE74341') {
    groups = as.character(tolower(gset@phenoData@data$title))
    groups = as.character(lapply(groups, function(x) ifelse(length(grep('eo|lo', x)) == 1, 'Case', x)))
    groups = as.character(lapply(groups, function(x) ifelse(length(grep('preterm|at', x)) == 1, 'Control', x)))
  }

  # if(ID =='GSE4707') {
  #   groups = as.character(tolower(gset@phenoData@data$title))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('eo', x)) == 1, 'EO', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('lo', x)) == 1, 'LO', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('_n', x)) == 1, 'Control', x)))
  # }
  
  if(ID =='GSE4707') {
    groups = as.character(tolower(gset@phenoData@data$title))
    groups = as.character(lapply(groups, function(x) ifelse(length(grep('eo|lo', x)) == 1, 'Case', x)))
    groups = as.character(lapply(groups, function(x) ifelse(length(grep('_n', x)) == 1, 'Control', x)))
  }
  
  # if(ID =='GSE24129') {
  #   groups = as.character(tolower(gset@phenoData@data$title))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('norm', x)) == 1, 'Control', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('fetal', x)) == 1, 'FGR', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('pre', x)) == 1, 'Case', x)))
  # }
  
  if(ID =='GSE24129') {
    groups = as.character(tolower(gset@phenoData@data$title))
    groups = as.character(lapply(groups, function(x) ifelse(length(grep('norm', x)) == 1, 'Control', x)))
    groups = as.character(lapply(groups, function(x) ifelse(length(grep('fetal|pre', x)) == 1, 'Case', x)))
  }
  
  if(ID == 'GSE25906') {
    groups = as.character(lapply(as.character(tolower(gset@phenoData@data$characteristics_ch1.2)), function(x) ifelse(length(grep('control', x)) == 1, 'Control', 'Case'))) 
  }
  
  if(ID == 'GSE60438') {
    groups = as.character(lapply(as.character(tolower(gset@phenoData@data$characteristics_ch1)), function(x) ifelse(length(grep('norm', x)) == 1, 'Control', 'Case'))) 
  }
  
  # if(ID == "GSE35574") {
  #   groups = as.character(tolower(gset@phenoData@data$characteristics_ch1.2))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('control', x)) == 1, 'Control', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('iugr', x)) == 1, 'IUGR', x)))
  #   groups = as.character(lapply(groups, function(x) ifelse(length(grep('pe', x)) == 1, 'Case', x)))
  # }

  if(ID == "GSE35574") {
    groups = as.character(tolower(gset@phenoData@data$characteristics_ch1.2))
    groups = as.character(lapply(groups, function(x) ifelse(length(grep('control', x)) == 1, 'Control', x)))
    groups = as.character(lapply(groups, function(x) ifelse(length(grep('iugr|pe', x)) == 1, 'Case', x)))
  }
  
  return(groups)
}

############################################
# Differential Expression Analysis
GEO_DEA <- function(targets, gset_pre, genes_annot, ID) {
  
  # create design matrix
  design <- stats::model.matrix(~ targets + 0, data.frame(gset_pre))
  colnames(design) <- gsub('targets', '', colnames(design))
  
  # make contrasts
  # if(ID == 'GSE4707') {
  #   
  #   contr.matrix <- limma::makeContrasts(
  #     EOvsC = EO-Control,
  #     LOvsC = LO-Control,
  #     levels = colnames(design))
  #   contr.matrix
  # }
  # 
  # if(ID == 'GSE74341') {
  #   
  #   contr.matrix <- limma::makeContrasts(
  #     EOvsPT = EO-PT,
  #     LOvsAT = LO-AT,
  #     levels = colnames(design))
  #   contr.matrix
  # }
  # 
  # if(ID == 'GSE35574') {
  #   
  #   contr.matrix <- limma::makeContrasts(
  #     IUGRvsC = IUGR-Control,
  #     CasevsControl = Case-Control,
  #     levels = colnames(design))
  #   contr.matrix
  # }
  # 
  # if(ID == 'GSE24129') {
  #   
  #   contr.matrix <- limma::makeContrasts(
  #     FGRvsC = FGR-Control,
  #     CasevsControl = Case-Control,
  #     levels = colnames(design))
  #   contr.matrix
  # }
  # 
  if(length(unique(targets)) == 2) {
    
    contr.matrix <- limma::makeContrasts(
      CasevsControl = Case-Control,
      levels = colnames(design))
    contr.matrix
  }
  
  # run limma
  fit <- limma::lmFit(data.frame(gset_pre), design)
  fit <- limma::contrasts.fit(fit, contrasts=contr.matrix)
  efit <- limma::eBayes(fit)
  
  # print number of up and down regulated genes between the samples
  top_dec_test <- decideTests(efit)
  print(summary(top_dec_test))
  
  # add gene information to results
  genes <- genes_annot[!duplicated(genes_annot$Probe),]
  genes$Probe <- sapply(genes$Probe, paste, collapse=",")
  rownames(genes) = genes$Probe
  genes <- genes[rownames(gset_pre),]
  
  # add gene information to efit
  efit$genes <- genes
  
  ## Process Results
  limma::write.fit(efit, top_dec_test, adjust='BH', file=paste('DE_Results/', ID ,'_DEgenes_studies.txt', sep = ''), sep = '\t')
  
}

############################################
# process all GSE ids in a list
gse_ids =c('GSE4707', 'GSE24129', 'GSE35574', 'GSE74341')

# GEO Main
# load libraries
pacman::p_load(Biobase)

# loop over ids
for (i in 1:length(gse_ids)) {
  # set ID
  ID = gse_ids[i]
  print(ID)
  
  # pull data from GEO
  gset = GEOData(ID)
  
  # get annotations 
  genes = Gene_Annotation(ID, gset)
  genes_annot = Annotate(genes)
  
  # remove control probes and collapse probes annotated to same gene
  exp <- Probe_collapse(genes, gset)
  dim(exp)[1]
  
  # Verify log-transformation & filter genes with exp below 1 quantile
  gset_pre <- GEO_preprocess(exp)
  dim(gset_pre)[1]
  
  # get targets
  targets <- Targets(ID, gset)
  print(table(targets))
  
  # differential expression analysis
  DE <- GEO_DEA(targets, gset_pre, genes_annot, ID)
}



############################################
############################################
## PROCESS RESULTS
# read in each of the DE results files - this includes the seprarate comparisons within each of the files
# set directory to location containing files
setwd("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_fullcomp")
# setwd("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_studies")
file_list <- list.files()

for (file in file_list){
  print(file)
  # if the merged dataset doesn't exist, create it
  if (!exists("pairwise_dataset")){
    pairwise_dataset <- read.table(file, header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("pairwise_dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
    pairwise_dataset<-rbind(pairwise_dataset, temp_dataset)
    rm(temp_dataset)
  }
}

# re-set directory to original project location
setwd("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R")

# loop over column and get first annotation
dataset <- pairwise_dataset
for (i in 1:length(dataset$symbol)) {
  print(i)
  gene = dataset$symbol[i]
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

# annotate symbols with Entrez ids
# load libraries
pacman::p_load(annotate, org.Hs.eg.db)

# annotate genes
genes_annot <- AnnotationDbi::select(org.Hs.eg.db, 
                                     keys=as.character(dataset$symbol),
                                     columns=c('ENTREZID',
                                               'ALIAS'), 
                                     keytype='SYMBOL')

# aggregate many:many mappings into single list
alias = aggregate(ALIAS ~ SYMBOL, 
                   data = genes_annot, 
                   paste, collapse = ' /// ')
names(alias) <- c("symbol", "alias")

# merge lists together
genes = merge(dataset, alias, by = 'symbol', all = TRUE)
genes = genes[,c(2,3,4,5,1,6)]

# autofill missing alias with gene symbols
genes$alias[is.na(genes$alias)] <- genes$symbol[is.na(genes$alias)]

# remove symbols with NA (removes 18,845 rows/14,948)
genes$symbol[genes$symbol == ""] <- NA
gene_data <- genes[complete.cases(genes$symbol), ]

# remove duplicates (removed 11,573 genes/9713 genes)
gene_data <- gene_data[!duplicated(gene_data[,c('alias', 'group')]),]

# re-label data and save
gene_data_studies <- gene_data
gene_data_pairwise <- gene_data

# Save the data objects
saveRDS(gene_data_pairwise, "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_fullcomp/pairwise_comparison_results.rds")
saveRDS(gene_data_studies, "~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_studies/studies_results.rds")

# Load the data objects
pairwise <- readRDS("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_fullcomp/pairwise_comparison_results.rds")
studies <- readRDS("~/Dropbox/Papers-Conferences-Projects/Hackathons/BioHackathon 2017/ignorenet/R/DE_Results/Final_studies/studies_results.rds")

# get counts by study
grp_symb = as.data.frame(table(grp = gene_data$group, symb = gene_data$alias))
tot_count <- aggregate(grp_symb$Freq, by=list(symbol=grp_symb$symb), FUN=sum)
table(tot_count$x)


