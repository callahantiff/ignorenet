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
# 
# ############################################
# ## Differential Expression Analyses
# # Analyze GSE Studies
# # http://wiki.bits.vib.be/index.php/Analyse_GEO2R_data_with_R_and_Bioconductor
# 
# ## USE CASE
# # gse_ids = c("GSE4707", "GSE35574", "GSE10588", "GSE25906", "GSE14722","GSE74341", "GSE44711", "GSE60438")
# 
# Get GEO data
GEOData <- function(ID) {
  # load libraries
  pacman::p_load(Biobase, GEOquery, limma)
  
  # gse file
  dirs <- grep("_out",list.dirs("ignorenet",recursive=FALSE),value=TRUE)
  gse_file <- list.files(dirs, ID, recursive=TRUE, full.names=TRUE)
  
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

  # remove file holders
  rm(gse_file)
  
  return(gset)
}


# Annotations - also needs to be fixed
Gene_Annotation <- function(ID, gset){
  
  # get gene annotation information
  fset <- Biobase::fData(gset)
  
  # identify probes to remove (control probes)
  control = grep("control", names(fset), value = TRUE, ignore.case = TRUE)
  if(length(control) != 0) {
    fset <- fset[which(fset[, control]=='FALSE'),]
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


# Access group assignments - harded coded for now, will need to be updated
Targets <- function(ID, gset){

  if(ID == 'GSE25906') {
    groups = as.character(lapply(as.character(gset@phenoData@data$characteristics_ch1.2), function(x) ifelse(length(grep('pree', x)) == 1, 'Case', 'Control')))
  }

  if(ID == 'GSE73374' | ID =='GSE10588' | ID =='GSE24129') {
    groups = as.character(lapply(as.character(tolower(gset@phenoData@data$title)), function(x) ifelse(length(grep('pre', x)) == 1, 'Case', 'Control')))
  }
  if(ID == 'GSE14722') {
    groups = as.character(lapply(as.character(tolower(gset@phenoData@data$title)), function(x) ifelse(length(grep('pree', x)) == 1, 'Case', 'Control')))
  }
  
  if(ID =='GSE4707') {
    groups = as.character(lapply(as.character(gset@phenoData@data$title), function(x) ifelse(length(grep('PE', x)) == 1, 'Case', 'Control')))
  }
  
  if(ID == 'GSE74341') {
    groups = as.character(lapply(as.character(tolower(gset@phenoData@data$title)), function(x) ifelse(length(grep('term', x)) == 1, 'Control', 'Case')))
  }

  if(ID == 'GSE44711') {
    groups = as.character(lapply(as.character(gset@phenoData@data$title), function(x) ifelse(length(grep('CON', x)) == 1, 'Control', 'Case')))
  }

  if(ID == 'GSE60438') {
    groups = as.character(lapply(as.character(gset@phenoData@data$source_name_ch1), function(x) ifelse(length(grep('pre-', x)) == 1, 'Case', 'Control')))
  }
  
  if(ID == 'GSE35574') {
    groups = as.character(lapply(as.character(gset@phenoData@data$characteristics_ch1.2), function(x) ifelse(length(grep('CON', x)) == 1, 'Control', 'Case')))
  }

  return(groups)
}

# Differential Expression Analysis
GEO_DEA <- function(targets, gset_pre, genes) {

  # create design matrix
  design <- stats::model.matrix(~ targets + 0, data.frame(gset_pre))
  colnames(design) <- gsub('targets', '', colnames(design))

  # make contrasts
  contr.matrix <- limma::makeContrasts(
    CasevsControl = Case-Control,
    levels = colnames(design))

  # run limma
  fit <- limma::lmFit(data.frame(gset_pre), design)
  fit <- limma::contrasts.fit(fit, contrasts=contr.matrix)
  efit <- limma::eBayes(fit)

  # get results
  topTab <- topTable(efit, adjust="BH", number=Inf)
  
  # add gene information to results
  genes <- genes[!duplicated(genes$Probe),]
  topTab <- data.frame(rownames(topTab), topTab, stringsAsFactors = FALSE)
  names(topTab) <- c('Probe', names(topTab[,-1]))
  results = merge(topTab, genes, by = 'Probe')
  
  return(results)
}


# GEO Main
GEO2R <- function(gse_ids){
  # load libraries
  pacman::p_load(Biobase)

  # list to write output
  geo_res <- list()

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
    
    # Verify log-transformation & filter genes with exp below 1 quantile
    gset_pre <- GEO_preprocess(exp)

    # get targets
    targets <- Targets(ID, gset)
    print(table(targets))

    # differential expression analysis
    DE <- GEO_DEA(targets, gset_pre, genes_annot)

    # append results with GSE
    geo_res[[i]] <- cbind(rep(ID, nrow(DE)), DE)
  }

  # aggregate data in list
  merged_res = do.call(rbind, geo_res)
  names(merged_res) <- c('gse', tolower(names(DE)))
  
  # remove ID column
  merged_res$id <- NULL

  # write outdata
  utils::write.table(merged_res,
                     'ignorenet/GEO_DE_results.txt',
                     row.names = FALSE,
                     col.names = TRUE,
                     quote = FALSE,
                     sep = '\t')
  
  return(merged_res)

}

# final_res <- GEO2R(gse_ids)

# # get a list of genes greater than or equal to 0.05
# filtered = merged_res[(merged_res$p.value) <= 0.05,]
# filtered_ordered <- filtered[order(filtered$p.value),]
# 
# # remove NAs
# filtered_clean <- filtered_ordered[complete.cases(filtered_ordered),]
# 
# jdk <- filtered_clean[filtered_clean$genename!= 'NA',]
# jdk <- jdk[c(1:2,9,3:5,8,6:7)]
# 
# jdk <- jdk[c(11,9)]
# 
# # jdk <- jdk[c("symbol", "genename")]
# utils::write.table(jdk,
#                    'ignorenet/jdk_results.tsv',
#                    row.names = FALSE,
#                    col.names = FALSE,
#                    quote = FALSE,
#                    sep = '\t')
# 




