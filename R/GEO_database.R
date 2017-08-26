## Load needed libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite("GEOmetadb")
library(GEOmetadb)

## Download and unzip database
# temporarily set working directory to sub-folder
current_wd = getwd()
setwd(paste(current_wd, "/GEO_Database", sep=""))

# download and compress the data using SQLite - checks for file before automatically downloading
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()

# reset working directory
setwd(current_wd)

# view SQlite file 
file.info('GEO_Database/GEOmetadb.sqlite')

## Connect to database
con <- dbConnect(SQLite(),'GEO_Database/GEOmetadb.sqlite')

# view tables
geo_tables <- dbListTables(con)

# view fields associated with each table
dbListFields(con,'gse')

# get table schema
dbGetQuery(con,'PRAGMA TABLE_INFO(gse)')

# query a table
rs <- dbGetQuery(con,'SELECT * from gse
                      WHERE gse = "GSE14722"')

# close db connection
dbDisconnect(con)

# remove db file in order to download new one
dbDisconnect(file.remove('GEO_Database/GEOmetadb.sqlite'))


