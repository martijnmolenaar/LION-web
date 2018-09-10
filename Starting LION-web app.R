
##install topOnto package and RSQLite v 0.11.4 with:
library(devtools)
install_github("hxin/topOnto")

source("https://bioconductor.org/biocLite.R")
biocLite()

install.packages('https://cran.r-project.org/src/contrib/Archive/RSQLite/RSQLite_0.11.4.tar.gz', repos=NULL, type='source')

##install topOnto.LION.db package with:

install_github("martijnmolenaar/topOnto.LION2.db/topOnto.LION.db") 

## load packages
library(RSQLite)            ## for databases
library(topOnto)            ## for ontology enrichment analysis
library(topOnto.LION.db)    ## database with lipid ontologies

#available via CRAN:

library(igraph)             ## for network visualization, use vs 1.0.1
library(visNetwork)
library(GMD)

library(shiny)
library(shinythemes)
library(shinyTree)
library(shinyWidgets)

library(data.table)
library(ggplot2)
library(ggthemes)

library(visNetwork)
library(shinyTree)
library(shinyWidgets)


## LION/web

library(shiny)

runApp("OntologyApp")



