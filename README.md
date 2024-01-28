# LION/web: LION enrichment analysis

App available online via [lipidontology.com](http://www.lipidontology.com)

<img src="https://raw.githubusercontent.com/martijnmolenaar/LION-web/master/OntologyApp/www/LIONicon%20enrich.png" alt="LION logo">

Please [cite](https://academic.oup.com/gigascience/article/8/6/giz061/5505544):
> **LION/web: a web-based ontology enrichment tool for lipidomic data analysis.**
> Martijn R Molenaar,  Aike Jeucken,  Tsjerk A Wassenaar,  Chris H A van de Lest, Jos F Brouwers,  J Bernd Helms. 
> *GigaScience, Volume 8, Issue 6, June 2019, giz061, https://doi.org/10.1093/gigascience/giz061*



## installation of R-packages required for LION/web
```R
if("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools",  repos = c(CRAN = "http://cran.rstudio.com"))}
library(devtools)

if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("topOnto" %in% rownames(installed.packages()) == FALSE) {install_github("martijnmolenaar/topOnto")}

## install RQRLite 0.11.4
if("RSQLite" %in% rownames(installed.packages()) == FALSE) {install.packages('https://cran.r-project.org/src/contrib/Archive/RSQLite/RSQLite_0.11.4.tar.gz', repos=NULL, type='source')}

##install topOnto.LION.db package:
if("topOnto.LION.db" %in% rownames(installed.packages()) == FALSE) {install_github("martijnmolenaar/topOnto.LION2.db/topOnto.LION.db")}

source("https://bioconductor.org/biocLite.R")
biocLite()

if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("ggthemes" %in% rownames(installed.packages()) == FALSE) {install.packages("ggthemes",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("httr" %in% rownames(installed.packages()) == FALSE) {install.packages("httr",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("igraph" %in% rownames(installed.packages()) == FALSE) {install.packages("igraph", repos = c(CRAN = "http://cran.rstudio.com"))}
if("shiny" %in% rownames(installed.packages()) == FALSE) {install.packages("shiny",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinyBS" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyBS",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinythemes" %in% rownames(installed.packages()) == FALSE) {install.packages("shinythemes",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinyTree" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyTree",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinyWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyWidgets",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("visNetwork" %in% rownames(installed.packages()) == FALSE) {install.packages("visNetwork",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("formattable" %in% rownames(installed.packages()) == FALSE) {install.packages("formattable",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("ggrepel" %in% rownames(installed.packages()) == FALSE) {install.packages("ggrepel",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("shinycssloaders" %in% rownames(installed.packages()) == FALSE) {install.packages("shinycssloaders",  repos = c(CRAN = "http://cran.rstudio.com"))}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {install.packages("jsonlite",  repos = c(CRAN = "http://cran.rstudio.com"))}

```

## Load all necessary libraries to check installation

```R
library(data.table)
library(ggplot2)
library(ggthemes)
library(httr)
library(igraph)
library(RSQLite)
library(shiny)
library(shinyBS)
library(shinythemes)
library(shinyTree)
library(shinyWidgets)
library(topOnto)
library('topOnto.LION.db')
library(visNetwork)
library(formattable)
library(jsonlite)
library(shinycssloaders)
library(ggrepel)
```

## run LION/web app
```R
runApp("OntologyApp")
```
