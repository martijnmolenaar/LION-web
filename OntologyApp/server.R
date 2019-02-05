options(stringsAsFactors = FALSE,shiny.sanitize.errors = TRUE)

 
## examples

PM <- read.csv(file = "data/PM-ER_pvalues.csv")
MT <- read.csv(file = "data/Mito-pvalues.csv")
ER <- read.csv(file = "data/ER_KLA-ER_CON_pvalues.csv")

backgroundlist_example <- paste(MT$lipids, "\n", collapse = '', sep = "")
sublist_example1 <- paste(PM$lipids[PM$pvalues < 0.20], "\n", collapse = '', sep = "")
sublist_example2 <- paste(MT$lipids[MT$pvalues < 0.20], "\n", collapse = '', sep = "")
sublist_example3 <- paste(ER$lipids[ER$pvalues < 0.20], "\n", collapse = '', sep = "")

pvalueExample1 <-  paste( paste(PM$lipids, "\t", PM$pvalues, sep = ""),    "\n",    collapse = '',    sep = ""  )
pvalueExample2 <-  paste( paste(MT$lipids, "\t", MT$pvalues, sep = ""),    "\n",    collapse = '',    sep = ""  )
pvalueExample3 <-  paste( paste(ER$lipids, "\t", ER$pvalues, sep = ""),    "\n",    collapse = '',    sep = ""  )

source(file  = "data/20190109 LION_tree_structure.R")

#### libraries

require(shiny)
require(visNetwork)
require(data.table)
require(GMD)
require(igraph)
require(ggplot2)
require(ggthemes)
library(shinyTree)
library(shinyWidgets)
library(shinyBS)
library(httr)

## loading lipid ontology data
require(RSQLite)
require(topOnto)
require('topOnto.LION.db')
topOnto::initONT('LION')



associationFile  <-  "data/20180614 LION_association.txt"


## LION term list for selection
LIONTermList <- read.table(file  = "data/20180209 LIONterms.txt", sep="\t")
LIONTermList$list <- as.list(LIONTermList[[2]])
names(LIONTermList$list) <- LIONTermList[[1]]
LIONTermList <- LIONTermList$list

## define functions

convertLipidNames <- function(name){
  
  options(stringsAsFactors = FALSE)
  
  sapply(name, function(input) {
    processed_input <- input
    if (processed_input==""){processed_input <- " "}
    
    #aliases 
    {
      processed_input <- gsub("DAG|diacylglycerol","DG",processed_input)
      processed_input <- gsub("TAG|triacylglycerol","TG",processed_input)
      
      processed_input <- gsub("GPA","PA",processed_input)
      processed_input <- gsub("GPEtn|phosphatidylethanolamine|Phosphatidylethanolamine","PE",processed_input)
      processed_input <- gsub("GPCho|phosphatidylcholine|Phosphatidylcholine","PC",processed_input)
      processed_input <- gsub("GPSer|phosphatidylserine|Phosphatidylserine","PS",processed_input)
      processed_input <- gsub("GPIns|phosphatidylinositol|Phosphatidylinositol","PI",processed_input)
      processed_input <- gsub("GPGro|phosphatidylglycerol|Phosphatidylglycerol","PG",processed_input)
      
      
      processed_input <- gsub("LBPA","BMP",processed_input)
      processed_input <- gsub("lyso|Lyso","L",processed_input)
      
      processed_input <- gsub("^FA","FFA",processed_input)
      processed_input <- gsub("ChE|CholE","CE",processed_input)
      
      processed_input <- gsub("\\Q(3'-sulfo)Galβ-Cer\\E|\\Q(3'-sulfo)GalCer\\E|Sulfogalactosyl ceramide[ ]*|[Ss]+(Gal|Hex)Cer|Sulfatide",
                              "SHexCer",processed_input)
      
      
      processed_input <- gsub("LacCer","Hex2Cer",processed_input)
      processed_input <- gsub("Gl[u]*c[β-]*Cer|Gal[β-]*Cer|GluCer","HexCer",processed_input)
    }    ### end aliases
    
    ### special cases:  'DG(0:0/18:1/16:0)' >> DG(18:1/16:0)
    if(grepl('DG',processed_input) & grepl('\\D0:0',processed_input)){
      processed_input <- gsub('/0:0',"",processed_input)
      processed_input <- gsub('0:0/',"",processed_input)
      
    }
    
    ### end special cases
    if(grepl("\\d+:\\d+;\\d+",processed_input) & !grepl("\\d+:\\d+;\\d+:\\d+",processed_input)){    ### 32:1;1 format >> 32:1
      processed_input <- gsub(";\\d+","",processed_input)
    }
    
    ###  >> remove (FA 16:0), not intensively tested
    processed_input <- gsub("\\(FA( )*", "(",processed_input)   
    
    ### 'CL(32:0)(34:1)' or 'TAG 32:0(FA 16:0)' >> sum FAs between (), not intensively tested
    FAs <- unlist(regmatches(processed_input, gregexpr("\\d+:\\d+", processed_input)))
    if(length(FAs)>0){
      FAs <- FAs[sapply(regmatches(FAs, gregexpr("\\d+", FAs)), function(i){  as.numeric(i)[1] > as.numeric(i)[2]  })]  ## C should be higher than DB  
    }
    if((grepl("CL", processed_input ) & (length(FAs) == 2 | length(FAs) == 3)) |
       (grepl("TAG|TG", processed_input ) & (length(FAs) == 2)))     {
      parts <- unlist(regmatches(processed_input, gregexpr("[\\(]*\\d+:\\d+[\\)]*", processed_input)) )
      prefix <- unlist(strsplit(processed_input, "[\\(]*\\d+:\\d+[\\)]*"))[1]
      suffix <- unlist(strsplit(processed_input, "[\\(]*\\d+:\\d+[\\)]*"))[2]
      processed_input <- paste(prefix, "(",
                               sum(sapply(regmatches(
                                 parts, gregexpr("\\d+", parts)
                               ), function(i) {
                                 as.numeric(i)[1]
                               })), ## number of Cs
                               ":",
                               sum(sapply(regmatches(
                                 parts, gregexpr("\\d+", parts)
                               ), function(i) {
                                 as.numeric(i)[2]
                               })), ## number of DBs
                               ")", suffix, sep = "")
    }
    ### end 'CL(32:0)(34:1)'
    
    
    if(length(FAs)==2 & grepl("^SM",    processed_input)){                ### voor 'SM 18/16:0)'
      if(!grepl("[dt]\\d+",processed_input)){                             ### double check: no 'd'18...?
        processed_input <- gsub("SM 18","SM d18",processed_input)         ### add 'd'
        processed_input <- gsub("SM\\(18","SM(d18",processed_input)
      }
      
    }
    
    processed_input <- gsub("_\\d*\\.*\\d*/\\d*\\.*\\d*", "",processed_input)   ## remove amu/RT info
    
    processed_input <- gsub(";\\d$","",processed_input)   ### in case of SHexCer 30:2;1
    
    if (grepl("DHHexosylceramide|Hexosylceramide|DHSphingomyelin|Sphingomyelin|DHCeramide|Ceramide", processed_input)){   ## for notation 'DHSphingomyelin 16'
      if (grepl("C\\d+",processed_input)){   ## replace 'C20:3 to 20:3
        C_containing_notation <- unlist(regmatches(processed_input,  gregexpr("C\\d+", processed_input)))
        processed_input <- sub(C_containing_notation,  gsub('C',"",C_containing_notation), processed_input)
      }
      
      if (!grepl("\\d+:\\d+",processed_input)){   ## replace 'C16 to 16:0
        C_number <- regmatches(processed_input, regexpr("\\d+", processed_input)) 
        processed_input <- sub(C_number, paste(C_number,":0",sep=""), processed_input)
      }
      
      if (grepl("^DH",processed_input)){   ## DH >> d18:0, otherwise >> d18:1
        processed_input <- gsub("DH","",processed_input)
        FAs <- regmatches(processed_input, regexpr("\\d+.*\\d*", processed_input)) 
        processed_input <-  sub(FAs, paste("d18:0/",FAs,sep=""), processed_input)
      } else {
        #processed_input <- gsub("DH","",processed_input)
        FAs <- regmatches(processed_input, regexpr("\\d+.*\\d*", processed_input)) 
        processed_input <-  sub(FAs, paste("d18:1/",FAs,sep=""), processed_input)
      }
      
      processed_input <- gsub("Hexosylceramide","HexCer",processed_input)
      processed_input <- gsub("Sphingomyelin","SM",processed_input)
      processed_input <- gsub("Ceramide","Cer",processed_input)
    }
    
    processed_input <-
      gsub("_|;", "/", processed_input)                                ## _ and ; into /
    
    for(substitute in unlist(regmatches(processed_input, gregexpr("\\d+-\\d+", processed_input)) )){
      processed_input <- gsub(substitute, gsub("-","/",substitute), processed_input) 
    }   ## - into / if it's between numbers
    
    
    processed_input <-
      gsub("\\((\\d+[EZ]{1},*)+\\)", "", processed_input)   ## removing (5Z,8Z,11Z,14Z)-like info
    
    if (grepl("\\d+:\\d+[pe]{1}", processed_input)){             ### ether lipids: in case of 18:1p in stead of P-18:1
      
      ether_FAs <- unlist(regmatches(processed_input, gregexpr("\\d+:\\d+[pe]{1}", processed_input)))  ## extract 18:1
      
      for(i in 1:length(ether_FAs)){
        processed_input <-
          gsub(ether_FAs[i],              ## replace 18:1p for P-18:1
               paste(
                 ifelse(grepl("p", ether_FAs[i]), "P-", "O-"),
                 regmatches(ether_FAs[i], regexpr("\\d+:\\d+", ether_FAs[i])),
                 sep = ""
               ),
               processed_input)
      }
    }
    
    if (!grepl("\\(", processed_input)) {
      ## if they don't contain /
      ## process PC 18:1/20:4 format to PC(..)
      
      middlePart <-         ### ..O-18:1/16:0
        regmatches(processed_input,
                   regexpr("[dtAOP-]*(\\d+:\\d+/*)+", processed_input))
      otherParts <-         ## SM
        unlist(strsplit(processed_input, "[dtAOP-]*(\\d+:\\d+/*)+"))
      otherParts <- gsub(" ", "", otherParts)
      text <- ""
      
      for (i in 1:max(c(length(otherParts),   length(middlePart)))) {
        ## reconstruct
        if (!is.na(otherParts[i])) {
          text <- paste(text, otherParts[i], sep = "")
        }
        
        if (!is.na(middlePart[i])) {
          text <- paste(text, "(", middlePart[i], ")", sep = "")
        }
      }
      processed_input <- text
    }
    
    
    if(grepl("-DiHex",processed_input)){          # in case of d18:0/12:0-DiHex
      processed_input <- gsub("-DiHex","",processed_input)
      processed_input <- paste("Hex2Cer",processed_input,sep="")
    }
    
    if(grepl("-MonoHex",processed_input)){          # in case of d18:2/26:0-MonoHex
      processed_input <- gsub("-MonoHex","",processed_input)
      processed_input <- paste("HexCer",processed_input,sep="")
    }
    
    if(grepl("\\D+ [\\(|\\d+]", processed_input)){        ## in case of 'CL (70:2)'
      processed_input <- gsub("\\D+ [\\(|\\d+]", 
                              gsub(" ","", regmatches(processed_input, regexpr("\\D+ [\\(|\\d+]", processed_input)) ),
                              processed_input)
    }
    
    if (grepl("^LP.*\\((\\d+:\\d+/*)+", processed_input)) {
      ## if LPC format >> format to PC(xx:0:0)
      
      
      middlePart <-         ### ..18:1 + /0:0
        paste(regmatches(
          processed_input,
          #regexpr("(\\d+:\\d+/*)+", processed_input)
          regexpr("([OP-]*\\d+:\\d+[ep/]*)+", processed_input)
        ), sep = "")
      
      if(identical(middlePart,character(0))){middlePart <- ""}
      ### if no 0:0 yet, add this
      if(!grepl("\\D0:0|^0:0", middlePart)){middlePart <- paste(middlePart,"/0:0",sep="")}
      
      
      otherParts <-         ## LPC
        unlist(strsplit(processed_input, "([OP-]*\\d+:\\d+[ep/]*)+"))
      otherParts <- gsub("L", "", otherParts)
      if(identical(otherParts,character(0))){otherParts <- ""}
      
      text <- ""
      
      for (i in 1:max(c(length(otherParts),   length(middlePart)))) {
        ## reconstruct
        if (!is.na(otherParts[i])) {
          text <- paste(text, otherParts[i], sep = "")
        }
        
        if (!is.na(middlePart[i])) {
          text <- paste(text, middlePart[i], sep = "")
        }
      }
      processed_input <- text
      
      
    }
    

    if(grepl('^P\\D\\(A-*',processed_input) ){      ########## ether lipids 20180605
      processed_input <- gsub('\\(A-*',"(O-",processed_input)    ##gsub('A-*',"(O-",processed_input)
    }
    if(grepl('^P\\D\\([OP]\\d+',processed_input) ){      ########## ether lipids 20180219
      processed_input <- gsub('\\(O',"(O-",processed_input)
      processed_input <- gsub('\\(P',"(P-",processed_input)
    }
    
    if(grepl("[[:upper:]]{1}[[:lower:]]{2}.+(sterol|enone)", processed_input)){         ### in case of Lanosterol / Cholestenone
      UpperLower <- regmatches(processed_input, regexpr("[[:upper:]]{1}[[:lower:]]{2}", processed_input)) 
      processed_input <- gsub("[[:upper:]]{1}[[:lower:]]{2}", tolower(UpperLower),processed_input )
    }
    
    
    ## reorder FAs
    if (
      (length(unlist(regmatches(processed_input, gregexpr("\\d+:\\d+",        processed_input)))) > 1) &                         ## when more than 1 FA)
      (length(unlist(regmatches(processed_input, gregexpr("[dtm](\\d+:\\d+)|[PO]-(\\d+:\\d+)", processed_input)))) < 1) &        ## unless it's SM(d18:1...) or P-18:1
      (!grepl("0:0", processed_input))                           ## unless it's containing a 0:0 FA (=lyso, should be on second position))
    ){                                                               
      FAs <- t(as.data.frame(strsplit(unlist(
        regmatches(
          processed_input,
          gregexpr("\\d+:\\d+", processed_input)
        )
      ), ":")))
      FAs <- as.data.frame(sapply(as.data.frame(FAs), as.numeric))
      FAs <- FAs[order(FAs$V1, FAs$V2), ]
      FAs <- apply(FAs, 1, function(row) {
        paste(row[1], ":", row[2], sep = "")
      })
      
      parts <- unlist(strsplit(processed_input, "\\d+:\\d+"))
      text <- ""
      
      for (i in 1:length(parts)) {
        if (!is.na(parts[i])) {
          text <- paste(text, parts[i], sep = "")
        }
        
        if (!is.na(FAs[i])) {
          text <- paste(text, FAs[i], sep = "")
        }
      }
      text
    } else {
      processed_input
    }
    
  })
  
  
  
}      ## 20190125

sendEmail <- function(subject = "LION/web usage", mail_message = "empty", from = "system"){
  url <- "emailserver"
  api_key <- "xxxxxxx"
  the_body <-
    list(
      from="webmaster",
      to="email@email.com",
      subject=subject,
      text=paste(mail_message,"\n\n=====================\nfrom: ",from,"\n=====================\nLION/web: www.lipidontology.com", sep= "")
    )
  req <- httr::POST(url, httr::authenticate("api", api_key),  encode = "form",  body = the_body)
  httr::stop_for_status(req)
  TRUE
}   ## mail function

associatedTerms <- function(lipid, ontologyObject, all = FALSE, reformat = FALSE) {    ## function: which LION-terms are associated with a lipid?
  if (!reformat) {
    associations <- sapply(lipid, function(lipid_i) {
      terms <-
        ontologyObject@termName[names(ontologyObject@termName) %in% names(genesInTerm(ontologyObject))[sapply(genesInTerm(ontologyObject), function(term) {
          any(term == lipid_i)
        })]]   ## find terms
      
      if (!all) {
        terms <- terms[grepl("LION", names(terms))]
      }     ## is all == TRUE, than remove CAT and all  terms
      return(terms)
    })
    return(associations)
    
  } else if (reformat) {
    
    terms <- names(genesInTerm(ontologyObject))
    if (!all) {
      terms <- terms[grepl("LION", terms)]
    }     ## is all == TRUE, than remove CAT and all  terms
    
    lipid_term_table <- t(sapply(lipid, function(lipid_i) {
      terms_per_lipid <-
        names(ontologyObject@termName[names(ontologyObject@termName) %in% names(genesInTerm(ontologyObject))[sapply(genesInTerm(ontologyObject), function(term) {
          any(term == lipid_i)
        })]])
      return(ifelse(terms %in% terms_per_lipid, "x", ""))
    }))
    
    names <- ontologyObject@termName[match(terms, names(ontologyObject@termName))]
    
    lipid_term_table <- rbind(terms, names, lipid_term_table)
    lipid_term_table <- cbind(rownames(lipid_term_table), lipid_term_table)
    rownames(lipid_term_table) <- NULL
    colnames(lipid_term_table) <- NULL
    lipid_term_table <- as.data.frame(lipid_term_table)
    lipid_term_table[[1]][1:2] <- c("LION-term", "LION-name")
    return(lipid_term_table)
  }
  
} ## which LION-terms are associated with a lipid?



## read associations
lipidID2TERM <- readMappings(file = associationFile)   ## topOnto function, but removes spaces

# Define server
function(input, output, session) {
  
  
  ### hide tabs at start-up
  hideTab(inputId = "tabs", target = "LION input")
  hideTab(inputId = "tabs", target = "LION enrichment table")
  hideTab(inputId = "tabs", target = "LION enrichment graph")
  hideTab(inputId = "tabs", target = "LION network view")
  
  
  ## pre-processing with CSVs:
  
  input_data <- reactive({
    
    req(input$file1)
    
    
    df <- try(read.csv(input$file1$datapath,
                   header = FALSE,
                   sep = ",",
                   quote = ""))
    
    
    ## error handling
    
    errors <- NULL
    if (!(is.data.frame(df))) {
      df <- data.frame(errors = df[1], column = 0)
      errors <- c(errors, "ERROR: File type not supported")
    } else {
      if (dim(df)[2]  == 1) {
        errors <- c(errors, "ERROR: No commas found to seperate columns")
      }
      if (dim(df)[2]  > 1 & dim(df)[2]  < 3) {
        errors <- c(errors, "ERROR: Unkown error")
      }
      if (sum(is.na(df[, -1]))   > 0) {
        errors <- c(errors, "ERROR: There are missing values")
      }
      if (sum( df[-c(1,2),-1] == "" ) > 0) {
        errors <- c(errors, "ERROR: Dataset contains empty cells")
      }
      if (all(table(as.character(df[1, -1])) < 2)) {
        errors <- c(errors, "ERROR: Some or all conditions are n < 2")
      }
      if (any(duplicated(df[2, -1]))) {
        errors <-
          c(errors, "ERROR: One or more samples names are not unique")
      }
    }
    errors <- paste(errors, collapse = "<br>")
    
    ## end error handling
    
    
    if(!(is.null(errors))){
      input_data <- list(df = df,
         matrix = sapply(df[-c(1,2),-1], as.numeric),
         conditions = unique(as.character(df[1,-1])),
         samples = as.character(df[2,-1]),
         meta = df[1:2,-1],
         IDs = df[-c(1,2),1],
         errors = errors)
    } else {
      input_data <-list(df = NULL,
           matrix = NULL,
           conditions = NULL,
           samples = NULL,
           meta = NULL,
           IDs = NULL,
           errors = errors)
    }
    
    return(input_data)
    
  })
  
  output$selectLocalStatisticUI <- renderUI({
    input_data <- input_data()
    conditions <- input_data$conditions
    
    if (input_data$errors == "") {
      
      wellPanel(
        popify(placement = "bottom", title = "Info",
        radioButtons("local_statistics", label = h4("Select local statistics to rank input identifiers"),
                     choices = list("one-tailed T-test (2 conditions)" = 1, "2-LOG[fold change] (2 conditions)" = 2, "one-way ANOVA F-test (>2 conditions)" = 3), 
                     selected = 1), content = 'Select a local statistic to rank the metabolites based on the provided dataset. LION-terms associated with lipids higher ranked than expected by chance will be reported as enriched.'),
      
        br(),
        uiOutput("conditionsUI")
      )
    } else {     ### if there are errors...
      HTML(paste('<font color="red">',
                 input_data$errors,
                 "</font>",sep=""))
    }
    
  })
  
  output$conditionsUI <- renderUI({
    input_data <- input_data()
    conditions <- input_data$conditions
    
    if (input$local_statistics == 1){        ## t-test pvalues 
      
      
      wellPanel(
        
        h5(tags$i("a T-test p-value will be calculated for every input identifier")),
        
        fluidRow(
          column(
            width = 6,
            selectInput(
              "conditionA",
              h4("condition of interest"),
              conditions,
              selected = conditions[1]
            )
          ),
          column(
            width = 6,
            selectInput(
              "conditionB",
              h4("control condition"),
              conditions,
              selected = conditions[2]
            )
          )
        ),
        actionButton(
          inputId = "calculatePreprocessing",
          label = "  Calculate local statistics",
          width = NULL,
          icon = icon("line-chart", lib = "font-awesome")
        ),
        br(),
        br(),
        plotOutput("p_value_plot"),
        br(),
        uiOutput("conditionsUI_pB")
      )
    } else if(input$local_statistics == 2) {   ### 2log fold change
      
      wellPanel(
        
        h5(tags$i("a 2-LOG[fold change] value will be calculated for every input identifier")),
        
        fluidRow(
          column(
            width = 6,
            selectInput(
              "conditionA",
              h4("condition of interest"),
              conditions,
              selected = conditions[1]
            )
          ),
          column(
            width = 6,
            selectInput(
              "conditionB",
              h4("control condition"),
              conditions,
              selected = conditions[2]
            )
          )
        ),
        actionButton(
          inputId = "calculatePreprocessing",
          label = "  Calculate local statistics",
          width = NULL,
          icon = icon("line-chart", lib = "font-awesome")
        ),
        br(),
        br(),
        plotOutput("p_value_plot"),
        br(),
        uiOutput("conditionsUI_pB")
      )
      
    } else {     ### F-test
      
      wellPanel(
        
       h5(tags$i("an F-test p-value will be calculated for every input identifier")),
        
       selectizeInput(inputId = "selectedConditions",  h4("conditions of interest"), 
                       choices = input_data()$conditions, selected = input_data()$conditions,
                       multiple = TRUE,options = list(plugins= list('remove_button'))
        ),
       
        actionButton(
          inputId = "calculatePreprocessing",
          label = "  Calculate local statistics",
          width = NULL,
          icon = icon("line-chart", lib = "font-awesome")
        ),
        br(),
        br(),
        htmlOutput("FtestError"),
        br(),
        plotOutput("p_value_plot"),
        br(),
        uiOutput("conditionsUI_pB")
      )
     
    }
  })
  
  output$FtestError <- renderText({ 
    message <- ifelse(length(input$selectedConditions)<2, "ERROR: Number of conditions should be 2 or higher","")
    HTML(paste('<font color="red">',
               message,
               "</font>",sep=""))
  })
  
  
  output$conditionsUI_pB <- renderUI({
    data <- pre_processed_data()[[1]]   ### show when 'pre_processed_data()' is constructed
    if(!(any(is.na(data$pValues)))){
      actionButton(inputId = "submitPreprocessing",
                   label = "  Use values as local statistics", 
                   width = NULL,
                   icon = icon("share-square-o", lib = "font-awesome"))
    }
    
    
  })

  pre_processed_data <-     eventReactive(input$calculatePreprocessing, {
    
    
    if(any(input$local_statistics %in% c(1,2))){
    setA <- input_data()$matrix[,which(input_data()$meta[1,] == input$conditionA)]
    setB <- input_data()$matrix[,which(input_data()$meta[1,] == input$conditionB)]
    }
    
    if(input$local_statistics == 1){   ### when t-test p-values are selected
      
      pValues <- apply(cbind(setA,setB), 1, function(row){
        t.test(x=row[1:dim(setA)[2]],
               y=row[(dim(setA)[2]+1):(dim(setA)[2]+dim(setB)[2])],
               alternative = "greater")$p.value 
        
      })
      outputList <- list(data.frame(IDs = input_data()$IDs,
                                    pValues = pValues))
      
    }
      
    if(input$local_statistics == 2){   ### when 2-log FC values are selected
      
      FC_values <- apply(cbind(setA,setB), 1, function(row){
        x <- mean(row[1:dim(setA)[2]], na.rm = TRUE)
        y <- mean(row[(dim(setA)[2]+1):(dim(setA)[2]+dim(setB)[2])], na.rm = TRUE)
        log(x = x/y, base = 2)
      
      })
      outputList <- list(data.frame(IDs = input_data()$IDs,
                                    pValues = FC_values))
      
    }
    
    if(input$local_statistics == 3){   ### when 2-log FC values are selected
      
      input_data <- input_data()
      
      if (length(input$selectedConditions) > 1) {      ## minimum of 2 conditions needed
        pValues <- sapply(1:length(input_data$IDs), function(lipid_nr) {
          df <- melt(
            data = sapply(input$selectedConditions, function(condition) {
              input_data$matrix[lipid_nr,  unlist(input_data$meta[1,]) == condition]
            }),
            value.name = "signal",
            variable.name = "condition"
          )
          names(df) <- c("l", "condition", "signal")
          
          model.df <- lm(signal ~ condition, data = df)
          anova(model.df)['condition', 'Pr(>F)']
        })
      } else {pValues = NA}
    
          
       
      outputList <- list(data.frame(IDs = input_data()$IDs,
                                    pValues = pValues))
      
    }
    
    names(outputList) <- paste("set",isolate(input$submit),sep="")
    outputList
    
  })
  
  output$p_value_plot <- renderPlot({
    
    df <- pre_processed_data()[[1]]
    df <- df[,2][order(df[,2])]
   
    ylab <- subset(data.frame(selection = 1:3, ylab = c('t-test p-value','2-LOG[fold change]','F-test p-value (log scale)')), 
           selection == isolate(input$local_statistics))$ylab
    
    if(ylab == 'F-test p-value (log scale)'){
      if (length(input$selectedConditions) > 1) { 
      qplot(y=df, main = "Distribution local statistics", 
            xlab = "metabolites", 
            ylab = ylab,
            log = "y"
            
      )} else {
        ggplot(data.frame(x = 1:length(df), y=NA), aes(x = x, y = y))+
          labs(x = "metabolites", y = ylab, title = "Distribution local statistics")+
          geom_blank()
      }
      
    } else {
      qplot(y=df, main = "Distribution local statistics", 
            xlab = "metabolites", 
            ylab = ylab
      )
    }
    
   
    
  })
  
  ## end of pre-processing
  

  ## examples by p-value
  
  observeEvent(input$submitPreprocessing, {
    pre_processed_data <- pre_processed_data()[[1]]
    pre_processed_data <- paste( paste(pre_processed_data$IDs, "\t", pre_processed_data$pValues, sep = ""),    "\n",    collapse = '',    sep = ""  )
    
    updateTextAreaInput(session, "listwPvalues",  value = pre_processed_data)
    ### update dependent on statistics choice
    if(isolate(input$local_statistics) == 2){   ### 2 == log fold change 
      updateRadioButtons(session, "ranking", selected = "descending" )
    } else {     ## in case of p-values
      updateRadioButtons(session, "ranking", selected = "ascending" )
    }
    
    ###
    updateTabsetPanel(session, "sub_method",
                      selected = "(ii) analysis")
  })
  
  output$examplePre1 <- downloadHandler(
    filename <- function() {
      paste("lipidomics set - adapted from Andreyev AY et al 2010","csv",sep=".")
    },
    
    content <- function(file) {
      file.copy("data/lipidomics set - adapted from Andreyev AY et al 2010.csv", file)
    },
    contentType = "text/csv"
    )
  output$examplePre2 <- downloadHandler(
    filename <- function() {
      paste("lipidomics set FA incorporation","csv",sep=".")
    },
    
    content <- function(file) {
      file.copy("data/lipidomics set FA incorporation.csv", file)
    },
    contentType = "text/csv"
  )
  output$examplePre3 <- downloadHandler(
    filename <- function() {
      paste("lipidomics set membrane fluidity after AA incubation","csv",sep=".")
    },
    
    content <- function(file) {
      file.copy("data/lipidomics set membrane fluidity after AA incubation.csv", file)
    },
    contentType = "text/csv"
  )
  
  
  observeEvent(input$exampleB1, {
    updateTextAreaInput(session, "listwPvalues",  value = pvalueExample1)
  })
  observeEvent(input$exampleB2, {
    updateTextAreaInput(session, "listwPvalues",  value = pvalueExample2)
  })
  observeEvent(input$exampleB3, {
    updateTextAreaInput(session, "listwPvalues",  value = pvalueExample3)
  })
  observeEvent(input$plot_click, {           ### display term info in console
    termInfo <- dataBlock()$to_ggplot_display$data
    
  
  })
  
  ## examples by sublist
  observeEvent(input$exampleA1, {
    updateTextAreaInput(session, "sublist",  value = sublist_example1)
    updateTextAreaInput(session, "background",  value = backgroundlist_example)
  })
  observeEvent(input$exampleA2, {
    updateTextAreaInput(session, "sublist",  value = sublist_example2)
    updateTextAreaInput(session, "background",  value = backgroundlist_example)
  })
  observeEvent(input$exampleA3, {
    updateTextAreaInput(session, "sublist",  value = sublist_example3)
    updateTextAreaInput(session, "background",  value = backgroundlist_example)
  })
  
  observeEvent(input$submitB, {
    ## first unhide tabs
    showTab(inputId = "tabs", target = "LION input")
    showTab(inputId = "tabs", target = "LION enrichment table")
    showTab(inputId = "tabs", target = "LION enrichment graph")
    showTab(inputId = "tabs", target = "LION network view")  
    
    dataBlock()
    
    ## goto LION input
    updateTabsetPanel(session, "tabs",
                      selected = "LION input")
  })
  observeEvent(input$submitA, {
    ## first unhide tabs
    showTab(inputId = "tabs", target = "LION input")
    showTab(inputId = "tabs", target = "LION enrichment table")
    showTab(inputId = "tabs", target = "LION enrichment graph")
    showTab(inputId = "tabs", target = "LION network view")  
    
    ## goto LION input
    updateTabsetPanel(session, "tabs",
                      selected = "LION input")
  })

  
  
  ### enrichment analysis
  dataBlock <- reactive({
    
    errorhandling <- NULL
    
    input$submitB
    input$submitA
    
    input_listwPvalues <- isolate(input$listwPvalues)
    
    ### some lipids contain comma's, so only the last comma should be regarded as a separator
    numericPattern <- unlist(regmatches(input_listwPvalues, gregexpr(", *(\\d|\\.)+(\n|$)", input_listwPvalues)) )   ## extract patterns to replace last comma
    
    for(p in numericPattern){            ### replace last comma of line to \t; this will be the seperator
      input_listwPvalues <- gsub(p, gsub(",", "\t", p), input_listwPvalues)
    }
    
    
    ####  is input appropriate??
    
    if (isolate(input$sub_method) == "(ii) analysis") {
      if(!(grepl("\t", input_listwPvalues) & 
           grepl("\n", input_listwPvalues) &
           grepl("\\D",input_listwPvalues) &
           grepl("\\d",input_listwPvalues) )){           ### definition of wrong input
        
        hideTab(inputId = "tabs", target = "LION enrichment table")
        hideTab(inputId = "tabs", target = "LION enrichment graph")
        hideTab(inputId = "tabs", target = "LION network view")
        
        return(NULL)
      }
    }
    
    if (isolate(input$method) == "Target-list mode") {        
      if( isolate(input$sublist)==""  | isolate(input$background) == ""){   ### definition of wrong input
        hideTab(inputId = "tabs", target = "LION enrichment table")
        hideTab(inputId = "tabs", target = "LION enrichment graph")
        hideTab(inputId = "tabs", target = "LION network view")
        
        return(NULL)
      }
    }
    ####   
    
  
    
    if(isolate(input$ranking)=="ascending"){direction <- 1} else {direction <- -1}
    
    ###  table 1
    
    withProgress(message = 'progress:', value = 0, {
      incProgress(.1, detail = paste("mapping input data"))
      if (isolate(input$method) == "Target-list mode") {
      if ((is.null(isolate(input$sublist)) ||
           isolate(input$sublist) == "") == FALSE) {
        lipidExistance <-
          data.frame(input = unique(c(unlist(strsplit(isolate(input$sublist), "\n")),  unlist(strsplit(isolate(input$background), "\n"))))  )
        lipidExistance$'simplified input' <- convertLipidNames(lipidExistance$input)
        
        
      } else {
        lipidExistance <- data.frame("input")
      }
      
    }
    
      if (isolate(input$sub_method) == "(ii) analysis") {
      if ((is.null(input_listwPvalues) ||
           input_listwPvalues == "") == FALSE) {

        pValueList <- list()
        
        pValueList$input <-
          transpose(strsplit(unlist(strsplit(
            input_listwPvalues, "\n"
          )), "[\t\\|]"))[[1]]
        
        pValueList$pvalues <-
          transpose(strsplit(unlist(strsplit(
            input_listwPvalues, "\n"
          )), "[\t\\|]"))[[2]]
        
        if(any(is.null(pValueList$input),is.null(pValueList$pvalues))){
          errorhandling <- "unexpected input"
          pValueList$input <- NULL
          pValueList$pvalues <- NULL
        }
  
        pValueList <- as.data.frame(pValueList)
        lipidExistance <- data.frame(input = pValueList$input)
        lipidExistance$'simplified input' <- convertLipidNames(lipidExistance$input)
        
        
      } else {
        lipidExistance <- data.frame("input")
      }
    }
      if (length(lipidExistance$input) < 1) {
      
    } else {
      
      lipidExistance$'LION ID' <-  unlist(lipidID2TERM)[match(lipidExistance$'simplified input',names(lipidID2TERM))]
      lipidExistance$'LION ID'[is.na(lipidExistance$'LION ID')] <- "not found"
      
      
      lipidExistance$'LION name' <- sapply(lipidExistance$'LION ID', function(ID){     ### look for LIONname
        output <- names(lipidID2TERM)[lipidID2TERM == ID]
        if(identical(output, character(0))){output <- "not found"}
        output[!grepl("^SLM:|^LM..\\d+",output )]                  ## don't return names like LMGP04010883 or SLM:000042385
      })
      
      lipidExistance$'simplified input' <- NULL        ## remove this column, not interesting for user
      
      ## lipidExistance <- lipidExistance[order(lipidExistance$'LION ID'),]   ## reorder based on matching
      lipidExistance
    }
    
    })
    
    
    if (!any(dim(lipidExistance) == c(1, 1))) {
      ### perform when there is input (dim 1 by 1 >> no input)
      
      withProgress(message = 'progress:', value = 0, {
        # Generate an HTML table view of the ontology data
        
        incProgress(.5, detail = paste("submitting data"))
        
        if (isolate(input$method) == "Target-list mode") {
          if ((is.null(isolate(input$sublist)) ||
               isolate(input$sublist) == "") == FALSE) {
            lipidIDs <-
              list(sublist = convertLipidNames(  unlist(strsplit(isolate(input$sublist), "\n"))  ),
                   backgroundlist = convertLipidNames(   unlist(strsplit(isolate(input$background), "\n"))  )
                   
                   )
            
            
            lipidIDlogical <-
              factor(as.integer(lipidIDs$backgroundlist %in% lipidIDs$sublist))
            names(lipidIDlogical) <- lipidIDs$backgroundlist
            
            if(any(names(lipidID2TERM) %in% names(lipidIDlogical)) &
                                   sum(lipidIDlogical=="1")  > 1  ){   ### can at least one lipid be matched?
              ONTdata <- new(
                "topONTdata",
                ontology = "LION",
                allGenes = lipidIDlogical,
                annot = annFUN.gene2GO,
                gene2GO = lipidID2TERM) 
            } else {
              
              names(lipidIDlogical)[1] <- "PC(38:1)"                   ### make a fake object
              levels(lipidIDlogical) <- c(levels(lipidIDlogical),"1")
              lipidIDlogical[1] <- "1"
              
              ONTdata <- new(
                "topONTdata",
                ontology = "LION",
                allGenes = lipidIDlogical,
                annot = annFUN.gene2GO,
                gene2GO = lipidID2TERM)
            }
            
            
            resultFis <-
              runTest(ONTdata,
                      algorithm = "classic",
                      statistic = "fisher")
            to_display <-
              GenTable(ONTdata,
                       'p-value' = resultFis,
                       topNodes = 2000)
            
            to_display <-
              to_display[grep("LION", to_display$TERM.ID), ]
            
            LUT <- data.frame(ID = names(ONTdata@termName),
                              Discription = ONTdata@termName)
            
            to_display$Term <- LUT$Discription[match(to_display$TERM.ID, LUT$ID)]      
            to_display <- to_display[to_display$Annotated > 2, ]        ### remove terms with 1 or 2 lipids
            
            ####  limiting by LION-term selection
            if(isolate(input$LIONselection)){    ### switch for LION-selection
              TermsOfInterest <- get_selected(isolate(input$tree), format = "names") 
              TermsOfInterest <- unique(c(unlist(TermsOfInterest),
                                          unlist(sapply(TermsOfInterest, function(element){attr(element,"ancestry")}))
              )) ## vector with LION-names
              to_display <- to_display[to_display$Term %in% TermsOfInterest,]  
            }
            
            ## limiting redundant/similar parent terms
            
            if(isolate(input$RemoveRedundantTerms)){    ### switch for LION-selection
              ONT_DAG <- ONTdata@graph
              ONT_DAG_lev <- buildLevels(ONT_DAG)
              DAG.env <- ONT_DAG_lev$nodes2level
            
              DAGlevel <-   sapply(to_display$TERM.ID, function(id) {
                DAG.env[[id]]
              })
              
              TermContent <- 
              sapply(to_display$TERM.ID, function(id) {
                genesInTerm(ONTdata)[[id]]
              })
              
              test_similarity <-
                sapply(names(TermContent), function(term_i) {
                  sapply(names(TermContent), function(term_j) {
                    if (term_i != term_j) {
                      if (length(TermContent[[term_i]]) == length(TermContent[[term_j]])) {
                        same <- all(TermContent[[term_i]] %in% TermContent[[term_j]])
                      } else {
                        same <- FALSE
                      }
                      
                      
                      outputList <- list(
                        term_a = term_i,
                        term_b = term_j,
                        isSame = same,
                        term_a_level = DAGlevel[names(DAGlevel) == term_i],
                        term_b_level = DAGlevel[names(DAGlevel) == term_j]
                      )
                      
                      
                      outputList$remove <- ""
                      
                      if ((outputList$term_a_level - outputList$term_b_level)==1) {
                        outputList$remove <- outputList$term_b
                      } 
                      if ((outputList$term_a_level - outputList$term_b_level)==-1) {
                        outputList$remove <- outputList$term_a
                      }
                      
                      outputList
                      
                    } else {
                      list(
                        term_a = term_i,
                        term_b = term_j,
                        remove = "",
                        isSame = FALSE
                      )
                    }
                  }, simplify = FALSE)
                })
              
              test_similarity <-
                test_similarity[unlist(lapply(test_similarity, function(n) {
                  n[['isSame']]
                }))]
              
              test_similarity <- 
              unique(unlist(lapply(test_similarity, function(n) {
                n[['remove']]
              })))
              
              print(to_display[to_display$TERM.ID %in% test_similarity,]  )
              to_display <- to_display[!(to_display$TERM.ID %in% test_similarity),]  
              
            }
            
            
            ###
            
            to_display$`p-value` <- gsub("< ","",to_display$`p-value`)    ### '< 1e-30' cannot be understood
            
            to_display$'FDR q-value' <-
              format(p.adjust(as.numeric(to_display$`p-value`), "fdr"), digits = 3)
            
            colnames(to_display) <-
              c(
                "Term ID",
                "Discription",
                "Annotated",
                "Significant",
                "Expected",
                "p-value",
                "FDR q-value"
              )
            
        
            
          } else {
            to_display <- data.frame()
          }
        } ## end if "By target list"
        
        if (isolate(input$sub_method) == "(ii) analysis") {
          if ((is.null(input_listwPvalues) ||
               input_listwPvalues == "") == FALSE) {
            pValueList <- list()
            ### seperate by tab or | >> comma gives problems with input: PC(18:1(3Z)/20:4(3Z,6Z,9Z,12Z))
            pValueList$input <-   convertLipidNames(  transpose(strsplit(unlist(strsplit(input_listwPvalues, "\n")), "[\t\\|]"))[[1]]  )
            pValueList$pvalues <-   transpose(strsplit(unlist(strsplit(input_listwPvalues, "\n")), "[\t\\|]")) [[2]] 
            pValueList <- as.data.frame(pValueList)
            pValueList$pvalues <- as.numeric(pValueList$pvalues)
            
            lipidIDrank <- rank(pValueList$pvalues * direction)    ## direction by radiobutton input$ranking
            names(lipidIDrank) <- pValueList$input   
            
            mySel <- function(allScore) {
              return(rep(TRUE, length(lipidIDrank)))
            }
            
            if(any(names(lipidID2TERM) %in% names(lipidIDrank))){   ### can at least one lipid be matched?
              ONTdata <- new(
                "topONTdata",
                ontology = "LION",
                allGenes = lipidIDrank,
                annot = annFUN.gene2GO,
                gene2GO = lipidID2TERM,
                geneSelectionFun = mySel
              )
            } else {       ### make a fake object
              names(lipidIDrank)[1] <- 'PC(32:0)'
              ONTdata <- new(
                "topONTdata",
                ontology = "LION",
                allGenes = lipidIDrank,
                annot = annFUN.gene2GO,
                gene2GO = lipidID2TERM,
                geneSelectionFun = mySel
              )
            }
            
            
            resultFis <-
              runTest(ONTdata,
                      algorithm = "classic",
                      statistic = "ks")
            
            incProgress(.7, detail = paste("enrichment statistics"))
            
            to_display <-
              GenTable(ONTdata,
                       'p-value' = resultFis,
                       topNodes = 2000)
            
            to_display <-
              to_display[grep("LION", to_display$TERM.ID), ]
            
    
            
            LUT <- data.frame(ID = names(ONTdata@termName),
                              Discription = ONTdata@termName)
            
            to_display$Term <- LUT$Discription[match(to_display$TERM.ID, LUT$ID)]        ### otherwise, some names are abbrev.
            to_display <- to_display[to_display$Annotated > 2, ]        ### remove terms with 1 or 2 lipids
            
            ####  limiting by LION-term selection
            if(isolate(input$LIONselection)){    ### switch for LION-selection
              TermsOfInterest <- get_selected(isolate(input$tree), format = "names") 
              TermsOfInterest <- unique(c(unlist(TermsOfInterest),
                                          unlist(sapply(TermsOfInterest, function(element){attr(element,"ancestry")}))
                                     )) ## vector with LION-names
              to_display <- to_display[to_display$Term %in% TermsOfInterest,]  
            }
            
            
            ## limiting redundant/similar parent terms
            
            ## limiting redundant/similar parent terms
            
            if(isolate(input$RemoveRedundantTerms)){    ### switch for LION-selection
              ONT_DAG <- ONTdata@graph
              ONT_DAG_lev <- buildLevels(ONT_DAG)
              DAG.env <- ONT_DAG_lev$nodes2level
              
              DAGlevel <-   sapply(to_display$TERM.ID, function(id) {
                DAG.env[[id]]
              })
              
              TermContent <- 
                sapply(to_display$TERM.ID, function(id) {
                  genesInTerm(ONTdata)[[id]]
                })
              
              test_similarity <-
                sapply(names(TermContent), function(term_i) {
                  sapply(names(TermContent), function(term_j) {
                    if (term_i != term_j) {
                      if (length(TermContent[[term_i]]) == length(TermContent[[term_j]])) {
                        same <- all(TermContent[[term_i]] %in% TermContent[[term_j]])
                      } else {
                        same <- FALSE
                      }
                      
                      
                      outputList <- list(
                        term_a = term_i,
                        term_b = term_j,
                        isSame = same,
                        term_a_level = DAGlevel[names(DAGlevel) == term_i],
                        term_b_level = DAGlevel[names(DAGlevel) == term_j]
                      )
                      
                      
                      outputList$remove <- ""
                      
                      if ((outputList$term_a_level - outputList$term_b_level)==1) {
                        outputList$remove <- outputList$term_b
                      } 
                      if ((outputList$term_a_level - outputList$term_b_level)==-1) {
                        outputList$remove <- outputList$term_a
                      }
                      
                      outputList
                      
                    } else {
                      list(
                        term_a = term_i,
                        term_b = term_j,
                        remove = "",
                        isSame = FALSE
                      )
                    }
                  }, simplify = FALSE)
                })
              
              test_similarity <-
                test_similarity[unlist(lapply(test_similarity, function(n) {
                  n[['isSame']]
                }))]
              
              test_similarity <- 
                unique(unlist(lapply(test_similarity, function(n) {
                  n[['remove']]
                })))
              
              print(to_display[to_display$TERM.ID %in% test_similarity,]  )
              to_display <- to_display[!(to_display$TERM.ID %in% test_similarity),]  
              
            }
            
            
            ###
             
            
            to_display$`p-value` <- gsub("< 1e","< 1.0e",to_display$`p-value`)    ### '< 1e-30' cannot be understood
            to_display$`p-value` <- gsub("< ","",to_display$`p-value`)    ### '< 1e-30' cannot be understood
            
            to_display$'FDR q-value' <-
              format(p.adjust(as.numeric(to_display$`p-value`), "fdr"), digits = 3)
            
            colnames(to_display) <-
              c(
                "Term ID",
                "Discription",
                "Annotated",
                "Significant",
                "Expected",
                "p-value",
                "FDR q-value"
              )
            to_display <- to_display[, c(1, 2, 3, 6, 7)]
            
        
          } else {
            to_display <- data.frame()
          }
        }
        
        
        lengthOfTable <- length(to_display$"Term ID")
        
        ### make detailed table with lipids per term
        to_display_detailed <- to_display
        
        
        identifiers <- lapply(to_display_detailed$'Term ID', function(ID){
          genesInTerm(ONTdata)[[ID]]
        })
        
  ####  
        
        if(length(to_display_detailed[,1]) > 0){     ### if this is zero, script goes wrong
          
          output <- lapply(1:length(to_display_detailed[,1]), function(row){   ### add rows with lipid info
            output <- as.data.frame(matrix(ncol = length(to_display_detailed[1,]),
                                           nrow = length(identifiers[[row]])+1,
                                           data=""))
            output[1,] <-   to_display_detailed[row,]
            output[-1,1] <- paste("   ",identifiers[[row]])
            
            output
          })
          
          to_display_detailed <- as.data.frame(rbindlist(output))
          colnames(to_display_detailed) <- colnames(to_display)
        } else {to_display_detailed <- to_display}
        
        
        ### lipid associations report
        
        lipidReport <- associatedTerms(lipid = ONTdata@allGenes, ontologyObject = ONTdata, reformat = TRUE)
        colnames(lipidReport) <- lipidReport[1,]
        lipidReport <- lipidReport[-1,]
        
        ### term associations report
        
        lipidsInTerms <- genesInTerm(ONTdata)
        
        lipidsInTerms <- data.frame(ID = names(lipidsInTerms),
                                 Discription = LUT$Discription[match(names(lipidsInTerms), LUT$ID)],
                                 nr_lipids = sapply(lipidsInTerms, function(term){length(term)}),
                                 lipids = sapply(lipidsInTerms, function(term){paste(term, collapse = '; ')}))
        
        lipidsInTerms <- subset(lipidsInTerms, Discription != lipids)
        colnames(lipidsInTerms) <- c("LION-term","LION-name", "nr_lipids", "lipid identifiers")
        
        
        ### ggplot of data
        
        to_ggplot_display <- to_display
        to_ggplot_display$color <-
          -log(as.numeric(to_ggplot_display$`FDR q-value`), base = 10)
        to_ggplot_display$color[to_ggplot_display$color > 6] <- 6
        
        to_ggplot_display <- ggplot(data = to_ggplot_display[1:min(lengthOfTable,40), ],
                                    aes(
                                      y = -log(as.numeric(`FDR q-value`), base = 10),
                                      x = reorder(Discription, -log(as.numeric(
                                        `FDR q-value`
                                      ) , base =  10))
                                    )) +
          geom_hline(yintercept = -log(0.05, base = 10),
                     alpha = .3) +
          geom_bar(stat = 'identity', aes(fill = color), width = .70) +
          scale_fill_gradient2(
            limits = c(0, 6),
            midpoint = -log(0.05, base = 10),
            low = "grey",
            mid = "grey",
            high = "red"
          ) +
          ggtitle("LION-term enrichment") +
          xlab("") +
          ylab("-LOG(FDR q-value)") +
          guides(fill = "none") +
          #scale_y_continuous(limits = c(0, 6))+
          coord_flip() +
          theme_pander()
        
        incProgress(.8, detail = paste("enrichment graph"))
        
        ### network
        
        
        if (isolate(input$method) == "Target-list mode") {
          resultTable <-
            GenTable(ONTdata,  'p-value' = resultFis, topNodes = 2000)
        }
        
        if (isolate(input$sub_method) == "(ii) analysis") {
          resultTable <-
            GenTable(ONTdata,  'p-value' = resultFis, topNodes = 2000)
          
        }
        
        if (exists("resultFis")) {
          
          resultTable$`p-value` <- gsub("< ","",resultTable$`p-value`)    ### '< 1e-30' cannot be understood
          
          resultTable$'FDR q-value' <-
            format(p.adjust(as.numeric(resultTable$`p-value`), "fdr"), digits = 3)
          nrOfSignNodes <-
            max(c(sum(
              as.numeric(resultTable$'FDR q-value') < 0.05
            ), 5))   ## with a minimum of 5
          
          
          #nr of nodes by FDR qvalue
          enrichmentGraph <-
            showSigOfNodes(
              ONTdata,
              score(resultFis),
              firstSigNodes = nrOfSignNodes,
              useInfo = 'all',
              swPlot = FALSE
            )
          
          enrichment_iGraph <-
            graph_from_graphnel(enrichmentGraph$dag)
          
          nodes <-
            data.frame(id = unique(c(
              get.data.frame(enrichment_iGraph)$from,
              get.data.frame(enrichment_iGraph)$to
            )))
          nodes$label <-
            resultTable$Term[match(nodes$id, resultTable$TERM.ID)]
          nodes$shape <- "dot"
          
          #val_col <- colorRampPalette(c("#BEBEBE","#BEBEBE","#E8AD43","#FFC100","#FFFF00","red"))(100)
          val_col <-
            colorRampPalette(c("#BEBEBE", "#BEBEBE", "yellow", "orange", "red"))(100)
          logPs <-
            -log(as.numeric(resultTable$`p-value`[match(nodes$id, resultTable$TERM.ID)]), 10)
          logPs[logPs < 1] <- 1
          logPs[logPs > 10] <- 10
          nodes$color <- val_col[logPs * 10]
          
          nodes$color[grep("CAT", nodes$id)] <- "#EEEEEE"
          nodes$shape[grep("CAT", nodes$id)] <- "box"
          nodes$color[grep("all", nodes$id)] <- "black"
          nodes$shape[grep("all", nodes$id)] <- "triangle"
          
          ## size as function of nr of annotions
          nodes$size <-
            0.2 * as.numeric(resultTable$Annotated[match(nodes$id, resultTable$TERM.ID)])
          ##
          
          nodes$size[nodes$size < 4] <- 4
          nodes$size[grep("all", nodes$id)] <- 20
          nodes$label[grep("all", nodes$id)] <- "ontology root"
          
          ## get levels
          ONT_DAG <- ONTdata@graph
          ONT_DAG <- buildLevels(ONT_DAG)
          DAG.env <- ONT_DAG$nodes2level
          DAG.env$`LION:0080986`
          nodes$level <-   sapply(nodes$id, function(id) {
            DAG.env[[id]]
          })
          
          edges <- get.data.frame(enrichment_iGraph)
          edges$color <- "grey"
          edges$arrows <- "from"
          
          ####  limiting by LION-term selection
          if(isolate(input$LIONselection)){    ### switch for LION-selection
            #TermsOfInterest ## vector with LION-names from previous section
            nodes <- nodes[nodes$label %in% c(TermsOfInterest,"ontology root"),]
            ## delete nodes not selected for analysis
          }
          ###
          
          ### delete dead ends of CATs w/o LION-terms
          CATs <-
            edges$from[paste(gsub(":\\d+", "", edges$from),
                             gsub(":\\d+", "", edges$to),
                             sep = "") %in% c("CATLION", "CATCAT")]
          CATs <- nodes$id[!(nodes$id %in% CATs)]
          CATs <- CATs[grepl("CAT", CATs)]
          
          nodes <- nodes[!nodes$id %in% CATs, ]
          
          to_display_network <- visNetwork(nodes, edges) %>%
            visPhysics(stabilization = TRUE) %>%
            visHierarchicalLayout(direction = "LR", levelSeparation = 250) %>%
            visInteraction(navigationButtons = TRUE)   %>% 
            visPhysics(
              solver = "forceAtlas2Based",
              forceAtlas2Based = list(gravitationalConstant = -20)
            )
          
          incProgress(.9, detail = paste("network"))
          
          ### return as list object
          
        }
        
        if(isolate(input$EmailMissingLipids)){   ## email missing lipids when option is set
          #sendEmail(subject = "missing annotations",mail_message = paste(subset(lipidExistance,`LION ID`=="not found")$input, collapse = "\n"))
          sendEmail(subject = "missing annotations", mail_message = paste(apply(lipidExistance, 1, function(row){paste(row, collapse = "\t")}), collapse = "\n"))
        }
        
        
        list(
          ONTdata = ONTdata,
          lipidExistance = lipidExistance,
          to_display = to_display,
          to_display_detailed = to_display_detailed,
          lipidReport = lipidReport,
          lipidsInTerms = lipidsInTerms,
          to_ggplot_display = to_ggplot_display,
          to_display_network = to_display_network
          
        )
        
      })
      } else {
        ### empty if there is no input
        list(
          ONTdata = NULL,
          lipidExistance = NULL,
          to_display = NULL,
          to_display_detailed = NULL,
          lipidReport = NULL,
          lipidsInTerms = NULL,
          to_ggplot_display = NULL,
          to_display_network = NULL
          
        )
        
      }  ### if perform when there is input (dim 1 by 1 >> no input)
    
    
  })
      
  ### mapping percentage
  output$mapping_percentage <- renderText({          
    
    if (is.null(dataBlock()$lipidExistance)){  
      "Unexpected input"  
    } else {
      mappings <- !(dataBlock()$lipidExistance$'LION ID' == 'not found')
      
      if (length(isolate(dataBlock()$to_display[,1])) == 0){
        comment <- "\nHowever, no upstream LION-terms were found for analysis."
        hideTab(inputId = "tabs", target = "LION enrichment table")
        hideTab(inputId = "tabs", target = "LION enrichment graph")
        hideTab(inputId = "tabs", target = "LION network view")
      } else comment <- ""
      
      paste(
        sum(mappings), " out of ",  length(mappings),  " ", "(", round(mean(mappings)  * 100, digits = 2), "%", ")",
        " identifiers could be matched to LION.", comment,
        sep = ""
      )
    }
  })    
      
      
      
      ###
  output$downloadInputCNTRL <- renderUI({
       if(is.null(dataBlock()$lipidExistance)){} else {
         downloadButton("downloadInput", "Download input table")
       }
    })
  
   
      output$value <- renderTable({
        table <- dataBlock()$lipidExistance
        #if(isolate(input$EmailMissingLipids)){   ## email missing lipids when option is set
        #  sendEmail(subject = "missing annotations",mail_message = paste(subset(table,`LION ID`=="not found")$input, collapse = "\n"))
        #}
        table
      })
      
      output$downloadInput <- downloadHandler(
        filename = function() {
          paste("LION-input-job",isolate(input$submitB)+isolate(input$submitA), ".csv", sep="")
          
        },
        content = function(file) {
          write.csv(dataBlock()$lipidExistance, 
                    file, row.names = FALSE,quote = TRUE)
        }
      )
      
      # Generate an HTML table view of the ontology data
      output$ontology.table <- renderTable({
        dataBlock()$to_display
        
      })

      output$downloadTable <- downloadHandler(
        filename = function() {
          paste("LION-enrichment-job",isolate(input$submitB)+isolate(input$submitA), ".csv", sep="")
          
        },
        content = function(file) {
          write.csv(dataBlock()$to_display, 
                    file, row.names = FALSE,quote = TRUE)
        }
      )
      
      output$downloadDetailTable <- downloadHandler(
        filename = function() {
          paste("LION-enrichment-report-job",isolate(input$submitB)+isolate(input$submitA), ".zip", sep="")
          
        },
        content = function(file) {
          write.csv(dataBlock()$to_display_detailed, 
                    file = "LION-enrichment-detailed-job.csv", row.names = FALSE,quote = TRUE)
          write.csv(dataBlock()$lipidReport, 
                    file = "LION-lipid-associations.csv", row.names = FALSE,quote = FALSE)
          write.csv(dataBlock()$lipidsInTerms, 
                    file = "LION-term-associations.csv", row.names = FALSE,quote = TRUE)
          
          
        
          zip(zipfile=file, files=c( "LION-enrichment-detailed-job.csv", "LION-lipid-associations.csv","LION-term-associations.csv"))
          
          
        }
      )
      

      output$downloadPlotPNG <- downloadHandler(
        filename = function() { 
          paste("LION-enrichment-plot-job",isolate(input$submitB)+isolate(input$submitA),  '.png', sep='') },
        content = function(file) {
          ggsave(file,dataBlock()$to_ggplot_display,   device = 'png', dpi=300)
        }
      )
      output$downloadPlotSVG <- downloadHandler(
        filename = function() { 
          paste("LION-enrichment-plot-job",isolate(input$submitB)+isolate(input$submitA),  '.svg', sep='') },
        content = function(file) {
          ggsave(file,dataBlock()$to_ggplot_display,   device = svg,   scale=1.5)
        }
      )
        
      output$ontology.graph <- renderPlot({
         dataBlock()$to_ggplot_display
               
      })
      
      
      ## ontology enrichment tree
      output$ontology.plot <- renderVisNetwork({
        dataBlock()$to_display_network
      })
      
      
      updateSelectizeInput(session, 'LIONterms', 
                           choices =  LIONTermList, 
                           selected = NULL )
      
      ## email
      observe({
        if(is.null(input$send) || input$send==0) return(NULL)
        from <- isolate(input$from)
        to <- "m.r.molenaar@uu.nl"
        subject <- isolate(input$subject)
        msg <- isolate(input$message)
        
        sendEmail(subject = subject, from = from, mail_message = msg)
        showModal(modalDialog(
          footer = modalButton("Ok"),
          size = "m",
          easyClose = TRUE,
          title = paste("Confirmation:",subject), 
          p("Thank you for contacting us. Your message was sent successfully.")
          
        ))
      })
      
      
      
      output$tree <- renderTree({        ### LION-tree
        LIONstructure
      })
      
      
}
