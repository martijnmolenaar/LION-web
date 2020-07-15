options(stringsAsFactors = FALSE,shiny.sanitize.errors = F)

 
## examples

cluster_example <- read.csv(file = "data/clusters_lipids.csv")
PM <- read.csv(file = "data/PM-ER_pvalues.csv")
MT <- read.csv(file = "data/Mito-pvalues.csv")
ER <- read.csv(file = "data/ER_KLA-ER_CON_pvalues.csv")

LIONterms_rules <- read.csv(file = 'data/20191008 LIONterms_rules.csv', header = T)
LIONterms_rules$RULE1[LIONterms_rules$RULE1 == ""] <- "-"
LIONterms_FAs <- read.csv(file = 'data/20191008 LIONterms_FAs.csv', header = T)

backgroundlist_example <- paste(cluster_example$lipids, "\n", collapse = '', sep = "")

sublist_example1 <- paste(cluster_example$lipids[cluster_example$cluster == 6], "\n", collapse = '', sep = "")
sublist_example2 <- paste(cluster_example$lipids[cluster_example$cluster == 7], "\n", collapse = '', sep = "")
sublist_example3 <- paste(cluster_example$lipids[cluster_example$cluster == 8], "\n", collapse = '', sep = "")

pvalueExample1 <-  paste( paste(PM$lipids, "\t", PM$pvalues, sep = ""),    "\n",    collapse = '',    sep = ""  )
pvalueExample2 <-  paste( paste(MT$lipids, "\t", MT$pvalues, sep = ""),    "\n",    collapse = '',    sep = ""  )
pvalueExample3 <-  paste( paste(ER$lipids, "\t", ER$pvalues, sep = ""),    "\n",    collapse = '',    sep = ""  )

FA_composition_table <- read.csv(file = 'data/FA_composition_table.csv')

source(file  = "data/20191010 LION_tree_structure.R")

#### libraries

require(shiny)
require(visNetwork)
require(data.table)
require(igraph)
require(ggplot2)
require(ggthemes)
library(shinyTree)
library(shinyWidgets)
library(shinyBS)
library(httr)
library(formattable)
library(jsonlite)
library(ggrepel)
library(shinycssloaders)


## loading lipid ontology data
require(RSQLite)
require(topOnto)
require('topOnto.LION.db')
topOnto::initONT('LION')


associationFile  <-  "data/20190704 LION_association.txt"


## load define functions
source('data/20191008 LIONweb functions.R')


## read associations
lipidID2TERM <- readMappings(file = associationFile)   ## topOnto function, but removes spaces


# Define server logic for random distribution application
function(input, output, session) {
  
  showNotification(ui = "",
                   action =  p("By using this app you agree with the", a('Terms of Usage.',
                               href="https://martijnmolenaar.github.io/lipidontology.com/faq.html#basics", 
                               target="_blank")),
                   duration = 14, type = "default")
 
  ### hide tabs at start-up
  hideTab(inputId = "tabs", target = "LION input")
  hideTab(inputId = "tabs", target = "LION enrichment table")
  hideTab(inputId = "tabs", target = "LION enrichment graph")
  hideTab(inputId = "tabs", target = "LION network view")
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    
    if(!is.null(query[['studyid']])){
      updatePrettyRadioButtons(
        session,
        inputId = "file_input",
        selected = "load external dataset"
      )
      
    }
  })
  
 
  ## pre-processing with CSVs:
  
  input_data <- reactive({
    
        
    ifelse(input$file_input=="file input",
               req(input$file1),
               req(input$MWfile1))
    
    

 
    isolate(if(input$file_input == "file input"){
      file_location <- input$file1$datapath
      metadata <- input$file1$name
    } else if(input$file_input == "load external dataset"){
      file_location <- paste("https://www.metabolomicsworkbench.org/rest/study/study_id/",
                             input$MWfile1,
                             "/lion/", sep = "")
      
      metadata <- fromJSON(paste("https://www.metabolomicsworkbench.org/rest/study/study_id/",
                                           input$MWfile1,
                                           "/summary", sep = ""))$study_title
    })
    
  
    withProgress(message = 'loading data:', value = 0, {
      if(grepl("metabolomicsworkbench.org/rest", file_location)){
        incProgress(.1, detail = paste("from Metabolomis Workbench"))
        
        df <- try(read.csv(file_location,
                           header = FALSE,
                           sep = "," #, quote = ""
        ))
        
        incProgress(.9, detail = paste("..done"))
      } else {
        
        df <- try(read.csv(file_location,
                           header = FALSE,
                           sep = "," #, quote = ""
        ))
      }
    })
    
    
    ###
    
    ## error handling
    
    errors <- NULL
    if (!(is.data.frame(df))) {
      df <- data.frame(errors = df[1], column = 0)
      if(isolate(input$file_input) == "load external dataset"){
        errors <- c(errors, "ERROR: Metabolomics Workbench study ID not found")
      } else {
        errors <- c(errors, "ERROR: File type not supported")
      }
      
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
      if (sum(apply(df[-1,-c(1:2)],1,function(row){all(row == 0)})) > 0) {
        errors <-
          c(errors, "ERROR: Dataset contains rows with only zeros")
      }
    }
    errors <- paste(errors, collapse = "<br>")
    
    ## end error handling
    
    
    isolate(if(input$file_input == "load external dataset"){
      studyID <- paste(input$MWfile1,": ",sep = "")
    } else {studyID <- ""})
    
    #if(!(is.null(errors))){
    if(errors==""){
      
       input_data <- list(df = df,
         matrix = sapply(df[-c(1,2),-1], as.numeric),
         conditions = unique(as.character(df[1,-1])),
         samples = as.character(df[2,-1]),
         meta = df[1:2,-1],
         IDs = df[-c(1,2),1],
         studyID = studyID,
         metadata = metadata,
         errors = errors)
    } else {
      input_data <-list(df = NULL,
           matrix = NULL,
           conditions = NULL,
           samples = NULL,
           meta = NULL,
           IDs = NULL,
           studyID = NULL,
           metadata = NULL,
           errors = errors)
    }
    
    return(input_data)
    
  })
  
  output$load_datasetUI <- renderUI({
    file_input <- input$file_input
    
    if (file_input == "file input") {
      
      fluidRow(column(offset = .3, width = 12,
      
      br(),
      strong("Choose CSV File:"),
      fluidRow(
        
        column(offset=0,#style='margin-left:2%;' ,
               width = 10,
               popify(placement = "bottom", title = "File-input info",
                      fileInput("file1", label = NULL,
                                multiple = FALSE,
                                accept = c("text/csv",
                                           "text/comma-separated-values,text/plain",
                                           ".csv") ),
                      content = 'Format your dataset as comma seperated value files (.csv), with the first column reserved for metabolites and the other columns for numeric data (containing decimal points). Use double column headers; with row 1 containing condition identifiers and row 2 containing sample identifiers. Submit at least duplicates per condition. Dataset should be normalized before submission. Download a dataset below for an example.'
                      
               )),
        column(offset=0, width = 1,style = "margin-top: 5px;",align="center",
               popify(placement = "right", title = "Lipid nomenclature", options=list(container="body"),
                      el = icon(name = "question",lib = "font-awesome", "fa-2x"),
                      content = 'Format lipids in LIPIDMAPS notation style: a class-prefix followed by (summed) fatty acid(s) surrounded by parentheses. Examples are: PC(32:1); PE(18:1/16:0); SM(d18:1/18:0); TAG(54:2); etc. Check www.lipidmaps.org for more examples. LION will try to reformat alternative notation styles into LIPIDMAPS format.'))),
      
      downloadLink("examplePre1", "example set 1 (organelle fractions) [1]"),
      br(),
      downloadLink("examplePre2", "example set 2 (CHO-k1 incubated with several FFAs) [2]"),
      br(),
      downloadLink("examplePre3", "example set 3 (CHO-k1 incubated with AA) [2]"),
      br(),
      br(),
      em("[1] Andreyev AY, et al., 2010"),
      br(),
      em("[2] Molenaar MR, et al., 2019"),
      br(),br(),
      uiOutput("selectLocalStatisticUI")
      ))
  
      
    } else {     ### to load Metabolomics Workbench-file
      fluidRow(column(offset = .3, width = 12,
        
        br(),
        
        fluidRow(
          
          column(offset=0,#style='margin-left:2%;' ,
                 width = 10,
                 popify(placement = "bottom", title = "Metabolomics Workbench",
                        searchInput(
                          inputId = "MWfile1",
                          label = "Enter Metabolomics Workbench study ID", 
                          placeholder = "study ID", 
                          btnSearch = icon("search"), 
                          btnReset = icon("remove"),
                          width = "100%"
                        ), content = "In this field, a  study ID from the Metabolomics Workbench data repository can be entered—the dataset will be loaded directly to LION/web in the required format. To browse studies, please visit the Metabolomics Workbench website (see below). Please note that LION/web does only process lipids."
                        
                 ))
          ),
        
        br(),
        em("example: "),br(),
        "ST001140: ",
        a('Changes in the Canine Plasma Lipidome after Short- and Long-Term Excess Glucocorticoid Exposure',  
          href="https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001140", target="_blank"),
        br(),br(),
        em("For more information, please visit the",  a('Metabolomics Workbench data repository.',  href="https://www.metabolomicsworkbench.org/data/browse.php", target="_blank")),
        br(),
        br(),
        uiOutput("selectLocalStatisticUI")
      ))
    }
    
  })
  
  
  observeEvent(input$file_input, {
    
    query <- parseQueryString(session$clientData$url_search)
    
    if(!is.null(query[['studyid']])){
      updateSearchInput(
        session,
        inputId = "MWfile1",
        value = query[['studyid']],
        trigger = TRUE,
        label = "Enter Metabolomics Workbench study ID"
        
      )
  
    }
    })

  output$selectLocalStatisticUI <- renderUI({
    
    input_data <- input_data()
    conditions <- input_data$conditions
    
    
    if (input_data$errors == "") {
      
      wellPanel(
        tags$b(input_data$metadata),
        br(),br(),
       
        popify(placement = "bottom", title = "Info",
        radioButtons("local_statistics", label = h5("Select local statistics to rank input identifiers"),
                     choices = list("one-tailed T-test (2 conditions)" = 1, "2-LOG[fold change] (2 conditions)" = 2, "one-way ANOVA F-test (>2 conditions)" = 3), 
                     selected = 1), content = 'Select a local statistic to rank the metabolites based on the provided dataset. LION-terms associated with lipids higher ranked than expected by chance will be reported as enriched.'),
        popify(
          placement = "bottom",
          title = "Info",
          checkboxInput(
            inputId = "normalization",
            label = "normalize signals as percentage",
            value = FALSE
          ),
          content = 'If checked, lipid signals are expressed as percentage of the total signal per sample. This is useful when input data is not normalized.'
        ), 
        br(),
        uiOutput("conditionsUI")
      )
    } else {     ### if there are errors...
      HTML(paste('<font color="red">',
                 input_data$errors,
                 "</font>",sep=""))
    }
    
  })
  
  output$KS2_UI<- renderUI({
    
    if(input$ks_sided == "ks2"){
      
      fluidRow(
        column(width = 1),
        column(width = 11,
               selectizeInput("split_direction_in_output", 
                              label = NULL, 
                              choices = c("Separate up- and downregulated terms in barchart" = "split", 
                                          "Combine up- and downregulated terms in barchart" = "combine",
                                          "Display data in volcano plot" = "volcano"), 
                              selected = "split"))
      )
      
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
              h5("condition of interest"),
              conditions,
              selected = conditions[1]
            )
          ),
          column(
            width = 6,
            selectInput(
              "conditionB",
              h5("control condition"),
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
      fluidRow(
      actionButton(inputId = "submitPreprocessing",
                   label = "  Use values as local statistics", 
                   width = NULL,
                   icon = icon("share-square-o", lib = "font-awesome")),
      downloadButton("download_pValues", "")
      )
    }
    
    
  })

  pre_processed_data <-     eventReactive(input$calculatePreprocessing, {

    input_data <- input_data()

    if(isolate(input$normalization)){  ## normalization option switched on
      input_data$matrix[] <- 
        apply(input_data$matrix,2,function(i){
          i / sum(i) * 100
        })
    }
    
    if(any(input$local_statistics %in% c(1,2))){
    setA <- input_data$matrix[,which(input_data$meta[1,] == input$conditionA)]
    setB <- input_data$matrix[,which(input_data$meta[1,] == input$conditionB)]
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
    
    if(input$local_statistics == 3){   ### when ANOVA is selected
      
      #input_data <- input_data()
      
      if (length(input$selectedConditions) > 1) {      ## minimum of 2 conditions needed
        pValues <- 
          sapply(1:length(input_data$IDs), function(lipid_nr){
                       df <- do.call("rbind",sapply(input_data$conditions, function(condition){
                         data.frame(l = NA,
                                    condition = condition, 
                                    signal = input_data$matrix[lipid_nr,  unlist(input_data$meta[1,]) == condition]
                         )}, simplify = FALSE))
                       
                       model.df <- lm(signal ~ condition, data = df)
                       anova(model.df)['condition','Pr(>F)']
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
  
  
  observe({
    input$download_network_bs  
      visNetworkProxy("ontology.plot") %>% visGetPositions()   
  })
  
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
    
    input_listwPvalues <- input$listwPvalues
     
    if(grepl("\t", input_listwPvalues) & 
         grepl("\n", input_listwPvalues) &
         grepl("\\D",input_listwPvalues) &
         grepl("\\d",input_listwPvalues) ){   ## input requirements
      
    #if(input$listwPvalues!=""){
      showTab(inputId = "tabs", target = "LION input")
      showTab(inputId = "tabs", target = "LION enrichment table")
      showTab(inputId = "tabs", target = "LION enrichment graph")
      showTab(inputId = "tabs", target = "LION network view") 
      
      dataBlock()
      
      ## goto LION input
      updateTabsetPanel(session, "tabs",
                      selected = "LION input")
    } else {
      
      hideTab(inputId = "tabs", target = "LION input")
      hideTab(inputId = "tabs", target = "LION enrichment table")
      hideTab(inputId = "tabs", target = "LION enrichment graph")
      hideTab(inputId = "tabs", target = "LION network view") 
      
      updateTabsetPanel(session, "tabs",
                        selected = "General information")
      showNotification(ui = "",
                       action =  "No or incorrect input found. Please reformat and submit data in the text box.",
                       duration = 3, type = 'error')
    }
  })
  observeEvent(input$submitA, {
    input_list <- c(input$sublist, input$background)
    
    if(all(grepl("\n", input_list) & grepl("\\D",input_list))){   ## input requirements
    
    #if(input$sublist!="" & input$background!=""){
      ## first unhide tabs
      showTab(inputId = "tabs", target = "LION input")
      showTab(inputId = "tabs", target = "LION enrichment table")
      showTab(inputId = "tabs", target = "LION enrichment graph")
      showTab(inputId = "tabs", target = "LION network view")  
    
      dataBlock()
    
    ## goto LION input
      updateTabsetPanel(session, "tabs",
                      selected = "LION input")
    } else {
      
      hideTab(inputId = "tabs", target = "LION input")
      hideTab(inputId = "tabs", target = "LION enrichment table")
      hideTab(inputId = "tabs", target = "LION enrichment graph")
      hideTab(inputId = "tabs", target = "LION network view") 
      
      updateTabsetPanel(session, "tabs",
                        selected = "General information")
      showNotification(ui = "",
                       action = "No or incorrect input found. Please reformat and submit data in the two text boxes.",
                       duration = 3, type = 'error')
    }
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
    
    if (isolate(input$sub_method) == "(ii) analysis" & input$method == "Ranking mode") {
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
        
        ### new 20191010
        cat("substracting input target-list \n")
        
        sublist <- strsplit(isolate(input$sublist), "\n"); names(sublist) <- "ID"
        background <- strsplit(isolate(input$background), "\n"); names(background) <- "ID"
        
        background = data.frame(ID = c(background$ID, sublist$ID[!sublist$ID %in% background$ID]))
        background$IDnr <-  1:length(background$ID)
        background$input <- paste(sprintf("[#%04d]", background$IDnr),background$ID)
        
        sublist$input <-
          do.call("rbind",
                  sapply(unique(sublist$ID), function(ID_i){
                    background[background$ID == ID_i,]
                  }, simplify = FALSE))$input
        
        background$input <-
          do.call("rbind",
                  sapply(unique(background$ID), function(ID_i){
                    background[background$ID == ID_i,]
                  }, simplify = FALSE))$input
        
        
        lipidExistance <-  data.frame(input = background$input )
        cat("converting input target-list \n")
        lipidExistance$'simplified input' <- convertLipidNames(gsub("^\\[#\\d+\\] ","",lipidExistance$input))
        ### end new 20191010
        
        
      } else {
        lipidExistance <- data.frame("input")
      }
      
    }
    
      if (isolate(input$method) == "Ranking mode") {                  ## (isolate(input$sub_method) == "(ii) analysis") {
      if ((is.null(input_listwPvalues) ||
           input_listwPvalues == "") == FALSE) {
        
        
        pValueList <- list()
        
        pValueList$input <-
          transpose(strsplit(unlist(strsplit(
            input_listwPvalues, "\n"
          )), "[\t\\|]"))[[1]]
        
        pValueList$input <- paste(sprintf("[#%04d]", 1:length(pValueList$input)),pValueList$input)
        
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
        
        lipidExistance$'simplified input' <- convertLipidNames(gsub("^\\[#\\d+\\] ","",lipidExistance$input))
        
        
      } else {
        lipidExistance <- data.frame("input")
      }
    }
      if (length(lipidExistance$input) < 1) {
      
    } else {
      
      ##### new 20191106
      lipidExistance_list <-
        lapply(1:dim(lipidExistance)[1], function(row_i) {
          
          lipid_i <- lipidExistance[row_i, 2]   ## lipid_i is lipid of this iteration
          
          LION_ID <-
            unlist(lipidID2TERM[lipid_i == names(lipidID2TERM)])
          
          if (is.null(LION_ID)) {
            ### if input is LION:xxxx
            LION_ID <-
              unlist(lipidID2TERM[lipid_i == unlist(lipidID2TERM)])
            LION_ID <-
              LION_ID[!grepl("^SLM:|^LM..\\d+", names(LION_ID))]   ## remove SwissLipids/LIPIDMAPS IDs
          }
          
          if (!is.null(LION_ID)) {
            output_df <-
              data.frame(
                input = lipidExistance[row_i, 1],
                'simplified input' = lipid_i,
                name = names(LION_ID),
                LION = LION_ID,
                match = "direct"
              )
          } else {
            ### no direct matching is found
            if (isolate(input$SmartMatching)) {
              ### 'smartmatching' is on
              lipid_index <-
                list(
                  generalized = simplifyLipidnames(lipid_i),
                  headgroup = getLipidHeadgroup(lipid_i),
                  linkage = getLipidLinkage(lipid_i),
                  FAs = getFAs(lipid_i)
                )
              
              generalized_LION_ID <-
                unlist(lipidID2TERM[lipid_index$generalized == names(lipidID2TERM)])
              if (!is.null(generalized_LION_ID)) {
                terms <-
                  data.frame(name = names(lipidID2TERM)[lipid_index$generalized == names(lipidID2TERM)],
                             LION = generalized_LION_ID)
              } else {
                ## no generalized lipid found
                LIONterms_rules_i <-
                  LIONterms_rules[lipid_index$headgroup == LIONterms_rules$RULE1 ,]
                
                if (all(LIONterms_rules_i$RULE2 == "")) {
                  terms <- LIONterms_rules_i[, 1:2]
                } else {
                  terms <-
                    LIONterms_rules_i[lipid_index$linkage == LIONterms_rules_i$RULE2,][, 1:2]
                }
              }
              
              
              terms <-
                rbind(terms, LIONterms_FAs[LIONterms_FAs$name %in% lipid_index$FAs,])
              if (dim(terms)[1] > 0) {
                output_df <-
                  data.frame(input = lipidExistance[row_i, 1],
                             'simplified input' = lipid_i,
                             terms,
                             match = "smart matching")
                
              } else {
                ## no match by smart matching
                output_df <-
                  data.frame(
                    input = lipidExistance[row_i, 1],
                    'simplified input' = lipid_i,
                    name = "not found",
                    LION = "not found",
                    match = ""
                  )
              }
              
              
            } else {
              ### 'smartmatching' is off, and no match found
              output_df <-
                data.frame(
                  input = lipidExistance[row_i, 1],
                  'simplified input' = lipid_i,
                  name = "not found",
                  LION = "not found",
                  match = ""
                )
            }
            
            
          }
          if(isolate(input$FAprediction)){    ## predict FAs if nescerrary
            
            lipid_index <-
              list(
                generalized = simplifyLipidnames(lipid_i),
                headgroup = getLipidHeadgroup(lipid_i),
                FAs = getFAs(lipid_i)
              )
            
            if(lipid_index$generalized == lipid_i &   length(lipid_index$FAs) == 1){   ## is FA prediction applicable?
              
              predicted_FAs <- predict_FAs(headgroup = lipid_index$headgroup, 
                                           summedFA = gsub("C","",lipid_index$FAs), 
                                           composition_table = FA_composition_table )
              
              if(dim(LIONterms_FAs[LIONterms_FAs$name %in% predicted_FAs,])[1]>0){    ## any result?
                output_df <-
                  rbind(output_df,
                        data.frame(input = lipidExistance[row_i, 1],
                                   'simplified input' = lipid_i,
                                   LIONterms_FAs[LIONterms_FAs$name %in% predicted_FAs,],
                                   match = "FA-prediction"))
              } 
              
              
            }
            
            
          }
          
          return(output_df)
        })
      
      matching_statistics <- data.frame(total = length(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})),
                                        matched = sum(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})),
                                        percent = mean(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})    ) * 100)
      
      
      lipidExistance <- do.call("rbind",lipidExistance_list)
     
      colnames(lipidExistance) <- c('input','simplified input','LION name','LION ID', "match")
      
      
      color_df <- data.frame(match = c("direct","smart matching", "FA-prediction",""),
                             color = c("#4f4d4d","#9c3c2d","#c4942b","#bdbbbf"))
      
      lipidExistance_toview <-
        do.call("rbind", lapply(lipidExistance_list, function(lipidExistance_i) {
          color_pattern <- color_df$color[match(lipidExistance_i$match, color_df$match)]
          
          data.frame(
            input = unique(lipidExistance_i$input),
            name = paste("<font color='",color_pattern,"'>",lipidExistance_i$name,"</font>" , sep ="", collapse = "<br>"),
            #LION = paste("<font color='",color_pattern,"'>",lipidExistance_i$LION,"</font>" , sep ="",  collapse = "<br>")
            LION = paste("<a href='",
                         'http://bioportal.bioontology.org/ontologies/LION/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F',
                         gsub(":","_",lipidExistance_i$LION),
                         "' style='color: ", color_pattern,
                         "' target='_blank'>",
                         lipidExistance_i$LION,
                         "</a>",
                         sep = "", collapse = "<br>")
            
            
          )
        }))
      
      
      colnames(lipidExistance_toview) <- c('input','LION name','LION ID')
      
      lipidExistance_toview[[3]] <- gsub(" href=.+Fnot found'","",lipidExistance_toview[[3]])
      
      ##### end new 20191006
      
      lipidExistance$'simplified input' <- NULL        ## remove this column, not interesting for user
      
      
      ## use matched lipids as assocation table 
    
      lipidExistance_feasable <- lipidExistance[lipidExistance$`LION ID` != "not found",]
      
      matched_lipidID2TERM <- 
        sapply(unique(lipidExistance_feasable$input), function(input){
          lipidExistance_feasable$`LION ID`[lipidExistance_feasable$input == input]
        }, simplify = FALSE)
      
      
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
              list(sublist =  sublist$input,                    #unlist(strsplit(isolate(input$sublist), "\n")),
                   backgroundlist =   background$input          #unlist(strsplit(isolate(input$background), "\n"))   
                   )
            
            
            lipidIDlogical <-
              factor(as.integer(lipidIDs$backgroundlist %in% lipidIDs$sublist), levels = c(0,1))
            names(lipidIDlogical) <- lipidIDs$backgroundlist
            
            if( any(names(lipidIDlogical) %in% names(matched_lipidID2TERM)) &
                                   sum(lipidIDlogical=="1")  > 1  ){   ### can at least one lipid be matched?
              ONTdata <- new(
                "topONTdata",
                ontology = "LION",
                allGenes = lipidIDlogical,
                annot = annFUN.gene2GO,
                gene2GO = matched_lipidID2TERM) 
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
              
              lipidsInTerms <- genesInTerm(ONTdata)
              
              TermContent <- 
              sapply(to_display$TERM.ID, function(id) {
                lipidsInTerms[[id]]
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
        
        if (isolate(input$sub_method == "(ii) analysis" & input$method == "Ranking mode")) {
          if ((is.null(input_listwPvalues) ||
               input_listwPvalues == "") == FALSE) {
            
            
            pValueList$pvalues <- as.numeric(pValueList$pvalues)
            
            lipidIDrank <- rank(pValueList$pvalues * direction)    ## direction by radiobutton input$ranking
            names(lipidIDrank) <- pValueList$input   
            
            mySel <- function(allScore) {
              return(rep(TRUE, length(lipidIDrank)))
            }
            
            
            if(any(names(lipidIDrank) %in% names(matched_lipidID2TERM))){   ### can at least one lipid be matched?
              ONTdata <- new(
                "topONTdata",
                ontology = "LION",
                allGenes = lipidIDrank,
                annot = annFUN.gene2GO,
                gene2GO = matched_lipidID2TERM,
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
                      statistic = isolate(input$ks_sided))
          
            resultFis@score <- abs(resultFis@score)

            ES <- runTest(ONTdata,  algorithm = "classic",
                            statistic = "ks.score")@score
            
            sign_i_df <- data.frame(LION = names(ES), 
                                    sign = sign(ES),
                                    ES = ES)
            
            
            
       
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
              
              lipidsInTerms <- genesInTerm(ONTdata)
              
              TermContent <- 
                sapply(to_display$TERM.ID, function(id) {
                  lipidsInTerms[[id]]
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
            
            if(!is.null(sign_i_df) & isolate(input$ks_sided) == "ks2"){   ### show scores when 2-sided
              to_display$ES <- sign_i_df$ES[match(to_display$`Term ID`, sign_i_df$LION)]
              to_display$Regulated <- as.character(factor(sign_i_df$sign[match(to_display$`Term ID`, sign_i_df$LION)], levels = c(1,-1), labels = c("UP","DOWN")))
              }
            
            
        
          } else {
            to_display <- data.frame()
          }
        }
        
        lengthOfTable <- length(to_display$"Term ID")
        
        ### make detailed table with lipids per term
        
        to_display_detailed <- to_display
        
        lipidsInTerms <- genesInTerm(ONTdata)
        
        identifiers <- lapply(to_display_detailed$'Term ID', function(ID){
          lipidsInTerms[[ID]]
          #genesInTerm(ONTdata)[[ID]]
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
        
        #lipidsInTerms <- genesInTerm(ONTdata)
        
        lipidsInTerms <- data.frame(ID = names(lipidsInTerms),
                                    Discription = LUT$Discription[match(names(lipidsInTerms), LUT$ID)],
                                    nr_lipids = sapply(lipidsInTerms, function(term){length(term)}),
                                    lipids = sapply(lipidsInTerms, function(term){ paste(gsub("^\\[#\\d+\\] ","",term), collapse = '; ')   })
                                   )
        
        lipidsInTerms <- subset(lipidsInTerms, Discription != lipids)
        colnames(lipidsInTerms) <- c("LION-term","LION-name", "nr_lipids", "lipid identifiers")
        
        ### ggplot of data
        to_ggplot_display <- to_display
        

        to_ggplot_display$color <-
          -log(as.numeric(to_ggplot_display$`FDR q-value`), base = 10)
        to_ggplot_display$color[to_ggplot_display$color > 6] <- 6
        
        
        if(!is.null(isolate(input$split_direction_in_output))){
          graph_option <- isolate(input$split_direction_in_output)
        } else {
          graph_option <-"combine"
        }
        
        
        main_title <- substitute(bold("LION enrichment analysis")~
                         italic(x), list(x = tolower(isolate(input$method))))
        
        sub_title <- 
          c(ID = ifelse(
            isolate(input$file_input == "load external dataset") & isolate(input$method) == "Ranking mode",
            paste(isolate(input_data()$studyID)," ", sep =""),
            ""
          ),
          A = ifelse(isolate(input$method) == "Target-list mode","",paste("",isolate(input$conditionA), sep = "")),
          B = ifelse(isolate(input$method) == "Target-list mode","",paste("",isolate(input$conditionB), sep = "")),
          vs  = ifelse(is.null(isolate(input$conditionA)) | is.null(isolate(input$conditionB)) | isolate(input$method) == "Target-list mode",
                       "", "vs."),
          mode = isolate(input$method))
        
        sub_title <- 
          substitute(bold(studyID) ~ A ~ italic(vs)~ B, 
                     list(studyID = sub_title["ID"],
                          A = sub_title["A"],
                          B = sub_title["B"],
                          vs  = sub_title["vs"])
          )
        
        if(!is.null(isolate(input$local_statistics))){
          if(isolate(input$local_statistics)=="3"){   ### in case of a F-test for local statistics
            sub_title <- ""
          }
        }
        
        
        
        if(isolate(input$method == "Ranking mode" & input$ks_sided == "ks2" & graph_option == "split")){
          to_ggplot_display$sign <- sign_i_df$sign[match(to_ggplot_display$`Term ID`,sign_i_df$LION)]
          to_ggplot_display$sign_factor <- factor(ifelse(to_ggplot_display$sign<0,"down","up"), levels = c("up","down"))
          
          to_ggplot_display <- 
            ggplot(data = to_ggplot_display[1:min(lengthOfTable,40), ],
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
            facet_grid(sign_factor~.,  space = "free", scales = "free")+
            labs(title = main_title,      
                 subtitle = sub_title)+
            xlab("") +
            ylab("-LOG(FDR q-value)") +
            guides(fill = "none") +
            coord_flip() +
            theme_pander()
        } else if(isolate(input$method == "Ranking mode" & input$ks_sided == "ks2" & graph_option == "volcano")) {
          ## display as volcano plot
         
          ES <- runTest(ONTdata,
                  algorithm = "classic",
                  statistic = "ks.score")@score
          

          to_ggplot_display$ES <- ES[match(to_ggplot_display$`Term ID`, names(ES))]
          
          to_ggplot_display <- 
          ggplot(data = to_ggplot_display,
                 aes(x = ES, 
                     y = -log(as.numeric(`FDR q-value`), base = 10),
                     size = Annotated)
          ) +
            geom_hline(yintercept = -log(0.05, base = 10),
                       alpha = .8, linetype = 2) +
            geom_vline(xintercept = 0,
                       alpha = .8, linetype = 2) +
            geom_point(shape = 21,aes(fill = color)) +
            geom_text_repel(data = to_ggplot_display[as.numeric(to_ggplot_display$`FDR q-value`) < 0.1,],
                            aes(label = Discription), size = 4)+
            scale_fill_gradient2(
              limits = c(0, 6),
              midpoint = -log(0.05, base = 10),
              low = "grey",
              mid = "grey",
              high = "red"
            ) +
            labs(x ="LION-enrichment score (ES)", y = "-LOG(FDR q-value)", 
                 size = "# of lipids", title =  main_title,
                 subtitle = sub_title)+
            guides(fill = "none") +
            theme_minimal()+
            theme(plot.title = element_text(hjust = 1), plot.subtitle = element_text(hjust = 1))
          
        } else {
          to_ggplot_display <- 
            ggplot(data = to_ggplot_display[1:min(lengthOfTable,40), ],
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
            labs(title = main_title,      
                 subtitle = sub_title)+
            xlab("") +
            ylab("-LOG(FDR q-value)") +
            guides(fill = "none") +
            coord_flip() +
            theme_pander()
        }
          
        
        
        incProgress(.8, detail = paste("enrichment graph"))
        
        ### network
        
        
        if (isolate(input$method) == "Target-list mode") {
          resultTable <-
            GenTable(ONTdata,  'p-value' = resultFis, topNodes = 2000)
        }
        
        if (isolate(input$method) == "Ranking mode") {
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
          
          ### score(resultFis) >>> 0 scores are not tolerated
          resultFis_score <- score(resultFis)
          resultFis_score[resultFis_score==0] <- 0.001
          
          #nr of nodes by FDR qvalue
          enrichmentGraph <-
            showSigOfNodes(
              ONTdata,
              resultFis_score, #k, #abs(resultFis_score), #score(resultFis),
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
            visOptions(highlightNearest = TRUE) %>%
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
        
        
        
        ### now LUT is available, match LION names correctly
        lipidExistance$'LION name' <- LUT$Discription[match(lipidExistance$'LION ID', LUT$ID)]
        lipidExistance$'LION name'[is.na(lipidExistance$'LION name')] <- "not found"
        
        
       
        list(
          ONTdata = ONTdata,
          lipidExistance = lipidExistance,
          lipidExistance_toview = lipidExistance_toview,
          matching_statistics = matching_statistics,
          to_display = to_display,
          to_display_detailed = to_display_detailed,
          lipidReport = lipidReport,
          lipidsInTerms = lipidsInTerms,
          to_ggplot_display = to_ggplot_display,
          edges_nodes = list(edges, nodes),
          to_display_network = to_display_network
          
        )
        
      })
      } else {
        ### empty if there is no input
        list(
          ONTdata = NULL,
          lipidExistance = NULL,
          lipidExistance_toview, NULL,
          matching_statistics = NULL,
          to_display = NULL,
          to_display_detailed = NULL,
          lipidReport = NULL,
          lipidsInTerms = NULL,
          to_ggplot_display = NULL,
          edges_nodes = NULL,
          to_display_network = NULL
          
        )
        
      }  ### if perform when there is input (dim 1 by 1 >> no input)
    
    
  })
  
  
   output$networkcoord <- downloadHandler(
    filename = paste("LION-network-job",isolate(input$submitB)+isolate(input$submitA), ".zip", sep=""),
    content = function(file) {
    
    visNetworkProxy("ontology.plot") %>% visGetPositions()   ### does actually only works the second time
    
    if (!is.null(input$ontology.plot_positions)){
      coords_base <-    do.call(rbind, input$ontology.plot_positions)
      
      edges <- isolate(dataBlock()$edges_nodes[[1]])
      nodes <- isolate(dataBlock()$edges_nodes[[2]])
      
      nodes$shape <- ifelse(nodes$shape == "box", "square", nodes$shape)
      nodes$shape <- ifelse(nodes$shape == "triangle", "square", nodes$shape)
      nodes$shape <- ifelse(nodes$shape == "dot", "circle", nodes$shape)
      
      nodes$size <- nodes$size / 2
      nodes$size[nodes$shape == "square"] <- 10
      
      nodes$label.dist <- sqrt(nodes$size)/3.5
      
      colnames(edges)[match(c("from","to"),colnames(edges))] <- c("to","from")
      edges <- edges[,colnames(edges)[c(2,1,3,4,5)]]
      
      edges <- edges[(edges$from %in% nodes$id) & (edges$to %in% nodes$id),]  
      
      net = graph_from_data_frame(d = edges ,vertices = nodes,directed = T)
      
      layout <- layout_as_tree(net, flip.y = F, root = 1 , rootlevel = nodes$level)[,c(2,1)]
      layout[,1] <- unlist(coords_base[,1]) / 200
      layout[,2] <- unlist(coords_base[,2]) / 200
      
      interactive_network <- isolate(dataBlock()$to_display_network)
      
      save(list = c("net","layout", "edges","nodes","interactive_network"), file = "network.Rdata")
           
      png(paste("LION-network-job",isolate(input$submitB)+isolate(input$submitA), ".png", sep=""), width = 6000, height = 6000, res = 600)
      plot(net, layout = layout,
           vertex.label.cex=.4,
           vertex.label.degree = 1.5*pi, edge.arrow.size = .40,
           vertex.label.color="black")
      dev.off()
      
      svg(paste("LION-network-job",isolate(input$submitB)+isolate(input$submitA), ".svg", sep=""), width = 12, height = 12)
      plot(net, layout = layout,
           vertex.label.cex=.5, 
           vertex.label.degree = 1.5*pi, edge.arrow.size = .40,
           vertex.label.color="black")
      dev.off()
      
      Rfile <- read.delim("data/display_network.R", header = F)
      write.table(file = "display_network.R",Rfile$V1, quote = F, row.names = F, col.names = F)
      
      zip(zipfile=file, files=c( paste("LION-network-job",isolate(input$submitB)+isolate(input$submitA), ".png", sep=""), 
                                 paste("LION-network-job",isolate(input$submitB)+isolate(input$submitA), ".svg", sep=""),
                                 "network.Rdata","display_network.R",include_directories=FALSE)
      )
      
      
      }
 
  })
  #####
  
  
  
  ### mapping percentage
  output$mapping_percentage <- renderText({          
    
    if (is.null(dataBlock()$lipidExistance)){  
      output_text <- "Unexpected input"  
    } else {
     
      matching_statistics <- dataBlock()$matching_statistics
      
      if (length(isolate(dataBlock()$to_display[,1])) == 0){
        comment <- "\nHowever, no upstream LION-terms were found for analysis."
        hideTab(inputId = "tabs", target = "LION enrichment table")
        hideTab(inputId = "tabs", target = "LION enrichment graph")
        hideTab(inputId = "tabs", target = "LION network view")
      } else comment <- ""
      
     
      color_df <- data.frame(match = c("direct matching","smart matching", "fatty acid-association prediction",""),
                             color = c("#4f4d4d","#9c3c2d","#c4942b","#bdbbbf"))
      
       output_text <- paste("<i>",
                           matching_statistics$matched, " out of ",  matching_statistics$total,  " ", "(", round(matching_statistics$percent, digits = 2), "%", ")",
                           " identifiers are matched to LION. </i><br><br>",   
                           paste("<font color='",color_df$color,"'>",color_df$match,"</font>" , sep ="", collapse = "<br>"),
                           sep = "")
    }
    cat("mapping % done\n")
    output_text
    
  })    
      
      
      
      ###
  output$downloadInputCNTRL <- renderUI({
       if(is.null(dataBlock()$lipidExistance)){} else {
         downloadButton("downloadInput", "Download input table")
       }
    })
  
   
  #output$value <- renderTable({
  output$value <- renderFormattable({
    table <- dataBlock()$lipidExistance_toview   #lipidExistance
    formattable(table,align = c("l","l","l"))
  })     #, sanitize.text.function=identity)
      
      output$downloadInput <- downloadHandler(
        
        filename = function() {
          paste("LION-input-job",isolate(input$submitB)+isolate(input$submitA), ".csv", sep="")
          
        },
        content = function(file) {
          write.csv(dataBlock()$lipidExistance, 
                    file, row.names = FALSE,quote = TRUE)
        }
      )
      
   
      output$ontology.table <- renderFormattable({
        table <- dataBlock()$to_display
        rownames(table) <- NULL
        ## add linkds
        table[[1]] <-
        paste("<a href='http://bioportal.bioontology.org/ontologies/LION/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F",
              gsub(":","_",table[[1]]), "' style='color: #4f4d4d' target='_blank'>",
              table[[1]],"</a>",      sep = "")
        
        formattable(table, align = "l")
        
        
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
      
      output$download_pValues <- downloadHandler(
        
        filename = function() {
          paste("LION-rankingvalues-job",isolate(input$submitB)+isolate(input$submitA), ".csv", sep="")
          
        },
        content = function(file) {
          pre_processed_data <- pre_processed_data()[[1]]
          colnames(pre_processed_data) <- c("ID", "value")
          write.csv(pre_processed_data, 
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
                    file = "LION-lipid-associations.csv", row.names = FALSE,quote = TRUE)
          write.csv(dataBlock()$lipidsInTerms, 
                    file = "LION-term-associations.csv", row.names = FALSE,quote = TRUE)
          
          
        
          zip(zipfile=file, files=c( "LION-enrichment-detailed-job.csv", "LION-lipid-associations.csv","LION-term-associations.csv"))
          
          
        }
      )
      

      output$downloadPlotPNG <- downloadHandler(
        filename = function() { 
          paste("LION-enrichment-plot-job",isolate(input$submitB)+isolate(input$submitA),  '.png', sep='') },
        content = function(file) {
          ggsave(file, width = 2.5*5.52, height = 2.5*3.66,
                 dataBlock()$to_ggplot_display+
                   theme(plot.subtitle = element_text(hjust = 1),
                         plot.title = element_text(hjust = 1)),   
                 device = 'png', dpi=300)
        }
      )
      output$downloadPlotSVG <- downloadHandler(
        filename = function() { 
          paste("LION-enrichment-plot-job",isolate(input$submitB)+isolate(input$submitA),  '.svg', sep='') },
        content = function(file) {
          ggsave(file, width = 1.5*5.52, height = 1.5*3.66,
                 dataBlock()$to_ggplot_display+
                   theme(plot.subtitle = element_text(hjust = 1),
                         plot.title = element_text(hjust = 1)),   
                 device = svg,   scale=1.5)
        }
      )
        
      output$ontology.graph <- renderPlot({
         dataBlock()$to_ggplot_display
               
      })
      
      
      ## ontology enrichment tree
      output$ontology.plot <- renderVisNetwork({
        dataBlock()$to_display_network
      })
      
   
      ## email
      observe({
        if(is.null(input$send) || input$send==0) return(NULL)
        from <- isolate(input$from)
        to <- "xxx@xxx.nl"
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
