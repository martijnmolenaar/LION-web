require(shiny)
require(shinythemes)
require(visNetwork)
library(shinyTree)
library(shinyWidgets)
library(shinyBS)

##


# Define UI for random distribution application 
fluidPage(theme = shinytheme("cosmo"),
          
          # Application title
          #titlePanel("LION/web | Lipid Ontology enrichment analysis for lipidomics"),
          br(),
          
          sidebarLayout(
            
            
            sidebarPanel( width = 4,
                          
                          tabsetPanel(id = "method", type = "tabs", selected = "Ranking mode",
                          
                                      tabPanel(id = "Settings", title = "", icon = icon("cog", lib = "font-awesome"),
                                               br(),
                                               br(),
                                               
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                 materialSwitch(
                                                   inputId = "RemoveRedundantTerms",
                                                   value = TRUE,
                                                   label = "Exclude parent terms with same associations as child",
                                                   status = 'primary',
                                                   right = TRUE
                                                 ),
                                                 content = 'When a LION-term contains the same lipids as its parent (for instance, glycerophosphocholines and diacylglycerophosphocholines), the parent term (here glycerophosphocholines) will be excluded in analysis. '
                                               ), 
                                               br(),
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                 materialSwitch(
                                                   inputId = "EmailMissingLipids",
                                                   label = "Automatically send unmatched lipids to LION-team",
                                                   status = 'primary',
                                                   right = TRUE
                                                 ),
                                                 content = "Names of lipids that could not be matched to LION are automatically send to the LION-team. This will help us to improve the coverage of the database. User information and associated data will NOT be sent."
                                               ), 
                                               br(),
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                 materialSwitch(
                                                   inputId = "LIONselection",
                                                   label = "Preselect LION-terms for analysis",
                                                   status = 'primary',
                                                   right = TRUE
                                                 ),
                                                 content = "By default, all LION-terms will be used in the analysis. With this option, you can pre-select LION-terms to use in the enrichment analysis."
                                               ), 
                                               em(""),
                                               shinyTree("tree",checkbox = TRUE,search = TRUE),
                                               br()
                                      ), 
                                      
                                      tabPanel("Ranking mode", 
                                               br(),
                                               tabsetPanel(id = "sub_method", type = "pills", selected = "(i) process input",
                                                           tabPanel("(i) process input", 
                                                                    br(),
                                                                    "Choose CSV File:",
                                                                    fluidRow(
                                                                      
                                                                      column(offset=0,#style='margin-left:2%;' ,
                                                                             width = 10,
                                                                             popify(placement = "bottom", title = "File-input info",
                                                                             fileInput("file1", label = NULL,
                                                                                       multiple = FALSE,
                                                                                       accept = c("text/csv",
                                                                                                  "text/comma-separated-values,text/plain",
                                                                                                  ".csv") ),
                                                                             content = 'Format your dataset as comma seperated value files (.csv), with the first row reserved for metabolites and the other columns for numeric data (containing decimal points). Use double column headers; with row 1 containing condition identifiers and row 2 containing sample identifiers. Submit at least duplicates per condition. Dataset should be normalized before submission. Download a dataset below for an example.'
                                                                      
                                                                    )),
                                                                    column(offset=0, width = 1,style = "margin-top: 5px;",align="center",
                                                                    popify(placement = "right", title = "Lipid nomenclature", options=list(container="body"),
                                                                           el = icon(name = "question",lib = "font-awesome", "fa-2x"),
                                                                           content = 'Format lipids in LIPIDMAPS notation style: a class-prefix followed by (summed) fatty acid(s) surrounded by parentheses. Examples are: PC(32:1); PE(18:1/16:0); SM(d18:1/18:0); TAG(54:2); etc. Check www.lipidmaps.org for more examples. LION will try to reformat alternative notation styles into LIPIDMAPS format.'))),
                                                                    
                                                                    downloadLink("examplePre1", "example set 1 (organelle fractions, adapted from Andreyev AY et al, 2010)"),
                                                                    
                                                                    br(),
                                                                    downloadLink("examplePre2", "example set 2 (CHO-k1 incubated with several FFAs)"),
                                                                    br(),
                                                                    downloadLink("examplePre3", "example set 3 (CHO-k1 incubated with AA)"),
                                                                    br(),
                                                                    br(),
                                                                    uiOutput("selectLocalStatisticUI")),
                                                           tabPanel("(ii) analysis", 
                                                                    br(),
                                                                    popify(
                                                                      placement = "bottom",
                                                                      title = "Info",
                                                                      textAreaInput(
                                                                        "listwPvalues",
                                                                        label = h5("lipids with values:"),
                                                                        value = "",
                                                                        placeholder = 'comma or tab-separated, i.e.;\nPC(16:0/18:1)\t0.7303\nPC(36:1)\t0.4325\nSM(d18:1/20:1),0.7654\nSLM:000088147\t0.1958\nLMGP01010007\t0.8450\n...\netc.',
                                                                        resize = "vertical",
                                                                        height = 400
                                                                      ),
                                                                      content = 'Provide a list of metabolites with numeric values (comma or tab seperated) to rank the input. Select the correct ranking directions with the option "ranking direction" below.'
                                                                    ),
                                                                    br(), 
                                                                    fluidRow(
                                                                      #column(width = 1,actionButton("KShelp", "",icon = icon(name = "question",lib = "font-awesome"))),
                                                                      popify(placement = "top", title = "Info",
                                                                      column(offset = 0,#style='margin-left:2%;' , 
                                                                             width = 12, 
                                                                             actionButton("submitB","  Submit",
                                                                                          icon = icon("lightbulb-o", lib = "font-awesome"))
                                                                      ),content = 'After clicking "submit", LION/web matches the input lipids to the LION-database. Only matched identifiers will be used in the enrichment analysis. Subsequently, LION/web ranks the input lipids by the provided numeric value. The distributions of all LION-terms associated with the dataset are compared to uniform distributions and evaluated by Kolmogorov-Smirnov-tests. A low p-value is returned when lipids associated with a certain LION-term are higher on top of the ranked list than one would expect by chance.')), 
                                                                    br(),
                                                                    br(),
                                                                    radioButtons("ranking", "ranking direction",
                                                                                 c("from low to high" = "ascending",
                                                                                   "from high to low" = "descending")),
                                                                    br(),
                                                                    actionLink("exampleB1", "example 1 (plasma membrane vs ER fraction)*"),
                                                                    br(),
                                                                    actionLink("exampleB2", "example 2 (mitochondrial fraction vs homogenate)*"),
                                                                    br(),
                                                                    actionLink("exampleB3", "example 3 (stimulated vs non-stimulated ER fractions)*"),
                                                                    br(),
                                                                    br(),
                                                                    "* adapted from Andreyev AY, et al., 2010"
                                                                    )
                                                           )),
                                          tabPanel("Target-list mode", 
                                               br(),
                                               popify(placement = "bottom", title = "Info",
                                               textAreaInput(
                                                 "sublist",
                                                 label = h5("lipid target list:"),
                                                 value = "",
                                                 placeholder = 'PC(16:0/18:1)\nSM(d18:1/20:1)\n...\netc.',
                                                 resize = "vertical",
                                                 height = 100
                                               ), content = 'Provide two lists of lipids; a "hitlist" containing lipids that are evaluated for enrichment (lipids associated with a condition of interest, a sublist obtained by data-clustering, etc.), and a background set containing all the lipids in the experiment.'),
                                               textAreaInput(
                                                 "background",
                                                 label = h5("lipid background list:"),
                                                 value = "",
                                                 placeholder = 'PC(16:0/18:1)\nPC(36:1)\nSM(d18:1/20:1)\nSLM:000088147\nLMGP01010007\n...\netc.',
                                                 resize = "vertical",
                                                 height = 250
                                               ),
                                               br(),
                                               fluidRow(
                                                 #column(width = 1,actionButton("FShelp", "",icon = icon(name = "question",lib = "font-awesome"))),
                                                 popify(placement = "top", title = "Info",
                                                 column(offset = 0,#style='margin-left:2%;' , 
                                                        width = 12,
                                                        actionButton("submitA","  Submit",
                                                                     icon = icon("lightbulb-o", lib = "font-awesome"))
                                                 ), content = 'After clicking "submit", LION/web matches the input lipids to the LION-database. Only matched identifiers will be used in the enrichment analysis. Subsequently, LION/web calculates a score of every LION-term that is associated with the dataset, for both the hitlist as the background set. These scores are compared by Fisher-exact tests. A low p-value is returned when a certain LION-term is higher represented in the hitlist than one would expect by chance.')),
                                               br(),
                                               br(),
                                               actionLink("exampleA1", "example 1 (plasma membrane vs ER fraction)*"),
                                               br(),
                                               actionLink("exampleA2", "example 2 (mitochondrial fraction vs homogenate)*"),
                                               br(),
                                               actionLink("exampleA3", "example 3 (stimulated vs non-stimulated ER fractions)*"),
                                               br(),
                                               br(),
                                               "* adapted from Andreyev AY, et al., 2010"
                                               ),
                                      tabPanel("Contact", 
                                               br(),
                                               br(),
                                               'Questions, feedback or a request to increase coverage? Contact us:',
                                               br(),
                                               br(),
                                               textAreaInput("from", label = "From:", value="", resize = "none", rows = 1, placeholder = 'your_email@gmail.com'),
                                               selectizeInput("subject", label = "Subject:", selected = "Question", options = list(create = TRUE),
                                                           choices = c("Question","Feedback","Request to add lipids to LION","")),
                                               
                                               textAreaInput("message", label = "Message:",value = "", resize = "vertical",
                                                               placeholder ="Write message here. Please include a list of identifiers if you want to request lipids to add."),
                                               br(),
                                               actionButton("send", "Send mail")
                                               
                                               
                                               
                                      )

                        
                          
                          

            )),
            
            mainPanel(
              tabsetPanel(id = "tabs", type = "tabs", 
                          tabPanel("General information", 
                                   #br(),
                                   #tableOutput("value"),
                                   #hr(),
                                   h3("About LION/web"),
                                   p('The Lipid Ontology (LION) enrichment analysis web application (LION/web) is a novel bioinformatics tool for lipidomics that',
                                     'enables users to search for enriched LION-terms in lipidomic subsets. LION-terms contain detailed lipid classification',
                                     'by LIPIDMAPS, biophysical data, lipid functions and organelle associations.'),
                                   br(),
                                  
                                   "LION/web is currently under review.",
                                   "Any comments, questions or suggestions? Do you want to add lipid annotations? Please use our contact form.",br(),
                                   br(),
                                   "A pre-print about LION/web is available on BioRxiv:",
                                   br(),
                                   a('LION/web: a web-based ontology enrichment tool for lipidomic data analysis.',
                                     href="https://doi.org/10.1101/398040", target="_blank"),br(),
                                   em("Martijn R. Molenaar, Aike Jeucken, Tsjerk A. Wassenaar, Chris H. A. van de Lest, Jos F. Brouwers, J. Bernd Helms."), "BioRxiv (2018)",
                                   br(),
                                   br(),
                                   em("Biochemistry & Cell Biology, Universiteit Utrecht, The Netherlands"),
                                   br(),
                                   em('LION/web v. 2019.01.10'),
                                   br(),
                                   
                                   br(),
                                   
                                   
                                   #h3('Example of the LION-structure for PS(34:2):'),
                                   br()
                                   #img(src="LION.png", align = "center", height = '927px', width = '1003px')
                                   
                                   
                                   ),
                          tabPanel("LION input", 
                                   br(),
                                   em(   textOutput("mapping_percentage")  ),
                                   br(),
                                   tableOutput("value"),
                                   uiOutput("downloadInputCNTRL"),
                                   
                                   br()
                                   
                                   
                          ),
                          tabPanel("LION enrichment table", 
                                   br(),
                                   tableOutput("ontology.table"),
                                   fluidRow(
                                     column(
                                       offset = 0,
                                       width = 1,
                                       style = "margin-top: 5px;",
                                       align = "right",
                                       popify(
                                         placement = "top",
                                         title = "LION enrichment",
                                         options = list(container = "body"),
                                         el = icon(name = "question", lib = "font-awesome", "fa-2x"),
                                         content = 'This table includes all LION-terms that are associated with the provided dataset and is sorted by enrichment. N.B.: enrichment of terms related to membrane fluidity and bilayer thickness is based on phospholipid composition. Other factors that are not taken into account (including protein composition and cholesterol concentration) might affect these properties.'
                                       )
                                     ),
                                     column(offset = 0,width = 10,
                                            downloadButton("downloadTable", "Download table"),
                                            downloadButton("downloadDetailTable", "Download report")
                                            )),
                                   br(),
                                   br()
                                   ),
                          
                          tabPanel("LION enrichment graph", 
                                   br(),
                                   plotOutput("ontology.graph", height = 650
                                               #click = "plot_click"
                                               ),
                                   fluidRow(
                                     column(
                                       offset = 0,
                                       width = 1,
                                       style = "margin-top: 5px;",
                                       align = "right",
                                       popify(
                                         placement = "top",
                                         title = "LION enrichment",
                                         options = list(container = "body"),
                                         el = icon(name = "question", lib = "font-awesome", "fa-2x"),
                                         content = 'This graph represents the 40 most enriched LION-terms. N.B.: enrichment of terms related to membrane fluidity and bilayer thickness is based on phospholipid composition. Other factors that are not taken into account (including protein composition and cholesterol concentration) might affect these properties.'
                                       )
                                     ),
                                     column(offset = 0,width = 10,
                                        downloadButton("downloadPlotSVG", "Download plot as SVG"), 
                                        downloadButton("downloadPlotPNG", "Download plot as PNG"))
                                    
                                   ),
                                   
                                   br()
                                   ),
                          
                          tabPanel("LION network view", 
                                   visNetworkOutput("ontology.plot", height = 650),
                                   img(src="legend.png", align = "right", height = '67px', width = '550px'))
                                   
                                   
              )
            )
          )
)
