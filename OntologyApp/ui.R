require(shiny)
require(shinythemes)
require(visNetwork)
library(shinyTree)
library(shinyWidgets)

##


# Define UI for random distribution application 
fluidPage(theme = shinytheme("cosmo"),
          
          #titlePanel("LION/web | Lipid Ontology enrichment analysis for lipidomics (v2018.02.01)"),
          br(),
          
          sidebarLayout(
            
            
            sidebarPanel( width = 4,
                          tabsetPanel(id = "method", type = "tabs", selected = "Pre-processing",
                                      tabPanel(id = "Settings", title = "", icon = icon("cog", lib = "font-awesome"),
                                               br(),
                                               br(),
                                               materialSwitch(inputId = "RemoveRedundantTerms", 
                                                              value = TRUE,
                                                              label = "Exclude parent terms with same associations as child",
                                                              status = 'primary', right = TRUE),
                                               br(),
                                               materialSwitch(inputId = "LIONselection", 
                                                              label = "Preselect LION-terms for analysis",
                                                              status = 'primary', right = TRUE),
                                               em("Only selected LION-terms will be included downstream enrichment analysis"),
                                               br(),
                                               shinyTree("tree",checkbox = TRUE,search = TRUE),
                                               br()
                                      ), 
                                      
                                      tabPanel("Pre-processing", 
                                               br(),
                                               "Choose CSV File:",
                                               fluidRow(
                                                 column(width = 1,actionButton("CSVhelp", "",icon = icon(name = "question",lib = "font-awesome"))),
                                                 column(offset=0,style='margin-left:2%;' , width = 10,
                                                        fileInput("file1", label = NULL,
                                                             multiple = FALSE,
                                                             accept = c("text/csv",
                                                                        "text/comma-separated-values,text/plain",
                                                                        ".csv")))
                                               ),
                                               
                                               downloadLink("examplePre1", "example set 1 (organelle fractions, adapted from Andreyev AY et al, 2010)"),
                                               br(),
                                               downloadLink("examplePre2", "example set 2 (CHO-k1 incubated with several FFAs)"),
                                               br(),
                                               downloadLink("examplePre3", "example set 3 (CHO-k1 incubated with AA)"),
                                               br(),
                                               br(),
                                               br(),
                                               uiOutput("conditionsUI")),
                                     
                                      
                                      tabPanel("By ranking", 
                                               br(),
                                               textAreaInput(
                                                 "listwPvalues",
                                                 label = h5("lipids with values:"),
                                                 value = "",
                                                 placeholder = 'comma or tab-separated, i.e.;\nPC(16:0/18:1)\t0.7303\nPC(36:1)\t0.4325\nSM(d18:1/20:1),0.7654\nSLM:000088147\t0.1958\nLMGP01010007\t0.8450\n...\netc.',
                                                 resize = "vertical",
                                                 height = 400
                                               ),    
                                               br(),
                                               fluidRow(
                                                 column(width = 1,actionButton("KShelp", "",icon = icon(name = "question",lib = "font-awesome"))),
                                                 column(offset = 0,style='margin-left:2%;' , width = 10,
                                                        actionButton("submitB","  Submit",
                                                                     icon = icon("lightbulb-o", lib = "font-awesome"))
                                                 )),
                                               br(),
                                               br(),
                                               radioButtons("ranking", "Kolmogorov-Smirnov ranking",
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
                                               
                                      ),
                                      
                                      tabPanel("By target list", 
                                               br(),
                                               textAreaInput(
                                                 "sublist",
                                                 label = h5("lipid target list:"),
                                                 value = "",
                                                 placeholder = 'PC(16:0/18:1)\nSM(d18:1/20:1)\n...\netc.',
                                                 resize = "vertical",
                                                 height = 100
                                               ),
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
                                                 column(width = 1,actionButton("FShelp", "",icon = icon(name = "question",lib = "font-awesome"))),
                                                 column(offset = 0,style='margin-left:2%;' , width = 10,
                                                        actionButton("submitA","  Submit",
                                                                     icon = icon("lightbulb-o", lib = "font-awesome"))
                                                 )),
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
                                               )

                        
                          
                          

            )),
            
            mainPanel(
              tabsetPanel(id = "tabs", type = "tabs", 
                          tabPanel("General information", 
                                   h3("About LION/web"),
                                   p('The Lipid Ontology (LION) enrichment analysis web-client (LION/web) is a novel bioinformatics tool for lipidomics that',
                                     'enables users to search for enriched LION-terms in lipidomic subsets. LION-terms contain detailed lipid classification',
                                     'by LIPIDMAPS, biophysical data, lipid functions and organelle associations.'),
                                   

                                   br(),
                            
                                   'To illustrate the workflow, we included a few examples adapted from a published dataset:',
                                 
                                   br(),
                                   em('Subcellular organelle lipidomics in TLR-4-activated macrophages Andreyev AY, Fahy E, Guan Z, Kelly S, Li X, McDonald JG, Milne S, Myers D, Park H, Ryan A, Thompson BM, Wang E, Zhao Y, Brown HA, Merrill AH, Raetz CR, Russell DW, Subramaniam S, Dennis EA. J Lipid Res 51, 2785-97 (2010)  '),
                                   br(),
                                   br(),
                                   "LION/web is currently under review.",
                                   "Any comments, questions or suggestions? Feel free to mail us: m.r.molenaar@uu.nl",br(),
                                   em("Biochemistry & Cell Biology, Universiteit Utrecht, The Netherlands"),
                                   br(),
                                   h3('Example of the LION-structure for PS(34:2):'),
                                   br(),
                                   img(src="LION.png", align = "center", height = '927px', width = '1003px')
                                   
                                   
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
                                   downloadButton("downloadTable", "Download table"),
                                   downloadButton("downloadDetailTable", "Download report"),
                                   br(),
                                   br()
                                   ),
                          
                          tabPanel("LION enrichment graph", 
                                   br(),
                                   plotOutput("ontology.graph", height = 650,
                                               click = "plot_click"),
                                   downloadButton("downloadPlotSVG", "Download plot as SVG"),
                                   downloadButton("downloadPlotPNG", "Download plot as PNG"),
                                   br()
                                   ),
                          
                          tabPanel("LION network view", 
                                
                                   visNetworkOutput("ontology.plot", height = 650),
                                 
                                   img(src="legend.png", align = "right", height = '67px', width = '550px'))
                                   
                                   
              )
            )
          )
)
