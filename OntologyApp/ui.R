require(shiny)
require(shinythemes)
require(visNetwork)
library(shinyTree)
library(shinyWidgets)
library(shinyBS)
library(formattable)
library(shinycssloaders)


##


# Define UI for random distribution application 
#fluidPage(theme = shinytheme("cosmo"),
fluidPage(theme = shinytheme("journal"),
          
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
                                                   inputId = "SmartMatching",
                                                   value = TRUE,
                                                   label = "Use 'smartmatching' (beta)",
                                                   status = 'primary',
                                                   right = TRUE
                                                 ),
                                               content = "If possible, unmatched lipids are associated with related LION-terms: `TG(O-16:0/18:2/22:6)`	is associated with `alkyldiacylglycerols`, `C16:0`, `C18:2`, and `C22:6`"
                                               ), 
                                               br(),
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                      materialSwitch(
                                                        inputId = "FAprediction",
                                                        value = FALSE,
                                                        label = "Predict fatty acid assocations (beta)",
                                                        status = 'primary',
                                                        right = TRUE
                                                      ),
                                                      content = "For sum-formatted phospholipids, e.g. PC(34:1), predict the most likely fatty acid assocations, based on reported fatty acid compositions. For example, PC(34:1) will be associated to C16:0 and C18:1. Only use for datasets of mammalian origin"
                                               ), 
                                               br(),
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                 materialSwitch(
                                                   inputId = "EmailMissingLipids",
                                                   label = "Automatically send unmatched lipids to LION-team",
                                                   status = 'primary',
                                                   right = TRUE
                                                 ),
                                                 content = "Names of lipids that cannot be matched to LION are automatically send to the LION-team. This will help us to improve the coverage of the database. User information and associated data will NOT be sent."
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
                                                                    #br(),
                                                                    prettyRadioButtons(
                                                                      inputId = "file_input",
                                                                      label = "",
                                                                      choices = c("file input", 
                                                                                  "load external dataset"),
                                                                      inline = TRUE, 
                                                                      status = "danger",
                                                                      fill = TRUE
                                                                    ),
                                                                    withSpinner(color = "#9c3c2d", size = .5,
                                                                      uiOutput("load_datasetUI"))
                                                                    ),
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
                                                                    
                                                                    bsCollapse(id = "KSsettings",
                                                                               bsCollapsePanel("K-S settings", 
                                                                                                      radioButtons("ranking", "ranking direction",
                                                                                                                     c("from low to high" = "ascending",
                                                                                                                       "from high to low" = "descending")),
                                                                                                      radioButtons("ks_sided", "alternative hypothesis",
                                                                                                                     c("one-tailed (ECDF is higher)" = "ks",
                                                                                                                       "two-tailed" = "ks2")),
                                                                                               uiOutput("KS2_UI")
                                                                               )
                                                                    ),
                                                                    fluidRow(
                                                                      popify(placement = "top", title = "Info",
                                                                             column(offset = 0,#style='margin-left:2%;' , 
                                                                                    width = 12, 
                                                                                    actionButton("submitB","  Submit",
                                                                                                 icon = icon("lightbulb-o", lib = "font-awesome"))
                                                                             ),content = 'After clicking "submit", LION/web matches the input lipids to the LION-database. Only matched identifiers will be used in the enrichment analysis. Subsequently, LION/web ranks the input lipids by the provided numeric value. The distributions of all LION-terms associated with the dataset are compared to uniform distributions and evaluated by Kolmogorov-Smirnov-tests. A low p-value is returned when lipids associated with a certain LION-term are higher on top of the ranked list than one would expect by chance.')), 
                                                                    br(),
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
                                               actionLink("exampleA1", "example 1 (cluster 6 in Fig. 2A*)"),
                                               br(),
                                               actionLink("exampleA2", "example 2 (cluster 7 in Fig. 2A*)"),
                                               br(),
                                               actionLink("exampleA3", "example 3 (cluster 8 in Fig. 2A*)"),
                                               br(),
                                               br(),
                                               "* from Molenaar MR, et al., 2019"
                                               ),
                                      tabPanel(id = "Contact", title = "", icon = icon("envelope-open-text", lib = "font-awesome"),
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
                                   
                                   h3("LION/web: LION enrichment analysis"),
                                   br(),
                                   img(src="LIONicon enrich.png", align = "left", height = '200px'), #, width = '550px'),
                                   br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                                   p('The Lipid Ontology (LION) enrichment analysis web application (LION/web) is a novel bioinformatics tool for lipidomics that',
                                     'enables users to search for enriched LION-terms in lipidomic subsets. LION-terms contain detailed lipid classification',
                                     'by LIPIDMAPS, biophysical data, lipid functions and organelle associations.'),
                                   p("Any comments, questions or suggestions? Do you want to add lipid annotations? Please use our contact form."),
                                   br(),
                                   p("Please cite:"),
                                   a('LION/web: a web-based ontology enrichment tool for lipidomic data analysis.',
                                     href="https://doi.org/10.1093/gigascience/giz061", target="_blank"),br(),
                                   tags$b("Gigascience. 2019 Jun 1;8(6). pii: giz061. doi: 10.1093/gigascience/giz061."),
                                   br(),
                                   em("Martijn R. Molenaar, Aike Jeucken, Tsjerk A. Wassenaar, Chris H. A. van de Lest, Jos F. Brouwers, J. Bernd Helms."), 
                                   br(),
                                   br(),
                                   br(),
                                   p("Division Cell Biology, Metabolism & Cancer, Department Biomolecular Health Sciences, Universiteit Utrecht, The Netherlands"),
                                   'v. 2020.07.14',
                                   br(),
                                   
                                   br(),
                                   
                                   br()
                                   
                                   ),
                          tabPanel("LION input", 
                                   br(),
                                   #em(   
                                   htmlOutput("mapping_percentage"),    #textOutput("mapping_percentage")  ),
                                   br(),
                                   formattableOutput("value"), 
                                   uiOutput("downloadInputCNTRL"),
                                   
                                   br()
                                   
                                   
                          ),
                          tabPanel("LION enrichment table", 
                                   br(),
                                   #tableOutput("ontology.table"),
                                   formattableOutput("ontology.table"),
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
                                   img(src="legend.png", align = "right", height = '67px', width = '550px'),
                                   br(),br(),br(),br(),
                                   bsCollapse(id = "download_network_bs",
                                              bsCollapsePanel(title = "Download network", 
                                                              column(2, downloadButton("networkcoord", "Download as ZIP")))
                                   ))
                                   
                                   
              )
            )
          )
)
