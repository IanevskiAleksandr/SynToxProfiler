library(shiny); library(plotly);  library(shinythemes); library(shinyjs); library(shinytoastr)
shinyUI(navbarPage(useShinyjs(), useToastr(),
  
  title=div(class="myheader", img(class = "imgLogo", src="logo2.png", height="60", width="260"),

            div(class="reference", tags$a(href="http://www.ncbi.nlm.nih.gov/pubmed/26949479", target="_blank",
                                          tags$p("Reference"))
             ),
            div(class="howto", actionButton('howtoin','How to use', icon = icon("file", lib = "glyphicon"), class = 'btn-primary', onclick="window.open('http://syntoxprofiler.fimm.fi/howto/')")
            ),
            div(class="author", icon("user", lib = "glyphicon"))
  ),
  
  tags$head(
    HTML('<link rel="shortcut icon" href="logo.png">'),
    HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Lora" />'),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"), 
    tags$script('$(document).on("shiny:connected",function(n){var e=screen.width;e2=screen.height;Shiny.onInputChange("GetScreenWidth",e);Shiny.onInputChange("GetScreenHeight",e2);});')
  ),
  br(),
  theme = shinytheme("journal"),
 
  sidebarLayout(
    sidebarPanel(
      br(),
      # Input file
      fluidRow(column(width = 8, fileInput('annotfile', 'Annotation file', accept = c('.csv', '.xlsx', '.txt'))),
               column(width = 4, br(),downloadButton('downlEx','Example data', class = 'btn-primary'))
      ),
      
      radioButtons("conrol_", "Toxicity calculations:",
                   list("From control data" = "Control",
                        "PrOCTOR based approximation" = "PrOCTOR",
                        "No toxicity (rank combinations only based on synergy and efficacy)" = "SynTox")),
      br(),
      
      fluidRow(column(width = 8, checkboxInput("synfin", "Calculate synergy scores using SynergyFinder", !0)),
               column(width = 4, conditionalPanel(
                                        condition = "input.synfin == 1",  selectInput("state", "Choose a synergy scoring model:",
                                                                                      list("ZIP","Bliss","Loewe","HSA")
                                  ))
                      )
      ),
      
      br(),
      
      helpText("SynToxProfiler web application helps to prioritize drug combinations based on integrated synergy, efficacy and toxicity profiles. For more information, on how to use, please read a technical documentation and watch a video tutorial."),
      br(),
      
      fluidRow(
        column(offset = 4, width = 4, actionButton("analysis_", "Start analysis", icon = icon("play"),class = "btn-primary"))
      )
      
      
    ),
    mainPanel(
      tabsetPanel(type = "tabs", id = 'tabs',
                  
                  tabPanel("Summary", 
                           br(), #tags$video(id="videoIntro", type = "video/mp4",src = "video.mp4", width="95%", controls = "controls"),
                           uiOutput("video"),
                           fluidRow(column(8,offset=1,
                                           plotOutput("barplot", hover = hoverOpts("plot_hover2", delay = 250, delayType = "throttle"), width = "95%"), uiOutput("hover_info2")),
                                    div(id = "col1id", column(2, offset = 0,
                                           br(),br(),br(),
                                           
                                           fluidRow(
                                             radioButtons("icons2", "Ranking (based on):",
                                                          choices =
                                                            list("synergy, toxicity, efficacy", "synergy, toxicity",
                                                                 "synergy, efficacy", "toxicity, efficacy", "synergy", "efficacy", "toxicity"),
                                                          selected = "synergy, toxicity, efficacy"
                                             )
                                           ),
                                           fluidRow(
                                             downloadButton("exp_excel", "Export summary table")
                                           ), br(),                                          
                                           fluidRow(
                                             downloadButton("exp_2d2", "Export figure")
                                           )
                                    ))
                            )
                           
                  ),
                  tabPanel("Combined landscape", uiOutput("plotlcont"), uiOutput("hover_info3")),
                  tabPanel("Individual plots",
                           
                           
                           fluidRow( br(),br(),br(),
                             column(7, offset = 1,
                                    plotOutput("scatterplot", 
                                               hover = hoverOpts("plot_hover", delay = 250, delayType = "throttle")), uiOutput("hover_info")),
                             column(2, offset = 1,
                                    br(),br(),br(),
                                    fluidRow(
                                      radioButtons("icons", "Plot type:",
                                                   choices =
                                                           list("synergy vs toxicity", "synergy vs inhibition",
                                                                "inhibition vs toxicity"),
                                                         selected = "synergy vs toxicity"
                                      )
                                    ),
                                    fluidRow(
                                      downloadButton("exp_2d", "Export figure")
                                    )
                             )
                            ),
                           br()
                           ), 
       selected = "Summary"
      ),
      tags$div(
        HTML("<img src = 'footer.png' class = 'fixFoot' ></img>"),
        HTML('<script type="text/javascript"> window._urq=window._urq||[];_urq.push(["initSite","ff2cea18-f1c2-4241-b463-184ef777939b"]);(function(){var a=document.createElement("script");a.type="text/javascript";a.async=!0;a.src="https:"==document.location.protocol?"https://cdn.userreport.com/userreport.js":"http://cdn.userreport.com/userreport.js";var b=document.getElementsByTagName("script")[0];b.parentNode.insertBefore(a,b)})();</script>')
      )
    )
  )
))
