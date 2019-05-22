#
#

library(shiny)
library(rphenoscape)
library(treeplyr)
library(RCurl)
library(rjson)
library(readr)
library(urltools)
library(pracma)
#library(httr)
library(devtools)
library(phytools)
library(phylogram)
library(geiger)

# import functions 
source("functions.R")

# Define UI ----
ui <- fluidPage(
  titlePanel("title"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("taxon", h4("Taxon"), 
                       value = "Characiformes"),
      textInput("trait", h4("Trait"), 
                value = "dermal bone"),
      # filter coverage for traits
      sliderInput("f_traits", h4("Traits filter coverage"),
                  min = 0, max = 0.5, value = 0),
      # filter coverage for taxa
      sliderInput("f_taxa", h4("Taxa filter coverage"),
                  min = 0, max = 0.5, value = 0),
      actionButton("plot", "Fetch plot")
    ),
    
    mainPanel(
      ### test
      textOutput("confirm"),
      
      textOutput("results"),
    
      plotOutput("graph")
    )
  )
  
)

# Define server logic ----
server <- function(input, output) {
  # confirm input values
  confirm_t <- eventReactive(input$plot, {
    paste("Fetching plot for ", input$taxon, "and ", input$trait, "\n;", 
          input$f_traits, " traits filtered and ", input$f_taxa, " taxa filtered")
  })
  
  output$confirm <- renderText({
    confirm_t()
  })
  
  # confirm # of results 
  confirm_r <- eventReactive(input$plot, {
    paste("# taxa and # traits returned")
  })
  
  output$results <- renderText({
    confirm_r()
  })
  
  # fetch plot based on user input
  fplot <- eventReactive(input$plot, {
     td <- Get_Tree_Data(input$taxon, input$trait)
     td <- filter_coverage(td, input$f_traits, input$f_taxa)
     njt <- makeTree(td)
     paste("# taxa and # traits returned")
     plotData(td, njt, show.tip.label=TRUE, cex=0.25)
  })
  
  output$graph <- renderPlot({
    fplot()
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)

