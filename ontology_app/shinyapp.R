#
# This app returns a plot, based on user input, of traits vs taxa data from the phenoscape database

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
  titlePanel("Ontology trait explorer"),
  
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
    paste("Fetching plot for ", input$taxon, "and ", input$trait)
  })
  
  output$confirm <- renderText({
    confirm_t()
  })
  
  # confirm # of results 
  confirm_r <- eventReactive(input$plot, {
    paste(length(td$phy), " taxa and ", length(td$dat), " traits returned")
  })
  
  output$results <- renderText({
    confirm_r()
  })
  
  # fetch plot based on user input
  fplot <- eventReactive(input$plot, {
    td <- Get_Tree_Data(input$taxon, input$trait)
    td <- filter_coverage(td, traits=input$f_traits, taxa=input$f_taxa)
    
    # issues a warning if we can't make a tree with the given # of traits and taxa
    if ((ncol(td$dat) <= 2) || length(td$phy$tip.label) <= 2){
      stop("Not enough traits or taxa to make tree")
    }
    
      njt <- makeTree(td)
      plotData(td, njt, show.tip.label=TRUE, cex=0.25)
    
  })
  
  output$graph <- renderPlot({
    fplot()
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)

