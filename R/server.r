library(shiny)

# this allows us to use the functions we developed in the functions.R file which came from the origional semantic similarity function mark up file
source("functions.R")


shinyServer
(
  #this is how we take in input and output
  function(input,output)
  {
    # this is me trying to make the button work
    observeEvent
    (
      input$Submit, 
      {
        textOutput("Thank you for clicking")
    
      }
    )
    #these are calls to our functions from the functions.R file they aren't set up to do anything yet but they are here 
    td <- Get_Tree_Data(input$taxonName, input$entityName)
    njt <- makeTree(td)
    plotData(td, njt)
    # me testing to see if the input text box was working
    output$value <- renderText({ input$taxonName })
  }
)