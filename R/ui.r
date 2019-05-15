
#this .r file is where we set up the visuals of our shiney app and connect our buttons and widgets to our functions in the server.r file. 
library(shiny)


shinyUI
(
  fluidPage
  (
    #this section is for the title of the webpage
    titlePanel( title = "Look Dr. Uyeda I did a thing"),
    #this is where things displayed in the sidebar should go 
    sidebarLayout
    (
       sidebarPanel(
         
         #this for the input boxes on the shiney app the first section is the variale name under which the input will be stored the second is
         #the text that will be displayed above the input box and the third is what text will be displayed in the box 
         textInput("taxonName:", "Taxon Name", "Enter Name Here"),
                    verbatimTextOutput("value"),
         
         #this is the same as the above section the only difference is this is for the entity name and not the taxon name
         textInput("entityName:", "Entity Name", "Enter Name Here"),
         verbatimTextOutput("value"), 
         
         #this is a button to make it do stuff, it does not in fact yet do anything....
         actionButton("Submit", "Click The Button To Make Plot")
         ), 
      # this is the main pannel and where the output will go
      mainPanel( plotOutput( plotData(td, njt)))
    )
  )
)


