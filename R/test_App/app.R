################################################################## 
### ignorenet: R Shiny App to generate disease-specific ignoromes
### version 1.0.0
### date: 07.11.17
################################################################## 

# load needed libraries
library(shiny)

#################################
## User Interface
ui <- fluidPage(
  
  # App title ----
  titlePanel("Crappy Appy"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
  
  #=================#
  # INPUTS
  # input1: selector for variable to plot against mpg ----
  # variables must match the columns of the data
  selectInput("variable", "Variable:", #user label
              c("Cylinders" = "cyl",
                "Transmission" = "am",
                "Gears" = "gear")),
  
  # input2: checkbox for whether outliers should be included ----
  checkboxInput("outliers", "Show Outliers", TRUE)
),
  
#===================#
  # Main panel for displaying outputs ----
  mainPanel(
    
    # OUTPUTS
    # output1: plot caption (h2 is the header size) ----
    h3(textOutput("caption")),
    
    # output2: plot ----
    plotOutput("mpgPlot")
    
  )
)

#################################
## Server Implementation

## data preprocessing
# Tweak the "am" variable to have nicer factor labels -- since this
# doesn't rely on any user inputs, we can do this once at startup
# and then use the value throughout the lifetime of the app
source("APP_Preprocessing/mtcars.R")

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  # Compute the formula text for plot caption + plot ----
  # This is in a reactive expression since it is shared by the
  # output$caption and output$mpgPlot functions
  formulaText <- reactive({
    paste("mpg ~", input$variable)
  })
  
  # plot caption ----
  output$caption <- renderText({
    paste("Plot: ", formulaText())
  })
  
  # Generate a plot of the requested variable against mpg ----
  # and only exclude outliers if requested
  output$mpgPlot <- renderPlot({
    boxplot(as.formula(formulaText()),
            data = mpgData,
            outline = input$outliers,
            col = "#75AADB", pch = 19)
  })
  
}


#################################
## Run the application 
shinyApp(ui = ui, server = server)

