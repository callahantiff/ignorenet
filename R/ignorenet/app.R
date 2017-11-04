################################################################## 
### ignorenet: R Shiny App to generate disease-specific ignoromes
### version 1.0.0
### date: 07.11.17
################################################################## 

## LIBRARIES
library(shiny)
library(shinydashboard)
library(rhandsontable)


## PREPROCESSING
# source preprocessing file
source('GEO_Preprocessing.R')

#################################
## User Interface
ui <- dashboardPage(
  
  # header
  dashboardHeader(title = "ignorenet"),
  
  # sidebar
  dashboardSidebar(
    
    # side menu
    sidebarMenu(
      # menuItem("Home", tabName = "home"),
      menuItem("Query GEO", tabName = "dashboard"),
      menuItem("Group Selection", tabName = "samp"),
      menuItem("Differential Expression Analysis", tabName = "dea"),
      menuItem("ignorenet", tabName = "net")
    )
    
),
  
  # body
  dashboardBody(
    tags$head(tags$style(HTML('
                          /* logo */
                              .skin-blue .main-header .logo {
                              background-color: white; color:        rgb(0,144,197);
                              font-weight: bold;font-size: 24px;text-align: Right;
                              }
                              
                              /* logo when hovered */
                              
                              .skin-blue .main-header .logo:hover {
                              background-color: white;
                              }
                              
                              
                              /* navbar (rest of the header) */
                              .skin-blue .main-header .navbar {
                              background-color: white;
                              }
                              
                              /* main sidebar */
                              .skin-blue .main-sidebar {
                              background-color: gray;;
                              }
                              
                              /* active selected tab in the sidebarmenu */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: rgb(107,194,0);
                              color: white;font-weight: bold;;
                              }
                              
                              /* other links in the sidebarmenu */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: gray;
                              color: white;font-weight: bold;
                              }
                              
                              /* other links in the sidebarmenu when hovered */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: rgb(232,245,251);color: rgb(0,144,197);font-weight: bold;
                              }
                              
                              /* toggle button color  */
                              .skin-blue .main-header .navbar .sidebar-toggle{
                              background-color: white;color:rgb(0,144,197);
                              }
                              
                              /* toggle button when hovered  */
                              .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: rgb(0,144,197);color:white;
                              }
                              
                              # '))),
    
    # study tabs
    tabItems(
    
    # home
    # tabItem("home",
    #         box(
    #           width = 10, status = "info", solidHeader = TRUE,
    #           title = "HOME PAGE NEEDS TEXT")),
    
    # GEO Study Search
    tabItem("dashboard",
            fluidRow(
              
              # selected identifiers
              box(
                width = 4, status = "info", solidHeader = TRUE,
                title = "GEO Query Keywords",
                shiny::textInput("disease", label="Disease"),
                
                shiny::selectizeInput(
                  inputId = "org", 
                  label = "Select Organism",
                  multiple  = FALSE,
                  choices = GEOChoices(),
                  options = list(placeholder = 'Select an Organism')),
                shiny::submitButton(text="Find Studies")),
              
              #run GEO query
              box(
                width = 4, status = "info", solidHeader = TRUE,
                title = "Selected GSE Accession Identifiers",
                tableOutput('rows_out'),
                submitButton(text="Select GSE Identifiers")),
              
              # GEO results
              box(
                width = 10, status = "info", solidHeader = TRUE,
                title = "GEO Results",
                dataTableOutput('tbTable'))
            )),
      
    # experimental sample selection
    tabItem("samp",
            box(
              shiny::selectizeInput(
              inputId = "gse", 
              label = "Select A GSE Accession Identifier",
              multiple  = FALSE,
              choices = 'rows_out',
              options = list(placeholder = 'Select a Study')),
            
            # shiny::submitButton(text="Next"),
            shiny::submitButton(text="Analyze")
            ),
            
            box(
              width = 10, status = "info", solidHeader = TRUE,
              title = "GSE Study Metadata",
              rHandsontableOutput("hot")
              
            )),
    
    # differential expression analysis
    tabItem("dea",
            box(
              submitButton(text="ignorenet"),
              br(),
              br(),
              width = 15, status = "info", solidHeader = TRUE,
              title = "Significantly Differentially Expressed Genes",
              dataTableOutput('DEQtable')
              ))
)))


#################################
## Server Implementation
server <- function(input, output, session) {
  values <- reactiveValues()
  
  ## TAB1: Querying GEO
  # drop-down menu
  observe({
    updateSelectizeInput(session, 'org', choices = GEOChoices())
    updateSelectizeInput(session, 'gse', choices = select_rows())
  })
  
  # generate query output from inital text
  data <- reactive({
    GEOQuery(input$disease, input$org)})
  
  output$tbTable <- renderDataTable(
    data(),
    options = list(pageLength = 10),
    escape = FALSE,
    callback = "function(table) {
    table.on('click.dt', 'tr', function() {
    $(this).toggleClass('selected');
    Shiny.onInputChange('rows',
    table.rows('.selected').indexes().toArray());
    });}")
  
  # generate list of choosen studies
  select_rows <- reactive({
    new_row = input$rows + 1
    data()[new_row, ][1]})

  output$rows_out <- renderTable(
    select_rows())
  
  ## TAB 2: Sample Selection
    DF <- reactive({
    DF <- read.table("example_data.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    DF$Group <- rep('enter group', nrow(DF))

   return(DF)
  })

  
  ## Handsontable
  observe({
    if (!is.null(input$hot)) {
      values[["previous"]] <- isolate(values[["DF"]])
      DF = hot_to_r(input$hot)
    } else {
      if (is.null(values[["DF"]]))
        DF <- DF()
      else
        DF <- values[["DF"]]
    }
    values[["DF"]] <- DF
  })
  
  output$hot <- renderRHandsontable({
    DF <- values[["DF"]]
    if (!is.null(DF))
      rhandsontable(DF, stretchH = "all")
  })
  
  ## Add column
  output$ui_newcolname <- renderUI({
    textInput("newcolumnname", "Name", sprintf("newcol%s", 1+ncol(values[["DF"]])))
  })
  
  observeEvent(input$addcolumn, {
    DF <- isolate(values[["DF"]])
    values[["previous"]] <- DF
    newcolumn <- eval(parse(text=sprintf('%s(nrow(DF))', isolate(input$newcolumntype))))
    values[["DF"]] <- setNames(cbind(DF, newcolumn, stringsAsFactors=FALSE), c(names(DF), isolate(input$newcolumnname)))
  })

  
  ## TAB3: Differential Expression Analysis
  # generate query output from inital text
  dea_res <- reactive({
    res <- read.table("DEA_results.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    return(res)
    })
  
  output$DEQtable <- renderDataTable(
    dea_res())
  
}

#################################
## Run the application 
shiny::shinyApp(ui = ui, server = server)

