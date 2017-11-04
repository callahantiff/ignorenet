library(rhandsontable)
library(shiny)
library(shinydashboard)


  ui <- dashboardPage(
    
    # header
    dashboardHeader(title = "ignornet"),
    
    # sidebar
    dashboardSidebar(
      
      # side menu
      sidebarMenu(
        # menuItem("Home", tabName = "home"),
        menuItem("Query GEO", tabName = "dashboard"),
        menuItem("Differential Expression Analysis", tabName = "dea")
      )
      
    ),
    
    dashboardBody(
      box(
        width = 6, status = "info", solidHeader = TRUE,
        title = "GSE Study Groups",
          div(class='row', 
              div(class="col-sm-5", 
                  uiOutput("ui_newcolname"),
                  actionButton("addcolumn", "Add")),
              div(class="col-sm-4", 
                  radioButtons("newcolumntype", "Type", c("integer", "double", "character"))),
              div(class="col-sm-3"))),
        
        box(
          width = 10, status = "info", solidHeader = TRUE,
          title = "GSE Study Metadata",
          actionButton("cancel", "Cancel last action"),
          rHandsontableOutput("hot")
        
      )
    )
  )
  
  server <- shinyServer(function(input, output) {
    values <- reactiveValues()
    
    DF <- reactive({
      DF <- data.frame(Value = 1:10, Status = TRUE, Name = LETTERS[1:10],
                       Date = seq(from = Sys.Date(), by = "days", length.out = 10),
                       stringsAsFactors = FALSE)
      
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
    
  })
  
  ## run app 
runApp(list(ui=ui, server=server))