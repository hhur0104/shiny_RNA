library(shiny)
library(shinyjqui)
library(shinydashboard)
library(plotly)
library(htmlwidgets)

header <- dashboardHeader(
  title = "RNA Differential Expression Analysis Dashboard",
  titleWidth = 350
)

body <- dashboardBody(
  tags$head(tags$style("#shiny-modal img { max-width: 100%; }")),
  uiOutput('multabs')
)

side <- dashboardSidebar(
  width = 350,
  sidebarMenu(
    radioButtons(
      "reference",
      h4("Reference Genome:"),
      choices = c("Human", "Mouse", "Other"),
      selected = "Human",
      inline = TRUE
    ),
    fileInput('file1', div(style="display: inline-block;",
                           fluidRow(style="vertical-align:top;",
                             column(9, h4("Choose count data (.tsv): ")), 
                             column(3, actionButton("file1tip", label = icon("info-circle")))
                           )),
              accept=c('text/tsv', 'text/tab-separated-values,text/plain')),
    
    fileInput('file2', div(style="display: inline-block;",
                           fluidRow(style="vertical-align:top;",
                                    column(9, h4("Choose design file (.csv): ")), 
                                    column(3, actionButton("file2tip", label = icon("info-circle")))
                           )),
                accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    uiOutput('comparisonBoxes'),
    uiOutput('batchCrrct'),
    numericInput("pval", "Adj.P.Value:", value = 0.05, min = 0, max = 100),
    numericInput("lgfch", "Log Fold Change:", value = 1.5),
    actionButton("runAnalysis", "Run"),
    uiOutput('export', class="shiny-input-container")
  )
)

dashboardPage(
  header,
  side,
  body
)