# app.R inside shiny-apps/mrna

# ---- Required packages ----
library(shiny)
library(DT)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(tidygraph)
library(pracma)

# ---- Load helper functions ----
# must be in the same directory as app.R
source("aux_functions_otegui_runtime.R")

# ---- Load data ----
# MUST be in the same directory as app.R
if (file.exists("net_data.Rdata")) {
  load("net_data.Rdata")   # loads Net and Module
} else {
  stop("net_data.Rdata not found in app directory.")
}

# ---- UI ----
ui <- fluidPage(
  titlePanel(title),

  sidebarLayout(
    sidebarPanel(
      textInput(
        "gene",
        label = "Gene (feature or common name):",
        value = default_gene
      ),
      actionButton("go", "Update"),
      helpText("Gene should match the 'feature' column in the mRNA network.")
    ),

    mainPanel(
      h3("Node summary"),
      htmlOutput("node_info"),

      tags$hr(),

      h3("Genes in module / local neighborhood"),
      DTOutput("gene_table")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  current_gene <- reactiveVal(default_gene)

  observeEvent(input$go, {
    g <- trimws(input$gene %||% "")
    if (nzchar(g)) current_gene(g)
  }, ignoreInit = TRUE)

  output$node_info <- renderUI({
    g <- current_gene()
    HTML(printNodeInfo(Net, g))
  })

  output$gene_table <- renderDT({
    g <- current_gene()
    genes <- searchForGene(Net, Module, g)

    if (length(genes) == 0 || all(is.na(genes))) {
      return(tibble(message = "No genes found for this query"))
    }

    tibble(feature = genes)
  },
  options = list(pageLength = 10),
  rownames = FALSE
  )
}

# ---- App object ----
shinyApp(ui = ui, server = server)

