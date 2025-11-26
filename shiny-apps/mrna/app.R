## =====================================================
## Shiny app for MERLIN mRNA network
## ======================================================

suppressPackageStartupMessages({
library(shiny)
library(DT)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(tidygraph)
library(igraph)
library(bslib)
})

## Load runtime helpers and net_data.Rdata
source("aux_functions_otegui_runtime.R")

## Objects for UI
node_df       <- Net %N>% as_tibble()
all_features  <- node_df$feature
regulator_ids <- regulators(Net)
go_terms      <- enriched_go_terms

## =====================================================
## UI
## =====================================================

app_theme <- bs_theme(
  version      = 4,
  bootswatch   = "flatly",
  primary      = "#1E88E5"
)

ui <- fluidPage(
  theme = app_theme,

  tags$head(
    tags$style(HTML("
      body {
        background-color: #f5f7fa;
      }
      .sidebar-panel-custom {
        background-color: #ffffff;
        border-radius: 10px;
        padding: 15px 18px;
        border: 1px solid #dee2e6;
        box-shadow: 0 2px 4px rgba(0,0,0,0.04);
      }
      .sidebar-panel-custom h4 {
        margin-top: 0;
        font-weight: 600;
      }
      .tabbable > .nav-tabs {
        margin-bottom: 15px;
        border-bottom: 1px solid #dee2e6;
      }
      .nav-tabs > li > a {
        padding: 8px 14px;
        font-weight: 500;
      }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        background-color: #1E88E5;
        color: #ffffff !important;
      }
      .nav-tabs > li > a:hover {
        background-color: #e3f2fd;
      }
      .btn-primary {
        font-weight: 500;
        border-radius: 999px;
      }
      .shiny-output-error-validation {
        color: #c62828;
        font-weight: 500;
      }
    "))
  ),

  titlePanel("mRNA Network Explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      class = "sidebar-panel-custom",

      h4("Query options"),

      selectInput(
        "query_type",
        "Query type",
        choices = c(
          "Gene       " = "gene",
          "Module ID"   = "module",
          "Regulator"   = "regulator",
          "GO term"     = "go",
          "Gene list"   = "gene_list"
        ),
        selected = "gene"
      ),

      ## Gene query
      conditionalPanel(
        condition = "input.query_type == 'gene'",
        textInput(
          "gene",
          "Gene (feature id)",
          value = default_gene,
          placeholder = "atg12"
        ),
        helpText("Gene should match gene name in mRNA network.")
      ),

      ## Module query
      conditionalPanel(
        condition = "input.query_type == 'module'",
        selectInput(
          "module_id",
          "Module ID",
          choices  = module_ids,
          selected = module_ids[1]
        )
      ),

      ## Regulator query
      conditionalPanel(
        condition = "input.query_type == 'regulator'",
        selectizeInput(
          "regulator",
          "Regulator",
          choices  = regulator_ids,
          options  = list(
            placeholder = "Type to search regulators",
            maxItems    = 1
          ),
          multiple = FALSE
        )
      ),

      ## GO term query
      conditionalPanel(
        condition = "input.query_type == 'go'",
        selectizeInput(
          "go_term",
          "GO term",
          choices  = go_terms,
          options  = list(
            placeholder = "Type or select a GO term",
            maxItems    = 1
          ),
          multiple = FALSE
        )
      ),

      ## Gene list query
      conditionalPanel(
        condition = "input.query_type == 'gene_list'",
        tags$label("Gene list (comma or whitespace separated)"),
        textAreaInput(
          "gene_list",
          label       = NULL,
          rows        = 5,
          placeholder = "atg12, atg8, atg5"
        ),
        checkboxGroupInput(
          "gene_list_options",
          "Gene list expansion:",
          choices = c(
            "Include modules of genes"   = "mod",
            "Include neighbors"          = "neigh"
          ),
          selected = c("mod", "neigh")
        )
      ),

      hr(),
      actionButton("run_query", "Run query", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Summary",
          br(),
          htmlOutput("summary_text"),
          br(),
          htmlOutput("node_text")
        ),
        tabPanel(
          "Network",
          plotOutput("network_plot", height = "600px")
        ),
        tabPanel(
          "Nodes",
          DTOutput("nodes_table")
        ),
        tabPanel(
          "Edges",
          DTOutput("edges_table")
        ),
        tabPanel(
          "Modules",
          DTOutput("modules_table")
        )
      )
    )
  )
)

## ======================================================
## SERVER
## ======================================================

server <- function(input, output, session) {

  ## Helper: parse gene list
  parse_gene_list <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(character(0))
    gl <- unlist(strsplit(txt, "[,\\s]+"))
    unique(gl[nzchar(gl)])
  }

  ## Core reactive: run requested query & build subgraph & tables
  query_result <- eventReactive(input$run_query, {
    mode       <- input$query_type
    gene_list  <- NULL
    subNet     <- NULL
    node_html  <- NULL

    ## 1) Build appropriate subgraph using helper functions
    if (mode == "gene") {
      g <- trimws(input$gene)
      validate(need(nzchar(g), "Enter a gene name."))

      if (!g %in% all_features) {
        validate(need(FALSE, paste0("Gene '", g, "' not found in network.")))
      }

      subNet    <- geneSubgraph(Net, Module, g)
      gene_list <- g
      node_html <- printNodeInfo(Net, g)

    } else if (mode == "module") {
      m <- input$module_id
      validate(need(!is.null(m) && !is.na(m), "Select a module ID."))

      subNet    <- moduleSubgraph(Net, Module, m)
      gene_list <- searchForModule(Module, m)

    } else if (mode == "regulator") {
      r <- trimws(input$regulator)
      validate(need(nzchar(r), "Select or enter a regulator id."))

      if (!r %in% all_features) {
        validate(need(FALSE, paste0("Regulator '", r, "' not found in network.")))
      }

      subNet    <- geneSubgraph(Net, Module, r)
      gene_list <- r
      node_html <- printNodeInfo(Net, r)

    } else if (mode == "go") {
      go <- input$go_term
      validate(need(!is.null(go) && nzchar(go), "Please select a GO term."))

      subNet    <- goSubgraph(Net, Module, enrich_2_module, go)
      gene_list <- NULL

    } else if (mode == "gene_list") {
      gl <- parse_gene_list(input$gene_list)
      validate(need(length(gl) > 0, "Could not parse any gene ids from the gene list."))

      opts      <- input$gene_list_options
      subNet    <- geneListSubgraph(Net, Module, gl, search_additional = opts)
      gene_list <- gl
    }

    validate(need(!is.null(subNet), "No subgraph could be constructed from this query."))

    if (gorder(subNet) == 0) {
      return(list(
        subNet       = subNet,
        node_html    = node_html,
        module_html  = "No nodes returned for this query.",
        nodes_tbl    = tibble(),
        edges_tbl    = tibble(),
        modules_tbl  = tibble()
      ))
    }

    ## 2) Build node + edge tables
    node_edge <- graph2NodeEdgeTables(subNet)
    nodes_tbl <- node_edge$nodes
    edges_tbl <- node_edge$edges

    ## 3) Module-level narrative + enrichment
    gl_for_enrich <- if (!is.null(gene_list)) gene_list else character(0)

    module_html <- printAllModuleInfo(
      SubNet    = subNet,
      Module    = Module,
      gene_list = gl_for_enrich,
      genes     = if (exists("genes")) genes else NULL
    )

    modules_tbl <- tibble()
    if (length(gl_for_enrich) > 0 && exists("genes")) {
      modules_tbl <- computeEnrichment(Module, gl_for_enrich, num_genes = length(genes)) %>%
        arrange(corrected_pval) %>%
        filter(!is.na(corrected_pval)) %>%
        filter(corrected_pval < 0.1)
    }

    list(
      subNet       = subNet,
      node_html    = node_html,
      module_html  = module_html,
      nodes_tbl    = nodes_tbl,
      edges_tbl    = edges_tbl,
      modules_tbl  = modules_tbl
    )
  })

  ## Outputs -------------------------------------------

  output$node_text <- renderUI({
    res <- query_result()
    if (is.null(res$node_html)) return(NULL)
    HTML(res$node_html)
  })

  output$summary_text <- renderUI({
    res <- query_result()
    HTML(res$module_html)
  })

  output$nodes_table <- renderDT({
    res <- query_result()
    req(nrow(res$nodes_tbl) > 0)
    datatable(
      res$nodes_tbl,
      escape  = FALSE,
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

  output$edges_table <- renderDT({
    res <- query_result()
    req(nrow(res$edges_tbl) > 0)
    datatable(
      res$edges_tbl,
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

output$modules_table <- renderDT({
  res <- query_result()
  req(nrow(res$modules_tbl) > 0)

  tbl <- res$modules_tbl %>%
    dplyr::rename(
      `Number of requests`  = module,
      `Module ID`           = module_id,
      `Number of genes from query list in module` = k,
      `Number of all genes in module` = K,
      `Length of gene list query`  = n
    )

  datatable(
  tbl,
  escape = FALSE,
  options = list(
    pageLength = 20,
    scrollX = TRUE,
    autoWidth = TRUE
  ),
  class = "display nowrap cell-border"
) %>%
  formatStyle(
    columns = names(tbl),
    whiteSpace = "normal",
    wordWrap = "break-word"
  )

})

  ## Network plot
output$network_plot <- renderPlot({
  res    <- query_result()
  subNet <- res$subNet
  req(gorder(subNet) > 0)

  ig <- as.igraph(subNet)

  # Use feature names as labels
  vertex_labels <- igraph::vertex_attr(ig, "feature")

  plot(
    ig,
	vertex.label        = vertex_labels,
    vertex.label.cex    = 1.7,
    vertex.label.color  = "black",
    vertex.color        = "#00FF0080",
	edge.width= 1.5,

    # Node size
    vertex.size         = 20,

    # Push label bottom-left
    vertex.label.dist   = 1.5,              # how far from node center
    vertex.label.degree = pi/4,   # node label position

    edge.arrow.size     = 0.2,
    layout              = layout_with_fr(ig)
  )
})

}

## ====================================================
## App object
## ====================================================

app <- shinyApp(ui = ui, server = server)
app

