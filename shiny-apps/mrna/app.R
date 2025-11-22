# =============================================================
# app.R - ShinyLive runtime app for maize autophagy mRNA network
# =============================================================

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
})

## ---- Helpers & precomputed data ----
source("aux_functions_otegui_runtime.R")
# Loads:
#   Net, Module, enrich_2_module, genes, enriched_go_terms, module_ids, default_gene

## ---- Global vectors for UI choices ----

node_df <- Net %>%
  activate(nodes) %>%
  as_tibble()

all_features <- node_df$feature

regulator_ids <- node_df %>%
  filter(regulator %in% c(TRUE, "scr")) %>%
  pull(feature)

if (!exists("module_ids")) {
  module_ids <- sort(unique(node_df$module[node_df$module != -9999]))
}

if (!exists("enriched_go_terms") && exists("enrich_2_module")) {
  enriched_go_terms <- sort(unique(enrich_2_module$go))
}

if (!exists("genes")) {
  genes <- all_features
}

## -------------------------------------------------------
## UI
## -------------------------------------------------------
ui <- fluidPage(
  titlePanel("mRNA Network Explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput(
        "query_type",
        "Query type",
        choices = c(
          "Gene (node)" = "gene",
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
        textInput("gene", "Gene (feature id)", value = default_gene)
      ),

      ## Module query
      conditionalPanel(
        condition = "input.query_type == 'module'",
        selectInput(
          "module_id",
          "Module ID",
          choices  = sort(unique(module_ids)),
          selected = sort(unique(module_ids))[1]
        )
      ),

      ## Regulator query
      conditionalPanel(
        condition = "input.query_type == 'regulator'",
        selectizeInput(
          "regulator",
          "Regulator (feature id)",
          choices  = sort(unique(regulator_ids)),
          options  = list(placeholder = "Type to search regulators"),
          multiple = FALSE
        )
      ),

      ## GO term query
      conditionalPanel(
        condition = "input.query_type == 'go'",
        selectizeInput(
          "go_term",
          "GO term",
          choices  = sort(unique(enriched_go_terms)),
          options  = list(placeholder = "Type or select a GO term"),
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
          placeholder = "e.g. gene1, gene2, gene3"
        ),
        checkboxGroupInput(
          "gene_list_options",
          "Gene list expansion:",
          choices = c(
            "Include modules of genes" = "mod",
            "Include neighbors"        = "neigh"
          ),
          selected = c("mod", "neigh")
        )
      ),

      hr(),
      actionButton("run_query", "Run query", class = "btn-primary")
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
          plotOutput("subgraph_plot", height = "600px")
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

## -------------------------------------------------------
## SERVER
## -------------------------------------------------------
server <- function(input, output, session) {

  ## Helper: parse gene list text
  parse_gene_list <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(character(0))
    gl <- unique(strsplit(txt, "[,\\s]+")[[1]])
    gl[nzchar(gl)]
  }

  ## Core reactive: run  query and make subgraph & tables
  query_result <- eventReactive(input$run_query, {

    mode       <- input$query_type
    gene_list  <- NULL
    subNet     <- NULL
    node_html  <- NULL

    ## 1) GENE
    if (mode == "gene") {
      g <- trimws(input$gene)
      validate(need(nzchar(g), "Please enter a gene."))

      if (!g %in% all_features) {
        validate(need(FALSE, paste0("Gene '", g, "' not found in network.")))
      }

      subNet    <- geneSubgraph(Net, Module, g)
      gene_list <- g
      node_html <- printNodeInfo(Net, g)

    ## 2) MODULE
    } else if (mode == "module") {
      m <- input$module_id
      validate(need(!is.null(m) && !is.na(m), "Please select a module ID."))

      subNet    <- moduleSubgraph(Net, Module, m)
      gene_list <- searchForModule(Module, m)

    ## 3) REGULATOR
    } else if (mode == "regulator") {
      r <- trimws(input$regulator)
      validate(need(nzchar(r), "Please enter/select a regulator."))

      if (!r %in% all_features) {
        validate(need(FALSE, paste0("Regulator '", r, "' not found in network.")))
      }

      subNet    <- geneSubgraph(Net, Module, r)
      gene_list <- r
      node_html <- printNodeInfo(Net, r)

    ## 4) GO TERM
    } else if (mode == "go") {
      go <- input$go_term
      validate(need(nzchar(go), "Please select a GO term."))

      subNet    <- goSubgraph(Net, Module, enrich_2_module, go)
      gene_list <- NULL

    ## 5) GENE LIST
    } else if (mode == "gene_list") {
      gl <- parse_gene_list(input$gene_list)
      validate(need(length(gl) > 0, "Could not parse any genes."))

      opts      <- input$gene_list_options
      subNet    <- geneListSubgraph(Net, Module, gl, search_additional = opts)
      gene_list <- gl
    }

    validate(need(!is.null(subNet), "No subgraph could be constructed."))

    if (gorder(subNet) == 0) {
      return(list(
        subNet      = subNet,
        node_html   = node_html,
        module_html = "No nodes returned for this query.",
        nodes_tbl   = tibble(),
        edges_tbl   = tibble(),
        modules_tbl = tibble()
      ))
    }

    ## Node + edge tables
    node_edge <- graph2NodeEdgeTables(subNet)
    nodes_tbl <- node_edge[[1]]
    edges_tbl <- node_edge[[2]]

    ## Module-level summary
    gl_for_enrich <- if (!is.null(gene_list)) gene_list else character(0)

    module_html <- printAllModuleInfo(
      SubNet    = subNet,
      Module    = Module,
      gene_list = gl_for_enrich,
      genes     = genes
    )

    ## Enrichment table
    modules_tbl <- tibble()
    if (length(gl_for_enrich) > 0) {
      modules_tbl <- computeEnrichment(Module, gl_for_enrich, num_genes = length(genes)) %>%
        arrange(corrected_pval) %>%
        filter(!is.na(corrected_pval)) %>%
        filter(corrected_pval < 0.1)
    }

    list(
      subNet      = subNet,
      node_html   = node_html,
      module_html = module_html,
      nodes_tbl   = nodes_tbl,
      edges_tbl   = edges_tbl,
      modules_tbl = modules_tbl
    )
  })

  ## Network plot
  output$subgraph_plot <- renderPlot({
    res <- query_result()
    g   <- res$subNet

    validate(
      need(!is.null(g) && gorder(g) > 0, "No nodes in subgraph to plot for this query.")
    )

    ig <- as.igraph(g)

    vlabels <- if ("feature" %in% igraph::vertex_attr_names(ig)) {
      igraph::vertex_attr(ig, "feature")
    } else {
      igraph::V(ig)$name
    }

    vcols <- if ("module" %in% igraph::vertex_attr_names(ig)) {
      as.factor(igraph::vertex_attr(ig, "module"))
    } else {
      "steelblue"
    }

    lay <- igraph::layout_with_fr(ig)

    plot(
      ig,
      layout              = lay,
      vertex.label        = vlabels,
      vertex.size         = 6,
      vertex.label.cex    = 0.6,
      vertex.label.family = "sans",
      vertex.color        = vcols,
      edge.arrow.size     = 0.3,
      main                = "Subgraph"
    )
  })

  ## Text outputs
  output$node_text <- renderUI({
    res <- query_result()
    if (is.null(res$node_html)) return(NULL)
    HTML(res$node_html)
  })

  output$summary_text <- renderUI({
    res <- query_result()
    HTML(res$module_html)
  })

  ## Tables
  output$nodes_table <- renderDT({
    res <- query_result()
    req(nrow(res$nodes_tbl) > 0)
    datatable(res$nodes_tbl, escape = FALSE, options = list(pageLength = 20))
  })

  output$edges_table <- renderDT({
    res <- query_result()
    req(nrow(res$edges_tbl) > 0)
    datatable(res$edges_tbl, options = list(pageLength = 20))
  })

  output$modules_table <- renderDT({
    res <- query_result()
    req(nrow(res$modules_tbl) > 0)
    datatable(res$modules_tbl, escape = FALSE, options = list(pageLength = 20))
  })
}

shinyApp(ui = ui, server = server)

