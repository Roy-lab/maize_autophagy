# aux_functions_otegui_runtime.R
# Runtime helpers for the ShinyLive app (no hard paths)

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(tidygraph)
library(pracma)
library(DT)

# ---- App defaults
title <- "Maize autophagy: mRNA network inference"
default_edge_color_pallette <- "RdBu"
default_node_color_pallette <- "Reds"
default_gene <- "glu2"
default_expression_heatmap <- "Reds"
default_expression_range <- c(0, 5)
default_expression_min <- 0
default_expression_max <- 10
default_tfa_palette_heatamp <- "PiYG"
default_tfa_range <- c(-2, 2)
default_tfa_min <- -10
default_tfa_max <- 10

# ---- Load precomputed network data ----
if (file.exists("net_data.Rdata")) {
  load("net_data.Rdata")  # loads Net, Module, genename_map, etc.
} else {
  stop("net_data.Rdata not found")
}

# ---- RUNTIME FUNCTIONS ONLY ----

searchForModule <- function(Module, moduleID) {
  modIdx <- which(Module$module == moduleID)
  gene_in_module <- Module$gene_list[[modIdx]]
  regulators <- Module$regulators[[modIdx]]$regulator
  mod_genes <- c(gene_in_module, regulators)
  return(mod_genes)
}

searchForGene <- function(Net, Module, gene) {
  moduleIDs <- Net %N>%
    filter(feature %in% gene) %>%
    pull(module)
  moduleIDs <- unique(moduleIDs[!is.na(moduleIDs)])
  genes <- list()
  for (ID in moduleIDs) {
    new_genes <- searchForModule(Module, ID)
    genes <- append(unlist(genes), unlist(new_genes))
  }

  neighbors <- Net %N>%
    filter(feature %in% gene) %>%
    pull(neighbors)

  genes <- append(unlist(genes), unlist(unlist(neighbors)))
  unique(genes)
}

printNodeInfo <- function(Net, node_name) {
  if (is.na(node_name)) {
    return(" ")
  } else {
    node_data <- Net %N>%
      filter(feature == node_name) %>%
      as_tibble()
    if (nrow(node_data) == 0) {
      return(sprintf("Node %s not found in network.", node_name))
    }

    node_id <- node_data %>% pull(id)

    neighbors <- setdiff(
      Net %>%
        convert(to_local_neighborhood, node_id) %N>%
        as_tibble() %>%
        pull(feature),
      node_name
    )

    text <- sprintf(
      "Node Name: %s<br/>Module: %d<br/>GO Terms: %s<br/>Neighbors: %s",
      node_name,
      node_data$module,
      paste(unlist(node_data$go), collapse = ", "),
      paste(unlist(neighbors), collapse = ", ")
    )
    return(text)
  }
}

