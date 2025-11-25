## =====================================================
## Runtime helpers for ShinyLive mRNA network app
## =====================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidygraph)
})

## ---- Load precomputed network data ----
## net_data.Rdata must define at least:
##   Net (tbl_graph), Module (module annotation table),
##   enrich_2_module (module–GO mapping),
##   genes (gene universe), enriched_go_terms, module_ids
if (file.exists("net_data.Rdata")) {
  load("net_data.Rdata")
} else {
  stop("net_data.Rdata not found in app directory.")
}

## ---- Single app default ----
default_gene <- "atg12"

## ====================================================
## 1. BASIC HELPERS
## ====================================================

# Noode features in the full network
get_net_features <- function(Net) {
  Net %N>% as_tibble() %>% pull(feature) %>% as.character()
}

# Table of module & feature pairs
module_gene_long <- function(Module) {
  Module %>%
    dplyr::select(module, gene_list) %>%
    tidyr::unnest_longer(gene_list, values_to = "feature") %>%
    dplyr::mutate(
      module  = as.integer(module),
      feature = as.character(feature)
    )
}

# All genes for one module ID
module_genes <- function(Module, module_id) {
  Module %>%
    dplyr::filter(module == !!module_id) %>%
    dplyr::pull(gene_list) %>%
    unlist() %>%
    as.character() %>%
    unique()
}


## Return all module IDs as a character vector
get_module_ids <- function(Module) {
  col_candidates <- c("module_id", "module", "Module", "mod_id")
  mod_col <- intersect(col_candidates, names(Module))
  if (length(mod_col) == 0L) {
    stop("Cannot determine module id column in Module data.frame.")
  }
  sort(unique(as.character(Module[[mod_col[1L]]])))
}

## Return module id column name (internally)
get_module_col <- function(Module) {
  col_candidates <- c("module_id", "module", "Module", "mod_id")
  mod_col <- intersect(col_candidates, names(Module))
  if (length(mod_col) == 0L) {
    stop("Cannot determine module id column in Module data.frame.")
  }
  mod_col[1L]
}

## Return feature / gene column name for Module
get_module_feature_col <- function(Module) {
  col_candidates <- c("feature", "gene", "gene_id", "name")
  fcol <- intersect(col_candidates, names(Module))
  if (length(fcol) == 0L) {
    stop("Cannot determine feature column in Module data.frame.")
  }
  fcol[1L]
}

## Return regulators as a character vector of feature IDs
regulators <- function(Net) {
  node_df <- Net %N>% as_tibble()
  if (!("regulator" %in% names(node_df))) {
    return(character(0))
  }
  node_df %>%
    filter(regulator == "scr" | regulator == TRUE) %>%
    pull(feature) %>%
    unique() %>%
    sort()
}

## Convenience vectors used by the app UI
module_ids <- get_module_ids(Module)

if (!exists("enriched_go_terms")) {
  ## If not provided, derive from enrich_2_module if available
  if (exists("enrich_2_module")) {
    if ("go_term" %in% names(enrich_2_module)) {
      enriched_go_terms <- sort(unique(enrich_2_module$go_term))
    } else if ("go_id" %in% names(enrich_2_module)) {
      enriched_go_terms <- sort(unique(enrich_2_module$go_id))
    } else {
      enriched_go_terms <- character(0)
    }
  } else {
    enriched_go_terms <- character(0)
  }
}

## ====================================================
## 2. SEARCH HELPERS
## ====================================================

## Return all features in a given module
searchForModule <- function(Module, module_id) {
  # Genes assigned to this module
  gene_in_module <- module_genes(Module, moduleID)

  # Regulators for this module (if any)
  regs_tbl <- Module$regulators[[which(Module$module == moduleID)]]
  regs <- if (!is.null(regs_tbl)) regs_tbl$regulator else character(0)

  unique(c(gene_in_module, regs))

}

## Return genes associated with a single query gene
searchForGene <- function(Net, Module, gene) {
  if (is.null(gene) || !nzchar(gene)) {
    return(character(0))
  }
  searchForGeneList(Net, Module, gene_list = gene)
}

## Main helper: starting from a gene list, expand to modules / neighbors
searchForGeneList <- function(Net, Module, gene_list, search_additional= c("mod", "neigh")) {

  # Normalize options
  search_additional <- intersect(search_additional,
                                 c("mod", "neigh"))

  # Start with genes actually present in the network
  result_list <- Net %N>%
    dplyr::filter(feature %in% gene_list) %>%
    dplyr::pull(feature)

  # Expand by modules
  if ("mod" %in% search_additional) {
    mod_ids <- Net %N>%
      dplyr::filter(feature %in% gene_list) %>%
      dplyr::pull(module) %>%
      unique()
    mod_ids <- mod_ids[!is.na(mod_ids)]

    for (ID in mod_ids) {
      result_list <- c(
        result_list,
        searchForModule(Module, ID)
      )
    }
  # Expand by neighbors
  if ("neigh" %in% search_additional) {
    more_neigh <- Net %N>%
      dplyr::filter(feature %in% gene_list) %>%
      dplyr::pull(neighbors)

    result_list <- c(result_list, unlist(more_neigh))
  }

  unique(result_list)

}



## ====================================================
## 3. SUBGRAPH CONSTRUCTION
## ====================================================

## Core safe inducer, to avoid "feature not found" errors
induceSubraph <- function(Net, features) {
  if (is.null(features) || length(features) == 0L) {
    return(Net %>% activate(nodes) %>% filter(FALSE))
  }

  features <- unique(features[nzchar(features)])

  node_df <- Net %N>% as_tibble()
  if (!("feature" %in% names(node_df))) {
    stop("Node data does not contain a 'feature' column.")
  }

  keep_features <- intersect(features, node_df$feature)
  if (length(keep_features) == 0L) {
    return(Net %>% activate(nodes) %>% filter(FALSE))
  }

  Net %>%
    activate(nodes) %>%
    filter(feature %in% keep_features) %>%
    induce_subgraph()
}

## Subgraph for a module
moduleSubgraph <- function(Net, Module, module_id) {
  genes <- searchForModule(Module, module_id)
  induceSubraph(Net, genes)
}

## Subgraph for a single gene
geneSubgraph <- function(Net, Module, gene) {
  genes <- searchForGene(Net, Module, gene)
  induceSubraph(Net, genes)
}

## Subgraph for a gene list
geneListSubgraph <- function(Net,
                             Module,
                             gene_list,
                             search_additional = c("mod", "neigh")) {
  genes <- searchForGeneList(
    Net           = Net,
    Module        = Module,
    gene_list     = gene_list,
    search_additional = search_additional
  )

  induceSubraph(Net, genes)
}

## Subgraph for a GO term: union of modules annotated with that term
goSubgraph <- function(Net, Module, enrich_2_module, go_term) {
  if (is.null(go_term) || !nzchar(go_term)) {
    return(Net %>% activate(nodes) %>% filter(FALSE))
  }
  if (!exists("enrich_2_module")) {
    warning("enrich_2_module not available; cannot build GO-based subgraph.")
    return(Net %>% activate(nodes) %>% filter(FALSE))
  }

  mcol <- get_module_col(Module)
  fcol <- get_module_feature_col(Module)

  ## Determine which column in enrich_2_module holds the GO term
  go_col <- if ("go_term" %in% names(enrich_2_module)) {
    "go_term"
  } else if ("go_id" %in% names(enrich_2_module)) {
    "go_id"
  } else {
    stop("Cannot determine GO column in enrich_2_module.")
  }

  mod_ids <- enrich_2_module %>%
    filter(.data[[go_col]] == go_term) %>%
    pull(.data[[mcol]]) %>%
    unique()

  if (length(mod_ids) == 0L) {
    return(Net %>% activate(nodes) %>% filter(FALSE))
  }

  genes <- Module %>%
    filter(.data[[mcol]] %in% mod_ids) %>%
    pull(.data[[fcol]]) %>%
    unique()

  induceSubraph(Net, genes)
}

## Placeholder for differential-score based subgraph
diffScoreSubgraph <- function(Net, diff_scores, threshold = 0) {
  if (is.null(diff_scores) || length(diff_scores) == 0L) {
    return(Net %>% activate(nodes) %>% filter(FALSE))
  }
  node_df <- Net %N>% as_tibble()
  if (!("feature" %in% names(node_df))) {
    stop("Node data does not contain a 'feature' column.")
  }

  df <- tibble(feature = names(diff_scores), score = as.numeric(diff_scores))
  keep_features <- df %>%
    filter(score >= threshold) %>%
    pull(feature)

  induceSubraph(Net, keep_features)
}

## =====================================================
## 4. NODE / EDGE TABLES
## =====================================================

graph2NodeEdgeTables <- function(SubNet) {
  node_tbl <- SubNet %N>% as_tibble()
  edge_tbl <- SubNet %E>% as_tibble()

  node_tbl <- node_tbl %>%
    mutate(across(where(is.logical), as.character))

  edge_tbl <- edge_tbl %>%
    mutate(across(where(is.logical), as.character))

  list(nodes = node_tbl, edges = edge_tbl)
}

## =====================================================
## 5. NARRATIVE HELPERS
## =====================================================

## Basic node-level HTML summary
printNodeInfo <- function(Net, gene) {
  node_df <- Net %N>% as_tibble()

  if (!("feature" %in% names(node_df))) {
    return(paste0("No 'feature' column found in node data.<br/>"))
  }

  row <- node_df %>%
    filter(feature == gene) %>%
    slice_head(n = 1)

  if (nrow(row) == 0L) {
    return(paste0("Gene '", gene, "' not found in network.<br/>"))
  }

  disp <- if ("display_name" %in% names(row)) row$display_name[1] else gene
  reg  <- if ("regulator" %in% names(row)) as.character(row$regulator[1]) else NA_character_
  mod  <- if ("module" %in% names(row)) as.character(row$module[1]) else
          if ("module_id" %in% names(row)) as.character(row$module_id[1]) else NA_character_

  out <- c()
  out <- c(out, sprintf("<b>Feature:</b> %s", gene))
  if (!is.na(disp)) out <- c(out, sprintf("<b>Display name:</b> %s", disp))
  if (!is.na(reg))  out <- c(out, sprintf("<b>Regulator:</b> %s", reg))
  if (!is.na(mod))  out <- c(out, sprintf("<b>Module:</b> %s", mod))

  paste(out, collapse = "<br/>")
}

## Map modules to GO terms, if available
module_go <- function(Module, enrich_2_module, module_id) {
  if (!exists("enrich_2_module")) {
    return(tibble(go = character(0)))
  }
  mcol <- get_module_col(Module)

  go_col <- if ("go_term" %in% names(enrich_2_module)) {
    "go_term"
  } else if ("go_id" %in% names(enrich_2_module)) {
    "go_id"
  } else {
    return(tibble(go = character(0)))
  }

  enrich_2_module %>%
    filter(.data[[mcol]] == module_id) %>%
    transmute(go = .data[[go_col]])
}

## Module-level summary
printModuleInfo <- function(Module,
                            enrich_2_module,
                            module_id,
                            gene_list = NULL) {
  fcol <- get_module_feature_col(Module)
  mcol <- get_module_col(Module)

  genes_in_mod <- Module %>%
    filter(.data[[mcol]] == module_id) %>%
    pull(.data[[fcol]]) %>%
    unique()

  n_genes      <- length(genes_in_mod)
  go_tbl       <- module_go(Module, enrich_2_module, module_id)

  txt <- sprintf("<b>Module %s</b><br/>", module_id)
  txt <- paste0(txt, sprintf("Size: %d genes<br/>", n_genes))

  if (!is.null(gene_list) && length(gene_list) > 0L) {
    overlap <- intersect(genes_in_mod, gene_list)
    txt <- paste0(
      txt,
      sprintf("Overlap with query gene list: %d genes<br/>", length(overlap))
    )
  }

  if (nrow(go_tbl) > 0L) {
    go_str <- paste(unique(go_tbl$go), collapse = "; ")
    txt <- paste0(txt, sprintf("GO terms: %s<br/>", go_str))
  }

  txt
}

## High-level summary over all modules touched by subgraph / gene list
printAllModuleInfo <- function(SubNet,
                               Module,
                               gene_list = NULL,
                               genes = NULL) {
  node_df <- SubNet %N>% as_tibble()
  if (!("feature" %in% names(node_df))) {
    return("No 'feature' column in node data.")
  }

  fcol <- get_module_feature_col(Module)
  mcol <- get_module_col(Module)

  present_features <- unique(node_df$feature)

  mod_ids <- Module %>%
    filter(.data[[fcol]] %in% present_features) %>%
    pull(.data[[mcol]]) %>%
    unique()

  if (length(mod_ids) == 0L) {
    return("No modules associated with this query.")
  }

  parts <- map_chr(
    mod_ids,
    ~ printModuleInfo(Module, enrich_2_module, module_id = .x, gene_list = gene_list)
  )

  paste(parts, collapse = "<br/><br/>")
}

## =====================================================
## 6. ENRICHMENT
## =====================================================

## Simple hypergeometric enrichment of query gene list against modules
computeEnrichment <- function(Module, gl, num_genes) {
  gl <- unique(gl[nzchar(gl)])
  if (length(gl) == 0L || num_genes <= 0) {
    mod <- Module
    mcol <- get_module_col(mod)
    tibble(
      module_id      = character(0),
      k              = integer(0),
      K              = integer(0),
      n              = integer(0),
      N              = integer(0),
      pval           = numeric(0),
      corrected_pval = numeric(0)
    )
  } else {
    fcol <- get_module_feature_col(Module)
    mcol <- get_module_col(Module)

    N <- num_genes
    n <- length(gl)

    Module %>%
      group_by(.data[[mcol]]) %>%
      summarise(
        module_id = as.character(first(.data[[mcol]])),
        K         = n_distinct(.data[[fcol]]),
        k         = length(intersect(gl, .data[[fcol]])),
        .groups   = "drop"
      ) %>%
      mutate(
        N    = N,
        n    = n,
        pval = phyper(k - 1, K, N - K, n, lower.tail = FALSE)
      ) %>%
      mutate(
        corrected_pval = p.adjust(pval, method = "BH")
      ) %>%
      arrange(corrected_pval)
  }
}

## =====================================================
## 7. SMALL UTILITIES
## =====================================================

## Get module id(s) for a given gene
getModuleID <- function(Module, gene) {
  fcol <- get_module_feature_col(Module)
  mcol <- get_module_col(Module)

  Module %>%
    filter(.data[[fcol]] == gene) %>%
    pull(.data[[mcol]]) %>%
    unique()
}

prepNodeTable <- function(SubNet) {
  SubNet %N>% as_tibble()
}

prepModuleTable <- function(Module, mod_ids) {
  fcol <- get_module_feature_col(Module)
  mcol <- get_module_col(Module)

  Module %>%
    filter(.data[[mcol]] %in% mod_ids) %>%
    arrange(.data[[mcol]], .data[[fcol]])
}

