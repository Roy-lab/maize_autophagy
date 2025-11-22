## ============================================================
## Runtime helpers for ShinyLive mRNA network app
## ============================================================

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

induce_subgraph <- function(graph, nodes) {
  # nodes can be logical, character, or integer indices – pass directly to igraph
  igraph::induced_subgraph(graph, vids = nodes)
}

## ==========================================================
## 1. SEARCH HELPERS
## ==========================================================

# Module table is assumed to have:
#  - module (ID)
#  - gene_list : list-column of genes per module
#  - regulators : list-column of regulator info (optional)
#  - GO : list-column of GO info (optional)

searchForModule <- function(Module, moduleID) {
  modIdx <- which(Module$module == moduleID)

  if (length(modIdx) == 0 || is.na(moduleID)) {
    return(character(0))
  }

  gene_in_module <- Module$gene_list[[modIdx]]

  regs_col <- Module$regulators[[modIdx]]
  regulators <- if (!is.null(regs_col)) {
    regs_col$regulator
  } else {
    character(0)
  }

  unique(c(unlist(gene_in_module), unlist(regulators)))
}

searchForGene <- function(Net, Module, gene) {
  # which modules does this gene belong to?
  moduleIDs <- Net %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(feature %in% gene) %>%
    pull(module)

  moduleIDs <- unique(moduleIDs[!is.na(moduleIDs)])
  genes_out <- character(0)

  # add module members
  if (length(moduleIDs) > 0) {
    for (ID in moduleIDs) {
      new_genes <- searchForModule(Module, ID)
      if (length(new_genes) > 0) {
        genes_out <- c(genes_out, new_genes)
      }
    }
  }

  # add neighbors from Net$neighbors list-column if present
  neighbors <- Net %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(feature %in% gene) %>%
    pull(neighbors)

  if (length(neighbors) > 0) {
    genes_out <- c(genes_out, unlist(neighbors))
  }

  unique(genes_out)
}

searchForGeneList <- function(Net, Module, gene_list, search_additional) {
  node_df <- Net %>%
    activate(nodes) %>%
    as_tibble()

  gene_list <- unique(gene_list[nzchar(gene_list)])
  gene_list <- intersect(gene_list, node_df$feature)

  if (length(gene_list) == 0) {
    return(character(0))
  }

  result_list <- node_df %>%
    filter(feature %in% gene_list) %>%
    pull(feature)

  ## Expand by modules
  if ("mod" %in% search_additional) {
    mod_genes <- character(0)

    moduleIDs <- unique(c(
      node_df %>%
        filter(feature %in% result_list, module != -9999) %>%
        pull(module),
      unlist(
        node_df %>%
          filter(feature %in% result_list) %>%
          pull(enriched_modules)
      )
    ))
    moduleIDs <- moduleIDs[!is.na(moduleIDs)]

    if (length(moduleIDs) > 0) {
      for (ID in moduleIDs) {
        new_genes <- searchForModule(Module, ID)
        if (length(new_genes) > 0) {
          mod_genes <- c(mod_genes, new_genes)
        }
      }
    }

    result_list <- c(result_list, mod_genes)
  }

  ## Expand by neighbors
  if ("neigh" %in% search_additional) {
    neighbors <- node_df %>%
      filter(feature %in% result_list) %>%
      pull(neighbors)

    if (length(neighbors) > 0) {
      result_list <- c(result_list, unlist(neighbors))
    }
  }

  unique(result_list)
}


## ===========================================================
## 2. ENRICHMENT
## ===========================================================

computeEnrichment <- function(Module, gl, num_genes) {
  gl <- unique(gl[nzchar(gl)])
  if (length(gl) == 0 || num_genes <= 0) {
    mod <- Module
    mod$enrich_pval    <- NA_real_
    mod$corrected_pval <- NA_real_
    mod[["Genes on List"]] <- vector("list", nrow(mod))
    return(mod)
  }

  n_mod <- nrow(Module)

  Module %>%
    rowwise() %>%
    mutate(
      m_size         = length(gene_list),
      intersect_size = length(intersect(gl, gene_list)),
      `Genes on List` = if (intersect_size > 0) list(intersect(gl, gene_list)) else list(character(0)),
      enrich_pval    = phyper(intersect_size, length(gl), num_genes, m_size, lower.tail = FALSE),
      corrected_pval = pmin(enrich_pval * n_mod, 1)
    ) %>%
    ungroup()
}


## ===========================================================
## 3. SUBGRAPH HELPERS
## ===========================================================

induceSubraph <- function(Net, genes) {
  genes <- unique(genes[nzchar(genes)])
  if (length(genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }
  Net %>%
    induce_subgraph(feature %in% genes)
}

moduleSubgraph <- function(Net, Module, module_id) {
  mod_genes <- searchForModule(Module, module_id)

  if (length(mod_genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  induceSubraph(Net, mod_genes) %>%
    activate(edges) %>%
    mutate(color_code = "#666")
}

geneSubgraph <- function(Net, Module, gene) {
  list_genes <- searchForGene(Net, Module, gene)

  if (length(list_genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  induceSubraph(Net, list_genes) %>%
    activate(edges) %>%
    mutate(color_code = "#666")
}

geneListSubgraph <- function(Net, Module, gene_list, search_additional) {
  list_genes <- searchForGeneList(Net, Module, gene_list, search_additional)

  if (length(list_genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  induceSubraph(Net, list_genes) %>%
    activate(edges) %>%
    mutate(color_code = "#666")
}

goSubgraph <- function(Net, Module, enrich_2_module, go_term) {
  modules_list <- enrich_2_module %>%
    filter(go == go_term) %>%
    pull(module) %>%
    unlist()

  modules_list <- unique(modules_list[!is.na(modules_list)])

  if (length(modules_list) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  mod_genes <- character(0)
  for (m in modules_list) {
    new_genes <- searchForModule(Module, m)
    if (length(new_genes) > 0) {
      mod_genes <- c(mod_genes, new_genes)
    }
  }

  mod_genes <- unique(mod_genes)

  if (length(mod_genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  induceSubraph(Net, mod_genes) %>%
    activate(edges) %>%
    mutate(color_code = "#666")
}

diffScoreSubgraph <- function(Net, min_targets, top_regs) {
  node_df <- Net %>%
    activate(nodes) %>%
    as_tibble()

  nodes <- node_df %>%
    arrange(desc(score)) %>%
    filter(regulator == "scr", degree >= min_targets)

  if (nrow(nodes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  n_take <- min(top_regs, nrow(nodes))
  topN   <- nodes %>% slice_head(n = n_take)

  genes <- unlist(topN$neighbors)
  genes <- unique(genes[nzchar(genes)])

  if (length(genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  induceSubraph(Net, genes) %>%
    activate(edges) %>%
    mutate(color_code = "#666")
}


## ============================================================
## 4. NODE/EDGE TABLES
## ============================================================

graph2NodeEdgeTables <- function(Net) {
  graph_nodes <- Net %>%
    activate(nodes) %>%
    as_tibble()

  graph_edges <- Net %>%
    activate(edges) %>%
    as_tibble()

  list(graph_nodes, graph_edges)
}


## ================================================================
## 5. PRINTING / MODULE INFO
## ================================================================

printNodeInfo <- function(Net, node_name) {
  if (is.na(node_name) || !nzchar(node_name)) {
    return(" ")
  }

  node_data <- Net %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(feature == node_name)

  if (nrow(node_data) == 0) {
    return(sprintf("Node %s not found in network.", node_name))
  }

  node_id <- node_data$id

  neighbors <- setdiff(
    Net %>%
      convert(to_local_neighborhood, node_id) %>%
      activate(nodes) %>%
      as_tibble() %>%
      pull(feature),
    node_name
  )

  sprintf(
    "Node Name: %s<br/>Module: %s<br/>GO Terms: %s<br/>Neighbors: %s",
    node_name,
    paste(node_data$module, collapse = ", "),
    paste(unlist(node_data$go), collapse = ", "),
    paste(unlist(neighbors), collapse = ", ")
  )
}

module_go <- function(enrich_2_module, module_id) {
  enrich_2_module %>%
    filter(module == module_id) %>%
    pull(go) %>%
    unique()
}

printModuleInfo <- function(Module, module_id, gene_list, genes) {
  if (is.na(module_id)) {
    return(" ")
  }

  gene_list <- unique(gene_list[nzchar(gene_list)])

  mod_local <- if (length(gene_list) > 0) {
    computeEnrichment(Module, gene_list, length(genes))
  } else {
    Module
  }

  module_info <- mod_local %>%
    filter(module == module_id)

  if (nrow(module_info) == 0) {
    return(sprintf("Module %s not found.<br/><br/>", module_id))
  }

  ## Regulators
  module_regulators <- module_info %>%
    select(regulators) %>%
    unnest(regulators, keep_empty = TRUE)

  if (nrow(module_regulators) == 0 || all(is.na(module_regulators$regulator))) {
    regulator_text <- "No enriched regulators<br/>"
  } else {
    regulator_text <- ""
    for (i in seq_len(nrow(module_regulators))) {
      r   <- module_regulators$regulator[i]
      pv  <- module_regulators$reg_correct_pvalue[i]
      tgt <- module_regulators$reg_target_genes[[i]]
      tgt <- paste(unlist(tgt), collapse = ", ")
      regulator_text <- sprintf(
        "%sRegulator Name: %s &emsp;P-Value: %.02e<br/>Target Genes: %s<br/>",
        regulator_text, r, pv, tgt
      )
    }
  }

  ## GO terms
  module_go_tbl <- module_info %>%
    select(GO) %>%
    unnest(GO, keep_empty = TRUE)

  if (nrow(module_go_tbl) == 0 || all(is.na(module_go_tbl$go))) {
    go_text <- "No enriched GO terms<br/>"
  } else {
    go_text <- ""
    for (i in seq_len(nrow(module_go_tbl))) {
      g   <- module_go_tbl$go[i]
      pv  <- module_go_tbl$go_correct_pvalue[i]
      gg  <- module_go_tbl$go_genes[[i]]
      gg  <- paste(unlist(gg), collapse = ", ")
      go_text <- sprintf(
        "%sGo Term: %s &emsp;P-Value: %.02e<br/>GO Genes: %s<br/>",
        go_text, g, pv, gg
      )
    }
  }

  gene_list_str <- ""
  enrich_p <- NA_real_
  corr_p   <- NA_real_

  if (length(gene_list) > 0) {
    gene_list_str <- paste(
      intersect(module_info$gene_list[[1]], gene_list),
      collapse = ", "
    )
    enrich_p <- module_info$enrich_pval
    corr_p   <- module_info$corrected_pval
  }

  if (length(gene_list) == 0) {
    sprintf(
      "Module Name: %s<br/>Module Genes: %s<br/><br/>Enriched Regulators:<br/>%s<br/>Enriched GO:<br/>%s<br/><br/>",
      module_id,
      paste(unlist(module_info$gene_list), collapse = ", "),
      regulator_text,
      go_text
    )
  } else {
    sprintf(
      paste0(
        "Module Name: %s<br/>",
        "module enrichment p-value: %.02e<br/>",
        "module corrected p-value: %.02e<br/>",
        "Genes from gene list: %s<br/>",
        "Module Genes: %s<br/><br/>",
        "Enriched Regulators:<br/>%s<br/><br/>",
        "Enriched GO:<br/>%s<br/><br/>"
      ),
      module_id,
      enrich_p,
      corr_p,
      gene_list_str,
      paste(unlist(module_info$gene_list), collapse = ", "),
      regulator_text,
      go_text
    )
  }
}

printAllModuleInfo <- function(SubNet, Module, gene_list, genes) {
  if (gorder(SubNet) == 0) {
    return("No nodes/modules in this subgraph.")
  }

  unique_modules <- SubNet %>%
    activate(nodes) %>%
    as_tibble() %>%
    pull(module) %>%
    unique()

  unique_modules <- unique_modules[!is.na(unique_modules)]

  if (length(unique_modules) == 0) {
    return("No module assignments for nodes in this subgraph.")
  }

  gene_list <- unique(gene_list[nzchar(gene_list)])

  texts <- vapply(
    unique_modules,
    function(id) printModuleInfo(Module, id, gene_list, genes),
    FUN.VALUE = character(1)
  )

  paste(texts, collapse = " ")
}


## ================================================================
## 6. MISC HELPERS
## ================================================================

getModuleID <- function(Net, node_name) {
  module_id <- Net %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(feature == node_name) %>%
    pull(module)

  if (length(module_id) == 0) {
    return(NA_integer_)
  }

  module_id
}

prepNodeTable <- function(Nodes_Table, disp_num = 5) {
  Nodes_Table
}

prepModuleTable <- function(Module_Table, method, disp_num = 5) {
  Module_Table
}

regulators <- function(Net) {
  Net %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(regulator == "scr" | regulator == TRUE)
}

