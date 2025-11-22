## ================================================================
## Runtime helper functions for ShinyLive mRNA app

## ---- Libraries ----
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(tidygraph)
library(pracma)

## ---- Load precomputed network data ----
if (file.exists("net_data.Rdata")) {
  load("net_data.Rdata")  # should define Net, Module, enrich_2_module, genes, etc.
} else {
  stop("net_data.Rdata not found in app directory. Build it on the cluster first.")
}

## ---- Set app default for a  default_gene ----
default_gene <- "atg12"

## -------------------------------------------------------------
## Search functions
## -------------------------------------------------------------

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
  moduleIDs <- Net %N%
    filter(feature %in% gene) %>%
    pull(module)

  moduleIDs <- unique(moduleIDs[!is.na(moduleIDs)])

  genes_out <- character(0)

  if (length(moduleIDs) > 0) {
    for (ID in moduleIDs) {
      new_genes <- searchForModule(Module, ID)
      if (length(new_genes) > 0) {
        genes_out <- c(genes_out, new_genes)
      }
    }
  }

  neighbors <- Net %N%
    filter(feature %in% gene) %>%
    pull(neighbors)

  if (length(neighbors) > 0) {
    genes_out <- c(genes_out, unlist(neighbors))
  }

  unique(genes_out)
}

searchForGeneList <- function(Net, Module, gene_list, search_additional) {
  gene_list <- unique(gene_list[nzchar(gene_list)])
  if (length(gene_list) == 0) {
    return(character(0))
  }

  result_list <- Net %N%
    filter(feature %in% gene_list) %>%
    pull(feature)

  ## Expand by modules
  if ("mod" %in% search_additional) {
    mod_genes <- character(0)

    moduleIDs <- unique(c(
      Net %N%
        filter(feature %in% result_list, module != -9999) %>%
        pull(module),
      unlist(
        Net %N%
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
    neighbors <- Net %N%
      filter(feature %in% result_list) %>%
      pull(neighbors)

    if (length(neighbors) > 0) {
      result_list <- c(result_list, unlist(neighbors))
    }
  }

  unique(result_list)
}

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

## -------------------------------------------------------------
## Diffusion-related helpers
## -------------------------------------------------------------

diffScoreSubgraph <- function(Net, min_targets, top_regs) {
  nodes <- Net %N%
    as_tibble() %>%
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

  induceSubraph(Net, genes) %E>%
    mutate(color_code = "#666")
}

## ------------------------------------------------------------
## Subgraph functions
## ------------------------------------------------------------

induceSubraph <- function(Net, list) {
  list <- unique(list[nzchar(list)])
  Net %>%
    convert(to_subgraph, feature %in% list)
}

moduleSubgraph <- function(Net, Module, module_id) {
  mod_genes <- searchForModule(Module, module_id)

  if (length(mod_genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  induceSubraph(Net, mod_genes) %E>%
    mutate(color_code = "#666")
}

geneSubgraph <- function(Net, Module, gene) {
  list_genes <- searchForGene(Net, Module, gene)

  if (length(list_genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  induceSubraph(Net, list_genes) %E>%
    mutate(color_code = "#666")
}

## --- Steiner tree helpers ---

getDistMatrix <- function(Net, gene_list) {
  gene_id <- Net %N%
    filter(feature %in% gene_list) %>%
    pull(id)

  gene_name <- Net %N%
    filter(feature %in% gene_list) %>%
    pull(feature)

  dist_matrix <- Net %N% as_tibble() %>% select(feature)

  for (idx in seq_along(gene_id)) {
    root <- gene_id[idx]
    name <- gene_name[idx]

    dist_matrix <- Net %N%
      mutate(dist = bfs_dist(root)) %>%
      as_tibble() %>%
      select(feature, dist) %>%
      rename(!!name := dist) %>%
      right_join(dist_matrix, by = "feature")
  }

  dist_matrix %>%
    mutate(across(-feature, ~ replace(.x, .x == 0, nrow(dist_matrix) + 1)))
}

buildSteinerTrees <- function(Net, gene_list) {
  if (length(gene_list) < 2) {
    return(Net %>% induce_subgraph(feature %in% gene_list))
  }

  Net_u <- Net %>% convert(to_undirected)
  dist_matrix <- getDistMatrix(Net_u, gene_list)
  gene_names  <- colnames(dist_matrix)[-1]

  dist2graph <- tibble(
    gene_names = gene_names,
    Dist       = Inf,
    Closest    = ""
  )

  steiner_tree <- tbl_graph()
  restrict_dist_matrix <- dist_matrix %>%
    filter(feature %in% gene_names)

  # initialize
  for (i in 2:ncol(restrict_dist_matrix)) {
    gene <- colnames(restrict_dist_matrix)[i]
    idx  <- which(dist2graph$gene_names == gene)

    col_vals <- restrict_dist_matrix[[i]]
    if (all(is.na(col_vals))) {
      dist2graph$Closest[idx] <- gene
      dist2graph$Dist[idx]    <- nrow(dist_matrix) + 1
    } else {
      dist <- min(col_vals, na.rm = TRUE)
      match_idx <- which(col_vals == dist)[1]
      closest   <- restrict_dist_matrix$feature[match_idx]
      dist2graph$Dist[idx]    <- dist
      dist2graph$Closest[idx] <- closest
    }
  }

  pick <- dist2graph %>%
    arrange(Dist) %>%
    slice(1)

  path_2_add <- c(pick$gene_names, pick$Closest)

  nodes_2_add <- Net_u %N%
    filter(feature %in% path_2_add) %>%
    pull(id)

  if (length(nodes_2_add) < 2) {
    steiner_tree <- Net_u %>% induce_subgraph(feature %in% path_2_add)
  } else {
    steiner_tree <- Net_u %>%
      convert(to_shortest_path, nodes_2_add[1], nodes_2_add[2])
  }

  nodes <- steiner_tree %N>% as_tibble() %>%  pull(feature)

  dist_matrix <- dist_matrix %>%
    select(-intersect(gene_list, nodes))
  dist2graph <- dist2graph %>%
    filter(gene_names %in% setdiff(gene_names, nodes))

  while (length(intersect(nodes, gene_names)) < length(gene_names) &&
         nrow(dist2graph) > 0) {

    restrict_dist_matrix <- dist_matrix %>%
      filter(feature %in% nodes)

    for (i in 2:ncol(restrict_dist_matrix)) {
      gene <- colnames(restrict_dist_matrix)[i]
      idx  <- which(dist2graph$gene_names == gene)

      col_vals <- restrict_dist_matrix[[i]]
      if (all(is.na(col_vals))) {
        dist2graph$Closest[idx] <- gene
        dist2graph$Dist[idx]    <- nrow(dist_matrix) + 1
      } else {
        dist <- min(col_vals, na.rm = TRUE)
        match_idx <- which(col_vals == dist)[1]
        closest   <- restrict_dist_matrix$feature[match_idx]
        dist2graph$Dist[idx]    <- dist
        dist2graph$Closest[idx] <- closest
      }
    }

    pick <- dist2graph %>%
      arrange(Dist) %>%
      slice(1)

    path_2_add <- c(pick$gene_names, pick$Closest)

    nodes_2_add <- Net_u %N%
      filter(feature %in% path_2_add) %>%
      pull(id)

    if (length(nodes_2_add) == 1) {
      node_info <- Net_u %N%
        filter(id %in% nodes_2_add) %>%
        as_tibble()
      steiner_tree <- steiner_tree %>% bind_nodes(node_info)
    } else {
      Path <- Net_u %>%
        convert(to_shortest_path, nodes_2_add[1], nodes_2_add[2])
      steiner_tree <- steiner_tree %>%
        graph_join(Path)
    }

    nodes <- steiner_tree %N>% pull(feature)

    dist_matrix <- dist_matrix %>%
      select(setdiff(colnames(dist_matrix), nodes))
    dist2graph <- dist2graph %>%
      filter(gene_names %in% setdiff(gene_names, nodes))
  }

  steiner_tree
}

geneListSubgraph <- function(Net, Module, gene_list, search_additional) {
  gene_list <- unique(gene_list[nzchar(gene_list)])

  if (length(gene_list) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  st <- NULL

  if ("stein" %in% search_additional && length(gene_list) > 1) {
    st <- buildSteinerTrees(Net, gene_list) %E>%
      mutate(is_steiner = TRUE)
    if (gorder(st) > 0) {
      gene_list <- st %N% pull(feature)
    }
  }

  list_genes <- searchForGeneList(Net, Module, gene_list, search_additional)

  if (length(list_genes) == 0) {
    return(Net %>% induce_subgraph(feature %in% character(0)))
  }

  sub_graph <- induceSubraph(Net, list_genes)

  if (!is.null(st) && gorder(st) > 0 && "stein" %in% search_additional) {
    sub_graph <- graph_join(sub_graph, st) %E>%
      mutate(
        is_steiner = replace_na(is_steiner, FALSE),
        color_code = if_else(is_steiner, "#fb8072", "#666")
      )
  } else {
    sub_graph <- sub_graph %E>%
      mutate(color_code = "#666")
  }

  sub_graph
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

  induceSubraph(Net, mod_genes) %E>%
    mutate(color_code = "#666")
}

graph2NodeEdgeTables <- function(Net) {
  graph_nodes <- Net %N%
    as_tibble() %>%
    mutate(id = id - 1)

  graph_edges <- Net %E%
    as_tibble() %>%
    mutate(from = from - 1, to = to - 1)

  list(graph_nodes, graph_edges)
}

## -------------------------------------------------------------
## Print / table helpers
## -------------------------------------------------------------

printNodeInfo <- function(Net, node_name) {
  if (is.na(node_name) || !nzchar(node_name)) {
    return(" ")
  }

  node_data <- Net %N%
    filter(feature == node_name) %>%
    as_tibble()

  if (nrow(node_data) == 0) {
    return(sprintf("Node %s not found in network.", node_name))
  }

  node_id <- node_data %>% pull(id)

  neighbors <- setdiff(
    Net %>%
      convert(to_local_neighborhood, node_id) %N%
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

printModuleInfo <- function(Module, module_id, gene_list, genes) {
  if (is.na(module_id)) {
    return(" ")
  }

  gene_list <- unique(gene_list[nzchar(gene_list)])

  mod_local <- if (!purrr::is_empty(gene_list)) {
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
  module_go <- module_info %>%
    select(GO) %>%
    unnest(GO, keep_empty = TRUE)

  if (nrow(module_go) == 0 || all(is.na(module_go$go))) {
    go_text <- "No enriched GO terms<br/>"
  } else {
    go_text <- ""
    for (i in seq_len(nrow(module_go))) {
      g   <- module_go$go[i]
      pv  <- module_go$go_correct_pvalue[i]
      gg  <- module_go$go_genes[[i]]
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

  if (!purrr::is_empty(gene_list)) {
    gene_list_str <- paste(
      intersect(module_info$gene_list[[1]], gene_list),
      collapse = ", "
    )
    enrich_p <- module_info$enrich_pval
    corr_p   <- module_info$corrected_pval
  }

  if (purrr::is_empty(gene_list)) {
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

  unique_modules <- SubNet %N%
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

getModuleID <- function(Net, node_name) {
  module_id <- Net %N%
    filter(feature == node_name) %>%
    as_tibble() %>%
    pull(module)

  if (length(module_id) == 0) {
    return(NA_integer_)
  }

  module_id
}

prepNodeTable <- function(Nodes_Table, disp_num) {
  Nodes_Table %>%
    mutate(
      .tidygraph_node_index = NULL,
      enriched_modules      = NULL
    ) %>%
    rowwise() %>%
    mutate(
      go = paste(
        unlist(setdiff(go[1:disp_num], NA)),
        collapse = " <br/>"
      ),
      neighbors = paste(unlist(neighbors), collapse = " | ")
    ) %>%
    ungroup() %>%
    select(-any_of(c("geneSuper", "expression", "Ortholog 1-1"))) %>%
    rename("Gene Name" = "feature") %>%
    mutate(
      id        = NULL,
      regulator = NULL
    ) %>%
    {
      out <- .
      out$module[out$module == -9999] <- NA
      out
    }
}

prepModuleTable <- function(Module_Table, method, disp_num = 5) {
  GO <- Module_Table %>%
    select(module, GO) %>%
    unnest(GO, keep_empty = TRUE) %>%
    group_by(module) %>%
    slice_head(n = disp_num)

  GO <- GO %>%
    rowwise() %>%
    mutate(
      GO = if (!is.na(go)) {
        sprintf(
          "GO term: %s<br/>module corrected p-value: %.02e<br/>Module Genes: %s<br/>",
          go,
          go_correct_pvalue,
          paste(unlist(go_genes), collapse = " | ")
        )
      } else {
        ""
      }
    ) %>%
    select(module, GO) %>%
    group_by(module) %>%
    nest(GO = "GO") %>%
    rowwise() %>%
    mutate(GO = paste(unlist(as.list(GO)), collapse = "<br/><br/>"))

  regulators_tbl <- Module_Table %>%
    select(module, regulators) %>%
    unnest(regulators, keep_empty = TRUE) %>%
    rowwise() %>%
    mutate(
      Regulators = if (!is.na(regulator)) {
        sprintf(
          "enriched regulator: %s<br/>regulator corrected p-value: %.02e<br/>regulator targets: %s<br/>",
          regulator,
          reg_correct_pvalue,
          paste(unlist(reg_target_genes), collapse = " | ")
        )
      } else {
        ""
      }
    ) %>%
    select(module, Regulators) %>%
    group_by(module) %>%
    nest(Regulators = "Regulators") %>%
    rowwise() %>%
    mutate(Regulators = paste(unlist(as.list(Regulators)), collapse = "<br/><br/>"))

  if (method == "list") {
    Module_Table <- Module_Table %>%
      select(module, `Genes on List`, gene_list, corrected_pval) %>%
      rowwise() %>%
      mutate(
        Genes = paste(sort(unlist(gene_list)), collapse = " | "),
        `Genes on List` = paste(sort(unlist(`Genes on List`)), collapse = " | ")
      ) %>%
      select(module, `Genes on List`, Genes, corrected_pval) %>%
      left_join(GO, by = "module") %>%
      left_join(regulators_tbl, by = "module") %>%
      rename(`Gene List enrichment p-value` = corrected_pval)
  } else {
    Module_Table <- Module_Table %>%
      select(module, gene_list) %>%
      rowwise() %>%
      mutate(Genes = paste(sort(unlist(gene_list)), collapse = " | ")) %>%
      select(module, Genes) %>%
      left_join(GO, by = "module") %>%
      left_join(regulators_tbl, by = "module")
  }

  Module_Table %>%
    mutate(
      module = sprintf(
        '<a href="#" onclick=Shiny.setInputValue("module_id_info", %s);">%s</a>',
        module, module
      )
    )
}

regulators <- function(Net) {
  Net %N%
    filter(regulator == "scr" | regulator == TRUE) %>%
    as_tibble()
}

