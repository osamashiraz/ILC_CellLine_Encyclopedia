# ==============================================================================
# Script 02: Fig 1, SupFig 2-4 - Molecular Subtyping, SET, and Fig 1B Sankey
# ==============================================================================
# Description: Runs the full Figure 1 pipeline via run_fig1_subtyping_and_figures():
#              1. Molecular subtyping (ims, ER/HER2, mRNA/RPPA/DNAm consensus clustering)
#              2. SupFig 2-4: similarity heatmaps, PCAs, Sankey diagrams
#              3. Fig 1C: SET signature analysis
#              4. Fig 1B: multiomics Sankey (mRNA -> RPPA -> DNAm)
#              Assigns fig1c_setheatmap and fig1b_sankey to .GlobalEnv.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - CL_Annots: Cell line annotations
#   - BRCA_CL_EXP: RNA expression data
#   - BRCA_CL_RPPA: Protein expression data
#   - BRCA_CL_DNAm: DNA methylation data
#
# Output:
#   - Molecular subtype assignments
#   - Consensus clustering results
#   - SET signature scores
#   - Multiple figure outputs (SupFig 2-4, Fig 1C)
#
# Author: Osama Shiraz Shah
# ==============================================================================

suppressPackageStartupMessages({
  library(BreastSubtypeR)
  library(SummarizedExperiment)
  library(ConsensusClusterPlus)
  library(ComplexHeatmap)
  library(ggplot2)
  library(circlize)
  library(ggsankey)
  library(forcats)
  library(dplyr)
  library(tidyr)
  library(ggthemes)
  library(CancerSubtypes)
})

# ==============================================================================
# FUNCTION: Calculate ER Status
# ==============================================================================
#' @param BRCA_CL_EXP RNA expression matrix
#' @param BRCA_CL_RPPA RPPA protein expression matrix
#' @param CL_Annots Cell line annotations
#' @return Data frame with ER status calls
calculate_er_status <- function(BRCA_CL_EXP, BRCA_CL_RPPA, CL_Annots) {
  
  # message("  Calculating ER status...")
  
  # Intersect samples
  cellids <- intersect(CL_Annots$Name, colnames(BRCA_CL_EXP))
  mat <- as.matrix(BRCA_CL_EXP)
  
  # RNA-based ER status
  ER_status <- data.frame(
    ID = CL_Annots[cellids, ]$Name,
    ER_rna = round(mat["ESR1", cellids], 2),
    ER_rna_flag = ifelse(
      mat["ESR1", cellids] > summary(mat["ESR1", cellids], na.rm = TRUE)[5],
      "ER+", "ER-"
    )
  )
  
  
  ER_status <- ER_status[!is.na(ER_status$ID), ]
  rownames(ER_status) <- ER_status$ID
  
  # Protein-based ER status
  ER_status$ER_protein <- round(
    BRCA_CL_RPPA["ER-A", match(ER_status$ID, colnames(BRCA_CL_RPPA))], 2
  )
  ER_status$ER_protein_flag <- ifelse(
    ER_status$ER_protein > summary(ER_status$ER_protein, na.rm = TRUE)[5],
    "ER+", "ER-"
  )
  
  # Integrated ER status (consensus of RNA and protein)
  ER_status$ER <- ifelse(
    rowSums(cbind(
      ER_status$ER_rna_flag == "ER+",
      ER_status$ER_protein_flag == "ER+"
    ), na.rm = TRUE) >= 1,
    "ER+", "ER-"
  )
  
  # ER_status$ER_rna_z <- scale(ER_status$ER_rna); ER_status$ER_protein_z <- scale(ER_status$ER_protein)
  # ER_status$ER_rna_level <- ifelse(ER_status$ER_rna_z >= 0.5, "High", ifelse(ER_status$ER_rna_z <= -0.5, "Low", "Mid"))
  # ER_status$ER_protein_level <- ifelse(ER_status$ER_protein_z >= 0.5, "High", ifelse(ER_status$ER_protein_z <= -0.5, "Low", "Mid"))
  # ER_status_ICLE <- ER_status[grep("-I", CL_Annots$Name, value = T),]
  # ER_status_ICLE[order(ER_status_ICLE$ID),]
  
  message("  ✓ ER status calculated for ", nrow(ER_status), " samples")
  
  return(ER_status)
}


calculate_pr_status <- function(BRCA_CL_EXP, BRCA_CL_RPPA, CL_Annots) {
  
  # message("  Calculating PR status...")
  
  # Intersect samples
  cellids <- intersect(CL_Annots$Name, colnames(BRCA_CL_EXP))
  mat <- as.matrix(BRCA_CL_EXP)
  
  # RNA-based PR status
  PR_status <- data.frame(
    ID = CL_Annots[cellids, ]$Name,
    PR_rna = round(mat["PGR", cellids], 2),
    PR_rna_flag = ifelse(
      mat["PGR", cellids] > summary(mat["PGR", cellids], na.rm = TRUE)[5],
      "PR+", "PR-"
    )
  )
  
  PR_status <- PR_status[!is.na(PR_status$ID), ]
  rownames(PR_status) <- PR_status$ID
  
  # Protein-based ER status
  PR_status$PR_protein <- round(
    BRCA_CL_RPPA["PR", match(PR_status$ID, colnames(BRCA_CL_RPPA))], 2
  )
  PR_status$PR_protein_flag <- ifelse(
    PR_status$PR_protein > summary(PR_status$PR_protein, na.rm = TRUE)[5],
    "PR+", "PR-"
  )
  
  # Integrated ER status (consensus of RNA and protein)
  PR_status$PR <- ifelse(
    rowSums(cbind(
      PR_status$PR_rna_flag == "PR+",
      PR_status$PR_protein_flag == "PR+"
    ), na.rm = TRUE) >= 1,
    "PR+", "PR-"
  )
  
  # PR_status$PR_rna_z <- scale(PR_status$PR_rna); PR_status$PR_protein_z <- scale(PR_status$PR_protein)
  # PR_status$PR_rna_level <- ifelse(PR_status$PR_rna_z >= 0.5, "High", ifelse(PR_status$PR_rna_z <= -0.5, "Low", "Mid"))
  # PR_status$PR_protein_level <- ifelse(PR_status$PR_protein_z >= 0.5, "High", ifelse(PR_status$PR_protein_z <= -0.5, "Low", "Mid"))
  # PR_status_ICLE <- PR_status[grep("-I", CL_Annots$Name, value = T),]
  # PR_status_ICLE[order(PR_status_ICLE$ID),]
  
  message("  ✓ PR status calculated for ", nrow(PR_status), " samples")
  
  return(PR_status)
}

# ==============================================================================
# FUNCTION: Calculate HER2 Status
# ==============================================================================
#' @param BRCA_CL_EXP RNA expression matrix
#' @param BRCA_CL_RPPA RPPA protein expression matrix
#' @param BRCA_CL_GISTIC GISTIC copy number calls
#' @param CL_Annots Cell line annotations
#' @return Data frame with HER2 status calls
calculate_her2_status <- function(BRCA_CL_EXP, BRCA_CL_RPPA, BRCA_CL_GISTIC, CL_Annots) {
  
  # message("\n  Calculating HER2 status...")
  
  cellids <- intersect(CL_Annots$Name, colnames(BRCA_CL_EXP))
  mat <- as.matrix(BRCA_CL_EXP)
  
  # RNA-based HER2 status
  HER2_status <- data.frame(
    ID = CL_Annots[cellids, ]$Name,
    HER2_rna = round(mat["ERBB2", cellids], 2),
    HER2_rna_flag = ifelse(
      mat["ERBB2", cellids] > summary(mat["ERBB2", cellids], na.rm = TRUE)[5],
      "HER2+", "HER2-"
    )
  )
  
  # Protein-based HER2 status
  HER2_status$HER2_protein <- round(
    BRCA_CL_RPPA["HER2-PY1248", match(HER2_status$ID, colnames(BRCA_CL_RPPA))], 2
  )
  HER2_status$HER2_protein_flag <- ifelse(
    HER2_status$HER2_protein > summary(HER2_status$HER2_protein, na.rm = TRUE)[5],
    "HER2+", "HER2-"
  )
  # Copy number-based HER2 status
  HER2_status$HER2_CN <- BRCA_CL_GISTIC[
    "ERBB2", match(HER2_status$ID, colnames(BRCA_CL_GISTIC))
  ]
  
  HER2_status <- HER2_status[!is.na(HER2_status$ID), ]
  rownames(HER2_status) <- HER2_status$ID
  
  # Integrated HER2 status (requires 2+ evidence types)
  HER2_status$HER2 <- ifelse(
    rowSums(cbind(
      HER2_status$HER2_rna_flag == "HER2+",
      HER2_status$HER2_protein_flag == "HER2+",
      HER2_status$HER2_CN >= 2
    ), na.rm = TRUE) >= 2,
    "HER2+", "HER2-"
  )
  
  # Manual corrections for validated cases
  her2_fix_pos <- c("JIMT1-C")
  HER2_status$HER2[HER2_status$ID %in% her2_fix_pos] <- "HER2+"
  
  # HER2_status$HER2_rna_z <- scale(HER2_status$HER2_rna); HER2_status$HER2_protein_z <- scale(HER2_status$HER2_protein)
  # HER2_status$HER2_rna_level <- ifelse(HER2_status$HER2_rna_z >= 0.5, "High", ifelse(HER2_status$HER2_rna_z <= -0.5, "Low", "Mid"))
  # HER2_status$HER2_protein_level <- ifelse(HER2_status$HER2_protein_z >= 0.5, "High", ifelse(HER2_status$HER2_protein_z <= -0.5, "Low", "Mid"))
  # HER2_status_ICLE <- HER2_status[grep("-I", CL_Annots$Name, value = T),]
  # HER2_status_ICLE[order(HER2_status_ICLE$ID),]
  
  message("  ✓ HER2 status calculated for ", nrow(HER2_status), " samples")
  
  return(HER2_status)
}

# ==============================================================================
# FUNCTION: Run Intrinsic Molecular Subtyping with BreastSubtypeR
# ==============================================================================
#' @param BRCA_CL_EXP RNA expression matrix
#' @param ER_status ER status data frame
#' @param HER2_status HER2 status data frame
#' @param CL_Annots Cell line annotations
#' @return List with Intrinsic Molecular results and summary
run_ims <- function(BRCA_CL_EXP, ER_status, HER2_status, CL_Annots) {
  
  # Gene name conversion
  ens_map <- suppressMessages(suppressWarnings(gene_name_conversion(
    geneNames = rownames(BRCA_CL_EXP),
    fromType = "SYMBOL",
    toType = c("ENTREZID", "ENSEMBL")
  )))
  
  ens_map <- ens_map[!is.na(ens_map$SYMBOL), ]
  ens_map <- ens_map[!duplicated(ens_map$SYMBOL), ]
  ens_map <- ens_map[!is.na(ens_map$ENTREZID), ]
  ens_map <- ens_map[!duplicated(ens_map$ENTREZID), ]
  
  mat <- as.matrix(BRCA_CL_EXP)
  ens_map <- subset(ens_map, SYMBOL %in% rownames(mat))
  
  mat <- mat[ens_map$SYMBOL, ]
  mat <- mat[!is.na(rownames(mat)), ]
  mat <- mat[!duplicated(rownames(mat)), ]
  
  # Gene metadata
  gene_metadata <- data.frame(
    ENTREZID = ens_map$ENTREZID,
    probe = ens_map$SYMBOL,
    stringsAsFactors = FALSE,
    row.names = ens_map$SYMBOL
  )
  
  # Cell metadata
  cell_metadata <- data.frame(
    CellLine = colnames(mat),
    PatientID = colnames(mat),
    Histology = CL_Annots[colnames(mat), "Histology"],
    ER = ER_status[colnames(mat), "ER"],
    HER2 = HER2_status[colnames(mat), "HER2"],
    stringsAsFactors = FALSE,
    row.names = colnames(mat)
  )
  
  cell_metadata <- cell_metadata[
    !is.na(cell_metadata$Histology) | 
      !is.na(cell_metadata$ER) | 
      !is.na(cell_metadata$HER2), 
  ]
  
  # Construct SummarizedExperiment
  se <- SummarizedExperiment(
    assays = list(data = mat[, cell_metadata$CellLine]),
    rowData = gene_metadata,
    colData = cell_metadata
  )
  
  # Run BreastSubtypeR
  data_input <- Mapping(se, method = "max", impute = FALSE, verbose = FALSE)
  ims_result <- suppressMessages(suppressWarnings(BS_Multi(
    data_input = data_input,
    methods = "AUTO",
    Subtype = FALSE,
    hasClinical = FALSE
  )))
  ims_raw_calls <- ims_result$res_subtypes
  
  # Get consensus call
  get_max_call <- function(x) {
    x <- unlist(x)
    tab <- table(x)
    max_val <- max(tab)
    winners <- names(tab[tab == max_val])
    list(
      TopCall = paste(winners, collapse = ";"),
      nVotes = max_val,
      Tie = length(winners) > 1
    )
  }
  BRCA_CL_Subtyping_Summary <- do.call(
    rbind,
    apply(ims_raw_calls[, -11], 1, get_max_call, simplify = FALSE)
  ) |> as.data.frame()
  BRCA_CL_Subtyping_Summary <- as.data.frame(
    t(apply(BRCA_CL_Subtyping_Summary, 1, unlist))
  )
  BRCA_CL_Subtyping_Summary <- cbind(
    cell_metadata[rownames(ims_raw_calls), -2],
    ims_raw_calls,
    BRCA_CL_Subtyping_Summary
  )
  BRCA_CL_Subtyping_Summary$Histology <- factor(
    BRCA_CL_Subtyping_Summary$Histology,
    levels = c("ILC", "ILC-like", "NST")
  )
  BRCA_CL_Subtyping_Summary$ER <- factor(
    BRCA_CL_Subtyping_Summary$ER,
    levels = c("ER+", "ER-")
  )
  BRCA_CL_Subtyping_Summary$HER2 <- factor(
    BRCA_CL_Subtyping_Summary$HER2,
    levels = c("HER2+", "HER2-")
  )
  
  BRCA_CL_Subtyping_Summary$TopCall <- as.character(BRCA_CL_Subtyping_Summary$TopCall)
  BRCA_CL_Subtyping_Summary$TopCall[
    BRCA_CL_Subtyping_Summary$TopCall == "Her2;LumA"
  ] <- "LumA;Her2"
  
  if (exists("annot_cols")) {
    BRCA_CL_Subtyping_Summary$TopCall <- factor(
      BRCA_CL_Subtyping_Summary$TopCall,
      levels = names(annot_cols$PAM50)
    )
  }
  
  message("  ✓ Intrinsic Molecular Subtyping subtyping complete")
  
  return(list(
    results = ims_raw_calls,
    summary = BRCA_CL_Subtyping_Summary
  ))
}

# ==============================================================================
# FUNCTION: Run mRNA Consensus Clustering
# ==============================================================================
#' @param BRCA_CL_EXP RNA expression matrix
#' @param CL_Annots Cell line annotations
#' @param n_features Number of top variable features to use
#' @param K Number of clusters
#' @return List with cluster assignments and consensus matrix
run_mrna_consensus_clustering <- function(BRCA_CL_EXP, CL_Annots, n_features = 6000, K = 8) {
  
  cellids <- intersect(CL_Annots$Name, colnames(BRCA_CL_EXP))
  mat <- as.matrix(BRCA_CL_EXP[, cellids])
  
  # Feature selection (with conditional plot suppression)
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    features <- rownames(quiet_run(CancerSubtypes::FSbyMAD(mat, cut.type = "topk", value = n_features)))
  } else {
    features <- rownames(CancerSubtypes::FSbyMAD(mat, cut.type = "topk", value = n_features))
  }
  
  # Scale data
  d <- pheatmap:::scale_rows(as.matrix(mat[features, ]))
  d <- as.matrix(na.omit(d))
  
  # Run consensus clustering (with conditional plot suppression)
  ccp_dir <- file.path(tempdir(), "ccp_mrna")
  if (exists("SUPPRESS_CONSENSUS_CLUSTER_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_CONSENSUS_CLUSTER_PLOTS", envir = .GlobalEnv)) {
    results <- suppressMessages(quiet_run(ConsensusClusterPlus::ConsensusClusterPlus(
      d, maxK = 10, reps = 2000, pItem = 0.8, pFeature = 0.9,
      clusterAlg = "hc", distance = "pearson", seed = 123,
      title = ccp_dir, plot = NULL
    )))
  } else {
    results <- ConsensusClusterPlus::ConsensusClusterPlus(
      d, maxK = 10, reps = 2000, pItem = 0.8, pFeature = 0.9,
      clusterAlg = "hc", distance = "pearson", seed = 123,
      title = ccp_dir
    )
  }
  
  # Clean up temporary directory
  if (dir.exists(ccp_dir)) {
    unlink(ccp_dir, recursive = TRUE)
  }
  
  # Extract clusters at K
  mrna_clusters <- cutree(results[[K]]$consensusTree, k = K)
  
  # Map clusters to subtypes (based on original analysis)
  mrna_clusters[mrna_clusters %in% 1] <- "Lum"
  mrna_clusters[mrna_clusters %in% 2] <- "Lum/HER2"
  mrna_clusters[mrna_clusters %in% 3] <- "Basal"
  mrna_clusters[mrna_clusters %in% 4] <- "Lum/HER2"
  mrna_clusters[mrna_clusters %in% 5] <- "HER2"
  mrna_clusters[mrna_clusters %in% 6] <- "Lum/HER2"
  mrna_clusters[mrna_clusters %in% 7] <- "Lum/HER2"
  mrna_clusters[mrna_clusters %in% 8] <- "Lum"
  
  mrna_clusters <- factor(mrna_clusters, levels = c("Lum", "Lum/HER2", "HER2", "Basal"))
  names(mrna_clusters) <- names(results[[K]]$consensusClass)
  
  message("  ✓ mRNA clustering complete")
  
  return(list(
    clusters = mrna_clusters,
    consensus_matrix = results[[K]]$consensusMatrix,
    results = results
  ))
}

# ==============================================================================
# FUNCTION: Run RPPA Consensus Clustering
# ==============================================================================
#' @param BRCA_CL_RPPA RPPA expression matrix
#' @param CL_Annots Cell line annotations
#' @param n_features Number of top variable features to use
#' @param K Number of clusters
#' @return List with cluster assignments and consensus matrix
run_rppa_consensus_clustering <- function(BRCA_CL_RPPA, CL_Annots, n_features = 50, K = 9) {
  
  rppa_cellids <- intersect(colnames(BRCA_CL_RPPA), CL_Annots$Name)  
  
  # Feature selection (with conditional plot suppression)
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    features <- rownames(quiet_run(CancerSubtypes::FSbyMAD(
      as.matrix(BRCA_CL_RPPA[, rppa_cellids]),
      cut.type = "topk", value = n_features
    )))
  } else {
    features <- rownames(CancerSubtypes::FSbyMAD(
      as.matrix(BRCA_CL_RPPA[, rppa_cellids]),
      cut.type = "topk", value = n_features
    ))
  }
  
  # Scale data
  d <- pheatmap:::scale_rows(as.matrix(
    t(scale(t(BRCA_CL_RPPA)))[features, rppa_cellids]
  ))
  d <- as.matrix(na.omit(d))
  
  # Run consensus clustering (with conditional plot suppression)
  ccp_dir <- file.path(tempdir(), "ccp_rppa")
  if (exists("SUPPRESS_CONSENSUS_CLUSTER_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_CONSENSUS_CLUSTER_PLOTS", envir = .GlobalEnv)) {
    results <- suppressMessages(quiet_run(ConsensusClusterPlus::ConsensusClusterPlus(
      d, maxK = 10, reps = 2000, pItem = 0.8, pFeature = 0.9,
      clusterAlg = "hc", distance = "pearson", seed = 123,
      title = ccp_dir, plot = NULL
    )))
  } else {
    results <- ConsensusClusterPlus::ConsensusClusterPlus(
      d, maxK = 10, reps = 2000, pItem = 0.8, pFeature = 0.9,
      clusterAlg = "hc", distance = "pearson", seed = 123,
      title = ccp_dir
    )
  }
  
  # Clean up temporary directory
  if (dir.exists(ccp_dir)) {
    unlink(ccp_dir, recursive = TRUE)
  }
  
  # Extract clusters at K
  rppa_clusters <- cutree(results[[K]]$consensusTree, k = K)
  
  # Map clusters to subtypes
  rppa_clusters[rppa_clusters %in% 1] <- "Lum"
  rppa_clusters[rppa_clusters %in% 2] <- "Basal"
  rppa_clusters[rppa_clusters %in% 3] <- "HER2"
  rppa_clusters[rppa_clusters %in% 4] <- "Lum/HER2"
  rppa_clusters[rppa_clusters %in% 5] <- "Basal"
  rppa_clusters[rppa_clusters %in% 6] <- "Lum"
  rppa_clusters[rppa_clusters %in% 7] <- "Lum/HER2"
  rppa_clusters[rppa_clusters %in% 8] <- "Basal"
  rppa_clusters[rppa_clusters %in% 9] <- "Basal"
  
  rppa_clusters <- factor(rppa_clusters, levels = c("Lum", "Lum/HER2", "HER2", "Basal"))
  names(rppa_clusters) <- names(results[[K]]$consensusClass)
  
  message("  ✓ RPPA clustering complete")
  
  return(list(
    clusters = rppa_clusters,
    consensus_matrix = results[[K]]$consensusMatrix,
    results = results
  ))
}

# ==============================================================================
# FUNCTION: Run DNAm Consensus Clustering
# ==============================================================================
#' @param BRCA_CL_DNAm DNAm beta value matrix
#' @param CL_Annots Cell line annotations
#' @param n_features Number of top variable features to use
#' @param K Number of clusters
#' @return List with cluster assignments and consensus matrix
run_dnam_consensus_clustering <- function(BRCA_CL_DNAm, CL_Annots, n_features = 10000, K = 7) {
  
  DNAm_cellids <- intersect(colnames(BRCA_CL_DNAm), CL_Annots$Name)
  
  # Load probe annotation
  if (exists("FILES")) {
    load(FILES$hm450k_probeset)
  } else {
    stop("HM450K probe set not loaded. Please load config.R")
  }
  
  allProbes <- intersect(HM450K_ProbeSet$probeID, rownames(BRCA_CL_DNAm))
  
  # Filter probes with no missing values
  BRCA_CL_DNAm_sub <- BRCA_CL_DNAm[
    rowSums(is.na(BRCA_CL_DNAm[allProbes, DNAm_cellids])) == 0,
    DNAm_cellids
  ]
  
  # Feature selection by variance (with conditional plot suppression)
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    features <- rownames(quiet_run(CancerSubtypes::FSbyVar(
      as.matrix(BRCA_CL_DNAm_sub),
      cut.type = "topk", value = n_features
    )))
  } else {
    features <- rownames(CancerSubtypes::FSbyVar(
      as.matrix(BRCA_CL_DNAm_sub),
      cut.type = "topk", value = n_features
    ))
  }
  features <- features[!is.na(rowMeans(BRCA_CL_DNAm[features, DNAm_cellids]))]
  
  # Scale data
  d <- pheatmap:::scale_rows(as.matrix(
    t(scale(t(BRCA_CL_DNAm)))[features, DNAm_cellids]
  ))
  d <- as.matrix(na.omit(d))
  
  # Run consensus clustering (with conditional plot suppression)
  ccp_dir <- file.path(tempdir(), "ccp_dnam")
  if (exists("SUPPRESS_CONSENSUS_CLUSTER_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_CONSENSUS_CLUSTER_PLOTS", envir = .GlobalEnv)) {
    results <- suppressMessages(quiet_run(ConsensusClusterPlus::ConsensusClusterPlus(
      d, maxK = 10, reps = 2000, pItem = 0.8, pFeature = 0.9,
      clusterAlg = "hc", distance = "pearson", seed = 123,
      title = ccp_dir, plot = NULL
    )))
  } else {
    results <- ConsensusClusterPlus::ConsensusClusterPlus(
      d, maxK = 10, reps = 2000, pItem = 0.8, pFeature = 0.9,
      clusterAlg = "hc", distance = "pearson", seed = 123,
      title = ccp_dir
    )
  }
  
  # Clean up temporary directory
  if (dir.exists(ccp_dir)) {
    unlink(ccp_dir, recursive = TRUE)
  }
  
  # Extract clusters at K
  dnam_clusters <- cutree(results[[K]]$consensusTree, k = K)
  
  # Map clusters to subtypes
  dnam_clusters[dnam_clusters %in% 1] <- "Lum/HER2"
  dnam_clusters[dnam_clusters %in% 2] <- "Lum/HER2"
  dnam_clusters[dnam_clusters %in% 3] <- "Basal"
  dnam_clusters[dnam_clusters %in% 4] <- "Lum/HER2"
  dnam_clusters[dnam_clusters %in% 5] <- "Basal"
  dnam_clusters[dnam_clusters %in% 6] <- "Basal"
  dnam_clusters[dnam_clusters %in% 7] <- "Lum/HER2"
  dnam_clusters[dnam_clusters %in% 8] <- "Basal"
  
  dnam_clusters <- factor(dnam_clusters, levels = c("Lum", "Lum/HER2", "HER2", "Basal"))
  names(dnam_clusters) <- names(results[[K]]$consensusClass)
  
  message("  ✓ DNAm clustering complete")
  
  return(list(
    clusters = dnam_clusters,
    consensus_matrix = results[[K]]$consensusMatrix,
    results = results
  ))
}

# ==============================================================================
# MAIN FUNCTION: Run Complete Molecular Subtyping Pipeline
# ==============================================================================
#' @param BRCA_CL_EXP RNA expression matrix
#' @param BRCA_CL_RPPA RPPA expression matrix
#' @param BRCA_CL_DNAm DNAm beta value matrix
#' @param BRCA_CL_GISTIC GISTIC copy number calls
#' @param CL_Annots Cell line annotations
#' @param output_dir Directory for saving results
#' @return Updated CL_Annots with all subtyping information
run_molecular_subtyping_pipeline <- function(BRCA_CL_EXP, BRCA_CL_RPPA, BRCA_CL_DNAm,
                                             BRCA_CL_GISTIC, CL_Annots,
                                             output_dir = NULL) {
  
  message("═══════════════════════════════════════════════════════")
  message("  Molecular Subtyping Pipeline")
  message("═══════════════════════════════════════════════════════\n")
  
  if (is.null(output_dir)) {
    output_dir <- DIRS$results_sub$molecular_subtyping
  }
  ensure_dir(output_dir)
  
  # Step 1: Calculate ER status
  message("  Step 1/7: Calculating ER status...")
  ER_status <- calculate_er_status(BRCA_CL_EXP, BRCA_CL_RPPA, CL_Annots)
  
  # Step 2: Calculate HER2 status
  message("\n  Step 2/7: Calculating HER2 status...")
  HER2_status <- calculate_her2_status(BRCA_CL_EXP, BRCA_CL_RPPA, BRCA_CL_GISTIC, CL_Annots)
  
  # Step 3: Run ims subtyping
  message("\n  Step 3/7: Running Intrinsic Molecular Subtyping subtyping with BreastSubtypeR...")
  ims_out <- run_ims(BRCA_CL_EXP, ER_status, HER2_status, CL_Annots)
  BRCA_CL_Subtyping_Summary <- ims_out$summary
  ims_raw_calls <- ims_out$results
  
  # Step 4: Run mRNA consensus clustering
  message("\n  Step 4/7: Running mRNA consensus clustering...")
  mrna_results <- run_mrna_consensus_clustering(BRCA_CL_EXP, CL_Annots, n_features = 6000, K = 8)
  BRCA_CL_Subtyping_Summary[, "mRNA Subtypes"] <- mrna_results$clusters[
    BRCA_CL_Subtyping_Summary$CellLine
  ]

  # Step 5: Run RPPA consensus clustering
  message("\n  Step 5/7: Running RPPA consensus clustering...")
  rppa_results <- run_rppa_consensus_clustering(BRCA_CL_RPPA, CL_Annots, n_features = 50, K = 9)
  BRCA_CL_Subtyping_Summary[, "RPPA Subtypes"] <- rppa_results$clusters[
    BRCA_CL_Subtyping_Summary$CellLine
  ]
  
  # Step 6: Run DNAm consensus clustering
  message("\n  Step 6/7: Running DNAm consensus clustering...")
  dnam_results <- run_dnam_consensus_clustering(BRCA_CL_DNAm, CL_Annots, n_features = 6000, K = 7)
  BRCA_CL_Subtyping_Summary[, "DNAm Subtypes"] <- dnam_results$clusters[
    BRCA_CL_Subtyping_Summary$CellLine
  ]
  
  # Step 7: Update CL_Annots with all new information
  message("\n  Step 7/7: Updating annotations and saving results...")
  # Add ER status
  if (!all(colnames(ER_status)[-1] %in% colnames(CL_Annots))) {
    CL_Annots <- cbind(
      CL_Annots,
      ER_status[CL_Annots$Name, setdiff(colnames(ER_status)[-1], colnames(CL_Annots)), drop = FALSE]
    )
  }
  
  # Add HER2 status
  if (!all(colnames(HER2_status)[-1] %in% colnames(CL_Annots))) {
    CL_Annots <- cbind(
      CL_Annots,
      HER2_status[CL_Annots$Name, setdiff(colnames(HER2_status)[-1], colnames(CL_Annots)), drop = FALSE]
    )
  }
  
  # Add subtyping results
  if (!all(colnames(BRCA_CL_Subtyping_Summary)[-c(1, 2, 3, 4)] %in% colnames(CL_Annots))) {
    CL_Annots <- cbind(
      CL_Annots,
      BRCA_CL_Subtyping_Summary[
        CL_Annots$Name,
        setdiff(colnames(BRCA_CL_Subtyping_Summary)[-c(1, 2, 3, 4)], colnames(CL_Annots)),
        drop = FALSE
      ]
    )
  }
  
  # Step 8: Save results
  write.table(
    BRCA_CL_Subtyping_Summary,
    file = file.path(output_dir, "BreastSubtypeR_Results.tsv"),
    row.names = FALSE, sep = "\t"
  )
  
  write.table(
    CL_Annots,
    file = file.path(output_dir, "Molecular_Subtyping_Summary.tsv"),
    row.names = FALSE, sep = "\t"
  )
  
  # Save individual subtype calls
  write.table(
    data.frame(
      ID = names(mrna_results$clusters),
      `mRNA Subtypes` = as.character(mrna_results$clusters),
      check.names = FALSE
    ),
    file = file.path(output_dir, "mRNA_subtypes.tsv"),
    row.names = FALSE, sep = "\t"
  )
  write.table(
    data.frame(
      ID = names(rppa_results$clusters),
      `RPPA Subtypes` = as.character(rppa_results$clusters),
      check.names = FALSE
    ),
    file = file.path(output_dir, "RPPA_subtypes.tsv"),
    row.names = FALSE, sep = "\t"
  )
  write.table(
    data.frame(
      ID = names(dnam_results$clusters),
      `DNAm Subtypes` = as.character(dnam_results$clusters),
      check.names = FALSE
    ),
    file = file.path(output_dir, "DNAm_subtypes.tsv"),
    row.names = FALSE, sep = "\t"
  )
  
  message("  ✓ Annotations updated and results saved\n")
  
  message("═══════════════════════════════════════════════════════")
  message("  Molecular Subtyping Complete")
  message("  Results saved to: ", output_dir)
  message("═══════════════════════════════════════════════════════\n")
  
  if ("mRNA Subtypes" %in% colnames(CL_Annots)) {
    CL_Annots$`mRNA Subtypes` <- factor(
      CL_Annots$`mRNA Subtypes`, 
      levels = SAMPLE_GROUPS$subtype_order
    )
  }
  
  if ("RPPA Subtypes" %in% colnames(CL_Annots)) {
    CL_Annots$`RPPA Subtypes` <- factor(
      CL_Annots$`RPPA Subtypes`, 
      levels = SAMPLE_GROUPS$subtype_order
    )
  }
  
  if ("DNAm Subtypes" %in% colnames(CL_Annots)) {
    CL_Annots$`DNAm Subtypes` <- factor(
      CL_Annots$`DNAm Subtypes`, 
      levels = SAMPLE_GROUPS$subtype_order
    )
  }
  
  # Assign to global environment
  assign("CL_Annots", CL_Annots, envir = .GlobalEnv)
  assign("BRCA_CL_Subtyping_Summary", BRCA_CL_Subtyping_Summary, envir = .GlobalEnv)
  assign("ER_status", ER_status, envir = .GlobalEnv)
  assign("HER2_status", HER2_status, envir = .GlobalEnv)
  assign("mRNA_subtypes", mrna_results$clusters, envir = .GlobalEnv)
  assign("RPPA_subtypes", rppa_results$clusters, envir = .GlobalEnv)
  assign("DNAm_subtypes", dnam_results$clusters, envir = .GlobalEnv)

  return(list(
    CL_Annots = CL_Annots,
    subtyping_summary = BRCA_CL_Subtyping_Summary,
    ims_raw_calls = ims_raw_calls,
    mrna_results = mrna_results,
    rppa_results = rppa_results,
    dnam_results = dnam_results,
    ER_status = ER_status,
    HER2_status = HER2_status
  ))
}

# ==============================================================================
# MAIN: Run full Fig 1 pipeline (subtyping, SupFig 2-4, Fig 1C, Fig 1B)
# ==============================================================================
#' Run molecular subtyping, SupFig 2-4, Fig 1C (SET), and Fig 1B (multiomics Sankey).
#' Assigns fig1c_setheatmap and fig1b_sankey to .GlobalEnv.
#' @param output_dir Directory for outputs (default: DIRS$results_sub$molecular_subtyping)
run_subtyping_and_set_analysis <- function(output_dir = NULL) {
  if (is.null(output_dir)) output_dir <- DIRS$results_sub$molecular_subtyping
  out_dir <- output_dir
  ensure_dir(out_dir)
  
  message("\n========================================")
  message("1. Molecular subtyping")
  message("========================================\n")
  
  pl <- run_molecular_subtyping_pipeline(
    BRCA_CL_EXP, BRCA_CL_RPPA, BRCA_CL_DNAm, BRCA_CL_GISTIC, CL_Annots,
    output_dir = out_dir
  )
  mrna_results <- pl$mrna_results
  rppa_results <- pl$rppa_results
  dnam_results <- pl$dnam_results
  ims_raw_calls <- pl$ims_raw_calls
  BRCA_CL_Subtyping_Summary <- pl$subtyping_summary
  CL_Annots <- pl$CL_Annots
  
  
  message("\n========================================")
  message("2. ERpos/neg and SET signature scores")
  message("========================================\n")
  
  source(file.path(DIRS$scripts$helpers, "03_Fig1C_SET_Signature.R"))
  SET_sig <- load_set_signature()
  SET_scores <- calculate_set_scores(BRCA_CL_EXP, SET_sig)
  ERpos <- setNames(SET_scores$ERpos, rownames(SET_scores))
  
  message("\n========================================")
  message("3. SupFig 2-4: Multi-omic subtypes (heatmaps, PCAs, Sankey)")
  message("========================================\n")
  
  colFun_similarity <- circlize::colorRamp2(c(0, 1), c("white", "blue"), space = "RGB")
  
  # ---- ims entropy heatmap ----
  set.seed(123)
  pdf(file.path(out_dir, "ims_entropy_heatmap.pdf"), width = 10, height = 6)
  draw(Heatmap(t(ims_raw_calls[, -9]), col = annot_cols$PAM50,
               column_order = order(ims_raw_calls$entropy),
               top_annotation = HeatmapAnnotation(Entropy = anno_barplot(ims_raw_calls$entropy))))
  dev.off()
  message("  ✓ ims entropy heatmap saved\n")
  
  # ---- mRNA (SupFig 2A, 2B) ----
  K_mrna <- 8
  htMat <- mrna_results$results[[K_mrna]]$consensusMatrix
  colnames(htMat) <- names(mrna_results$results[[K_mrna]]$consensusClass)
  rownames(htMat) <- colnames(htMat)
  
  mrna_clusters <- mrna_results$clusters
  cellids <- names(mrna_clusters)
  
  col_labs <- colnames(htMat)
  col_labs[grep("-I", col_labs, invert = TRUE)] <- ""
  
  set.seed(123)
  mRNA_sim_HT <- Heatmap(htMat, name = "Similarity\nScore", col = colFun_similarity,
                         row_split = mrna_clusters[colnames(htMat)], column_split = mrna_clusters[colnames(htMat)],
                         cluster_rows = TRUE, cluster_columns = TRUE, cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                         border = TRUE, border_gp = gpar(col = "black", lwd = 0.8), show_row_names = FALSE, column_labels = col_labs,
                         column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
                         column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
                         row_title = " ", gap = unit(0, "mm"), heatmap_legend_param = heatmap_legend_param,
                         top_annotation = HeatmapAnnotation(
                           annotation_name_side = "right", annotation_legend_param = heatmap_legend_param,
                           annotation_name_gp = gpar(fontsize = 12, fontface = "bold", col = "grey25", fontfamily = helv),
                           gap = unit(0, "mm"), gp = gpar(col = "black", lwd = 0.6), border = FALSE,
                           ims = BRCA_CL_Subtyping_Summary[colnames(htMat), "TopCall"],
                           `ER Pos Score` = scale(ERpos)[colnames(htMat), ],
                           `HER2 mRNA` = scale(BRCA_CL_EXP["ERBB2", match(colnames(htMat), colnames(BRCA_CL_EXP))]),
                           col = list(Histology = annot_cols$Histology, ims = annot_cols$PAM50, ER = annot_cols$ER,
                                      ERa = annot_cols$RPPA, HER2_phos = annot_cols$RPPA,
                                      `ER Pos Score` = annot_cols$RNA_zscore, `ER neg` = annot_cols$RNA_zscore,
                                      `HER2 CN` = annot_cols$CN, HER2 = annot_cols$HER2, `HER2 mRNA` = annot_cols$RNA_zscore)
                         ), use_raster = TRUE, raster_quality = 4
  )
  pdf(file.path(out_dir, "SupFig2A_mRNA_subtypes.pdf"), width = 10, height = 8)
  draw(mRNA_sim_HT, merge_legends = TRUE, heatmap_legend_side = "right")
  dev.off()
  message("  ✓ SupFig 2A: mRNA heatmap saved")
  
  # Feature selection with conditional plot suppression
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    feat_mrna <- rownames(quiet_run(CancerSubtypes::FSbyMAD(as.matrix(BRCA_CL_EXP[, cellids]), cut.type = "topk", value = 6000)))
  } else {
    feat_mrna <- rownames(CancerSubtypes::FSbyMAD(as.matrix(BRCA_CL_EXP[, cellids]), cut.type = "topk", value = 6000))
  }
  pca_mrna <- prcomp(t(BRCA_CL_EXP[feat_mrna, cellids]))
  prop_mrna <- as.data.frame(summary(pca_mrna)[[6]])
  Propf_mrna <- round(prop_mrna[2, ], 3) * 100
  plt_df <- data.frame(PC1 = pca_mrna$x[, 1], PC2 = pca_mrna$x[, 2], labels = rownames(pca_mrna$x),
                       Study = ifelse(grepl("-I", rownames(pca_mrna$x)), "ICLE", "CCLE"),
                       mRNA_Subtypes = mrna_clusters[rownames(pca_mrna$x)])
  plt <- ggplot(plt_df, aes(x = PC1, y = PC2, color = mRNA_Subtypes, shape = Study)) +
    geom_point(size = 5, alpha = 1) + scale_color_manual("Subtype", values = annot_cols$Subtypes) +
    scale_shape_manual(values = c(8, 18)) + theme_bw(20) +
    xlab(paste0("PC1 (", Propf_mrna[1], "%)")) + ylab(paste0("PC2 (", Propf_mrna[2], "%)"))
  ggsave(file.path(out_dir, "SupFig2B_mRNA_pca_subtypes.pdf"), plt, width = 7, height = 5)
  message("  ✓ SupFig 2B: mRNA PCA saved\n")
  
  # ---- RPPA (SupFig 3A, 3B) ----
  K_rppa <- 9
  rppa_cellids <- names(rppa_results$clusters)
  htMat <- rppa_results$results[[K_rppa]]$consensusMatrix
  colnames(htMat) <- names(rppa_results$results[[K_rppa]]$consensusClass)
  rownames(htMat) <- colnames(htMat)
  rppa_clusters <- rppa_results$clusters
  col_labs <- colnames(htMat)
  col_labs[grep("-I", col_labs, invert = TRUE)] <- ""
  
  set.seed(123)
  rppa_sim_HT <- Heatmap(htMat, name = "Similarity\nScore", col = colFun_similarity,
                         row_split = rppa_clusters[colnames(htMat)], column_split = rppa_clusters[colnames(htMat)],
                         cluster_rows = TRUE, cluster_columns = TRUE, cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                         border = TRUE, border_gp = gpar(col = "black", lwd = 0.8), show_row_names = FALSE, column_labels = col_labs,
                         column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
                         column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
                         row_title = " ", gap = unit(0, "mm"), heatmap_legend_param = heatmap_legend_param,
                         top_annotation = HeatmapAnnotation(
                           annotation_name_side = "right", annotation_legend_param = heatmap_legend_param,
                           annotation_name_gp = gpar(fontsize = 12, fontface = "bold", col = "grey25", fontfamily = helv),
                           gap = unit(0, "mm"), gp = gpar(col = "black", lwd = 0.6), border = FALSE,
                           ims = BRCA_CL_Subtyping_Summary[colnames(htMat), "TopCall"],
                           `ER Pos Score` = scale(ERpos)[colnames(htMat), ],
                           `HER2 mRNA` = scale(BRCA_CL_EXP["ERBB2", match(colnames(htMat), colnames(BRCA_CL_EXP))]),
                           col = list(Histology = annot_cols$Histology, ims = annot_cols$PAM50, ER = annot_cols$ER,
                                      ERa = annot_cols$RPPA, HER2_phos = annot_cols$RPPA,
                                      `ER Pos Score` = annot_cols$RNA_zscore, `ER neg` = annot_cols$RNA_zscore,
                                      `HER2 CN` = annot_cols$CN, HER2 = annot_cols$HER2, `HER2 mRNA` = annot_cols$RNA_zscore)
                         ), use_raster = TRUE, raster_quality = 4
  )
  pdf(file.path(out_dir, "SupFig3A_RPPA_subtypes.pdf"), width = 10, height = 8)
  draw(rppa_sim_HT, merge_legends = TRUE, heatmap_legend_side = "right")
  dev.off()
  message("  ✓ SupFig 3A: RPPA heatmap saved")
  
  # Feature selection with conditional plot suppression
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    feat_rppa <- rownames(quiet_run(CancerSubtypes::FSbyMAD(as.matrix(BRCA_CL_RPPA[, rppa_cellids]), cut.type = "topk", value = 50)))
  } else {
    feat_rppa <- rownames(CancerSubtypes::FSbyMAD(as.matrix(BRCA_CL_RPPA[, rppa_cellids]), cut.type = "topk", value = 50))
  }
  pca_rppa <- prcomp(t(BRCA_CL_RPPA[feat_rppa, rppa_cellids]))
  prop_rppa <- as.data.frame(summary(pca_rppa)[[6]])
  Propf_rppa <- round(prop_rppa[2, ], 3) * 100
  plt_df <- data.frame(PC1 = pca_rppa$x[, 1], PC2 = pca_rppa$x[, 2], labels = rownames(pca_rppa$x),
                       Study = ifelse(grepl("-I", rownames(pca_rppa$x)), "ICLE", "Marcotte"),
                       RPPA_Subtypes = rppa_clusters[rownames(pca_rppa$x)])
  plt <- ggplot(plt_df, aes(x = PC1, y = PC2, color = RPPA_Subtypes, shape = Study)) +
    geom_point(size = 5, alpha = 1) + scale_color_manual("Subtype", values = annot_cols$Subtypes) +
    scale_shape_manual(values = c(8, 18)) + theme_bw(20) +
    xlab(paste0("PC1 (", Propf_rppa[1], "%)")) + ylab(paste0("PC2 (", Propf_rppa[2], "%)"))
  ggsave(file.path(out_dir, "SupFig3B_RPPA_pca_subtypes.pdf"), plt, width = 7, height = 5)
  message("  ✓ SupFig 3B: RPPA PCA saved\n")
  
  # ---- DNAm (SupFig 4A, 4B) ----
  K_dnam <- 7
  DNAm_cellids <- names(dnam_results$clusters)
  htMat <- dnam_results$results[[K_dnam]]$consensusMatrix
  colnames(htMat) <- names(dnam_results$results[[K_dnam]]$consensusClass)
  rownames(htMat) <- colnames(htMat)
  dnam_clusters <- dnam_results$clusters
  col_labs <- colnames(htMat)
  col_labs[grep("-I", col_labs, invert = TRUE)] <- ""
  
  set.seed(123)
  dnam_sim_HT <- Heatmap(htMat, name = "Similarity\nScore", col = colFun_similarity,
                         row_split = dnam_clusters[colnames(htMat)], column_split = dnam_clusters[colnames(htMat)],
                         cluster_rows = TRUE, cluster_columns = TRUE, cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                         border = TRUE, border_gp = gpar(col = "black", lwd = 0.8), show_row_names = FALSE, column_labels = col_labs,
                         column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
                         column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
                         row_title = " ", gap = unit(0, "mm"), heatmap_legend_param = heatmap_legend_param,
                         top_annotation = HeatmapAnnotation(
                           annotation_name_side = "right", annotation_legend_param = heatmap_legend_param,
                           annotation_name_gp = gpar(fontsize = 12, fontface = "bold", col = "grey25", fontfamily = helv),
                           gap = unit(0, "mm"), gp = gpar(col = "black", lwd = 0.6), border = FALSE,
                           ims = BRCA_CL_Subtyping_Summary[colnames(htMat), "TopCall"],
                           `ER Pos Score` = scale(ERpos)[colnames(htMat), ],
                           `HER2 mRNA` = scale(BRCA_CL_EXP["ERBB2", match(colnames(htMat), colnames(BRCA_CL_EXP))]),
                           col = list(Histology = annot_cols$Histology, ims = annot_cols$PAM50, ER = annot_cols$ER,
                                      ERa = annot_cols$RPPA, HER2_phos = annot_cols$RPPA,
                                      `ER Pos Score` = annot_cols$RNA_zscore, `ER neg` = annot_cols$RNA_zscore,
                                      `HER2 CN` = annot_cols$CN, HER2 = annot_cols$HER2, `HER2 mRNA` = annot_cols$RNA_zscore)
                         ), use_raster = TRUE, raster_quality = 4
  )
  pdf(file.path(out_dir, "SupFig4A_DNAm_subtypes.pdf"), width = 10, height = 8)
  draw(dnam_sim_HT, merge_legends = TRUE, heatmap_legend_side = "right")
  dev.off()
  message("  ✓ SupFig 4A: DNAm heatmap saved")
  
  load(FILES$hm450k_probeset)
  allProbes <- intersect(HM450K_ProbeSet$probeID, rownames(BRCA_CL_DNAm))
  BRCA_CL_DNAm_sub <- BRCA_CL_DNAm[rowSums(is.na(BRCA_CL_DNAm[allProbes, DNAm_cellids])) == 0, DNAm_cellids]
  # Feature selection with conditional plot suppression
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    feat_dnam <- rownames(quiet_run(CancerSubtypes::FSbyVar(as.matrix(BRCA_CL_DNAm_sub), cut.type = "topk", value = 10000)))
  } else {
    feat_dnam <- rownames(CancerSubtypes::FSbyVar(as.matrix(BRCA_CL_DNAm_sub), cut.type = "topk", value = 10000))
  }
  feat_dnam <- feat_dnam[!is.na(rowMeans(BRCA_CL_DNAm[feat_dnam, DNAm_cellids]))]
  pca_dnam <- prcomp(t(BRCA_CL_DNAm[feat_dnam, DNAm_cellids]))
  prop_dnam <- as.data.frame(summary(pca_dnam)[[6]])
  Propf_dnam <- round(prop_dnam[2, ], 3) * 100
  plt_df <- data.frame(PC1 = pca_dnam$x[, 1], PC2 = pca_dnam$x[, 2], labels = rownames(pca_dnam$x),
                       Study = ifelse(grepl("-I", rownames(pca_dnam$x)), "ICLE", "Sanger"),
                       DNAm_Subtypes = dnam_clusters[rownames(pca_dnam$x)])
  plt <- ggplot(plt_df, aes(x = PC1, y = PC2, color = DNAm_Subtypes, shape = Study)) +
    geom_point(size = 5, alpha = 1) + scale_color_manual("Subtype", values = annot_cols$Subtypes) +
    scale_shape_manual(values = c(8, 18)) + theme_bw(20) +
    xlab(paste0("PC1 (", Propf_dnam[1], "%)")) + ylab(paste0("PC2 (", Propf_dnam[2], "%)"))
  ggsave(file.path(out_dir, "SupFig4B_DNAm_pca_subtypes.pdf"), plt, width = 7, height = 5)
  message("  ✓ SupFig 4B: DNAm PCA saved\n")
  
  # ---- Sankey (SupFig 2C, 3C, 4C) ----
  cl_levels <- rev(c("MDAMB134VI", "MDAMB330", "SUM44PE", "CAMA1", "BCK4", "HCC2185", "IPH926", "UACC3133",
                     "HCC2218", "MDAMB453", "ZR7530", "WCRC25", "OCUBM", "600MPE", "SKBR3", "MDAMB468", "HCC1187",
                     "Lum", "LumA", "LumB", "LumA;LumB", "LumA;Normal", "LumA;Her2", "Lum/HER2", "Her2", "HER2", "Normal", "Basal"))
  
  df_rna <- CL_Annots %>% filter(grepl("-I", CL_Annots$Name)) %>% dplyr::select(`mRNA Subtypes`, TopCall) %>%
    make_long(`mRNA Subtypes`, TopCall) %>% drop_na(node) %>%
    mutate(node = factor(node, levels = cl_levels), next_node = factor(next_node, levels = cl_levels))
  
  rna_ims <- ggplot(df_rna, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) +
    geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) +
    scale_fill_manual("Molecular Subtypes", values = c(annot_cols$PAM50, annot_cols$Subtypes)) +
    theme_sankey(base_size = 15) + xlab("")
  
  ggsave(file.path(out_dir, "SupFig2C_RNA_ims.pdf"), rna_ims, width = 7, height = 5)
  message("  ✓ SupFig 2C: mRNA–ims Sankey saved")
  
  df_rppa <- CL_Annots %>% filter(Study == "ICLE") %>% dplyr::select(`RPPA Subtypes`, TopCall) %>%
    make_long(`RPPA Subtypes`, TopCall) %>% drop_na(node) %>%
    mutate(node = factor(node, levels = cl_levels), next_node = factor(next_node, levels = cl_levels))
  rppa_sankey <- ggplot(df_rppa, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) +
    geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) +
    scale_fill_manual("Molecular Subtypes", values = c(annot_cols$PAM50, annot_cols$Subtypes)) +
    theme_sankey(base_size = 15) + xlab("")
  ggsave(file.path(out_dir, "SupFig3C_RPPA_ims.pdf"), rppa_sankey, width = 7, height = 5)
  message("  ✓ SupFig 3C: RPPA–ims Sankey saved")
  
  df_dnam <- CL_Annots %>% filter(Study == "ICLE") %>% dplyr::select(`DNAm Subtypes`, TopCall) %>%
    make_long(`DNAm Subtypes`, TopCall) %>% drop_na(node) %>%
    mutate(node = factor(node, levels = cl_levels), next_node = factor(next_node, levels = cl_levels))
  dnam_sankey <- ggplot(df_dnam, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) +
    geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) +
    scale_fill_manual("Molecular Subtypes", values = c(annot_cols$PAM50, annot_cols$Subtypes)) +
    theme_sankey(base_size = 15) + xlab("")
  ggsave(file.path(out_dir, "SupFig4C_DNAm_ims.pdf"), dnam_sankey, width = 7, height = 5)
  message("  ✓ SupFig 4C: DNAm–ims Sankey saved\n")
  
  
  
  # Fig 1B: Multiomics subtypes Sankey (mRNA -> RPPA -> DNAm)
  message("\n========================================")
  message("4. Fig 1B: Multiomics Sankey")
  message("========================================\n")
  cl_levels_fig1b <- rev(c("MDAMB134VI", "MDAMB330", "SUM44PE", "CAMA1", "BCK4",
                           "HCC2185", "IPH926", "UACC3133", "HCC2218", "MDAMB453", "ZR7530", "WCRC25",
                           "OCUBM", "600MPE", "SKBR3", "MDAMB468", "HCC1187", "Lum", "LumA", "LumB",
                           "LumA;LumB", "LumA;Normal", "LumA;Her2", "Lum/HER2", "Her2", "HER2", "Normal", "Basal"))
  df_fig1b <- CL_Annots %>%
    filter(Study == "ICLE") %>%
    dplyr::select(Sample, `mRNA Subtypes`, `RPPA Subtypes`, `DNAm Subtypes`) %>%
    make_long(Sample, `mRNA Subtypes`, `RPPA Subtypes`, `DNAm Subtypes`) %>%
    drop_na(node) %>%
    mutate(
      node = factor(node, levels = cl_levels_fig1b),
      next_node = factor(next_node, levels = cl_levels_fig1b)
    )
  fig1b_sankey <- ggplot(df_fig1b, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                                       fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) +
    geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) +
    scale_fill_manual("Molecular Subtypes", values = c(annot_cols$Subtypes)) +
    theme_sankey(base_size = 15) + xlab("")
  
  
  message("\n========================================")
  message("5. Fig 1C: SET heatmap")
  message("========================================\n")
  
  sample_ids <- subset(CL_Annots, Study == "ICLE")$Name
  set_res <- generate_SET_score_heatmap(BRCA_CL_EXP, SET_sig = NULL, sample_ids, CL_Annots, BRCA_CL_RPPA, annot_cols, random_seed = 123)
  fig1c_setheatmap <- set_res$SET_heatmap
  
  write.table(set_res$SET_scores,
              file = file.path(out_dir, "SET_Scores.tsv"),
              sep = "\t", col.names = NA
  )
  
  assign("fig1b_sankey", fig1b_sankey, envir = .GlobalEnv)
  assign("fig1c_setheatmap", fig1c_setheatmap, envir = .GlobalEnv)

  message("═══════════════════════════════════════════════════════")
  message("  MOLECULAR SUBTYPING COMPLETE")
  message("═══════════════════════════════════════════════════════\n")
}

# ==============================================================================
# Setup and run (when script is sourced)
# ==============================================================================
if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R")
  else if (file.exists("../config.R")) source("../config.R")
  else stop("config.R not found. Run from 2-Analysis or set wd accordingly.")
}
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}
if (!exists("BRCA_CL_GAM", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = FALSE)
}
required_data <- c("BRCA_CL_EXP", "BRCA_CL_RPPA", "BRCA_CL_DNAm", "BRCA_CL_GISTIC", "CL_Annots")
missing_data <- required_data[!sapply(required_data, exists, envir = .GlobalEnv)]
if (length(missing_data) > 0) {
  message("  Loading missing data objects: ", paste(missing_data, collapse = ", "))
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = FALSE)
}
required_helpers <- c("annot_cols", "helv", "heatmap_legend_param")
missing_helpers <- required_helpers[!sapply(required_helpers, exists, envir = .GlobalEnv)]
if (length(missing_helpers) > 0) {
  stop("Missing required helper objects: ", paste(missing_helpers, collapse = ", "),
       ". Ensure Helper_Functions.R is loaded.")
}
suppressWarnings(run_subtyping_and_set_analysis())
