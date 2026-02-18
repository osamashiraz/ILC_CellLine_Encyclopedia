# ==============================================================================
# Script 23: Figure 5 - RNAi Differential Dependencies (ILC vs NST)
# ==============================================================================
# Description: Performs gene dependency analysis comparing ILC vs NST cell lines
#              using D2 scores and SIMEM. Identifies consensus dependencies,
#              performs KEGG pathway enrichment, and identifies druggable targets.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - RNAi dependency data (via config paths)
#   - CL_Annots: Cell line annotations
#
# Output:
#   - Differential dependency analysis results (TSV files)
#   - KEGG pathway enrichment results
#   - Druggable dependency lists
#   - Figure 5 visualizations (fig5b_consensus_dep_plt, supfig12_dep_ht, 
#     fig5c_pathway_plt, fig5d_pathway_ht, fig5e_drug_ht) assigned to .GlobalEnv
#
# Author: Osama Shiraz Shah
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: Load Configuration and Dependencies
# ------------------------------------------------------------------------------

# Check for config and load if needed
if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) {
    source("config.R")
  } else {
    source("../config.R")
  }
}

# Load helper functions if not already loaded
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}

# Ensure %notin% operator exists
if (!exists("%notin%", envir = .GlobalEnv)) {
  `%notin%` <- Negate(`%in%`)
}

# Setup output directory
ensure_dir(DIRS$results_sub$dependencies)
dep_dir <- DIRS$results_sub$dependencies

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(limma)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(hypeR)
  library(msigdbr)
  library(reshape2)
})

# ------------------------------------------------------------------------------
# SECTION 1: Helper Functions
# ------------------------------------------------------------------------------

# ==============================================================================
# FUNCTION: Run SIMEM Analysis
# ==============================================================================
#' Run SIMEM (Statistical Inference for Mixed Effects Models) analysis
#' 
#' Performs differential dependency analysis comparing ILC vs NST cell lines
#' using SIMEM on Marcotte RNAi screens. Only runs if results file doesn't exist.
#' 
#' @param output_file Path to output TSV file (default: FILES$rnai_simem)
#' @param output_rdata Path to output Rdata file (default: FILES$rnai_simem_rdata)
#' @param parallel_nodes Number of parallel nodes (default: 1). Note: On mac studio parallel runs did not work. Run sequentially
#' @return Invisible NULL (results saved to files)
#' 
#' @details
#'   - Loads breast screens, hairpin annotations, and RNAi annotations
#'   - Runs SIMEM analysis with specified parameters
#'   - Saves results to TSV and Rdata files
run_simem_analysis <- function(output_file = NULL, output_rdata = NULL, parallel_nodes = 1) {
  if (is.null(output_file)) output_file <- FILES$rnai_simem
  if (is.null(output_rdata)) output_rdata <- FILES$rnai_simem_rdata
  
  # Load required SIMEM libraries
  helpers <- DIRS$scripts$helpers
  source(file.path(helpers, "simem", "data_format_lib.R"))
  source(file.path(helpers, "simem", "model_lib.R"))
  source(file.path(helpers, "simem", "simem_lib.R"))
  source(file.path(helpers, "simem", "plot_lib.R"))
  
  # Load additional required packages
  suppressPackageStartupMessages({
    library(Biobase)
    library(blme)
    library(doParallel)
    library(genefilter)
    library(locfit)
    library(MASS)
    library(plyr)
    library(preprocessCore)
    library(reshape)
  })
  
  # Load breast screens data
  load(FILES$rnai_breast_screens)
  
  # Fix sample name inconsistencies
  sampleNames(breast_screens) <- gsub("sum44", "sum44pe", sampleNames(breast_screens))
  breast_screens@phenoData@data$cell_line <- gsub("sum44", "sum44pe", breast_screens@phenoData@data$cell_line)
  breast_screens@phenoData@data$sample_id <- gsub("sum44", "sum44pe", breast_screens@phenoData@data$sample_id)
  breast_screens@phenoData@data$replicate_group <- gsub("sum44", "sum44pe", breast_screens@phenoData@data$replicate_group)
  
  # Create feature map
  featureMap <- breast_screens@featureData@data[, c("gene_id", "symbol")]
  featureMap <- featureMap[!duplicated(featureMap$gene_id), ]
  
  # Load hairpin annotations
  hairpin_annots <- read.delim(FILES$rnai_hairpin)
  hairpin_weights <- hairpin_annots[, c("trcn_id", "gene_id", "weight")]
  
  # Load RNAi annotations and prepare status
  RNAi_CL_Annots <- read.delim(FILES$rnai_annots)
  status <- RNAi_CL_Annots[, c("Sample", "Group")]
  status$cell_line <- tolower(status$Sample)
  rownames(status) <- status$cell_line
  status <- status[, c("cell_line", "Group")]
  
  # Find common cells
  commonCells <- intersect(status$cell_line, names(table(breast_screens$cell_line)))
  filter <- colnames(breast_screens)[breast_screens@phenoData@data$cell_line %in% commonCells]
  
  # Select features
  allGenes <- unique(hairpin_annots$symbol)
  selectedFeatures <- unique(subset(hairpin_annots, symbol %in% allGenes)$gene_id)
  
  rownames(status) <- status$cell_line
  status <- status[commonCells, ]
  
  # Subset breast screens
  breast_screens_sub <- breast_screens[, filter]
  
  # Run SIMEM
  message("  Running SIMEM (this may take a while)...")
  results <- simem(
    screens = breast_screens_sub,
    geneIds = selectedFeatures,
    covariate = "Group",
    covariateFactorOrder = c("NST", "ICLE"),
    reagentWeights = hairpin_weights,
    annotationsPerCellLine = status,
    inverseVarianceWeights = TRUE,
    signalProbWeights = TRUE,
    analyzeReagents = TRUE,
    parallelNodes = parallel_nodes
  )
  
  # Save results
  ensure_dir(dirname(output_file))
  save(results, file = output_rdata)
  write.table(results$gene, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")
  message("  ✓ SIMEM results saved to: ", dirname(output_file))
  
  return(invisible(NULL))
}

# ==============================================================================
# FUNCTION: Load and Prepare RNAi Data
# ==============================================================================
#' Load and prepare RNAi dependency data
#' 
#' @return List containing RNAi_CL_Annots, RNAi_D2, RNAi_D2_discrete, RNAi_D2_zcore
load_rnai_data <- function() {
  message("Loading RNAi data...")
  
  # Load RNAi annotations
  RNAi_CL_Annots <- read.delim(FILES$rnai_annots)
  rownames(RNAi_CL_Annots) <- RNAi_CL_Annots$Sample
  
  RNAi_CL_Annots$Group_factor <- factor(RNAi_CL_Annots$Group, levels = c("ICLE", "NST"), labels = c("ICLE", "NST"))
  
  # Load RNAi D2 scores
  RNAi_D2 <- as.data.frame(vroom::vroom(FILES$rnai_d2, show_col_types = FALSE)) 
  
  rownames(RNAi_D2) <- RNAi_D2$...1
  RNAi_D2 <- RNAi_D2[, grep(x = colnames(RNAi_D2), pattern = "_BREAST", value = TRUE)]
  
  # Clean row and column names
  rownames(RNAi_D2) <- gsub(x = rownames(RNAi_D2), pattern = " [(][[:digit:]|&]+[)]", replacement = "")
  colnames(RNAi_D2) <- gsub(x = colnames(RNAi_D2), pattern = "_BREAST", replacement = "")
  colnames(RNAi_D2) <- gsub(x = colnames(RNAi_D2), pattern = "SUM185PE", replacement = "SUM185")
  
  # Find common samples between RNAi_D2 columns and annotations
  common_samples <- intersect(colnames(RNAi_D2), RNAi_CL_Annots$Sample)
  
  if (length(common_samples) == 0) {
    stop("No matching samples found between RNAi_D2 columns and RNAi_CL_Annots$Sample.\n",
         "RNAi_D2 columns (first 10): ", paste(head(colnames(RNAi_D2), 10), collapse = ", "), "\n",
         "RNAi_CL_Annots$Sample (first 10): ", paste(head(RNAi_CL_Annots$Sample, 10), collapse = ", "))
  }
  
  if (length(common_samples) < length(RNAi_CL_Annots$Sample)) {
    missing_samples <- setdiff(RNAi_CL_Annots$Sample, colnames(RNAi_D2))
    warning("Some samples in annotations are missing from RNAi_D2: ", 
            paste(missing_samples, collapse = ", "))
  }
  
  # Filter to common samples only
  RNAi_D2 <- RNAi_D2[, common_samples]
  
  # Also filter annotations to only include samples that exist in RNAi_D2
  RNAi_CL_Annots <- RNAi_CL_Annots[intersect(RNAi_CL_Annots$Sample, common_samples), ]
  
  # Filter genes with too many NAs
  naMat <- t(apply(RNAi_D2, 1, is.na))
  RNAi_D2 <- RNAi_D2[rowSums(naMat) < 5, ]
  message("  ✓ Loaded ", nrow(RNAi_D2), " genes across ", ncol(RNAi_D2), " cell lines")
  
  # Create discrete dependency categories
  RNAi_D2_discrete <- RNAi_D2
  RNAi_D2_discrete[RNAi_D2 <= -4] <- "VERY STRONG"
  RNAi_D2_discrete[RNAi_D2 <= -2] <- "STRONG"
  RNAi_D2_discrete[RNAi_D2 <= -1] <- "MODERATE"
  RNAi_D2_discrete[RNAi_D2 <= -0.5] <- "WEAK"
  RNAi_D2_discrete[RNAi_D2 >= 1] <- "NONE"
  RNAi_D2_discrete[RNAi_D2 < 1 & RNAi_D2 > -0.5] <- ""
  RNAi_D2_discrete[is.na(RNAi_D2_discrete)] <- ""
  
  # Create z-scored matrix
  RNAi_D2_zcore <- t(scale(t(RNAi_D2)))
  
  return(list(
    RNAi_CL_Annots = RNAi_CL_Annots,
    RNAi_D2 = RNAi_D2,
    RNAi_D2_discrete = RNAi_D2_discrete,
    RNAi_D2_zcore = RNAi_D2_zcore
  ))
}

# ==============================================================================
# FUNCTION: Run D2 Differential Analysis
# ==============================================================================
#' Run D2 differential dependency analysis using limma
#' 
#' @param RNAi_D2 D2 dependency score matrix (genes x samples)
#' @param ICLE_cells Vector of ICLE cell line IDs
#' @param NST_cells Vector of NST cell line IDs
#' @param selectedFeatures Vector of gene IDs to analyze
#' @return Data frame with differential analysis results
run_d2_differential_analysis <- function(RNAi_D2, ICLE_cells, NST_cells, selectedFeatures) {
  # Build data matrix
  D <- RNAi_D2[selectedFeatures, c(ICLE_cells, NST_cells)]
  
  # Build design matrix
  cell_ids <- colnames(D)
  group <- ifelse(cell_ids %in% ICLE_cells, "ICLE", "NST")
  design <- model.matrix(~ 0 + factor(group))
  colnames(design) <- c("ICLE", "NST")
  
  # Fit linear model and set contrast
  fit <- lmFit(D, design = design)
  contr <- makeContrasts(NSTvsICLE = ICLE - NST, levels = design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2, robust = TRUE)
  
  # Get results
  tab <- topTable(fit2, coef = "NSTvsICLE", number = Inf, sort.by = "none")
  
  # Calculate group means
  ICLE_avgD2 <- rowMeans(RNAi_D2[, ICLE_cells, drop = FALSE], na.rm = TRUE)
  NST_avgD2 <- rowMeans(RNAi_D2[, NST_cells, drop = FALSE], na.rm = TRUE)
  
  # Create results data frame
  rnai_res <- data.frame(
    gene = selectedFeatures,
    ICLE_D2_avg = ICLE_avgD2[selectedFeatures],
    NST_D2_avg = NST_avgD2[selectedFeatures],
    D2_DIFF = (ICLE_avgD2 - NST_avgD2)[selectedFeatures],
    D2_pval = tab[selectedFeatures, ]$P.Value,
    D2_fdr = tab[selectedFeatures, ]$adj.P.Val
  )
  rownames(rnai_res) <- rnai_res$gene
  
  message("  ✓ Analyzed ", nrow(rnai_res), " genes")
  return(rnai_res)
}

# ==============================================================================
# FUNCTION: Compute Consensus Dependencies
# ==============================================================================
#' Compute consensus dependencies from D2 and SIMEM results
#' 
#' @param rnai_res D2 differential analysis results
#' @param RNAi_SIMEM_DA SIMEM differential analysis results
#' @param alpha P-value cutoff (default: 0.05)
#' @param eff Effect size cutoff (default: 0.2)
#' @return Updated rnai_res data frame with consensus flags
compute_consensus_dependencies <- function(rnai_res, RNAi_SIMEM_DA, alpha = 0.05, eff = 0.2) {
  # Add SIMEM results
  all_genes <- intersect(rownames(rnai_res), rownames(RNAi_SIMEM_DA))
  rnai_res[all_genes, "simem_DIFF"] <- RNAi_SIMEM_DA[all_genes, "difference_icle"]
  rnai_res[all_genes, "simem_pval"] <- RNAi_SIMEM_DA[all_genes, "pvalue_icle"]
  rnai_res[all_genes, "simem_fdr"] <- RNAi_SIMEM_DA[all_genes, "fdr_icle"]
  
  # Helper function for sign
  sgn <- function(x) ifelse(x > 0, 1L, ifelse(x < 0, -1L, NA_integer_))
  
  # Classify dependencies
  rnai_res <- rnai_res %>%
    mutate(
      sig_d2 = !is.na(D2_pval) & D2_pval <= alpha,
      sig_si = !is.na(simem_pval) & simem_pval <= alpha,
      d2_big = abs(D2_DIFF) >= eff,
      si_big = abs(simem_DIFF) >= eff,
      dir_agree = (sgn(D2_DIFF) == sgn(simem_DIFF)) & d2_big & si_big,
      type = case_when(
        sig_d2 & sig_si & dir_agree ~ "BOTH",
        sig_si & !sig_d2 & dir_agree ~ "SIMEM ONLY",
        sig_d2 & !sig_si & dir_agree ~ "D2 ONLY",
        TRUE ~ "NS"
      ),
      type = factor(type, levels = c("BOTH", "SIMEM ONLY", "D2 ONLY", "NS"))
    )
  
  # Add group labels
  rnai_res$Group <- ifelse(
    rnai_res$type != "NS" & rnai_res$simem_DIFF > 0, "NST",
    ifelse(rnai_res$type != "NS" & rnai_res$simem_DIFF < 0, "ILC", "")
  )
  
  rnai_res$Group_type <- ifelse(
    rnai_res$type != "NS", 
    paste0(rnai_res$Group, " - ", rnai_res$type), 
    rnai_res$type
  )
  
  # Filter to genes with SIMEM results
  rnai_res <- rnai_res[!is.na(rnai_res$simem_DIFF), ]
  rownames(rnai_res) <- rnai_res$gene
  
  # Add consensus flag
  rnai_res$Consensus_Dependency <- ifelse(rnai_res$type == "BOTH", "Y", "N")
  
  message("  ✓ Consensus dependencies: ", sum(rnai_res$Consensus_Dependency == "Y"))
  return(rnai_res)
}

# ==============================================================================
# FUNCTION: Add Drug-Gene Interaction Information
# ==============================================================================
#' Add drug-gene interaction information from DGIdb
#' 
#' @param rnai_res Differential dependency results
#' @return Updated rnai_res with DGI information
add_drug_gene_interactions <- function(rnai_res) {
  message("Adding drug-gene interaction information...")
  
  # Load DGI database
  DGI_db <- read.delim(FILES$dgidb)
  DGI_db <- subset(DGI_db, interaction_score != "NULL")
  DGI_db <- subset(DGI_db, approved %in% TRUE & gene_name %in% rnai_res$gene)
  DGI_db$interaction_score <- as.numeric(DGI_db$interaction_score)
  
  # Create DGI matrix
  DGI_db_mat <- reshape2::dcast(DGI_db, gene_name ~ gene_name, 
                                value.var = "interaction_score", 
                                fun.aggregate = max, na.rm = TRUE)
  rownames(DGI_db_mat) <- DGI_db_mat[, 1]
  DGI_db_mat <- DGI_db_mat[, -1]
  DGI_db_mat <- as.matrix(DGI_db_mat)
  mode(DGI_db_mat) <- "numeric"
  DGI_db_mat[is.infinite(DGI_db_mat)] <- NA
  
  # Create unique DGI data frame
  DGI_db_df_unique <- reshape2::melt(DGI_db_mat)
  DGI_db_df_unique <- DGI_db_df_unique[!is.na(DGI_db_df_unique$value), ]
  DGI_db_df_unique$Var1 <- as.character(DGI_db_df_unique$Var1)
  DGI_db_df_unique$Var2 <- as.character(DGI_db_df_unique$Var2)
  names(DGI_db_df_unique) <- c("gene", "gene", "score")
  rownames(DGI_db_df_unique) <- DGI_db_df_unique$gene
  
  # Add DGI information to results
  rnai_res$DGI <- ifelse(rnai_res$gene %in% rownames(DGI_db_df_unique), "Y", "N")
  rnai_res$DGI_SCORE <- DGI_db_df_unique[rnai_res$gene, "score"]
  rnai_res$DRUG <- DGI_db[match(rnai_res$DGI_SCORE, DGI_db$interaction_score), "drug_claim_name"]
  
  message("  ✓ Added DGI information for ", sum(rnai_res$DGI == "Y"), " genes")
  return(rnai_res)
}

# ------------------------------------------------------------------------------
# SECTION 2: Load Data
# ------------------------------------------------------------------------------

message("=", strrep("=", 78))
message("Section 2: Loading RNAi Data")
message("=", strrep("=", 78))

# Load RNAi data
rnai_data <- load_rnai_data()
RNAi_CL_Annots <- rnai_data$RNAi_CL_Annots
RNAi_D2 <- rnai_data$RNAi_D2
RNAi_D2_discrete <- rnai_data$RNAi_D2_discrete
RNAi_D2_zcore <- rnai_data$RNAi_D2_zcore

# Define cell line groups
ICLE_cells <- subset(RNAi_CL_Annots, Group %in% "ICLE")$Sample
NST_cells <- subset(RNAi_CL_Annots, Group %in% "NST")$Sample

# ------------------------------------------------------------------------------
# SECTION 3: SIMEM Differential Analysis
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("Section 3: SIMEM Differential Analysis")
message("=", strrep("=", 78))

# Check if SIMEM results exist, otherwise run analysis
if (!file.exists(FILES$rnai_simem)) {
  message("SIMEM results not found. Running SIMEM analysis...")
  run_simem_analysis()
} else {
  message("Loading existing SIMEM results...")
}

# Load SIMEM results
RNAi_SIMEM_DA <- read.delim(FILES$rnai_simem)

rownames(RNAi_SIMEM_DA) <- RNAi_SIMEM_DA$symbol
RNAi_SIMEM_DA <- RNAi_SIMEM_DA[!duplicated(RNAi_SIMEM_DA$symbol), ]
message("  ✓ Loaded SIMEM results for ", nrow(RNAi_SIMEM_DA), " genes")

# ------------------------------------------------------------------------------
# SECTION 4: D2 Differential Analysis
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("Section 4: D2 Differential Analysis")
message("=", strrep("=", 78))

# Select features
all_genes <- intersect(rownames(RNAi_D2), rownames(RNAi_SIMEM_DA))
selectedFeatures <- all_genes

# Run D2 differential analysis
rnai_res <- run_d2_differential_analysis(RNAi_D2, ICLE_cells, NST_cells, selectedFeatures)

# ------------------------------------------------------------------------------
# SECTION 5: Compute Consensus Dependencies
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("Section 5: Computing Consensus Dependencies")
message("=", strrep("=", 78))

# Compute consensus
rnai_res <- compute_consensus_dependencies(rnai_res, RNAi_SIMEM_DA, alpha = 0.05, eff = 0.2)

# Add drug-gene interactions
rnai_res <- add_drug_gene_interactions(rnai_res)

# Save results
rnai_res_col_formatted <- rnai_res[, setdiff(colnames(rnai_res), 
                                             c("sig_d2", "sig_si", "d2_big", "si_big", "dir_agree"))]

columns <- c("Gene", "ICLE D2 Score Average", "NST D2 Score Average",
             "D2 Effect Size", "D2 P-value", "D2 FDR", 
             "SIMEM Effect Size", "SIMEM P-value", "SIMEM FDR",
             "Significance", "Group", "Group & Significance", "Consensus Dependency", 
             "Present in DGI", "DGI Score", "Drug")
colnames(rnai_res_col_formatted) <- columns
rnai_res_col_formatted <- rnai_res_col_formatted[, columns]

write.table(rnai_res_col_formatted, 
            file = file.path(dep_dir, "RNAi_ILC_Consensus_Dependencies.tsv"), 
            row.names = FALSE, sep = "\t")
message("  ✓ Saved results to: ", file.path(dep_dir, "RNAi_ILC_Consensus_Dependencies.tsv"))

# Extract consensus significant dependencies
consensus_SIG <- subset(rnai_res, Consensus_Dependency == "Y")
consensus_SIG <- consensus_SIG[order(consensus_SIG$simem_pval, decreasing = FALSE), ]

# ------------------------------------------------------------------------------
# SECTION 6: Generate Fig 5B - Consensus Dependencies Scatter Plot
# ------------------------------------------------------------------------------

message("\n========================================")
message("Figure 5B: Consensus dependencies scatter plot")
message("========================================")

fig5b_consensus_dep_plt <- ggplot(rnai_res, aes(simem_DIFF, D2_DIFF)) +
  geom_point(data = subset(rnai_res, type == "NS") %>% mutate(Group_type = "NS"), 
             aes(color = Group_type), size = 0.1) + 
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.3) + 
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.3) + 
  geom_point(data = subset(rnai_res, type %in% c("BOTH", "D2 ONLY", "SIMEM ONLY")), 
             aes(color = Group_type), size = 1.5) + 
  scale_color_manual("Group", values = c(
    "ILC - BOTH" = "red", 
    "NST - BOTH" = "#395c9d", 
    "ILC - D2 ONLY" = "#F8BBD0", 
    "ILC - SIMEM ONLY" = "#BA68C8", 
    "NST - D2 ONLY" = "#7cc242", 
    "NST - SIMEM ONLY" = "#00c1f3", 
    "NS" = "gray80"
  )) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  ggrepel::geom_text_repel(
    data = rbind(
      head(subset(consensus_SIG[order(consensus_SIG$simem_DIFF + consensus_SIG$D2_DIFF, decreasing = TRUE), ], 
                  simem_DIFF > 0.2), 10),
      head(subset(consensus_SIG[order(consensus_SIG$simem_DIFF + consensus_SIG$D2_DIFF, decreasing = FALSE), ], 
                  simem_DIFF < -0.2), 10)
    ),
    size = 3, fontface = "bold", alpha = 0.7, 
    segment.color = "black",
    aes(label = gene), color = "black", box.padding = unit(3, "mm"),
    max.overlaps = 100000
  ) +
  theme_bw(20) + 
  theme(panel.grid.major = element_line(linetype = "dotted", size = 0.8), 
        panel.grid.minor = element_blank()) + 
  xlab("SIMEM Effect Size") + 
  ylab("D2 Effect Size") + 
  scale_x_continuous(breaks = seq(-1, 1, 0.2), labels = seq(-1, 1, 0.2)) + 
  scale_y_continuous(breaks = seq(-1, 1, 0.2), labels = seq(-1, 1, 0.2))

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("fig5b_consensus_dep_plt", fig5b_consensus_dep_plt, envir = .GlobalEnv)
message("  ✓ Fig 5B complete (fig5b_consensus_dep_plt assigned).")

# ------------------------------------------------------------------------------
# SECTION 7: Generate SupFig 12 - Consensus Dependencies Heatmap
# ------------------------------------------------------------------------------

message("\n========================================")
message("SupFig 12: Consensus dependencies heatmap")
message("========================================")

# Define color function
colFun <- circlize::colorRamp2(
  c(-2, -1.5, -1, -0.5, -0.25, 0, 0.5, 0.75, 1),
  colors = c('#e10900', "#ec6a00", "#fab700", "#ffea84", "white", "white", 
             "#90CAF9", '#57a1d8', "#3b4796"), 
  space = 'RGB'
)

# Create consensus dependency matrix
consensus_dep_mat <- RNAi_D2[consensus_SIG$gene, c(ICLE_cells, NST_cells)]

# Create heatmap
heatmap_legend_param <- list(
  labels_gp = gpar(fontsize = 12), 
  title_gp = gpar(fontsize = 12), 
  direction = "horizontal"
)

set.seed(123)
supfig12_dep_ht <- Heatmap(
  as.matrix(consensus_dep_mat), 
  col = colFun, 
  name = "Dependency", 
  show_row_names = TRUE, 
  cluster_column_slices = FALSE, 
  cluster_row_slices = FALSE, 
  border = TRUE, 
  border_gp = gpar(col = "black", lwd = 0.8), 
  column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  row_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
  use_raster = TRUE, 
  raster_quality = 4, 
  width = unit(7.5, "cm"), 
  row_split = rnai_res[rownames(consensus_dep_mat), "Group"], 
  row_title = paste0(c("N = ", "N = "), table(rnai_res[rownames(consensus_dep_mat), "Group"])), 
  column_title = " ", 
  show_row_dend = TRUE, 
  show_column_dend = TRUE,
  column_split = RNAi_CL_Annots[colnames(consensus_dep_mat), "Group_factor"], 
  heatmap_legend_param = heatmap_legend_param
)

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("supfig12_dep_ht", supfig12_dep_ht, envir = .GlobalEnv)
message("  ✓ SupFig 12 complete (supfig12_dep_ht assigned).")

# ------------------------------------------------------------------------------
# SECTION 8: Generate Fig 5C - KEGG Pathway Enrichment
# ------------------------------------------------------------------------------

message("\n========================================")
message("Figure 5C: KEGG pathway enrichment")
message("========================================")

# Load KEGG gene sets
KEGG <- msigdbr::msigdbr(species = "Homo sapiens", db_species = "HS", 
                         collection = "C2", subcollection = "CP:KEGG_LEGACY")
KEGG$gs_name <- gsub(x = KEGG$gs_name, pattern = "KEGG_", replacement = "")
KEGG$gs_name <- gsub(x = KEGG$gs_name, pattern = "_", replacement = " ")
KEGG$gs_name <- tools::toTitleCase(tolower(KEGG$gs_name))
KEGG <- split(KEGG$gene_symbol, KEGG$gs_name)

# Get ILC and NST genes
rnai_res_path <- subset(rnai_res, type %in% c("D2 ONLY", "SIMEM ONLY", "BOTH"))
ICLE_genes <- subset(rnai_res_path, grepl(x = Group, pattern = "ILC"))$gene
NST_genes <- subset(rnai_res_path, grepl(x = Group, pattern = "NST"))$gene

# Run GSEA
KEGG_ALL <- doGSEA(gset_db = "KEGG", KEGG, list(NST = NST_genes, ILC = ICLE_genes), plot_top_n = 7)

# Create plots
plt1 <- KEGG_ALL$lobular$plt + 
  theme(text = element_text(size = 20), axis.title.x = element_text(size = 20, color = "black")) + 
  xlab(expression(-log[10](P-value))) + 
  ggtitle("")

plt2 <- KEGG_ALL$ductal$plt + 
  theme(text = element_text(size = 20), axis.title.x = element_text(size = 20, color = "black")) + 
  xlab(expression(-log[10](P-value))) + 
  ggtitle("")

fig5c_pathway_plt <- ggpubr::ggarrange(plotlist = list(plt1, plt2), nrow = 1, widths = c(1, 1))

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("fig5c_pathway_plt", fig5c_pathway_plt, envir = .GlobalEnv)
message("  ✓ Fig 5C complete (fig5c_pathway_plt assigned).")

# ------------------------------------------------------------------------------
# SECTION 9: Generate Fig 5D - Pathway Level Dependencies Heatmap
# ------------------------------------------------------------------------------

message("\n========================================")
message("Figure 5D: Pathway level dependencies heatmap")
message("========================================")

# Extract pathway genes
lob_paths <- split(KEGG_ALL$lobular$df$hits, KEGG_ALL$lobular$df$label)[
  (KEGG_ALL$lobular$df)[c(1:7), 'label']
]
duct_paths <- split(KEGG_ALL$ductal$df$hits, KEGG_ALL$ductal$df$label)[
  (KEGG_ALL$ductal$df)[c(1, 2, 5), 'label']
]

# Calculate pathway scores
path_scores <- list()
gene_list <- list()

for (x in names(c(lob_paths, duct_paths))) {
  genes <- trimws(unlist(strsplit(x = c(lob_paths, duct_paths)[[x]], split = " ,")))
  gene_list[[x]] <- genes
  path_scores[[x]] <- colMeans(as.matrix(RNAi_D2[genes, c(ICLE_cells, NST_cells)]))
}

path_scores <- t(as.data.frame(path_scores))
rownames(path_scores) <- gsub(x = rownames(path_scores), pattern = '[.]', replacement = " ")

# Create heatmap
set.seed(123)
fig5d_pathway_ht <- Heatmap(
  as.matrix(path_scores), 
  col = colFun, 
  name = "Dependency", 
  show_row_names = TRUE, 
  cluster_column_slices = FALSE, 
  cluster_row_slices = FALSE, 
  border = TRUE, 
  border_gp = gpar(col = "black", lwd = 0.8), 
  width = unit(7.5, "cm"),
  column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  row_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
  use_raster = TRUE, 
  raster_quality = 4,
  row_split = ifelse(rownames(path_scores) %in% names(lob_paths), "ICLE", "NST"), 
  column_title = " ", 
  show_row_dend = TRUE, 
  show_column_dend = TRUE,
  column_split = RNAi_CL_Annots[colnames(path_scores), "Group_factor"], 
  heatmap_legend_param = heatmap_legend_param
)

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("fig5d_pathway_ht", fig5d_pathway_ht, envir = .GlobalEnv)
message("  ✓ Fig 5D complete (fig5d_pathway_ht assigned).")

# ------------------------------------------------------------------------------
# SECTION 10: Generate Fig 5E - Druggable Dependencies Heatmap
# ------------------------------------------------------------------------------

message("\n========================================")
message("Figure 5E: Druggable dependencies heatmap")
message("========================================")

# Identify druggable genes
druggable_genes_ILC <- subset(rnai_res, Group == "ILC" & Consensus_Dependency == "Y" & DGI %in% "Y" & DGI_SCORE > 1)
druggable_genes_ILC$label <- paste0(druggable_genes_ILC$DRUG, " (", druggable_genes_ILC$gene, ")")

# Create dependency matrix
drug_dep_mat <- RNAi_D2[druggable_genes_ILC$gene, c(ICLE_cells, NST_cells)]
drug_dep_mat[is.na(drug_dep_mat)] <- 0

# Create heatmap
set.seed(123)
fig5e_drug_ht <- Heatmap(
  as.matrix(drug_dep_mat), 
  col = colFun, 
  name = "Dependency", 
  show_row_names = TRUE, 
  cluster_column_slices = FALSE, 
  cluster_row_slices = FALSE, 
  row_order = order(druggable_genes_ILC[rownames(drug_dep_mat), "DGI_SCORE"], decreasing = TRUE),
  row_labels = druggable_genes_ILC[rownames(drug_dep_mat), "label"], 
  width = unit(7, "cm"), height = unit(5, "cm"),
  border = TRUE, 
  border_gp = gpar(col = "black", lwd = 0.8), 
  column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  row_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
  use_raster = TRUE, 
  raster_quality = 4, 
  row_split = rnai_res[rownames(drug_dep_mat), "Group"], 
  column_title = " ", 
  show_row_dend = TRUE, 
  show_column_dend = TRUE,
  column_split = RNAi_CL_Annots[colnames(drug_dep_mat), "Group_factor"], 
  heatmap_legend_param = heatmap_legend_param
)

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("fig5e_drug_ht", fig5e_drug_ht, envir = .GlobalEnv)
message("  ✓ Fig 5E complete (fig5e_drug_ht assigned).")

message("\n========================================")
message("  Figure 5 (5B–5E) and SupFig 12 complete")
message("========================================")
message("  - Fig 5B: fig5b_consensus_dep_plt")
message("  - SupFig 12: supfig12_dep_ht")
message("  - Fig 5C: fig5c_pathway_plt")
message("  - Fig 5D: fig5d_pathway_ht")
message("  - Fig 5E: fig5e_drug_ht")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/22_Fig5_RNAi_Differential_Dependencies.R")
#
# load_all_icle_data(load_external = TRUE)
# # Script runs on source. Uses FILES$rnai_annots, rnai_d2, rnai_simem, dgidb, etc.
#
# Outputs: 
#   - RNAi_ILC_Consensus_Dependencies.tsv in DIRS$results_sub$dependencies
#   - All figures assigned to global environment for saving in Main_Data_Analysis.Rmd
