# ==============================================================================
# Script 01: Supplementary Figure 1 - Genotypic Similarity Analysis
# ==============================================================================
# Description: Performs genotypic similarity analysis between ICLE and Marcotte
#              cell line cohorts using SNP array data. Generates similarity
#              heatmap and clustering visualizations.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - CL_Annots: Cell line annotations
#   - SNP array data (loaded via data loading functions)
#
# Output:
#   - supFig_1: Genotypic similarity heatmap (assigned to .GlobalEnv)
#   - Clustering results
#
# Author: Osama Shiraz Shah
# ==============================================================================

########################################
# Libraries
########################################
library(magrittr)
library(ComplexHeatmap)
library(viridis)
library(grid)

########################################
# Check for Required Data Objects
########################################
# Ensure config and required data objects are loaded
if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}

# Check if CL_Annots exists, load if missing
if (!exists("CL_Annots", envir = .GlobalEnv)) {
  message("CL_Annots not found. Loading annotations...")
  if (!exists("load_all_icle_data", envir = .GlobalEnv)) {
    source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  }
  # Load at least annotations (step 0)
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "00_load_annotations.R"))
  annot_data <- load_cell_line_annotations()
  assign("CL_Annots", annot_data$CL_Annots, envir = .GlobalEnv)
  assign("CL_Annots_simple", annot_data$CL_Annots_simple, envir = .GlobalEnv)
}

# Check for helper functions
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}


########################################
# 1. Load Data 
########################################

load_and_merge_genotype_data <- function(ICLE_path, Marcotte_path) {
  # Load raw genotype call matrices
  if (!exists("ICLE_Genotype_Mat")) load(ICLE_path)        # Loads ICLE_Genotype_Mat
  if (!exists("Marcotte_Genotype_Mat")) load(Marcotte_path)    # Loads Marcotte_Genotype_Mat
  
  # Filter to shared SNPs (rows), relabel columns by cohort
  common_snps <- intersect(rownames(ICLE_Genotype_Mat), rownames(Marcotte_Genotype_Mat))
  if (length(common_snps) == 0) stop("No common SNPs found between datasets.")
  ICLE_Genotype_Mat     <- ICLE_Genotype_Mat[common_snps, ]
  Marcotte_Genotype_Mat <- Marcotte_Genotype_Mat[common_snps, ]
  
  # Label columns with suffixes to track cohort
  colnames(ICLE_Genotype_Mat)     <- paste0(colnames(ICLE_Genotype_Mat), "-I")
  colnames(Marcotte_Genotype_Mat) <- paste0(colnames(Marcotte_Genotype_Mat), "-M")
  
  # Merge and remove any rows with missing values
  merged <- cbind(Marcotte_Genotype_Mat, ICLE_Genotype_Mat)
  merged <- na.omit(merged)
  
  return(merged)
}


########################################
# 2. Convert Genotype Calls to Integers
########################################

convert_genotype_calls <- function(mat) {
  # Convert AA, AB, BB to integers: 0, 1, 2
  mat <- as.matrix(mat)
  temp <- mat
  mat[temp == "AA"] <- 0
  mat[temp == "AB"] <- 1
  mat[temp == "BB"] <- 2
  mode(mat) <- "integer"
  return(mat)
}


########################################
# 3. Compute Pairwise Genotypic Similarity
########################################

compute_similarity_matrix <- function(geno_mat, random_sample = NULL) {
  # Compute similarity matrix (percentage of matching SNPs between each sample pair)
  
  if (!is.null(random_sample)) {
    geno_mat <- geno_mat[sample(rownames(geno_mat), random_sample), ]
  }
  
  storage.mode(geno_mat) <- "character"
  
  similarity_matrix <- outer(
    1:ncol(geno_mat),
    1:ncol(geno_mat),
    Vectorize(function(i, j) {
      mean(geno_mat[, i] == geno_mat[, j], na.rm = TRUE) * 100
    })
  )
  
  colnames(similarity_matrix) <- colnames(geno_mat)
  rownames(similarity_matrix) <- colnames(geno_mat)
  return(similarity_matrix)
}


########################################
# 4. Generate Genotypic Heatmap
########################################

plot_genotype_heatmap <- function(similarity_matrix, CL_Annots, annot_cols, heatmap_legend_param) {
  # Build study annotation
  sample_labels <- sapply(
    colnames(similarity_matrix),
    function(x) strsplit(x, "-", fixed = TRUE)[[1]][1]
  ) %>% as.character()
  
  study_annot <- c(
    rep("Marcotte", sum(grepl("-M$", colnames(similarity_matrix)))),
    rep("ICLE",     sum(grepl("-I$", colnames(similarity_matrix))))
  )
  names(study_annot) <- colnames(similarity_matrix)
  
  # Select ILC-like/ILC samples from metadata
  selected_samples <- sample_labels %in% subset(CL_Annots, Histology %in% c("ILC", "ILC-like"))$Sample

  sim_mat <- similarity_matrix[selected_samples, selected_samples, drop = FALSE]

  # Highlight samples with high similarity to at least one other
  # highly_similar <- names(which(rowSums(sim_mat > 90) > 1))
  # col_highlight <- ifelse(colnames(sim_mat) %in% highly_similar, "#33691E", "#A1887F")
  col_highlight <- "black"

  
  # Generate and draw heatmap
  ht <- Heatmap(
    matrix = t(sim_mat),
    name = "Genotypic Similarity",
    col = viridis(100),
    border = TRUE, use_raster = T,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
    show_row_dend = TRUE, show_column_dend = TRUE,
    column_dend_side = "bottom",
    show_row_names = TRUE, show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = gpar(fontsize = 12, col = col_highlight),
    width = unit(10, "cm"), height = unit(10, "cm"),
    heatmap_legend_param = heatmap_legend_param,
    top_annotation = HeatmapAnnotation(
      Study = study_annot[colnames(sim_mat)],
      show_annotation_name = TRUE,
      annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
      annotation_legend_param = heatmap_legend_param,
      col = list(
        Study = annot_cols$Study,
        Histology = annot_cols$Histology
      )
    )
  )
  
  return(ht)
}


########################################
# 5. Run Complete Pipeline
########################################

run_genotype_similarity_pipeline <- function(ICLE_path = NULL, Marcotte_path = NULL, 
                                               CL_Annots = NULL, annot_cols = NULL, 
                                               random_sample = NULL, 
                                               heatmap_legend_param = NULL) {
  
  message("═══════════════════════════════════════════════════════")
  message("  Supplementary Figure 1: Genotype Similarity Analysis")
  message("═══════════════════════════════════════════════════════\n")
  
  # Load CL_Annots if not provided
  if (is.null(CL_Annots)) {
    if (exists("CL_Annots", envir = .GlobalEnv)) {
      CL_Annots <- get("CL_Annots", envir = .GlobalEnv)
    } else {
      message("  CL_Annots not found. Loading annotations...")
      source(file.path(DIRS$scripts$helpers, "Data_Loading", "00_load_annotations.R"))
      annot_data <- load_cell_line_annotations()
      CL_Annots <- annot_data$CL_Annots
      assign("CL_Annots", CL_Annots, envir = .GlobalEnv)
    }
  }
  
  # Load annot_cols if not provided
  if (is.null(annot_cols)) {
    if (exists("annot_cols", envir = .GlobalEnv)) {
      annot_cols <- get("annot_cols", envir = .GlobalEnv)
    } else {
      source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
      annot_cols <- get("annot_cols", envir = .GlobalEnv)
    }
  }
  
  # Load heatmap_legend_param if not provided
  if (is.null(heatmap_legend_param)) {
    if (exists("heatmap_legend_param", envir = .GlobalEnv)) {
      heatmap_legend_param <- get("heatmap_legend_param", envir = .GlobalEnv)
    } else {
      source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
      heatmap_legend_param <- get("heatmap_legend_param", envir = .GlobalEnv)
    }
  }
  
  # Use config paths if not provided
  if (is.null(ICLE_path) && exists("FILES", envir = .GlobalEnv)) {
    ICLE_path <- FILES$icle_genotype
  }
  if (is.null(Marcotte_path) && exists("FILES", envir = .GlobalEnv)) {
    Marcotte_path <- FILES$marcotte_genotype
  }
  
  # Load and Merge genotype data
  message("  Step 1/4: Loading and merging genotype data...")
  merged_data <- load_and_merge_genotype_data(ICLE_path, Marcotte_path)
  message("  ✓ Genotype data merged: ", nrow(merged_data), " SNPs x ", ncol(merged_data), " samples\n")
  
  # Convert calls
  message("  Step 2/4: Converting genotype calls to integers...")
  merged_data <- convert_genotype_calls(merged_data)
  message("  ✓ Genotype calls converted\n")
  
  # Compute pairwise similarity matrix
  message("  Step 3/4: Computing genotypic similarity matrix...")
  similarity_matrix <- compute_similarity_matrix(merged_data, random_sample)
  message("  ✓ Similarity matrix computed\n")
  
  # Plot heatmap for selected samples
  message("  Step 4/4: Generating genotypic similarity heatmap...")
  sim_ht <- plot_genotype_heatmap(similarity_matrix, CL_Annots, annot_cols, heatmap_legend_param)
  message("  ✓ Heatmap generated\n")
  
  message("═══════════════════════════════════════════════════════")
  message("  Genotype Similarity Analysis Complete")
  message("═══════════════════════════════════════════════════════")
  return(sim_ht)
}


supFig_1 <- run_genotype_similarity_pipeline(
  CL_Annots = if(exists("CL_Annots")) CL_Annots else NULL,
  annot_cols = if(exists("annot_cols")) annot_cols else NULL,
  random_sample = 10000,
  heatmap_legend_param = if(exists("heatmap_legend_param")) heatmap_legend_param else NULL
)

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/00_load_annotations.R")
# source("Helper_Scripts/01_SupFig1_Genotype_Similarity.R")
#
# annot_data <- load_cell_line_annotations()
# CL_Annots <- annot_data$CL_Annots
#
# supFig_1 <- run_genotype_similarity_pipeline(
#   CL_Annots = CL_Annots,
#   annot_cols = annot_cols,
#   random_sample = 10000,
#   heatmap_legend_param = heatmap_legend_param
# )
#
# pdf(file = file.path(DIRS$results, "SupFig1_Genotype_Similarity.pdf"), width = 10, height = 10)
# draw(supFig_1, merge_legends = TRUE)
# dev.off()