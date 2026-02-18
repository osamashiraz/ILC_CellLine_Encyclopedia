# ==============================================================================
# Script 03 (Fig 1C): SET Signature Analysis
# ==============================================================================
# Description: Calculates SET scores and generates visualization heatmaps
#              showing ER-positive and ER-negative gene expression signatures
#              across ICLE cell lines.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - CL_Annots: Cell line annotations
#   - BRCA_CL_EXP: RNA expression data
#
# Output:
#   - SET signature scores
#   - Figure 1C heatmap visualizations
#
# Author: Osama Shiraz Shah
# ==============================================================================

library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)
library(dplyr)

# ==============================================================================
# FUNCTION: Load SET Signature Genes
# ==============================================================================
#' Load SET (Sensitivity to Endocrine Therapy) gene signatures
#'
#' @param set_file Path to SET gene signature file (optional, uses config if NULL)
#' @return List with ER positive and negative gene sets
#' @export
load_set_signature <- function(set_file = NULL) {
  
  # Check if SET_sig already exists in global environment
  if (exists("SET_sig", envir = .GlobalEnv)) {
    return(get("SET_sig", envir = .GlobalEnv))
  }
  
  if (is.null(set_file)) {
    if (!exists("FILES")) {
      stop("config.R must be loaded or set_file path must be provided")
    }
    set_file <- FILES$set_genes
  }
  
  if (!file.exists(set_file)) {
    stop("SET signature file not found: ", set_file)
  }
  
  SET_sig <- read_excel(set_file)
  SET_sig <- split(SET_sig$`Gene Symbol`, SET_sig$Type)
  
  return(SET_sig)
}

# ==============================================================================
# FUNCTION: Calculate SET Scores
# ==============================================================================
#' Calculate SET (Sensitivity to Endocrine Therapy) scores
#'
#' @param BRCA_CL_EXP RNA expression matrix (log2 transformed)
#' @param SET_sig List with Pos and Neg gene sets
#' @return Data frame with ERpos, ERneg, and SET scores
#' @export
calculate_set_scores <- function(BRCA_CL_EXP, SET_sig) {
  
  # Convert log2 to linear scale and calculate module scores
  ERpos_mat <- 2^(BRCA_CL_EXP[match(SET_sig$Pos, rownames(BRCA_CL_EXP)), , drop = FALSE]) - 1
  ERneg_mat <- 2^(BRCA_CL_EXP[match(SET_sig$Neg, rownames(BRCA_CL_EXP)), , drop = FALSE]) - 1
  
  # Calculate mean expression across genes (per sample)
  ERpos <- colMeans(ERpos_mat, na.rm = TRUE)
  ERneg <- colMeans(ERneg_mat, na.rm = TRUE)
  
  # Calculate SET score (difference scaled)
  SET <- as.numeric(scale(ERpos - ERneg))
  names(SET) <- colnames(BRCA_CL_EXP)
  
  # Calculate scaled version for visualization
  # SET_scaled <- 10 * ((SET - min(SET)) / (max(SET) - min(SET)) + 0.01)
  
  return(data.frame(
    Sample = names(SET),
    ERpos = ERpos,
    ERpos_z = scale(ERpos),
    ERneg = ERneg,
    ERneg_z = scale(ERneg),
    SET = SET,
    # SET_scaled = SET_scaled,
    row.names = names(SET)
  ))
}

# ==============================================================================
# FUNCTION: Create SET Signature Heatmap
# ==============================================================================
#' Generate heatmap for SET signature expression
#'
#' @param ICLE_SET_scores Expression matrix (signatures x samples)
#' @param CL_Annots Data frame with histology annotations
#' @param BRCA_CL_RPPA RPPA data matrix
#' @param annot_cols Named list of annotation colors
#' @param col_fun Color function for expression values (optional)
#' @param rppa_col_fun Color function for RPPA values (optional)
#' @param random_seed Random seed for reproducibility (default: 123)
#' @return ComplexHeatmap object
#' @export
create_set_signature_heatmap <- function(ICLE_SET_scores, 
                                         CL_Annots,
                                         BRCA_CL_RPPA,
                                         annot_cols,
                                         random_seed = 123) {
  
  # message("Creating SET signature heatmap...")
  
  # Get sample IDs
  sample_ids <- colnames(ICLE_SET_scores)
  
  # Create heatmap
  set.seed(random_seed)
  
  ht <- Heatmap(
    ICLE_SET_scores,
    name = "Mean Exp",
    col = annot_cols$RNA_zscore,
    
    # Dimensions
    height = unit(1, "cm"),
    width = unit(6, "cm"),
    
    # Borders
    border = TRUE,
    border_gp = gpar(col = "black"),
    
    # Row parameters
    row_names_gp = gpar(fontsize = 9, fontface = "bold", col = "black"),
    row_title_gp = gpar(fontsize = 9, col = "#757575", fontface = "bold"),
    show_row_dend = FALSE,
    
    # Column parameters
    column_names_gp = gpar(fontsize = 10, col = "#757575"),
    column_title = " ",
    show_column_names = TRUE,
    show_column_dend = FALSE,
    column_split = CL_Annots[sample_ids, "mRNA Subtypes"],
    cluster_column_slices = FALSE,
    
    # Legend
    heatmap_legend_param = list(direction = "vertical"),
    na_col = "gray",
    
    # Top annotations
    top_annotation = HeatmapAnnotation(
      Subtype = CL_Annots[sample_ids, "mRNA Subtypes"],
      # Histology = CL_Annots[sample_ids, "Histology"],
      `ER-A` = as.numeric(BRCA_CL_RPPA["ER-A", match(sample_ids, colnames(BRCA_CL_RPPA))]),
      PR = as.numeric(BRCA_CL_RPPA["PR", match(sample_ids, colnames(BRCA_CL_RPPA))]),
      `HER2-PY1248` = as.numeric(BRCA_CL_RPPA["HER2-PY1248", match(sample_ids, colnames(BRCA_CL_RPPA))]),
      
      col = c(
        annot_cols,
        list(
          Subtype = annot_cols$Subtypes,
          `ER-A` = annot_cols$RPPA,
          PR = annot_cols$RPPA,
          
          `HER2-PY1248` = annot_cols$RPPA
        )
      ),
      
      border = FALSE,
      annotation_legend_param = list(
        legend_direction = "vertical",
        legend_position = "bottom"
      ),
      annotation_name_gp = gpar(fontface = "bold", fontsize = 9, col = "#757575"),
      simple_anno_size = unit(3, "mm"),
      annotation_name_side = "right"
    )
  )
  
  # message("  ✓ SET signature heatmap created")
  
  return(ht)
}

# ==============================================================================
# FUNCTION: Create SET Score Heatmap
# ==============================================================================
#' Generate heatmap for SET scores
#'
#' @param SET_scores Data frame with SET scores (must have SET_scaled column)
#' @param sample_ids Sample IDs to include
#' @param col_fun Color function for SET scores (optional)
#' @param random_seed Random seed for reproducibility (default: 123)
#' @return ComplexHeatmap object
#' @export
create_set_score_heatmap <- function(SET_scores,
                                     sample_ids,
                                     random_seed = 123) {
  
  # message("Creating SET score heatmap...")
  
  # Prepare SET score matrix
  set_mat <- SET_scores[sample_ids, "SET", drop = F]
  
  # Create heatmap
  set.seed(random_seed)
  
  ht <- Heatmap(
    t(set_mat),
    name = "Sensitivity to Endocrine Therapy",
    col = annot_cols$SET,
    
    # Dimensions
    height = unit(1, "cm"),
    
    # Borders
    border = TRUE,
    
    # Row parameters
    row_names_gp = gpar(fontsize = 9, fontface = "bold", col = "black"),
    row_title_rot = 0,
    
    # Column parameters
    column_names_gp = gpar(fontsize = 10, col = "#757575"),
    show_column_names = TRUE,
    
    # Legend
    heatmap_legend_param = list(
      legend_direction = "vertical",
      legend_position = "bottom"
    )
  )
  
  message("  ✓ SET score heatmap created")
  
  return(ht)
}

# ==============================================================================
# FUNCTION: Generate Complete Figure 1C
# ==============================================================================
#' Generate complete Figure 1C with SET signature and scores
#'
#' @param BRCA_CL_EXP RNA expression matrix
#' @param SET_sig SET gene signatures (or NULL to load from config)
#' @param sample_ids Sample IDs to include
#' @param CL_Annots Histology annotations
#' @param BRCA_CL_RPPA RPPA data
#' @param annot_cols Annotation color scheme
#' @param random_seed Random seed for reproducibility (default: 123)
#' @return List with \code{SET_scores} (data frame) and \code{SET_heatmap} (ComplexHeatmap)
#' @export
generate_SET_score_heatmap <- function(BRCA_CL_EXP,
                               SET_sig = NULL,
                               sample_ids,
                               CL_Annots,
                               BRCA_CL_RPPA,
                               annot_cols,
                               random_seed = 123) {

  # Load SET signatures if not provided
  message("  Step 1/4: Loading SET signature genes...")
  SET_sig <- load_set_signature()
  message("  ✓ Loaded ", length(SET_sig$Pos), " ER+ genes and ", length(SET_sig$Neg), " ER- genes\n")
  
  # Calculate SET scores
  message("  Step 2/4: Calculating SET scores...")
  SET_scores <- calculate_set_scores(BRCA_CL_EXP, SET_sig)
  rownames(SET_scores) <- SET_scores$Sample
  SET_scores <- SET_scores[,-1]
  message("  ✓ SET scores calculated for ", nrow(SET_scores), " samples\n")
  
  # Create SET signature heatmap
  message("  Step 4/4: Creating SET signature heatmaps...")
  sig_heatmap <- create_set_signature_heatmap(
    t(SET_scores[sample_ids, c("ERpos_z", "ERneg_z")]),
    CL_Annots,
    BRCA_CL_RPPA,
    annot_cols,
    random_seed = random_seed
  )
  
  # Create SET score heatmap
  score_heatmap <- create_set_score_heatmap(
    SET_scores, sample_ids,
    random_seed = random_seed
  )
  
  # Combine heatmaps
  SET_heatmap <- sig_heatmap %v% score_heatmap
  message("  ✓ Heatmaps created and combined\n")
  
  return(list(
    SET_scores = SET_scores,
    SET_heatmap = SET_heatmap
  ))
}



# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/02_Fig1_SupFig2_3_4_Molecular_Subtyping.R")
# 
# load_all_icle_data(load_external = FALSE)
# # run 01_Molecular_Subtyping if CL_Annots etc. not yet available
# 
# set_res <- generate_SET_score_heatmap(
#   BRCA_CL_EXP,
#   SET_sig = NULL,
#   sample_ids = subset(CL_Annots, Study == "ICLE")$Name,
#   CL_Annots,
#   BRCA_CL_RPPA,
#   annot_cols,
#   random_seed = 123
# )
# set_heatmap <- set_res$SET_heatmap
# draw(set_heatmap, merge_legends = TRUE)
#
# pdf(file.path(DIRS$results_sub$molecular_subtyping, "Fig1C_SET_Signature.pdf"), width = 8, height = 6)
# draw(set_heatmap, merge_legends = TRUE)
# dev.off()

