# ==============================================================================
# Load RNA-seq Data
# ==============================================================================
# This script loads and prepares RNA-seq expression data
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load RNA-seq Data
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing RNA-seq data and metadata
#' @export
load_rna_data <- function(config_loaded = TRUE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  # Check if data already exists
  if (exists("BRCA_CL_EXP", envir = .GlobalEnv)) {
    message("RNA-seq data already loaded (skipping)")
    return(invisible(NULL))
  }
  
  message("Loading RNA-seq data...")
  
  # ------------------------------------------------------------------------------
  # Load RNA-seq expression matrix (from preprocessed file; no 1-Datasets script)
  # ------------------------------------------------------------------------------
  rna_file <- FILES$rnaseq_data
  
  if (!file.exists(rna_file)) {
    stop("RNA-seq data file not found: ", rna_file)
  }
  
  # Load the RData file (contains BRCA_CL_EXP)
  load(rna_file, envir = .GlobalEnv)
  
  message("  ✓ RNA-seq data loaded successfully")
  
  # ------------------------------------------------------------------------------
  # Verify data was loaded
  # ------------------------------------------------------------------------------
  if (exists("BRCA_CL_EXP", envir = .GlobalEnv)) {
    BRCA_CL_EXP <- get("BRCA_CL_EXP", envir = .GlobalEnv)
    message("  ✓ Loaded: ", nrow(BRCA_CL_EXP), " genes x ", 
            ncol(BRCA_CL_EXP), " samples")
    
    # Check if data is log-transformed
    max_val <- max(BRCA_CL_EXP, na.rm = TRUE)
    if (max_val > 50) {
      message("  ℹ Data appears to be raw counts (max value: ", round(max_val), ")")
    } else {
      message("  ℹ Data appears to be log2-transformed (max value: ", round(max_val, 2), ")")
    }
  } else {
    stop("BRCA_CL_EXP object not found after loading RNA-seq data")
  }
  
  return(invisible(NULL))
}

#' Calculate ER Positive Score
#' 
#' @param BRCA_CL_EXP RNA expression matrix
#' @param er_genes Vector of ER-positive signature genes
#' @return Numeric vector of ER-positive scores per sample
#' @export
calculate_er_positive_score <- function(BRCA_CL_EXP, er_genes = NULL) {
  
  # Default ER-positive signature genes
  if (is.null(er_genes)) {
    er_genes <- c("ESR1", "PGR", "GATA3", "FOXA1", "XBP1", "AGR2", "TFF1", "TFF3")
  }
  
  # Intersect with available genes
  available_genes <- intersect(er_genes, rownames(BRCA_CL_EXP))
  
  if (length(available_genes) == 0) {
    warning("No ER-positive signature genes found in expression data")
    return(NULL)
  }
  
  message("  ✓ Calculating ER-positive score using ", length(available_genes), " genes")
  
  # Calculate mean expression across ER genes
  ERpos_score <- colMeans(BRCA_CL_EXP[available_genes, , drop = FALSE], na.rm = TRUE)
  
  return(ERpos_score)
}

#' Calculate ER Negative Score
#' 
#' @param BRCA_CL_EXP RNA expression matrix
#' @param er_neg_genes Vector of ER-negative signature genes
#' @return Numeric vector of ER-negative scores per sample
#' @export
calculate_er_negative_score <- function(BRCA_CL_EXP, er_neg_genes = NULL) {
  
  # Default ER-negative signature genes (basal markers)
  if (is.null(er_neg_genes)) {
    er_neg_genes <- c("KRT5", "KRT14", "KRT17", "EGFR", "CAV1", "CAV2")
  }
  
  # Intersect with available genes
  available_genes <- intersect(er_neg_genes, rownames(BRCA_CL_EXP))
  
  if (length(available_genes) == 0) {
    warning("No ER-negative signature genes found in expression data")
    return(NULL)
  }
  
  message("  ✓ Calculating ER-negative score using ", length(available_genes), " genes")
  
  # Calculate mean expression across ER-negative genes
  ERneg_score <- colMeans(BRCA_CL_EXP[available_genes, , drop = FALSE], na.rm = TRUE)
  
  return(ERneg_score)
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/03_load_rna_data.R")
# 
# load_rna_data()
# ERpos <- calculate_er_positive_score(BRCA_CL_EXP)
# ERneg <- calculate_er_negative_score(BRCA_CL_EXP)
