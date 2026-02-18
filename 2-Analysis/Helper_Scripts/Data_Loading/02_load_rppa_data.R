# ==============================================================================
# Load RPPA (Reverse Phase Protein Array) Data
# ==============================================================================
# This script loads and prepares RPPA protein expression data
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load RPPA Data
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing RPPA data and metadata
#' @export
load_rppa_data <- function(config_loaded = TRUE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  # Check if data already exists
  if (exists("BRCA_CL_RPPA", envir = .GlobalEnv)) {
    message("RPPA data already loaded (skipping)")
    return(invisible(NULL))
  }
  
  message("Loading RPPA data...")
  
  # ------------------------------------------------------------------------------
  # Source the RPPA preparation script
  # Sources: 1-Datasets/ICLE/RPPA/Prepare_RPPA_Data.R
  # ------------------------------------------------------------------------------
  rppa_script <- file.path(DIRS$icle$rppa, "Prepare_RPPA_Data.R")
  
  if (!file.exists(rppa_script)) {
    stop("RPPA preparation script not found: ", rppa_script)
  }
  
  # Source the script (sets working directory temporarily)
  source(rppa_script, chdir = TRUE)
  
  message("  ✓ RPPA data loaded successfully")
  
  # ------------------------------------------------------------------------------
  # Verify objects were created
  # ------------------------------------------------------------------------------
  # BRCA_CL_RPPA is required
  if (!exists("BRCA_CL_RPPA", envir = .GlobalEnv)) {
    stop("BRCA_CL_RPPA not found after loading RPPA data")
  }
  
  # gene_antibody_map is optional (used for visualizations)
  if (!exists("gene_antibody_map", envir = .GlobalEnv)) {
    warning("gene_antibody_map not found. RPPA visualizations may not work correctly.")
  }
  
  # Return summary
  if (exists("BRCA_CL_RPPA", envir = .GlobalEnv)) {
    BRCA_CL_RPPA <- get("BRCA_CL_RPPA", envir = .GlobalEnv)
    message("  ✓ Loaded: ", nrow(BRCA_CL_RPPA), " proteins x ", 
            ncol(BRCA_CL_RPPA), " samples")
  }
  
  return(invisible(NULL))
}

#' Generate RPPA Visualizations
#' 
#' @param BRCA_CL_RPPA RPPA expression matrix
#' @param CL_Annots Cell line annotations
#' @param CL_Annots_simple Cell line annotations (simplified version)
#' @param gene_antibody_map Gene to antibody mapping
#' @param pathway_dir Directory containing pathway definitions
#' @param target_gene Gene to highlight (default: "CDH1")
#' @return List of visualization objects
#' @export
generate_rppa_visualizations <- function(BRCA_CL_RPPA, CL_Annots, CL_Annots_simple, 
                                         gene_antibody_map, pathway_dir = NULL, 
                                         target_gene = "CDH1") {
  
  if (is.null(pathway_dir) && exists("FILES")) {
    pathway_dir <- FILES$pathway_dir
  }
  
  # Check if rppa_data_visualizations function exists
  if (!exists("rppa_data_visualizations")) {
    warning("rppa_data_visualizations function not found. Skipping visualizations.")
    return(NULL)
  }
  
  message("Generating RPPA visualizations...")
  
  rppa_plots <- rppa_data_visualizations(
    BRCA_CL_RPPA, 
    CL_Annots, 
    CL_Annots_simple, 
    gene_antibody_map, 
    pathway_dir, 
    target_gene
  )
  
  message("  ✓ RPPA visualizations generated")
  
  return(rppa_plots)
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/02_load_rppa_data.R")
# 
# load_rppa_data()
# # Note: The function checks if data is already loaded and skips if present
#
# rppa_plots <- generate_rppa_visualizations(BRCA_CL_RPPA, CL_Annots, CL_Annots_simple, gene_antibody_map, target_gene = "CDH1")
