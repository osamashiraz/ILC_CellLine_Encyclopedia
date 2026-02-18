# ==============================================================================
# Load Single Nucleotide Variation (SNV) Data
# ==============================================================================
# This script loads and prepares SNV/mutation data from WES
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load SNV Data
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing SNV data and metadata
#' @export
load_snv_data <- function(config_loaded = TRUE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  # Check if data already exists
  if (exists("BRCA_CL_SNV", envir = .GlobalEnv) && exists("BRCA_CL_MAF", envir = .GlobalEnv)) {
    message("SNV data already loaded (skipping)")
    return(invisible(NULL))
  }
  
  message("Loading SNV/mutation data...")
  
  # ------------------------------------------------------------------------------
  # Source the SNV preparation script
  # Sources: 1-Datasets/ICLE/WES/Prepare_SNV_Data.R
  # ------------------------------------------------------------------------------
  snv_script <- file.path(DIRS$icle$wes, "Prepare_SNV_Data.R")
  
  if (!file.exists(snv_script)) {
    stop("SNV preparation script not found: ", snv_script)
  }
  
  # Source the script (sets working directory temporarily)
  source(snv_script, chdir = TRUE)
  
  # message("  ✓ SNV data loaded successfully")
  
  # ------------------------------------------------------------------------------
  # Verify objects were created
  # ------------------------------------------------------------------------------
  expected_objects <- c("BRCA_CL_SNV", "BRCA_CL_MAF")
  
  loaded_objects <- c()
  for (obj in expected_objects) {
    if (exists(obj, envir = .GlobalEnv)) {
      loaded_objects <- c(loaded_objects, obj)
      obj_data <- get(obj, envir = .GlobalEnv)
      if (is.data.frame(obj_data) || is.matrix(obj_data)) {
        message("  ✓ Loaded ", obj, ": ", nrow(obj_data), " entries")
      }
    }
  }
  
  if (length(loaded_objects) == 0) {
    warning("No expected SNV objects found in global environment")
  }
  
  return(invisible(NULL))
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/05_load_snv_data.R")
# 
# load_snv_data()
