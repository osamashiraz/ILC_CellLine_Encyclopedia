# ==============================================================================
# Load DNA Methylation (DNAm) Data
# ==============================================================================
# This script loads and prepares DNAm data from Illumina 850K arrays
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load DNAm Data
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing DNAm data and metadata
#' @export
load_dnam_data <- function(config_loaded = TRUE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  # Check if data already exists
  if (exists("BRCA_CL_DNAm", envir = .GlobalEnv)) {
    message("DNAm data already loaded (skipping)")
    return(invisible(NULL))
  }
  
  message("Loading DNAm data...")
  
  # ------------------------------------------------------------------------------
  # Source the DNAm preparation script
  # Sources: 1-Datasets/ICLE/DNAm/Prepare_DNAm_Data.R
  # ------------------------------------------------------------------------------
  dnam_script <- file.path(DIRS$icle$dnam, "Prepare_DNAm_Data.R")
  
  if (!file.exists(dnam_script)) {
    stop("DNAm preparation script not found: ", dnam_script)
  }
  
  # Source the script (sets working directory temporarily)
  source(dnam_script, chdir = TRUE)
  
  message("  ✓ DNAm data loaded successfully")
  
  # ------------------------------------------------------------------------------
  # Verify objects were created
  # ------------------------------------------------------------------------------
  expected_objects <- c("BRCA_CL_DNAm")
  
  loaded_objects <- c()
  for (obj in expected_objects) {
    if (exists(obj, envir = .GlobalEnv)) {
      loaded_objects <- c(loaded_objects, obj)
      obj_data <- get(obj, envir = .GlobalEnv)
      message("  ✓ Loaded ", obj, ": ", nrow(obj_data), " probes x ", 
              ncol(obj_data), " samples")
    }
  }
  
  if (length(loaded_objects) == 0) {
    warning("No expected DNAm objects found in global environment")
  }
  
  return(invisible(NULL))
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/06_load_dnam_data.R")
# 
# load_dnam_data()
