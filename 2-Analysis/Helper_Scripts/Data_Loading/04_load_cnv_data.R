# ==============================================================================
# Load Copy Number Variation (CNV) Data
# ==============================================================================
# This script loads and prepares CNV data from CytoSNP arrays
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load CNV Data
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing CNV data and metadata
#' @export
load_cnv_data <- function(config_loaded = TRUE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  # # ------------------------------------------------------------------------------
  # # Check if required downstream objects already exist in GlobalEnv
  # # ------------------------------------------------------------------------------
  # required_objects <- c("BRCA_CL_GISTIC")
  # missing_objects <- required_objects[!sapply(required_objects, exists, envir = .GlobalEnv)]
  # 
  # # If all required objects exist, skip loading
  # if (length(missing_objects) == 0) {
  #   message("CNV data already loaded (all required objects present in GlobalEnv)")
  #   message("  ✓ BRCA_CL_GISTIC: present")
  #   return(invisible(NULL))
  # }
  # 
  # # Report which objects are missing
  # if (length(missing_objects) > 0) {
  #   message("Loading CNV data (missing objects: ", paste(missing_objects, collapse = ", "), ")...")
  # }
  
  # ------------------------------------------------------------------------------
  # Source the CNV preparation script
  # Sources: 1-Datasets/ICLE/CytoSNP/2_GenomeStudio/Prepare_SNP_Data.R
  # ------------------------------------------------------------------------------
  cnv_script <- file.path(DIRS$icle$cytosnp, "2_GenomeStudio", "Prepare_SNP_Data.R")
  
  if (!file.exists(cnv_script)) {
    stop("CNV preparation script not found: ", cnv_script)
  }
  
  # Source the script (sets working directory temporarily)
  source(cnv_script, chdir = TRUE)
  
  message("  ✓ CNV data loaded successfully")
  
  # ------------------------------------------------------------------------------
  # Verify objects were created
  # ------------------------------------------------------------------------------
  loaded_objects <- c()
  for (obj in required_objects) {
    if (exists(obj, envir = .GlobalEnv)) {
      loaded_objects <- c(loaded_objects, obj)
      obj_data <- get(obj, envir = .GlobalEnv)
      if (is.data.frame(obj_data) || is.matrix(obj_data)) {
        message("  ✓ Loaded ", obj, ": ", nrow(obj_data), " genes x ", 
                ncol(obj_data), " samples")
      } else {
        message("  ✓ Loaded ", obj)
      }
    } else {
      warning("Expected object not found after loading: ", obj)
    }
  }
  
  if (length(loaded_objects) == 0) {
    stop("No expected CNV objects found in global environment after loading")
  }
  
  # Check if we successfully loaded all required objects
  still_missing <- required_objects[!required_objects %in% loaded_objects]
  if (length(still_missing) > 0) {
    warning("Some required CNV objects are still missing: ", paste(still_missing, collapse = ", "))
  }
  
  return(invisible(NULL))
}

#' Load CNV Segmentation Data
#' 
#' @param seg_file Path to segmentation file (optional, uses config if NULL)
#' @return Data frame with segmentation data
#' @export
load_cnv_segments <- function(seg_file = NULL) {
  
  if (is.null(seg_file)) {
    if (!exists("FILES")) {
      stop("config.R must be loaded or seg_file path must be provided")
    }
    seg_file <- FILES$cnv_seg
  }
  
  if (!file.exists(seg_file)) {
    stop("Segmentation file not found: ", seg_file)
  }
  
  message("Loading CNV segmentation data...")
  
  # Read segmentation file
  seg_data <- read.delim(seg_file, stringsAsFactors = FALSE)
  
  message("  ✓ Loaded ", nrow(seg_data), " segments for ", 
          length(unique(seg_data$ID)), " samples")
  
  return(seg_data)
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/04_load_cnv_data.R")
# 
# load_cnv_data()
# seg_data <- load_cnv_segments()
