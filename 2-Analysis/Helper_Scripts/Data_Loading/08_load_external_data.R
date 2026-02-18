# ==============================================================================
# Load External Datasets
# ==============================================================================
# This script loads external datasets (TCGA, MSK, CCLE, etc.)
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load TCGA Data
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing TCGA data
#' @export
load_tcga_data <- function(config_loaded = TRUE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  # ------------------------------------------------------------------------------
  # Check if required downstream objects already exist in GlobalEnv
  # ------------------------------------------------------------------------------
  # Core required objects (annotations are critical)
  required_tcga_objects <- c("TCGA_Annots")
  
  # Optional but commonly used objects
  optional_tcga_objects <- c("TCGA_BRCA_GAM", "TCGA_BRCA_DNAm", "TCGA_BRCA_EXP_log2CPM", 
                             "TCGA_BRCA_RPPA", "TCGA_BRCA_CN")
  
  # Check for required objects
  missing_required <- required_tcga_objects[!sapply(required_tcga_objects, exists, envir = .GlobalEnv)]
  existing_optional <- optional_tcga_objects[sapply(optional_tcga_objects, exists, envir = .GlobalEnv)]
  
  # If all required objects exist, skip loading
  if (length(missing_required) == 0) {
    message("TCGA data already loaded (all required objects present in GlobalEnv)")
    message("  ✓ TCGA_Annots: present")
    if (length(existing_optional) > 0) {
      message("  ✓ Also present: ", paste(existing_optional, collapse = ", "))
    }
    return(invisible(NULL))
  }
  
  # Report which objects are missing
  if (length(missing_required) > 0) {
    message("Loading TCGA data (missing required objects: ", paste(missing_required, collapse = ", "), ")...")
  }
  
  # ------------------------------------------------------------------------------
  # Source the TCGA loading script
  # Sources: 1-Datasets/External/TCGA/Load_TCGA_Data.R
  # ------------------------------------------------------------------------------
  tcga_script <- file.path(DIRS$external$tcga, "Load_TCGA_Data.R")
  
  if (!file.exists(tcga_script)) {
    warning("TCGA loading script not found: ", tcga_script)
    return(invisible(NULL))
  }
  
  # Source the script (sets working directory temporarily)
  source(tcga_script, chdir = TRUE)
  
  # ------------------------------------------------------------------------------
  # Verify objects were created and assign to global environment
  # ------------------------------------------------------------------------------
  # Ensure TCGA_Annots is in global environment (critical)
  if (!exists("TCGA_Annots", envir = .GlobalEnv)) {
    stop("TCGA_Annots not found after loading TCGA data. This is a required object.")
  } else {
    # Ensure it's assigned to global environment
    assign("TCGA_Annots", get("TCGA_Annots", envir = .GlobalEnv), envir = .GlobalEnv)
    tcga_annots <- get("TCGA_Annots", envir = .GlobalEnv)
    message("  ✓ TCGA_Annots loaded: ", nrow(tcga_annots), " samples")
  }
  
  # Report on other loaded objects
  loaded_optional <- c()
  for (obj in optional_tcga_objects) {
    if (exists(obj, envir = .GlobalEnv)) {
      loaded_optional <- c(loaded_optional, obj)
      obj_data <- get(obj, envir = .GlobalEnv)
      if (is.data.frame(obj_data) || is.matrix(obj_data)) {
        message("  ✓ Loaded ", obj, ": ", nrow(obj_data), " x ", ncol(obj_data))
      } else {
        message("  ✓ Loaded ", obj)
      }
    }
  }
  
  message("  ✓ TCGA data loaded successfully")
  
  return(invisible(NULL))
}

#' Load MSK Data
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing MSK data
#' @export
load_msk_data <- function(config_loaded = TRUE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  # ------------------------------------------------------------------------------
  # Check if required downstream objects already exist in GlobalEnv
  # ------------------------------------------------------------------------------
  # Core required objects (annotations are critical)
  required_msk_objects <- c("MSK_Annots")
  
  # Optional but commonly used objects
  optional_msk_objects <- c("MSK_BRCA_GAM", "MSK_BRCA_MET_GAM")
  
  # Check for required objects
  missing_required <- required_msk_objects[!sapply(required_msk_objects, exists, envir = .GlobalEnv)]
  existing_optional <- optional_msk_objects[sapply(optional_msk_objects, exists, envir = .GlobalEnv)]
  
  # If all required objects exist, skip loading
  if (length(missing_required) == 0) {
    message("MSK data already loaded (all required objects present in GlobalEnv)")
    message("  ✓ MSK_Annots: present")
    if (length(existing_optional) > 0) {
      message("  ✓ Also present: ", paste(existing_optional, collapse = ", "))
    }
    return(invisible(NULL))
  }
  
  # Report which objects are missing
  if (length(missing_required) > 0) {
    message("Loading MSK data (missing required objects: ", paste(missing_required, collapse = ", "), ")...")
  }
  
  # ------------------------------------------------------------------------------
  # Source the MSK loading script
  # Sources: 1-Datasets/External/MSK/Load_MSK_Data.R
  # ------------------------------------------------------------------------------
  msk_script <- file.path(DIRS$external$msk, "Load_MSK_Data.R")
  
  if (!file.exists(msk_script)) {
    warning("MSK loading script not found: ", msk_script)
    return(invisible(NULL))
  }
  
  # Source the script (sets working directory temporarily)
  source(msk_script, chdir = TRUE)
  
  # ------------------------------------------------------------------------------
  # Verify objects were created and assign to global environment
  # ------------------------------------------------------------------------------
  # Ensure MSK_Annots is in global environment (critical)
  if (!exists("MSK_Annots", envir = .GlobalEnv)) {
    stop("MSK_Annots not found after loading MSK data. This is a required object.")
  } else {
    # Ensure it's assigned to global environment
    assign("MSK_Annots", get("MSK_Annots", envir = .GlobalEnv), envir = .GlobalEnv)
    msk_annots <- get("MSK_Annots", envir = .GlobalEnv)
    message("  ✓ MSK_Annots loaded: ", nrow(msk_annots), " samples")
  }
  
  # Report on other loaded objects
  loaded_optional <- c()
  for (obj in optional_msk_objects) {
    if (exists(obj, envir = .GlobalEnv)) {
      loaded_optional <- c(loaded_optional, obj)
      obj_data <- get(obj, envir = .GlobalEnv)
      if (is.data.frame(obj_data) || is.matrix(obj_data)) {
        message("  ✓ Loaded ", obj, ": ", nrow(obj_data), " x ", ncol(obj_data))
      } else {
        message("  ✓ Loaded ", obj)
      }
    }
  }
  
  message("  ✓ MSK data loaded successfully")
  
  return(invisible(NULL))
}

#' Load All External Datasets
#' 
#' @param datasets Character vector of datasets to load (default: all)
#' @return List containing all external data
#' @export
load_external_data <- function(datasets = c("tcga", "msk")) {
  
  message("Loading external datasets...")
  
  results <- list()
  
  if ("tcga" %in% tolower(datasets)) {
    results$tcga <- load_tcga_data()
  }
  
  if ("msk" %in% tolower(datasets)) {
    results$msk <- load_msk_data()
  }
  
  message("✓ External data loading complete\n")
  
  return(invisible(results))
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/08_load_external_data.R")
# 
# load_external_data()
# # Or load specific datasets:
# load_tcga_data()
# load_msk_data()
