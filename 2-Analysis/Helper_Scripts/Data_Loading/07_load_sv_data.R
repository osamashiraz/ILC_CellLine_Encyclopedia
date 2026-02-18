# ==============================================================================
# Load Structural Variation (SV) Data from Bionano Optical Genome Mapping
# ==============================================================================
# This script loads and prepares structural variation data including:
#   - Bionano SV calls (deletions, duplications, translocations, inversions)
#   - Putative gene fusions
#   - SV filtering and standardization
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load Structural Variation Data
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @param recalculate Logical, whether to recalculate from raw data (default: FALSE)
#' @return List containing ICLE_SV data frame and summary statistics
#' @export
load_sv_data <- function(config_loaded = TRUE, recalculate = FALSE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  # Check if data already exists (unless recalculating)
  if (!recalculate && exists("ICLE_SV", envir = .GlobalEnv)) {
    message("SV data already loaded (skipping)")
    return(list(
      ICLE_SV = get("ICLE_SV", envir = .GlobalEnv),
      metadata = list(load_time = Sys.time(), already_loaded = TRUE)
    ))
  }
  
  message("Loading structural variation data...")
  
  # ------------------------------------------------------------------------------
  # 1. Check for preprocessed data
  # ------------------------------------------------------------------------------
  sv_script_path <- file.path(DIRS$icle$bionano, "2_Structural_Variations", "Prepare_SV_Data.R")
  output_file <- file.path(DIRS$icle$bionano, "2_Structural_Variations", "ICLE_SV_Filtered.csv")
  
  if (!file.exists(sv_script_path)) {
    stop("SV preparation script not found: ", sv_script_path)
  }
  
  # ------------------------------------------------------------------------------
  # 2. Load or Generate SV Data
  # Sources (when generating): 1-Datasets/ICLE/Bionano/2_Structural_Variations/Prepare_SV_Data.R
  # ------------------------------------------------------------------------------
  if (file.exists(output_file) && !recalculate) {
    message("  ✓ Loading preprocessed SV data from: ", basename(output_file))
    ICLE_SV <- read.csv(output_file, stringsAsFactors = FALSE)
  } else {
    message("  ⟳ Generating SV data from raw files...")
    message("    This may take a few minutes...")
    
    # Source the preparation script
    current_dir <- getwd()
    setwd(file.path(DIRS$icle$bionano, "2_Structural_Variations"))
    source("Prepare_SV_Data.R", chdir = TRUE)
    setwd(current_dir)
    
    if (!exists("ICLE_SV")) {
      stop("Failed to load ICLE_SV data from preparation script")
    }
    
    message("  ✓ SV data generated and saved")
  }
  
  # ------------------------------------------------------------------------------
  # 3. Ensure correct data types
  # ------------------------------------------------------------------------------
  ICLE_SV$RefStartPos <- as.integer(ICLE_SV$RefStartPos)
  ICLE_SV$RefEndPos <- as.integer(ICLE_SV$RefEndPos)
  
  # ------------------------------------------------------------------------------
  # 4. Validate SV Data
  # ------------------------------------------------------------------------------
  required_cols <- c("Sample", "Type", "RefcontigID1", "RefcontigID2", 
                     "RefStartPos", "RefEndPos", "SVsize")
  missing_cols <- setdiff(required_cols, colnames(ICLE_SV))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in SV data: ", paste(missing_cols, collapse = ", "))
  }
  
  message("  ✓ Loaded ", nrow(ICLE_SV), " structural variations")
  
  # ------------------------------------------------------------------------------
  # 5. Summary Statistics
  # ------------------------------------------------------------------------------
  sv_counts <- table(ICLE_SV$Sample)
  # message("\nSV counts per sample:")
  # print(data.frame(
  #   Sample = names(sv_counts), 
  #   SV_Count = as.numeric(sv_counts),
  #   row.names = NULL
  # ))
  
  sv_type_counts <- table(ICLE_SV$Type)
  # message("\nSV type distribution:")
  # print(sv_type_counts)
  
  # Check for fusions
  if ("PutativeGeneFusion" %in% colnames(ICLE_SV)) {
    n_fusions <- sum(ICLE_SV$PutativeGeneFusion != "-")
    message("  ✓ Loaded ", n_fusions, " nPutative gene fusions")
  }
  
  # ------------------------------------------------------------------------------
  # 6. Load SV visualization function
  # ------------------------------------------------------------------------------
  # Source the visualization function from the preparation script
  if (!exists("plot_sv_type_distribution", envir = .GlobalEnv)) {
    # message("\n  ⟳ Loading SV visualization functions...")
    source(sv_script_path, chdir = TRUE)
  }
  
  # ------------------------------------------------------------------------------
  # 7. Return Results
  # ------------------------------------------------------------------------------
  results <- list(
    ICLE_SV = ICLE_SV,
    metadata = list(
      n_sv_total = nrow(ICLE_SV),
      n_samples = length(unique(ICLE_SV$Sample)),
      sv_types = names(sv_type_counts),
      sv_counts_by_sample = sv_counts,
      load_time = Sys.time()
    )
  )
  
  return(results)
}

#' Get SV counts for specific samples
#' 
#' @param ICLE_SV SV data frame
#' @param samples Character vector of sample names
#' @return Data frame with SV counts by type
#' @export
get_sv_counts <- function(ICLE_SV, samples = NULL) {
  
  if (!is.null(samples)) {
    ICLE_SV <- ICLE_SV[ICLE_SV$Sample %in% samples, ]
  }
  
  sv_summary <- aggregate(
    list(Count = rep(1, nrow(ICLE_SV))),
    by = list(Sample = ICLE_SV$Sample, Type = ICLE_SV$Type),
    FUN = sum
  )
  
  return(sv_summary)
}

#' Get putative fusions
#' 
#' @param ICLE_SV SV data frame
#' @param samples Character vector of sample names (optional)
#' @return Data frame with fusion information
#' @export
get_fusions <- function(ICLE_SV, samples = NULL) {
  
  if (!"PutativeGeneFusion" %in% colnames(ICLE_SV)) {
    stop("PutativeGeneFusion column not found in SV data")
  }
  
  fusions <- ICLE_SV[ICLE_SV$PutativeGeneFusion != "-", ]
  
  if (!is.null(samples)) {
    fusions <- fusions[fusions$Sample %in% samples, ]
  }
  
  return(fusions)
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/07_load_sv_data.R")
# 
# sv_data <- load_sv_data()
# ICLE_SV <- sv_data$ICLE_SV
# 
# # Get SV counts for specific samples
# icle_samples <- c("MDAMB134VI", "BCK4", "HCC2185")
# sv_counts <- get_sv_counts(ICLE_SV, icle_samples)
# 
# # Get fusions
# fusions <- get_fusions(ICLE_SV)
sv_data <- load_sv_data()
assign("ICLE_SV", sv_data$ICLE_SV, envir = .GlobalEnv)