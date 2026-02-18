# ==============================================================================
# Script 01: Master Data Loading Orchestrator
# ==============================================================================
# This script orchestrates all data loading in the correct order
# as defined in Main_Data_Analysis.Rmd
#
# Execution Order:
#   0. Annotations (required for all subsequent steps)
#   1. RPPA (Protein)
#   2. RNA-seq (mRNA Expression)
#   3. Structural Variations (Bionano)
#   4. Copy Number Variations (CytoSNP)
#   5. Single Nucleotide Variations (WES)
#   6. DNA Methylation (DNAm)
#   7. External Datasets (TCGA, MSK, etc.)
#   8. Genomic Alteration Matrices (GAMs)
#
# Data sources: ICLE data from 1-Datasets/ICLE/... (paths in DIRS$icle, FILES);
#               external data from 1-Datasets/External/... via DIRS$external, FILES.
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load All ICLE Data in Correct Order
#' 
#' This function loads all ICLE datasets into the global environment.
#' Individual datasets are placed directly in global space (no list wrappers).
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @param load_external Logical, whether to load external datasets (default: TRUE)
#' @param verbose Logical, whether to print detailed messages (default: TRUE)
#' @return NULL (invisible) - all data objects are loaded into global environment
#' @export
load_all_icle_data <- function(config_loaded = TRUE, 
                                load_external = TRUE, 
                                verbose = TRUE) {
  
  # Verify config is loaded
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }

  load_status <- list()
  start_time <- Sys.time()
  
  # ------------------------------------------------------------------------------
  # STEP 0: Load Cell Line Annotations (REQUIRED)
  # ------------------------------------------------------------------------------
  if (exists("CL_Annots", envir = .GlobalEnv) && exists("CL_Annots_simple", envir = .GlobalEnv)) {
    if (verbose) message("STEP 0/8: Cell Line Annotations already loaded (skipping)\n")
    load_status$annotations <- "ALREADY_LOADED"
  } else {
    if (verbose) message("STEP 0/8: Loading Cell Line Annotations...")
    tryCatch({
      source(file.path(DIRS$scripts$helpers, "Data_Loading", "00_load_annotations.R"))
      annot_data <- load_cell_line_annotations()
      assign("CL_Annots", annot_data$CL_Annots, envir = .GlobalEnv)
      assign("CL_Annots_simple", annot_data$CL_Annots_simple, envir = .GlobalEnv)
      load_status$annotations <- "SUCCESS"
      if (verbose) message("  ✓ Annotations loaded\n")
    }, error = function(e) {
      load_status$annotations <- paste("FAILED:", e$message)
      if (verbose) message("  ✗ Error loading annotations: ", e$message, "\n")
    })
  }
  
  # ------------------------------------------------------------------------------
  # STEP 1: Load RPPA Data
  # ------------------------------------------------------------------------------
  if (exists("BRCA_CL_RPPA", envir = .GlobalEnv)) {
    if (verbose) message("STEP 1/8: RPPA data already loaded (skipping)\n")
    load_status$rppa <- "ALREADY_LOADED"
  } else {
    if (verbose) message("STEP 1/8: Loading RPPA (Protein) Data...")
    tryCatch({
      source(file.path(DIRS$scripts$helpers, "Data_Loading", "02_load_rppa_data.R"))
      load_rppa_data()
      load_status$rppa <- "SUCCESS"
      if (verbose) message("  ✓ RPPA data loaded\n")
    }, error = function(e) {
      load_status$rppa <- paste("FAILED:", e$message)
      if (verbose) message("  ✗ Error loading RPPA: ", e$message, "\n")
    })
  }
  
  # ------------------------------------------------------------------------------
  # STEP 2: Load RNA-seq Data
  # ------------------------------------------------------------------------------
  if (exists("BRCA_CL_EXP", envir = .GlobalEnv)) {
    if (verbose) message("STEP 2/8: RNA-seq data already loaded (skipping)\n")
    load_status$rnaseq <- "ALREADY_LOADED"
  } else {
    if (verbose) message("STEP 2/8: Loading RNA-seq (mRNA) Data...")
    tryCatch({
      source(file.path(DIRS$scripts$helpers, "Data_Loading", "03_load_rna_data.R"))
      load_rna_data()
      load_status$rnaseq <- "SUCCESS"
      if (verbose) message("  ✓ RNA-seq data loaded\n")
    }, error = function(e) {
      load_status$rnaseq <- paste("FAILED:", e$message)
      if (verbose) message("  ✗ Error loading RNA-seq: ", e$message, "\n")
    })
  }
  
  
  # ------------------------------------------------------------------------------
  # STEP 3: Load Copy Number Variation Data
  # ------------------------------------------------------------------------------
  # Check if all required CNV objects exist in GlobalEnv
  required_cnv_objects <- c("BRCA_CL_GISTIC")
  missing_cnv_objects <- required_cnv_objects[!sapply(required_cnv_objects, exists, envir = .GlobalEnv)]
  
  if (length(missing_cnv_objects) == 0) {
    if (verbose) message("STEP 3/8: CNV data already loaded (all required objects present in GlobalEnv)\n")
    load_status$cnv <- "ALREADY_LOADED"
  } else {
    if (verbose) {
      message("STEP 3/8: Loading Copy Number Variation (CytoSNP) Data...")
    }
    tryCatch({
      source(file.path(DIRS$scripts$helpers, "Data_Loading", "04_load_cnv_data.R"))
      load_cnv_data()
      load_status$cnv <- "SUCCESS"
      message("  Missing objects: ", paste(missing_cnv_objects, collapse = ", "))
      if (verbose) message("  ✓ CNV data loaded\n")
    }, error = function(e) {
      load_status$cnv <- paste("FAILED:", e$message)
      if (verbose) message("  ✗ Error loading CNV: ", e$message, "\n")
    })
  }
  
  # ------------------------------------------------------------------------------
  # STEP 4: Load Single Nucleotide Variation Data
  # ------------------------------------------------------------------------------
  if (exists("BRCA_CL_SNV", envir = .GlobalEnv) && exists("BRCA_CL_MAF", envir = .GlobalEnv)) {
    if (verbose) message("STEP 4/8: SNV data already loaded (skipping)\n")
    load_status$snv <- "ALREADY_LOADED"
  } else {
    if (verbose) message("STEP 4/8: Loading Single Nucleotide Variation (WES) Data...")
    tryCatch({
      source(file.path(DIRS$scripts$helpers, "Data_Loading", "05_load_snv_data.R"))
      load_snv_data()
      load_status$snv <- "SUCCESS"
      if (verbose) message("  ✓ SNV data loaded\n")
    }, error = function(e) {
      load_status$snv <- paste("FAILED:", e$message)
      if (verbose) message("  ✗ Error loading SNV: ", e$message, "\n")
    })
  }
  
  # ------------------------------------------------------------------------------
  # STEP 5: Load DNA Methylation Data
  # ------------------------------------------------------------------------------
  if (exists("BRCA_CL_DNAm", envir = .GlobalEnv)) {
    if (verbose) message("STEP 5/8: DNAm data already loaded (skipping)\n")
    load_status$dnam <- "ALREADY_LOADED"
  } else {
    if (verbose) message("STEP 5/8: Loading DNA Methylation (DNAm) Data...")
    tryCatch({
      source(file.path(DIRS$scripts$helpers, "Data_Loading", "06_load_dnam_data.R"))
      load_dnam_data()
      load_status$dnam <- "SUCCESS"
      if (verbose) message("  ✓ DNAm data loaded\n")
    }, error = function(e) {
      load_status$dnam <- paste("FAILED:", e$message)
      if (verbose) message("  ✗ Error loading DNAm: ", e$message, "\n")
    })
  }
  
  
  # ------------------------------------------------------------------------------
  # STEP 6: Load Structural Variation Data
  # ------------------------------------------------------------------------------
  if (exists("ICLE_SV", envir = .GlobalEnv)) {
    if (verbose) message("STEP 6/8: SV data already loaded (skipping)\n")
    load_status$sv <- "ALREADY_LOADED"
  } else {
    if (verbose) message("STEP 6/8: Loading Structural Variation (Bionano) Data...")
    tryCatch({
      source(file.path(DIRS$scripts$helpers, "Data_Loading", "07_load_sv_data.R"))
      sv_data <- load_sv_data()
      assign("ICLE_SV", sv_data$ICLE_SV, envir = .GlobalEnv)
      load_status$sv <- "SUCCESS"
      if (verbose) message("  ✓ SV data loaded\n")
    }, error = function(e) {
      load_status$sv <- paste("FAILED:", e$message)
      if (verbose) message("  ✗ Error loading SV: ", e$message, "\n")
    })
  }
  
  
  # ------------------------------------------------------------------------------
  # STEP 7: Load External Datasets (Optional)
  # ------------------------------------------------------------------------------
  if (load_external) {
    # Check if required external data objects already exist in GlobalEnv
    # TCGA_Annots and MSK_Annots are critical required objects
    required_tcga_objects <- c("TCGA_Annots")
    required_msk_objects <- c("MSK_Annots")
    
    missing_tcga <- required_tcga_objects[!sapply(required_tcga_objects, exists, envir = .GlobalEnv)]
    missing_msk <- required_msk_objects[!sapply(required_msk_objects, exists, envir = .GlobalEnv)]
    
    # If all required objects exist, skip loading
    if (length(missing_tcga) == 0 && length(missing_msk) == 0) {
      if (verbose) {
        message("STEP 7/8: External data already loaded (all required objects present in GlobalEnv)")
        message("  ✓ TCGA_Annots: present")
        message("  ✓ MSK_Annots: present\n")
      }
      load_status$external <- "ALREADY_LOADED"
    } else {
      if (verbose) {
        message("STEP 7/8: Loading External Datasets (TCGA, MSK)...")
        if (length(missing_tcga) > 0) {
          message("  Missing TCGA objects: ", paste(missing_tcga, collapse = ", "))
        }
        if (length(missing_msk) > 0) {
          message("  Missing MSK objects: ", paste(missing_msk, collapse = ", "))
        }
      }
      tryCatch({
        source(file.path(DIRS$scripts$helpers, "Data_Loading", "08_load_external_data.R"))
        load_external_data()
        load_status$external <- "SUCCESS"
        if (verbose) message("  ✓ External data loaded\n")
      }, error = function(e) {
        load_status$external <- paste("FAILED:", e$message)
        if (verbose) message("  ✗ Error loading external data: ", e$message, "\n")
      })
    }
  } else {
    load_status$external <- "SKIPPED"
    if (verbose) message("STEP 7/8: Skipping external datasets\n")
  }
  
  # ------------------------------------------------------------------------------
  # STEP 8: Generate Genomic Alteration Matrices
  # ------------------------------------------------------------------------------
  if (exists("BRCA_CL_GAM", envir = .GlobalEnv)) {
    if (verbose) message("STEP 8/8: GAM data already loaded (skipping)\n")
    load_status$gam <- "ALREADY_LOADED"
  } else {
    if (verbose) message("STEP 8/8: Generating Genomic Alteration Matrices...")
    tryCatch({
      gam_script <- file.path(DIRS$scripts$helpers, "Data_Loading", "09_generate_gams.R")
      if (file.exists(gam_script)) {
        source(gam_script, chdir = TRUE)
        load_status$gam <- "SUCCESS"
        
      } else {
        load_status$gam <- "SCRIPT NOT FOUND"
        if (verbose) message("  ℹ GAM generation script not found\n")
      }
    }, error = function(e) {
      load_status$gam <- paste("FAILED:", e$message)
      if (verbose) message("  ✗ Error generating GAM: ", e$message, "\n")
    })
  }
  
  # ------------------------------------------------------------------------------
  # Summary
  # ------------------------------------------------------------------------------
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  
  if (verbose) {
    message("")
    message("═══════════════════════════════════════════════════════")
    message("  Data Loading Complete")
    message("═══════════════════════════════════════════════════════")
    message("")
    message("Loading Status Summary:")
    for (dataset in names(load_status)) {
      status_symbol <- if (load_status[[dataset]] == "SUCCESS") "✓" else 
        if (load_status[[dataset]] == "SKIPPED") "⊘" else 
        if (load_status[[dataset]] == "ALREADY_LOADED") "⊚" else "✗"
      message(sprintf("  %s %-20s: %s", status_symbol, toupper(dataset), load_status[[dataset]]))
    }
    message("")
    message(sprintf("Total Time: %.2f minutes", as.numeric(elapsed_time)))
    message("\nKey Objects Available in Global Environment:")
    for (obj in c("CL_Annots", "BRCA_CL_RPPA", "BRCA_CL_EXP", "BRCA_CL_GISTIC", 
                  "BRCA_CL_SNV", "BRCA_CL_DNAm", "ICLE_SV", "BRCA_CL_GAM")) {
      if (exists(obj, envir = .GlobalEnv)) {
        obj_data <- get(obj, envir = .GlobalEnv)
        if (is.data.frame(obj_data) || is.matrix(obj_data)) {
          message(sprintf("  • %-20s: %s", obj, paste(dim(obj_data), collapse = " x ")))
        } else {
          message(sprintf("  • %-20s: loaded", obj))
        }
      }
    }
    message("═══════════════════════════════════════════════════════")
  }
  
  # Return nothing - all data is already in global environment
  return(invisible(NULL))
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# 
# # Load all data (objects placed in global environment)
# load_all_icle_data()
# 
# # Load only ICLE data (skip external)
# load_all_icle_data(load_external = FALSE)
# 
# # Silent loading
# load_all_icle_data(verbose = FALSE)
#
# ==============================================================================
# File Naming Convention
# ==============================================================================
# All data loaders follow numeric prefix naming:
#   00_load_annotations.R   - Cell line annotations (REQUIRED FIRST)
#   01_load_all_data.R      - This orchestrator
#   02_load_rppa_data.R     - Protein expression
#   03_load_rna_data.R      - RNA-seq expression
#   04_load_cnv_data.R      - Copy number variations
#   05_load_snv_data.R      - Single nucleotide variations
#   06_load_dnam_data.R     - DNA methylation
#   07_load_sv_data.R       - Structural variations
#   08_load_external_data.R - External datasets
#   09_generate_gams.R      - Genomic alteration matrices
