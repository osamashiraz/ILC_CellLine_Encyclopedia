# ==============================================================================
# Script 00: Load Cell Line Annotations and Histology Data
# ==============================================================================
# This script loads and prepares cell line annotation data including:
#   - Histology classifications
#   - Molecular subtyping information
#   - Sample metadata
#
# Author: Osama Shiraz Shah
# ==============================================================================

#' Load Cell Line Annotations
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing CL_Annots_simple and CL_Annots data frames
#' @export
load_cell_line_annotations <- function(config_loaded = TRUE) {
  
  # Check if already loaded
  if (exists("CL_Annots", envir = .GlobalEnv) && exists("CL_Annots_simple", envir = .GlobalEnv)) {
    message("Cell line annotations already loaded. Returning existing objects.")
    return(list(
      CL_Annots_simple = get("CL_Annots_simple", envir = .GlobalEnv),
      CL_Annots = get("CL_Annots", envir = .GlobalEnv),
      metadata = list(
        n_samples = nrow(get("CL_Annots", envir = .GlobalEnv)),
        load_time = Sys.time(),
        already_loaded = TRUE
      )
    ))
  }
  
  # Verify config is loaded
  if (!config_loaded || !exists("FILES", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this function. Run: source('config.R')")
  }
  
  message("Loading cell line annotations...")
  
  # ------------------------------------------------------------------------------
  # 1. Load Histology Data
  # ------------------------------------------------------------------------------
  if (!file.exists(FILES$cl_annots_simple)) {
    stop("Cell line annotations simple file not found: ", FILES$cl_annots_simple)
  }
  
  CL_Annots_simple <- read.delim(FILES$cl_annots_simple, stringsAsFactors = FALSE)
  rownames(CL_Annots_simple) <- CL_Annots_simple$Name
  
  # Validate required columns
  required_cols <- c("Name", "Histology")
  missing_cols <- setdiff(required_cols, colnames(CL_Annots_simple))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in annotations simple data: ", paste(missing_cols, collapse = ", "))
  }
  
  # message("  ✓ Loaded annotations simple data: ", nrow(CL_Annots_simple), " cell lines")
  
  # ------------------------------------------------------------------------------
  # 2. Load Molecular Subtyping Annotations
  # ------------------------------------------------------------------------------
  if (!file.exists(FILES$cl_annotations)) {
    stop("Cell line annotations file not found: ", FILES$cl_annotations)
  }
  
  CL_Annots <- read.delim(FILES$cl_annotations, stringsAsFactors = FALSE)
  rownames(CL_Annots) <- CL_Annots$Name
  
  # Clean column names (replace dots with spaces)
  colnames(CL_Annots) <- gsub("[.]", " ", colnames(CL_Annots))
  
  message("  ✓ Loaded annotation data: ", nrow(CL_Annots), " cell lines")
  
  # ------------------------------------------------------------------------------
  # 3. Set Factor Levels
  # ------------------------------------------------------------------------------
  # This ensures consistent ordering in visualizations
  
  if ("mRNA Subtypes" %in% colnames(CL_Annots)) {
    CL_Annots$`mRNA Subtypes` <- factor(
      CL_Annots$`mRNA Subtypes`, 
      levels = SAMPLE_GROUPS$subtype_order
    )
  }
  
  if ("RPPA Subtypes" %in% colnames(CL_Annots)) {
    CL_Annots$`RPPA Subtypes` <- factor(
      CL_Annots$`RPPA Subtypes`, 
      levels = SAMPLE_GROUPS$subtype_order
    )
  }
  
  if ("DNAm Subtypes" %in% colnames(CL_Annots)) {
    CL_Annots$`DNAm Subtypes` <- factor(
      CL_Annots$`DNAm Subtypes`, 
      levels = SAMPLE_GROUPS$subtype_order
    )
  }
  
  if ("Histology" %in% colnames(CL_Annots)) {
    CL_Annots$Histology <- factor(
      CL_Annots$Histology, 
      levels = SAMPLE_GROUPS$histology_order
    )
  }
  
  # message("  ✓ Set factor levels for molecular subtypes and histology")
  
  # ------------------------------------------------------------------------------
  # 4. Data Validation
  # ------------------------------------------------------------------------------
  # Check for consistency between datasets
  common_samples <- intersect(CL_Annots_simple$Name, CL_Annots$Sample)
  
  # if (length(common_samples) == 0) {
  #   warning("No common samples found between histology and annotation data!")
  # } else {
  #   message("  ✓ ", length(common_samples), " samples found in both datasets")
  # }
  
  # Check for NA values in critical columns
  # critical_cols <- c("Name", "Histology", "Study")
  # for (col in critical_cols) {
  #   if (col %in% colnames(CL_Annots)) {
  #     na_count <- sum(is.na(CL_Annots[[col]]))
  #     if (na_count > 0) {
  #       warning("Found ", na_count, " NA values in column: ", col)
  #     }
  #   }
  # }
  
  # ------------------------------------------------------------------------------
  # 5. Summary Statistics
  # ------------------------------------------------------------------------------
  if ("Histology" %in% colnames(CL_Annots)) {
    message("\nHistology distribution:")
    print(table(CL_Annots$Histology, useNA = "ifany"))
  }
  
  if ("Study" %in% colnames(CL_Annots)) {
    message("\nStudy distribution:")
    print(table(CL_Annots$Study, useNA = "ifany"))
  }
  
  # ------------------------------------------------------------------------------
  # 6. Return Results
  # ------------------------------------------------------------------------------
  results <- list(
    CL_Annots_simple = CL_Annots_simple,
    CL_Annots = CL_Annots,
    metadata = list(
      n_samples = nrow(CL_Annots),
      n_common = length(common_samples),
      load_time = Sys.time()
    )
  )
  
  # message("\n✓ Annotation loading complete\n")
  
  return(results)
}

#' Get ICLE-specific cell lines
#' 
#' @param CL_Annots Cell line annotations data frame
#' @param exclude_overlap Exclude cell lines overlapping with CCLE
#' @return Character vector of ICLE cell line names
#' @export
get_icle_samples <- function(CL_Annots, exclude_overlap = TRUE) {
  
  # Filter for ICLE study
  icle_filter <- CL_Annots$Study == "ICLE"
  
  # Optionally exclude CCLE overlaps
  if (exclude_overlap && "overlapWCCLE" %in% colnames(CL_Annots)) {
    icle_filter <- icle_filter & (CL_Annots$overlapWCCLE != "Y")
  }
  
  # Further filter by histology if needed
  if ("Histology" %in% colnames(CL_Annots)) {
    icle_samples <- CL_Annots[icle_filter & 
                              CL_Annots$Histology %in% c("ILC", "ILC-like"), 
                              "Name"]
  } else {
    icle_samples <- CL_Annots[icle_filter, "Name"]
  }
  
  return(icle_samples)
}

#' Get NST cell lines
#' 
#' @param CL_Annots Cell line annotations data frame
#' @param exclude_overlap Exclude cell lines overlapping with CCLE
#' @return Character vector of NST cell line names
#' @export
get_nst_samples <- function(CL_Annots, exclude_overlap = TRUE) {
  
  nst_filter <- CL_Annots$Histology == "NST"
  
  if (exclude_overlap && "overlapWCCLE" %in% colnames(CL_Annots)) {
    nst_filter <- nst_filter & (CL_Annots$overlapWCCLE != "Y")
  }
  
  nst_samples <- CL_Annots[nst_filter, "Name"]
  
  return(nst_samples)
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/00_load_annotations.R")
# 
# annot_data <- load_cell_line_annotations()
# CL_Annots_simple <- annot_data$CL_Annots_simple
# CL_Annots <- annot_data$CL_Annots
# 
# icle_samples <- get_icle_samples(CL_Annots)
# nst_samples <- get_nst_samples(CL_Annots)
