# ==============================================================================
# Script 08: Figure 2C - CDH1 Exonic Deletions Analysis
# ==============================================================================
# Description: Analyzes CDH1 exonic deletions from WES data and generates
#              a heatmap visualization showing deletion patterns across cell lines.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - CDH1 GTF annotation file (via config paths)
#   - WES pileup files for CDH1 reads (via config paths)
#   - CL_Annots: Cell line annotations
#
# Output:
#   - fig2c_cdh1_exonic_del_heatmap: CDH1 exonic deletion heatmap (assigned to .GlobalEnv)
#   - PDF saved via Main_Data_Analysis.Rmd
#
# Author: Osama Shiraz Shah
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ComplexHeatmap)
  library(stringr)
})

# %notin% operator is defined in Helper_Functions.R
if (!exists("%notin%", envir = .GlobalEnv)) {
  `%notin%` <- Negate(`%in%`)
}

# ==============================================================================
# Helper Functions
# ==============================================================================

# ' Load and Process CDH1 Annotation Data
# ' 
# ' Loads CDH1 GTF annotation file and extracts exon information for the
# ' canonical transcript NM_004360 (longest isoform).
# ' 
# ' @param annotation_file Path to CDH1 GTF bed file
# ' @return Data frame with CDH1 exon annotations
# ' @export
# load_cdh1_annotation <- function(annotation_file) {
  
#   message("Loading CDH1 annotation data...")
  
#   if (!file.exists(annotation_file)) {
#     stop("CDH1 annotation file not found: ", annotation_file)
#   }
  
#   # Read annotation file
#   CDH1_annotation <- read.delim(annotation_file, header = FALSE, 
#                                 stringsAsFactors = FALSE)
  
#   # Extract exon numbers
#   CDH1_annotation$exon <- as.character(
#     stringr::str_extract_all(
#       string = CDH1_annotation$V10, 
#       pattern = "exon_number [0-9]+", 
#       simplify = TRUE
#     )
#   )
#   CDH1_annotation$exon <- as.numeric(
#     gsub(x = CDH1_annotation$exon, 
#          pattern = "exon_number ", 
#          replacement = "")
#   )
  
#   # Extract transcript IDs
#   CDH1_annotation$transcript <- as.character(
#     stringr::str_extract_all(
#       string = CDH1_annotation$V10, 
#       pattern = "transcript_id NM_[0-9]+", 
#       simplify = TRUE
#     )
#   )
#   CDH1_annotation$transcript <- gsub(
#     x = CDH1_annotation$transcript, 
#     pattern = "transcript_id ", 
#     replacement = ""
#   )
  
#   # Select and rename relevant columns
#   CDH1_exon_annot <- CDH1_annotation[, c(1, 2, 3, 8, 11, 12)]
#   CDH1_exon_annot$exon[is.na(CDH1_exon_annot$exon)] <- 0
#   names(CDH1_exon_annot) <- c("chr", "start", "end", "feature", 
#                                "exon_number", "transcript")
  
#   # Filter for NM_004360 (canonical/longest isoform)
#   # Reference: https://www.ncbi.nlm.nih.gov/gene/999
#   CDH1_exon_annot <- subset(
#     CDH1_exon_annot, 
#     transcript == "NM_004360" & exon_number != 0
#   )# [, -6]  # Remove transcript column
  
#   # Add additional metadata
#   CDH1_exon_annot$width <- CDH1_exon_annot$end - CDH1_exon_annot$start
#   CDH1_exon_annot$symbol <- "CDH1"
  
#   # Filter for exon features only
#   CDH1_exon_annot <- subset(
#     CDH1_exon_annot, 
#     exon_number != 0 & grepl(pattern = "exon", x = feature)
#   )
  
#   message("  ✓ Loaded ", nrow(CDH1_exon_annot), " exons for CDH1 (NM_004360)")
  
#   return(CDH1_exon_annot)
# }


#' Load WES CDH1 Read Counts
#' 
#' Loads WES pileup files containing CDH1 read counts across samples.
#' 
#' @param pileup_dir Directory containing pileup files
#' @return Named list of read count vectors
#' @export
load_wes_cdh1_counts <- function(pileup_dir) {
  
  message("  Loading WES CDH1 read counts...")
  
  if (!dir.exists(pileup_dir)) {
    stop("Pileup directory not found: ", pileup_dir)
  }
  
  # Find all pileup files
  WES_IDs <- list.files(
    path = pileup_dir, 
    pattern = "*dupsmarked*",
    full.names = FALSE
  )
  
  if (length(WES_IDs) == 0) {
    stop("No pileup files found in: ", pileup_dir)
  }
  
  message("  Found ", length(WES_IDs), " WES samples")

  annots <- read.delim(file.path(pileup_dir, WES_IDs[1]), header = FALSE, stringsAsFactors = FALSE)[,c(1,2,3,4,8,10)]
  names(annots) <- c("chr", "start", "end", "feature", "location", "additional_info")
  annots$transcript <- gsub(x = annots$additional_info, pattern = ".*(transcript_id NM_[0-9]+).*", replacement = "\\1")
  annots$transcript <- gsub(x = annots$transcript, pattern = "transcript_id ", replacement = "")
  annots$exon_number <- gsub(x = annots$additional_info, pattern = ".*(exon_number [0-9]+).*", replacement = "\\1")
  annots$exon_number <- gsub(x = annots$exon_number, pattern = "exon_number ", replacement = "")
  annots$exon_number <- as.numeric(annots$exon_number)

  annots <- annots[,c("chr", "start", "end", "feature", "location", "transcript", "exon_number")]

  WES_CDH1_counts <- as.data.frame(matrix(nrow = nrow(annots), ncol = ncol(annots)+length(WES_IDs), data = NA))
  WES_CDH1_counts[,1:ncol(annots)] <- as.data.frame(annots)
  colnames(WES_CDH1_counts)[1:ncol(annots)] <- colnames(annots)  

  simple_names <- gsub(x = WES_IDs, pattern = "[-]dupsmarked[.]bam_WES[.]tsv|CDH1_", replacement = "")
  simple_names <- gsub(x = simple_names, pattern = "MPE600", replacement = "600MPE")
  simple_names <- gsub(x = simple_names, pattern = "MM", replacement = "MDAMB")
  simple_names <- gsub(x = simple_names, pattern = "MDAMB134", replacement = "MDAMB134VI")
  colnames(WES_CDH1_counts)[ncol(annots)+1:length(WES_IDs)] <- simple_names

  # Load read counts from each file
  for (i in 1:length(WES_IDs)) {
    file_path <- file.path(pileup_dir, WES_IDs[i])
    tmp <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)[,c(11)]
    WES_CDH1_counts[,simple_names[i]] <- as.matrix(tmp)
  }
  message("  ✓ Loaded read counts for ", length(simple_names), " samples")
  
  return(WES_CDH1_counts)
}


#' Generate CDH1 Exonic Deletion Heatmap
#' 
#' Creates a heatmap showing CDH1 exonic deletion patterns across samples.
#' 
#' @param CDH1_cds_counts Combined CDH1 count data
#' @param samples_to_plot Character vector of sample names to include
#' @return ComplexHeatmap object
#' @export
generate_cdh1_exon_deletion_heatmap <- function(WES_CDH1_counts, 
                                           samples_to_plot = c("OCUBM", "MDAMB134VI", 
                                                              "HCC2218", "HCC2185", "SKBR3")) {
  
  message("  Generating CDH1 deletion heatmap...")
  
  # Filter for samples of interest
  missing_samples <- setdiff(samples_to_plot, colnames(WES_CDH1_counts))
  if (length(missing_samples) > 0) {
    warning("Missing samples: ", paste(missing_samples, collapse = ", "))
    samples_to_plot <- intersect(samples_to_plot, colnames(WES_CDH1_counts))
  }

  # Extract count matrix
  mat <- WES_CDH1_counts[, samples_to_plot, drop = FALSE] %>% as.matrix()
  
  # Convert to deletion matrix (0 = reads present, 1 = deletion)
  mat[mat > 0] <- 1
  mat <- 1 - mat
  
  # Create row split labels (genomic position - exon number)
  row_labels <- paste0("Exon ", WES_CDH1_counts$exon_number)
  explicit_order <- order(WES_CDH1_counts$exon_number)
  row_labels <- factor(row_labels, levels = unique(row_labels[explicit_order]))
  
  # Generate heatmap
  set.seed(123)
  cdh1_exon_del_ht <- Heatmap(
    matrix = mat, 
    show_row_names = FALSE, 
    row_names_gp = gpar(fontsize = 8), 
    name = "Deletion",
    col = c("0" = "gray", "1" = "#0d47a1"), 
    column_split = colnames(mat),
    width = unit(5, "cm"), 
    height = unit(8, "cm"),
    row_split = row_labels, 
    column_title = " ", 
    row_title_rot = 0,
    column_title_rot = 90, 
    show_column_dend = FALSE, 
    show_column_names = TRUE, 
    row_title_gp = gpar(fontsize = 10),
    cluster_row_slices = FALSE, 
    cluster_rows = FALSE, 
  )
  
  message("  ✓ Heatmap generated with ", nrow(mat), " exons and ", 
          ncol(mat), " samples")
  
  return(cdh1_exon_del_ht)
}


# ==============================================================================
# Main Analysis Function
# ==============================================================================

#' Run Complete CDH1 Exonic Deletion Analysis
#' 
#' Performs end-to-end analysis of CDH1 exonic deletions from WES data.
#' 
#' @param config_loaded Logical, whether config.R has been loaded
#' @return List containing annotation data, counts, and heatmap
#' @export
run_cdh1_exonic_deletion_analysis <- function(config_loaded = TRUE) {
  
  message("\n========================================")
  message("Figure 2C: CDH1 Exonic Deletion Analysis")
  message("========================================")
  
  if (!config_loaded || !exists("DIRS", envir = .GlobalEnv)) {
    stop("config.R must be loaded before running this analysis. Run: source('config.R')")
  }
  if (!exists("FILES", envir = .GlobalEnv)) stop("FILES not found; config.R may be incomplete.")
  annotation_file <- FILES$cdh1_gtf
  pileup_dir <- FILES$cdh1_pileups_dir
  
  # Step 1: Load CDH1 annotation & WES read counts
  # CDH1_exon_annot <- load_cdh1_annotation(annotation_file)
  WES_CDH1_counts <- load_wes_cdh1_counts(pileup_dir)
  WES_CDH1_counts = subset(WES_CDH1_counts, location == "exon" & transcript == "NM_004360")
  
  # Step 4: Generate heatmap
  cdh1_exon_del_ht <- generate_cdh1_exon_deletion_heatmap(
    WES_CDH1_counts,
    samples_to_plot = c("OCUBM", "MDAMB134VI", "HCC2218", "HCC2185", "SKBR3")
  )
  
  assign("fig2c_cdh1_exonic_del_heatmap", cdh1_exon_del_ht, envir = .GlobalEnv)
  message("  ✓ CDH1 exonic deletion analysis complete (fig2c_cdh1_exonic_del_heatmap assigned).\n")
  
  return(invisible(cdh1_exon_del_ht))
}

fig2c_cdh1_exonic_del_heatmap <-run_cdh1_exonic_deletion_analysis()


# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/08_Fig2C_CDH1_exonic_deletions.R")
#
# fig2c_cdh1_exonic_del_heatmap <- run_cdh1_exonic_deletion_analysis()#
# fig2c_cdh1_exonic_del_heatmap assigned to .GlobalEnv; save PDF inline in Rmd.





