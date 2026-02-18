# ==============================================================================
# Script 17: Figure 3C - Chromothripsis Identification (ShatterSeek)
# ==============================================================================
# Description: Runs ShatterSeek on ICLE CNV + SV data to call chromothripsis events.
#              Identifies complex chromosomal rearrangements characteristic of
#              chromothripsis across ICLE cell lines.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - ICLE_SV: Structural variation data
#   - CNV data (via config paths)
#   - CL_Annots: Cell line annotations
#
# Output:
#   - ShatterSeek analysis results
#   - Chromothripsis event calls
#
# Author: Osama Shiraz Shah
# ==============================================================================

# %notin% operator is defined in Helper_Functions.R
if (!exists("%notin%", envir = .GlobalEnv)) {
  `%notin%` <- Negate(`%in%`)
}

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(reshape2)
  library(GenomicRanges)
  library(ShatterSeek)
  library(gridExtra)
  library(ComplexHeatmap)
  library(grid)
})



# -----------------------
# 2. Load and Preprocess CNV Data
# -----------------------
load_cnv_data <- function(cnv_path, sample_patterns = c("-I|SKBR3-M")) {
  cnv_df <- vroom::vroom(cnv_path, show_col_types = FALSE) %>% as.data.frame()
  # Standardize sample names by removing '-I' or '-M' suffixes
  cnv_df$Sample <- gsub("-I|-M", "", cnv_df$ID)
  # Select samples matching '-I' in ID and 'SKBR3-M'
  selected_samples <- c(grep(sample_patterns, unique(cnv_df$ID), value = TRUE))
  cnv_df <- subset(cnv_df, ID %in% selected_samples)
  # Convert chromosome 23 to 'X'
  cnv_df$chrom <- gsub("23", "X", cnv_df$chrom)
  
  
  cnv_df$total_cn <- round(2^cnv_df$seg.mean * 2)
  cnv_df = cnv_df[, c("Sample", "chrom", "loc.start", "loc.end", "seg.mean", "total_cn")]
  colnames(cnv_df) = c("Sample", "chromosome", "start", "end", "seg.mean", "total_cn")
  
  return(cnv_df)
}

# -----------------------
# 3. Merge Adjacent CN Segments
# -----------------------
merge_adjacent_cn_segments <- function(cn_df) {
  cn_df$total_cn[cn_df$total_cn == 0] <- 150000
  cn_df$total_cn[is.na(cn_df$total_cn)] <- 0
  
  gr <- as(cn_df, "GRanges")
  cov <- coverage(gr, weight = gr$total_cn)
  merged_gr <- as(cov, "GRanges")
  merged_df <- merged_gr %>% as.data.frame() %>% filter(score != 0)
  merged_df <- merged_df[, c(1, 2, 3, 6)]
  names(merged_df) <- names(cn_df)[1:4]
  merged_df$total_cn[merged_df$total_cn == 150000] <- 0
  
  return(merged_df)
}


# -----------------------
# 4. Load and Preprocess SV Data
# -----------------------
load_sv_data <- function(sv_path) {
  sv_df <- read.csv(sv_path)
  sv_df$RefStartPos <- as.integer(sv_df$RefStartPos)
  sv_df$RefEndPos   <- as.integer(sv_df$RefEndPos)
  return(sv_df)
}


# -----------------------
# 5. Prepare ShatterSeek Inputs
# -----------------------
prepare_shatterseek_inputs <- function(sample_name, cn_df, sv_df) {
  cn_input <- CNVsegs(
    chrom    = as.character(cn_df$chromosome),
    start    = cn_df$start,
    end      = cn_df$end,
    total_cn = cn_df$total_cn
  )
  sample_sv <- subset(sv_df, Sample == sample_name)
  sample_sv$Type <- dplyr::recode(sample_sv$Type,
                                  "Del" = "DEL",
                                  "Tranloc-inter" = "TRA",
                                  "Tranloc-intra" = "TRA",
                                  "Dup" = "DUP",
                                  "INS" = "INS")
  sample_sv$RefcontigID1 <- gsub("chr", "", sample_sv$RefcontigID1)
  sample_sv$RefcontigID2 <- gsub("chr", "", sample_sv$RefcontigID2)
  sample_sv <- subset(sample_sv, RefcontigID2 != "-1")
  sv_input <- SVs(
    chrom1  = as.character(sample_sv$RefcontigID1),
    pos1    = as.numeric(sample_sv$RefStartPos),
    chrom2  = as.character(sample_sv$RefcontigID2),
    pos2    = as.numeric(sample_sv$RefEndPos),
    SVtype  = as.character(sample_sv$Type),
    strand1 = rep("+", nrow(sample_sv)),
    strand2 = rep("+", nrow(sample_sv))
  )
  return(list(SV_data = sv_input, CN_data = cn_input))
}

# -----------------------
# 6. Run ShatterSeek on a Single Sample
# -----------------------
run_shatterseek_sample <- function(sample_name, cnv_df, sv_df) {
  sample_cnv <- subset(cnv_df, Sample == sample_name)[, c("chromosome", "start", "end", "total_cn")]
  merged_cn <- merge_adjacent_cn_segments(sample_cnv)
  inputs <- prepare_shatterseek_inputs(sample_name, merged_cn, sv_df)
  ss_result <- shatterseek(SV.sample = inputs$SV_data,
                           seg.sample = inputs$CN_data,
                           genome = "hg19")
  summary_df <- ss_result@chromSummary[, c(
    "chrom", "start", "end", "number_DEL", "number_DUP", "number_TRA",
    "number_h2hINV", "number_t2tINV", "clusterSize_including_TRA",
    "max_number_oscillating_CN_segments_2_states",
    "chr_breakpoint_enrichment", "pval_exp_cluster", "pval_fragment_joins",
    "inter_pval_fragment_joins", "inter_other_chroms_coords_all"
  )]
  summary_df$total_intra_SVs <- rowSums(ss_result@chromSummary[, c("number_DEL", "number_DUP", "number_h2hINV", "number_t2tINV")])
  summary_df$confidence <- ifelse(
    summary_df$clusterSize_including_TRA >= 12 &
      summary_df$max_number_oscillating_CN_segments_2_states >= 7 &
      summary_df$pval_fragment_joins <= 0.001 &
      (summary_df$chr_breakpoint_enrichment <= 0.001 | summary_df$pval_exp_cluster <= 0.001),
    "High",
    ifelse(
      summary_df$total_intra_SVs >= 8 &
        summary_df$max_number_oscillating_CN_segments_2_states >= 4 &
        summary_df$pval_fragment_joins <= 0.05 &
        (summary_df$chr_breakpoint_enrichment <= 0.05 | summary_df$pval_exp_cluster <= 0.05),
      "Low", ""
    )
  )
  summary_df$Sample <- sample_name
  return(summary_df)
}

# -----------------------
# 7. Run ShatterSeek for All Samples
# -----------------------
run_shatterseek_all <- function(cnv_df, sv_df) {
  all_samples <- unique(sv_df$Sample)
  result_list <- lapply(all_samples, function(sample) {
    message("Processing: ", sample)
    run_shatterseek_sample(sample, cnv_df, sv_df)
  })
  names(result_list) <- all_samples
  return(do.call(rbind, result_list))
}



# -----------------------
# 8. Save Results
# -----------------------
save_results <- function(summary_df, confidence_matrix, out_dir) {
  save(summary_df, file = file.path(out_dir, "ShatterSeek_df.Rdata"))
  # write.csv(summary_df,
  #             file = file.path(out_dir, "ShatterSeek_df.csv"),
  #             quote = FALSE,row.names = FALSE)
  confidence_matrix = as.data.frame(confidence_matrix)
  confidence_matrix$chrom = rownames(confidence_matrix)
  write.table(confidence_matrix,
              file = file.path(out_dir, "ShatterSeek_matrix.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}




# =======================
# 9. Chromothripsis Visualization
# =======================

visualize_chromothripsis_heatmap <- function(thripsisMat, annot_cols = NULL, heatmap_legend_param) {
  if (is.null(annot_cols)) annot_cols <- get("annot_cols", envir = .GlobalEnv)
  
  # Prepare sample names to match annotation format
  sample_names <- paste0(colnames(thripsisMat), "-I")
  
  set.seed(123); thripsis_ht <- 
    Heatmap(thripsisMat, border = TRUE, col = c("High" = "#eded2a", "Low" = "#a1d99b", "ns" = "white"),
            column_order = order(colSums(thripsisMat != "ns"), decreasing = TRUE), 
            width = unit(6, "cm"), height = unit(8, "cm"), 
            na_col = "gray", heatmap_legend_param = heatmap_legend_param,
            column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv), 
            column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv), 
            row_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv), 
            row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
            name = "Event Confidence", row_title = "Chromosomes", row_names_side = "left",
            row_order = order(rowSums(thripsisMat != "ns"), decreasing = TRUE),
            top_annotation = HeatmapAnnotation(Subtypes = CL_Annots[sample_names, "mRNA Subtypes"],
                                               # Histology = CL_Annots[sample_names, "Histology"],
                                               annotation_legend_param = heatmap_legend_param,
                                               simple_anno_size = unit(4, "mm"),
                                               col = annot_cols, # border = TRUE, gp = gpar(col = NA, lwd = 0.6),
                                               annotation_name_gp = gpar(fontsize = 12, col = "black", border = T, 
                                                                         fontfamily = helv), gap = unit(0.5, "mm")
            ), border_gp = gpar(col = "black"),rect_gp = gpar(col = "gray", lwd = 0.5) )
  
  return(thripsis_ht)
}




# =======================
# Main Function to Run Shatter Seek
# =======================

run_shatterseek_main <- function(cnv_path, sv_df, output_dir) {
  
  # Load and preprocess data
  cnv_df <- load_cnv_data(cnv_path)
  
  # Run chromothripsis analysis
  print("Running Shatter Seek on All ICLE Samples to identify Chromothripsis Events")
  chromothripsis_df <- run_shatterseek_all(cnv_df, sv_df)
  
  confidence_matrix <- reshape2::dcast(chromothripsis_df, chrom ~ Sample, value.var = "confidence")
  rownames(confidence_matrix) <- confidence_matrix$chrom
  confidence_matrix <- confidence_matrix[, -1] %>% as.matrix()
  confidence_matrix[is.na(confidence_matrix) | confidence_matrix == ""] <- "ns"
  
  # Save results
  print(paste0("Saving Chromothripsis Analysis Results to: ", output_dir))
  save_results(chromothripsis_df, confidence_matrix, output_dir)
  
  return(list(chromothripsis_df = chromothripsis_df, confidence_matrix = confidence_matrix, cnv_df = cnv_df, sv_df = sv_df))
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/13_Fig3_Identify_Chromothripsis.R")
#
# load_all_icle_data(load_external = TRUE)
#
# shatterseek_outs <- run_shatterseek_main(
#   cnv_path = FILES$cnv_seg,
#   sv_df = ICLE_SV,
#   output_dir = DIRS$results_sub$ogm
# )
#
# Outputs: chromothripsis TSV/PDF in output_dir; used by Fig3C and 14_Fig3D.


message("\n========================================")
message("Figure 3C: Chromothripsis landscape")
message("========================================")

# Check if required data objects exist, load if needed
if (!exists("ICLE_SV", envir = .GlobalEnv)) {
  message("  Loading ICLE_SV data...")
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  suppressMessages(load_all_icle_data(load_external = TRUE))
}

if (!exists("run_shatterseek_main", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "17_Fig3C_Identify_Chromothripsis.R"), chdir = TRUE)
}

if (!exists("shatterseek_outs", envir = .GlobalEnv)) {
  cnv_path <- FILES$cnv_seg
  output_dir <- DIRS$results_sub$ogm
  shatterseek_outs <- run_shatterseek_main(cnv_path, ICLE_SV, output_dir)
}

if (!exists("visualize_chromothripsis_heatmap", envir = .GlobalEnv)) {
  message("  Warning: visualize_chromothripsis_heatmap not found. Fig 3C heatmap may fail.")
  fig3c_thripsis_ht <- NULL
} else {
  fig3c_thripsis_ht <- visualize_chromothripsis_heatmap(
    shatterseek_outs$confidence_matrix,
    annot_cols,
    heatmap_legend_param
  )
}

assign("shatterseek_outs", shatterseek_outs, envir = .GlobalEnv)
if (!is.null(fig3c_thripsis_ht)) assign("fig3c_thripsis_ht", fig3c_thripsis_ht, envir = .GlobalEnv)
message("  ✓ Fig 3C complete (shatterseek_outs, thripsis_ht assigned).\n")

