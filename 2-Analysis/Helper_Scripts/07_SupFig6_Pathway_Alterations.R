# ==============================================================================
# Script 07: Supplementary Figure 6 - Alterations in Key Pathways (ICLE)
# ==============================================================================
# Description: Generates pathway-level alteration matrix and heatmap for ICLE
#              cell lines. Aggregates gene-level alterations to pathway-level
#              summaries for visualization.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - BRCA_CL_GAM: Cell line genomic alteration matrix
#   - CL_Annots: Cell line annotations
#
# Output:
#   - SupFigS6: Pathway alteration heatmap (assigned to .GlobalEnv)
#   - pathway_alt_mat: Pathway alteration matrix (assigned to .GlobalEnv)
#
# Author: Osama Shiraz Shah
# ==============================================================================

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}

suppressPackageStartupMessages({ library(dplyr); library(ComplexHeatmap) })

message("\n========================================")
message("SupFig 6: Pathway alterations (ICLE)")
message("========================================")

# Check if required data objects exist, load if needed
if (!exists("BRCA_CL_GAM", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = TRUE)
}

# Load pathway databases if not already loaded
if (!exists("ddr_genes", envir = .GlobalEnv)) {
  ddr_genes <- read.delim(FILES$ddr_genes)
} else {
  ddr_genes <- get("ddr_genes", envir = .GlobalEnv)
}

if (!exists("pathway_db", envir = .GlobalEnv)) {
  pathway_db <- read.delim(FILES$maftools_pathways)
} else {
  pathway_db <- get("pathway_db", envir = .GlobalEnv)
}

if (!exists("oncoVar", envir = .GlobalEnv)) {
  if (exists("FILES", envir = .GlobalEnv) && "oncovar" %in% names(FILES)) {
    oncoVar <- read.delim(FILES$oncovar)
    rownames(oncoVar) <- oncoVar$Gene_symbol
  } else {
    stop("oncoVar not found and FILES$oncovar not available")
  }
} else {
  oncoVar <- get("oncoVar", envir = .GlobalEnv)
}

brca_drivers <- subset(oncoVar, Driver.Level %in% c(3, 4))$Gene_symbol

pathway_db <- rbind(pathway_db, subset(ddr_genes, Gene %notin% pathway_db$Gene)) %>%
  filter(Gene %in% brca_drivers)
rownames(pathway_db) <- pathway_db$Gene

pathway_alt_mat <- BRCA_CL_GAM[intersect(pathway_db$Gene, rownames(BRCA_CL_GAM)), ]
pathway_alt_mat[pathway_alt_mat == "LOH"] <- ""
pathway_alt_mat[pathway_alt_mat == "GAIN"] <- ""
pathway_alt_mat[pathway_alt_mat == ""] <- "WT"
pathway_alt_mat[pathway_alt_mat == "AMP"] <- "AMP"
pathway_alt_mat[pathway_alt_mat == "DEL"] <- "DEL"
pathway_alt_mat[pathway_alt_mat == "MUT"] <- "MUT"
pathway_alt_mat[pathway_alt_mat == "MUT;LOH"] <- "MUT"
pathway_alt_mat[pathway_alt_mat == "MUT;GAIN"] <- "MUT"
pathway_alt_mat[pathway_alt_mat == "MUT;DEL"] <- "MUT"
pathway_alt_mat[pathway_alt_mat == "MUT;AMP"] <- "MUT"

ICLE_cells <- grep(x = colnames(pathway_alt_mat), pattern = "-I", value = TRUE)
pathway_alt_mat <- pathway_alt_mat[, ICLE_cells]
pathway_alt_mat <- pathway_alt_mat[rowSums(pathway_alt_mat != "WT") > 0, ]

pathway_scores <- list()
for (p in unique(pathway_db[rownames(pathway_alt_mat), ]$Pathway)) {
  pathway_scores[[p]] <- sum(pathway_alt_mat[subset(pathway_db, Pathway == p & Gene %in% rownames(pathway_alt_mat))$Gene, ] != "WT")
}

row_split <- pathway_db[rownames(pathway_alt_mat), ]$Pathway
pathway_alt_freq <- 100 * sapply(split(pathway_db[rownames(pathway_alt_mat), "Gene"], pathway_db[rownames(pathway_alt_mat), "Pathway"]), length) /
  sapply(split(subset(pathway_db, Pathway %in% unique(row_split))$Gene, subset(pathway_db, Pathway %in% unique(row_split))$Pathway), length)
row_split <- factor(row_split, levels = names(pathway_alt_freq)[order(pathway_alt_freq, decreasing = TRUE)])

set.seed(123)
SupFigS6 <- Heatmap(
  pathway_alt_mat,
  col = annot_cols$Alt_Simple,
  name = "Alteration",
  width = unit(8, "cm"),
  height = unit(15, "cm"),
  row_names_gp = gpar(fontface = "bold", fontsize = 10, col = "black"),
  row_title_rot = 0,
  column_names_gp = gpar(fontsize = 14, col = "black"),
  show_row_dend = FALSE,
  column_labels = gsub("-I", "", ICLE_cells),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_column_dend = TRUE,
  row_order = order(rowSums(pathway_alt_mat != "WT"), decreasing = TRUE),
  row_split = row_split,
  column_split = CL_Annots[ICLE_cells, "mRNA Subtypes"],
  column_title = " ",
  row_title = gsub("_", " ", unique(paste0(levels(row_split), " (", round(pathway_alt_freq[levels(row_split)], 0), "%)"))),
  top_annotation = HeatmapAnnotation(
    Subtype = CL_Annots[ICLE_cells, "mRNA Subtypes"],
    col = list(Subtype = annot_cols$Subtypes),
    simple_anno_size = unit(4, "mm"),
    annotation_legend_param = heatmap_legend_param
  ),
  column_title_rot = 90,
  heatmap_legend_param = heatmap_legend_param,
  right_annotation = HeatmapAnnotation(
    Freq = anno_barplot(rowSums(pathway_alt_mat != "WT"), ylim = c(0, 15), gp = gpar(fill = "darkgray", col = "black"), bar_width = 0.7),
    annotation_legend_param = heatmap_legend_param,
    which = "row",
    simple_anno_size = unit(8, "mm")
  ),
  border = TRUE,
  border_gp = gpar(col = "black"),
  rect_gp = gpar(col = "gray", lwd = 1)
)

assign("SupFigS6", SupFigS6, envir = .GlobalEnv)
assign("pathway_alt_mat", pathway_alt_mat, envir = .GlobalEnv)
message("  ✓ SupFig 6 complete (SupFigS6, pathway_alt_mat assigned).\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# load_all_icle_data(load_external = TRUE)
# source("Helper_Scripts/07_SupFig6_Pathway_Alterations.R")
#
# pdf(file.path(DIRS$results, "SupFig6_Pathway_Alterations.pdf"), width = 10, height = 15)
# draw(SupFigS6, merge_legends = TRUE)
# dev.off()
# ensure_dir(DIRS$results_sub$molecular_resemblance)
# ensure_dir(DIRS$results_sub$molecular_resemblance)
# write.table(pathway_alt_mat, 
#             file.path(DIRS$results_sub$molecular_resemblance, "SupTable_ICL_pathway_alterations.tsv"), 
#             sep = "\t", quote = FALSE, col.names = NA)
