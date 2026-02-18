# ==============================================================================
# Script 16: Figure 3B - Chromosomal Topography of Translocation Breakpoints
# ==============================================================================
# Description: Generates translocation breakpoint heatmap showing chromosomal
#              distribution of structural variation breakpoints across ICLE
#              cell lines.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - ICLE_SV: Structural variation data
#   - CL_Annots: Cell line annotations
#   - alt_count_chr: Alteration counts per chromosome (from 14_SupFig8_TMB_SV_Preparation.R)
#
# Outputs (assigned to .GlobalEnv):
#   - fig3b_transloc_breakpoints_ht: Translocation breakpoint heatmap (Fig 3B)
#
# Author: Osama Shiraz Shah
# ==============================================================================

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

message("\n========================================")
message("Figure 3B: Translocation breakpoint topography")
message("========================================")

# Check if required data objects exist, load if needed
if (!exists("ICLE_SV", envir = .GlobalEnv)) {
  message("  Loading ICLE_SV data...")
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = TRUE)
}

# Filter translocation SVs
sv_df <- subset(ICLE_SV, Type %in% c("Tranloc-inter", "Tranloc-intra"))[, 
            c("RefcontigID1", "RefcontigID2")]
colnames(sv_df) <- c("chrA", "chrB")

# Define chromosome order
chroms <- paste0("chr", c(1:22, "X"))
sv_df$chrA <- factor(sv_df$chrA, levels = chroms)
sv_df$chrB <- factor(sv_df$chrB, levels = chroms)

# Create count matrix
mat <- as.matrix(table(sv_df$chrA, sv_df$chrB))

# Filter chromosomes with at least 1 translocation
if (exists("alt_count_chr", envir = .GlobalEnv)) {
  alt_count_chr <- get("alt_count_chr", envir = .GlobalEnv)
  selected_chr <- alt_count_chr$Chr[alt_count_chr$SV_transloc > 0]
} else {
  selected_chr <- chroms[rowSums(mat) > 0 | colSums(mat) > 0]
}
mat <- mat[intersect(chroms, selected_chr), intersect(chroms, selected_chr)]

# Color function
col_fun <- colorRamp2(c(0, max(mat)), c("white", "#F57F17"))

# Row annotation - barplot of translocation counts
row_ha <- rowAnnotation(
  "Count" = anno_barplot(rowSums(mat), 
                        border = FALSE,
                        gp = gpar(fill = "gray", col = "black"),
                        axis_param = list(side = "bottom")),
  width = unit(1, "cm")
)

# Clean chromosome labels
rownames(mat) <- gsub("chr", "", rownames(mat))
colnames(mat) <- gsub("chr", "", colnames(mat))

# Create heatmap
fig3b_transloc_breakpoints_ht <- Heatmap(
  mat,
  name = "Count",
  col = col_fun,
  heatmap_legend_param = heatmap_legend_param,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "gray", lwd = 0.5),
  right_annotation = row_ha,
  row_names_side = "left",
  column_names_side = "bottom",
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (mat[i, j] > 5) {
      grid.text(mat[i, j], x, y, gp = gpar(fontsize = 15))
    }
  }
)

assign("fig3b_transloc_breakpoints_ht", fig3b_transloc_breakpoints_ht, envir = .GlobalEnv)
message("  ✓ Fig 3B complete (fig3b_transloc_breakpoints_ht assigned).\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# load_all_icle_data(load_external = TRUE)
# source("Helper_Scripts/14_SupFig8_TMB_SV_Preparation.R")  # Provides alt_count_chr
# source("Helper_Scripts/16_Fig3B_Translocation_Breakpoints.R")
#
# pdf(file.path(DIRS$results_sub$ogm, "Fig3B_Translocation_Distribution.pdf"), width = 6.5, height = 5)
# draw(fig3b_transloc_breakpoints_ht, merge_legends = TRUE)
# dev.off()
