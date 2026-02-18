# ==============================================================================
# Script 15: Figure 3A - ICLE Genomic Instability Metrics and SV Overview
# ==============================================================================
# Description: Generates Fig 3A left (genomic instability metrics heatmap) and
#              Fig 3A right (SV type distribution plot). Visualizes genomic
#              instability patterns across ICLE cell lines.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - CL_Annots: Cell line annotations
#   - ICLE_SV: Structural variation data
#   - Genomic instability metrics (FGA, TMB, etc.)
#
# Outputs (assigned to .GlobalEnv):
#   - fig3a_genomic_instability: Genomic instability heatmap (Fig 3A left)
#   - fig3a_sv_distribution: SV type distribution ggplot (Fig 3A right)
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
  library(grid)
  library(ggplot2)
})

message("\n========================================")
message("Figure 3A: Genomic instability metrics and SV distribution")
message("========================================")

# Fig 3A Left: Genomic Instability Metrics
sv_cells <- paste0(names(sort(table(ICLE_SV$Sample))), '-I')

myAnnoBarplot_row <- function(x, col = "#B0BEC5", which = 'row') {
  return(anno_barplot(-x, gp = gpar(fill = col), axis = FALSE, which = which, 
                     border_gp = gpar(col = "gray"), width = unit(6, "mm"), 
                     height = unit(4, "mm"), bar_width = 0.4))
}

left_annotation <- HeatmapAnnotation(
  Subtype = CL_Annots[rev(sv_cells), "mRNA Subtypes"], 
  height = unit(7, "cm"), 
  width = unit(3, "cm"),
  FGA = myAnnoBarplot_row(CL_Annots[rev(gsub("SKBR3-I", "SKBR3-C", sv_cells)), "FGA"], "#F48FB1", 'row'),
  TMB = myAnnoBarplot_row(CL_Annots[rev(gsub("SKBR3-I", "SKBR3-C", sv_cells)), "TMB"], "#D1C4E9", 'row'),
  Fusions = myAnnoBarplot_row(CL_Annots[rev(sv_cells), "Fusions"], "#C8E6C9", 'row'),
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold", col = "black", border = TRUE, fontfamily = "Helvetica"),
  gap = unit(0.5, "mm"),
  gp = gpar(col = NA, lwd = 0.6),
  annotation_legend_param = heatmap_legend_param,
  simple_anno_size = unit(3, "mm"),
  annotation_name_side = "top",
  col = list(Subtype = annot_cols$Subtypes),
  which = "row"
)

set.seed(123)
fig3a_genomic_instability <- Heatmap(
  t(t(rev(sv_cells))),
  row_title = " ",
  name = "Cell Lines",
  col = setNames(rep("white", length(sv_cells)), sv_cells),
  show_heatmap_legend = FALSE,
  left_annotation = left_annotation,
  row_labels = gsub("-I", "", rev(sv_cells)),
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  border = TRUE,
  border_gp = gpar(col = "black"),
  width = unit(0, "cm"),
  height = unit(7, "cm"),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_column_names = TRUE,
  column_title = " ",
  na_col = "gray",
  heatmap_legend_param = heatmap_legend_param,
  column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = "Helvetica"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Helvetica"),
  row_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = "Helvetica"),
  row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Helvetica")
)





# Fig 3A Right: SV Distribution

#' Plot SV type distribution by sample
#' 
#' @param sv_data Data frame with SV data
#' @param clin_data Data frame with clinical data
#' @param sample_col Column name for sample IDs
#' @param type_col Column name for SV type
#' @param histology_col Column name for histology
#' @param ilc_label Label for ILC samples
#' @return ggplot object
plot_sv_type_distribution <- function(
    sv_data, clin_data, sample_col = "Sample",
    type_col = "Type", histology_col = "Histology", ilc_label = "ILC"
) {
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(scales)
  
  # Summarize SV counts and frequencies by sample and SV type
  sv_type_summary <- sv_data %>%
    group_by(SampleID = .data[[sample_col]], SV_Type = .data[[type_col]]) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(SampleID) %>%
    mutate(
      Frequency = Count / sum(Count),
      TotalSVs = sum(Count)
    ) %>%
    ungroup()
  
  sv_type_summary$SV_Type <- as.character(sv_type_summary$SV_Type)
  
  # Order samples by total SV count
  total_svs_per_sample <- sv_type_summary %>%
    group_by(SampleID) %>%
    summarise(TotalSVs = sum(Count), .groups = "drop")
  sample_order <- order(total_svs_per_sample$TotalSVs)
  
  # Assign colors by histology subtype
  ilc_samples <- clin_data %>%
    filter(.data[[histology_col]] == ilc_label) %>%
    pull(.data[[sample_col]])
  sample_colors <- annot_cols$Histology[ifelse(
    unique(sv_type_summary$SampleID) %in% ilc_samples, "ILC", "ILC-like")]
  
  # Assign distinct colors for each SV type
  sv_type_palette <- ggsci::pal_npg()(length(unique(sv_type_summary$SV_Type)))
  names(sv_type_palette) <- sort(unique(sv_type_summary$SV_Type))
  
  # Plot: Stacked bar chart of SV types per sample
  p <- ggplot(sv_type_summary, aes(x = reorder(SampleID, TotalSVs), y = Count, fill = SV_Type)) +
    geom_bar(position = "stack", stat = "identity", width = 0.70, color = "black", lwd = 0.2) + 
    xlab("") + ylab("Count") + labs(fill = "SV Type") + 
    ggpubr::theme_classic2(base_size = 14) + 
    theme(
      axis.text = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      # axis.text.y = element_text(colour = sample_colors[sample_order])
    ) +
    coord_flip() +
    scale_fill_manual(values = sv_type_palette)
  
  return(p)
}

fig3a_sv_distribution <- plot_sv_type_distribution(ICLE_SV, CL_Annots) +
    scale_y_continuous(breaks = seq(0, 550, 100)) +
    theme(panel.grid.major = element_line(color = "gray60", linetype = "dashed", linewidth = 0.1))

assign("fig3a_sv_distribution", fig3a_sv_distribution, envir = .GlobalEnv)
assign("fig3a_genomic_instability", fig3a_genomic_instability, envir = .GlobalEnv)
message("  ✓ Fig 3A complete\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# load_all_icle_data(load_external = TRUE)
# source("Helper_Scripts/14_Fig3A_Genomic_Instability.R")
#
# ensure_dir(DIRS$results_sub$ogm)
# pdf(file.path(DIRS$results_sub$ogm, "Fig3A_left_Metrics_of_GenomicInstability.pdf"), width = 8, height = 8)
# draw(fig3a_genomic_instability)
# dev.off()
# ggsave(file.path(DIRS$results_sub$ogm, "Fig3A_right_SV_Distribution.pdf"), fig3a_sv_distribution, width = 6, height = 5)
