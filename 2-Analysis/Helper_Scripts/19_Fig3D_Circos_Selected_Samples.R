# ==============================================================================
# Script 20: Figure 3D - Circos Plots for Selected Samples (BCK4, 600MPE, HCC2185, ZR7530)
# ==============================================================================
# Description: Generates circos plots for 4 selected ICLE cell lines showing
#              structural variations, copy number alterations, and chromothripsis
#              events. Requires shatterseek_outs from 19_Fig3C_Chromothripsis.R.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#   - 18_Fig3D_SupFig9_Prepare_Circos_Visualizations.R (sourced by this script; provides
#     draw_circos_plot_for_all, genome_tracks, prepare_sample_tracks, draw_circos_plot)
#   - 19_Fig3C_Chromothripsis.R (for shatterseek_outs)
#
# Required Data Objects:
#   - ICLE_SV: Structural variation data
#   - shatterseek_outs: ShatterSeek analysis results (from Fig3C)
#   - CL_Annots: Cell line annotations
#
# Output:
#   - Circos plot objects for selected samples (assigned to .GlobalEnv)
#
# Author: Osama Shiraz Shah
# ==============================================================================

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}

message("\n========================================")
message("Figure 3D: Circos plots for selected samples")
message("========================================")

if (!exists("shatterseek_outs", envir = .GlobalEnv)) {
  stop("shatterseek_outs not found. Source 19_Fig3C_Chromothripsis.R first.")
}

trk_bck4 <- c(prepare_sample_tracks("BCK4", ICLE_SV, shatterseek_outs$cnv_df, shatterseek_outs$chromothripsis_df), genome_tracks)
trk_600 <- c(prepare_sample_tracks("600MPE", ICLE_SV, shatterseek_outs$cnv_df, shatterseek_outs$chromothripsis_df), genome_tracks)
trk_hcc <- c(prepare_sample_tracks("HCC2185", ICLE_SV, shatterseek_outs$cnv_df, shatterseek_outs$chromothripsis_df), genome_tracks)
trk_zr <- c(prepare_sample_tracks("ZR7530", ICLE_SV, shatterseek_outs$cnv_df, shatterseek_outs$chromothripsis_df), genome_tracks)

assign("trk_bck4", trk_bck4, envir = .GlobalEnv)
assign("trk_600", trk_600, envir = .GlobalEnv)
assign("trk_hcc", trk_hcc, envir = .GlobalEnv)
assign("trk_zr", trk_zr, envir = .GlobalEnv)
message("  ✓ Fig 3D complete (trk_bck4, trk_600, trk_hcc, trk_zr assigned).\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# load_all_icle_data(load_external = TRUE)
# source("Helper_Scripts/19_Fig3C_Chromothripsis.R")  # provides shatterseek_outs
# source("Helper_Scripts/20_Fig3D_Circos_Selected_Samples.R")
#
# fig3d_dir <- DIRS$results_sub$ogm
# pdf(file.path(fig3d_dir, "Fig3D_Circos_BCK4.pdf"), width = 8, height = 8)
# draw_circos_plot(trk_bck4)
# dev.off()
# # Repeat for other samples...
