# ==============================================================================
# Script 14 (orchestrator): Figure 3 + SupFig 8–10 – run all Fig 3 / SV scripts in order
# ==============================================================================
# Description: Sources 14 (TMB/SV prep, SupFig8), 15 (Fig3A), 16 (Fig3B),
#              17 (Fig3C), 18 (SupFig9 circos prep), 19 (Fig3D), 20 (Fig3E/3F/SupFig10)
#              so the Rmd can call this single script.
#              Correct execution order: 14 first (provides alt_count_chr etc.),
#              then 15, 16, 17, 18, 19, 20.
#
# Author: Osama Shiraz Shah
# ==============================================================================

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}

helpers <- DIRS$scripts$helpers
message("═══════════════════════════════════════════════════════")
message("  Figure 3 + SupFig 8–10 – running all panels in order")
message("═══════════════════════════════════════════════════════\n")

suppressWarnings({
  source(file.path(helpers, "14_SupFig8_TMB_SV_Preparation.R"), chdir = TRUE)
  source(file.path(helpers, "15_Fig3A_Genomic_Instability.R"), chdir = TRUE)
  source(file.path(helpers, "16_Fig3B_Translocation_Breakpoints.R"), chdir = TRUE)
  source(file.path(helpers, "17_Fig3C_Identify_Chromothripsis.R"), chdir = TRUE)
  source(file.path(helpers, "18_Fig3D_SupFig9_Prepare_Circos_Visualizations.R"), chdir = TRUE)
  source(file.path(helpers, "19_Fig3D_Circos_Selected_Samples.R"), chdir = TRUE)
  source(file.path(helpers, "20_Fig3E_3F_SupFig10_SV_Fusions.R"), chdir = TRUE)
})
message("\n═══════════════════════════════════════════════════════")
message("  Figure 3 + SupFig 8–10 complete")
message("═══════════════════════════════════════════════════════\n")
