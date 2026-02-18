# ==============================================================================
# Script 08 (orchestrator): Figure 2 (CDH1) – run all Fig 2 scripts in order
# ==============================================================================
# Description: Sources 08_Fig2C through 13_Fig2H in order so the Rmd can call
#              this single script for the full CDH1 alteration landscape.
#              All figure objects (fig2c_*, fig2d_*, etc.) are assigned by
#              the child scripts to .GlobalEnv.
#
# Order: 08 (Fig2C exonic deletions) → 09 (Fig2D FMI) → 10 (Fig2E Protein Paint)
#        → 11 (Fig2F allele freq) → 12 (Fig2G cell/tumor barplots) → 13 (Fig2H landscape)
#
# Author: Osama Shiraz Shah
# ==============================================================================

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}

helpers <- DIRS$scripts$helpers
message("═══════════════════════════════════════════════════════")
message("  Figure 2C-H – running all panels in order")
message("═══════════════════════════════════════════════════════\n")

suppressWarnings({
  source(file.path(helpers, "08_Fig2C_CDH1_exonic_deletions.R"), chdir = TRUE)
  source(file.path(helpers, "09_Fig2D_Tumor_local_distant_Alteration_barplots.R"), chdir = TRUE)
  source(file.path(helpers, "10_Fig2E_CDH1_Protein_paint_plot.R"), chdir = TRUE)
  source(file.path(helpers, "11_Fig2F_CDH1_Allele_Frequency.R"), chdir = TRUE)
  source(file.path(helpers, "12_Fig2G_Cell_Tumor_Alteration_barplots.R"), chdir = TRUE)
  source(file.path(helpers, "13_Fig2H_CDH1_Alteration_Landscape.R"), chdir = TRUE)
})
message("═══════════════════════════════════════════════════════")
message("  Figure 2C-H complete")
message("═══════════════════════════════════════════════════════\n")
