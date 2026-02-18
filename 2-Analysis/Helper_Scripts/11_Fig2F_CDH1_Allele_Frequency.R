# ==============================================================================
# Script 11: Figure 2F - CDH1 Allele Frequency in Cell Lines
# ==============================================================================
# Description: Analyzes and visualizes CDH1 mutation allele frequencies across
#              ICLE cell lines by histology. Generates jitter plot showing
#              allele frequency distribution by mutation type.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - CL_Annots: Cell line annotations
#   - BRCA_CL_MAF: Cell line mutation data
#
# Output:
#   - fig2f_cdh1_af: CDH1 allele frequency plot (assigned to .GlobalEnv)
#
# Author: Osama Shiraz Shah
# ==============================================================================

# Check for config and DIRS
if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}

# Check for %notin% operator
if (!exists("%notin%", envir = .GlobalEnv)) {
  `%notin%` <- Negate(`%in%`)
}

# Check for required data objects
if (!exists("CL_Annots", envir = .GlobalEnv)) {
  message("  CL_Annots not found. Loading annotations...")
  if (!exists("load_all_icle_data", envir = .GlobalEnv)) {
    source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  }
  load_all_icle_data(load_external = FALSE)
}

if (!exists("BRCA_CL_MAF", envir = .GlobalEnv)) {
  message("  BRCA_CL_MAF not found. Loading mutation data...")
  if (!exists("load_all_icle_data", envir = .GlobalEnv)) {
    source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  }
  load_all_icle_data(load_external = FALSE)
}

if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})


message("\n========================================")
message("Figure 2F: CDH1 allele frequency plot")
message("========================================")

CL_af_ids <- unique(c(subset(CL_Annots, overlapWCCLE != "Y" & `mRNA Subtypes` %notin% c("Basal"))$Name, grep("-I", CL_Annots$Name, value = TRUE)))

cdh1_mutations <- BRCA_CL_MAF %>%
  filter(Hugo_Symbol == "CDH1", !is.na(t_AF), !is.na(Histology), Broad_Variant_Classification != "Other",
         Tumor_Sample_Barcode %in% subset(CL_Annots, overlapWCCLE != "Y")$Name)
cdh1_mutations$Broad_Variant_Classification <- factor(cdh1_mutations$Broad_Variant_Classification, levels = c("Truncating", "Missense", "Silent"))

fig2f_cdh1_af <- ggplot(cdh1_mutations, aes(x = Histology, y = t_AF, shape = Histology)) +
  geom_jitter(aes(fill = Broad_Variant_Classification, color = Broad_Variant_Classification), width = 0.2, size = 3, stroke = 1, alpha = 1) +
  scale_color_manual("Type", values = annot_cols$Mutation) + scale_fill_manual("Type", values = annot_cols$Mutation) +
  scale_y_continuous(name = "CDH1 Mutation Allele Frequency", breaks = seq(0, 1, 0.1)) +
  scale_shape_manual(values = c(21, 22, 2)) +
  xlab("") + theme_clean(15) + theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1))
assign("fig2f_cdh1_af", fig2f_cdh1_af, envir = .GlobalEnv)

message("  ✓ Fig 2F complete (fig2f_cdh1_af assigned).\n")


# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# load_all_icle_data(load_external = FALSE)
# source("Helper_Scripts/11_Fig2F_CDH1_Allele_Frequency.R")
#
# # Save figure
# ensure_dir(DIRS$results_sub$cdh1)
# ggsave(file.path(DIRS$results_sub$cdh1, "Fig2F_CDH1_Allele_Frequency.pdf"), 
#        fig2f_cdh1_af, width = 7, height = 5)