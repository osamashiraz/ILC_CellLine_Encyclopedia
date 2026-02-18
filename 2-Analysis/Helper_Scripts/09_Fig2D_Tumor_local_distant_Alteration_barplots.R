# ==============================================================================
# Script 09: Figure 2D - FMI CDH1 Alterations (Local vs Distant Metastasis)
# ==============================================================================
# Description: Foundation Medicine breast diagnostic cohort: CDH1 alteration
#              barplots by histology (ILC vs NST) and metastasis status
#              (Local vs Distant).
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - Foundation Medicine data files (via config paths)
#
# Output:
#   - Figure 2D barplot visualizations
#
# Author: Osama Shiraz Shah
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(readxl)
  library(ggplot2)
  library(ggthemes)
  library(gt)
})

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}
# %notin% operator is defined in Helper_Functions.R
if (!exists("%notin%", envir = .GlobalEnv)) {
  `%notin%` <- Negate(`%in%`)
}
ensure_dir(DIRS$results_sub$cdh1)
out_dir <- DIRS$results_sub$cdh1

message("\n========================================")
message("Figure 2D: FMI CDH1 alterations (Local vs Distant)")
message("========================================")

# Load FMI data if not already loaded
if (!exists("FMI_annot", envir = .GlobalEnv)) {
  FMI_annot <- readxl::read_excel(FILES$fmi_breast_diag, sheet = 3) %>% as.data.frame()
  FMI_annot <- FMI_annot[, c("xrn", "disease_ontology_term", "tissue", "local_met_status", "ER_positive")]
  rownames(FMI_annot) <- FMI_annot$xrn
  FMI_annot$disease_ontology_term[grep(x = FMI_annot$disease_ontology_term, pattern = "ilc")] <- "ILC"
  FMI_annot$disease_ontology_term[grep(x = FMI_annot$disease_ontology_term, pattern = "idc")] <- "NST"
  FMI_annot <- subset(FMI_annot, disease_ontology_term %in% c("NST", "ILC"))
} else {
  FMI_annot <- get("FMI_annot", envir = .GlobalEnv)
}

if (!exists("FMI", envir = .GlobalEnv)) {
  FMI <- readxl::read_excel(FILES$fmi_breast_diag, sheet = 2) %>% as.data.frame()
} else {
  FMI <- get("FMI", envir = .GlobalEnv)
}

collapse_values <- function(x) paste0(x, collapse = "-")



FMI_sub <- subset(FMI, xrn %in% FMI_annot$xrn)
FMI_sub$alt <- ifelse(FMI_sub$alteration_type == "CN" & FMI_sub$coding_type %in% c("deletion", "amplification"),
                      paste0(FMI_sub$alteration_type, "-", FMI_sub$coding_type), FMI_sub$alteration_type)

FMI_mat <- reshape2::dcast(FMI_sub, xrn ~ gene, fun.aggregate = collapse_values, value.var = "alt", fill = "WT")
FMI_mat <- FMI_mat %>% mutate(CDH1_Simple = case_when(
  CDH1 %in% c("CN-amplification", "CN-amplification-RE") ~ "GAIN",
  CDH1 %in% c("CN-deletion-RE", "RE") ~ "RE",
  CDH1 %in% c("CN-deletion") ~ "DEL",
  CDH1 %in% c("CN-deletion-SV", "RE-SV", "SV", "SV-SV", "SV-SV-SV") ~ "MUT"
))
rownames(FMI_mat) <- FMI_mat$xrn
FMI_mat[FMI_annot$xrn, "disease_ontology_term"] <- FMI_annot$disease_ontology_term
FMI_mat[FMI_annot$xrn, "tissue"] <- FMI_annot$tissue
FMI_mat[FMI_annot$xrn, "local_met_status"] <- FMI_annot$local_met_status
names(FMI_mat) <- c("ID", "CDH1_Alt", "CDH1_Alt_Simple", "Histology", "Tissue_Source", "Metastasis_Status")
FMI_mat[is.na(FMI_mat$CDH1_Alt_Simple), "CDH1_Alt_Simple"] <- "WT"
FMI_mat[is.na(FMI_mat$CDH1_Alt), "CDH1_Alt"] <- "WT"
FMI_mat <- subset(FMI_mat, Metastasis_Status %notin% c("unknown", "ambiguous", "ln"))
FMI_mat$Metastasis_Status <- ifelse(FMI_mat$Metastasis_Status == "local", "Local", "Distant")

FMI_mat_sub <- subset(FMI_mat, CDH1_Alt_Simple %notin% c("MUT", "GAIN"))

fig2d_fmi_alts_tbl <- suppressMessages(FMI_mat_sub %>% group_by(Histology, Metastasis_Status, CDH1_Alt_Simple) %>%
  summarise(N = n()) %>% mutate(prop = 100 * N / sum(N)))

fig2d_fmi_alts_tbl$Histology <- factor(fig2d_fmi_alts_tbl$Histology, levels = c("ILC", "NST"))
fig2d_fmi_alts_tbl$Metastasis_Status <- factor(fig2d_fmi_alts_tbl$Metastasis_Status, levels = c("Local", "Distant"))

fig2d_fmi_alts <- ggplot(subset(fig2d_fmi_alts_tbl, CDH1_Alt_Simple %notin% c("WT")),
                   aes(y = prop, x = Histology, fill = reorder(CDH1_Alt_Simple, prop))) +
  geom_bar(stat = "identity", show.legend = TRUE) + theme_clean(14) +
  guides(fill = guide_legend(title = "CDH1 Alteration Types")) +
  ylab("Fraction") + xlab("") + scale_fill_manual("Alterations", values = annot_cols$Alt_FMI) +
  facet_grid(. ~ Metastasis_Status) + scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))

FMI_mat_sub$Event <- ifelse(FMI_mat_sub$CDH1_Alt_Simple %notin% c("WT", "GAIN", "MUT"), "DEL_RE", "No_DEL_RE")
fmi_pri_chi_pval <- chisq.test(table(subset(FMI_mat_sub, Metastasis_Status == "Local")$Event,
                                     subset(FMI_mat_sub, Metastasis_Status == "Local")$Histology))
fmi_met_chi_pval <- chisq.test(table(subset(FMI_mat_sub, Metastasis_Status == "Distant")$Event,
                                     subset(FMI_mat_sub, Metastasis_Status == "Distant")$Histology))

fig2d_fmi_alts_tbl_pval <- FMI_mat_sub %>% group_by(Histology, Metastasis_Status, Event) %>%
  summarise(N = n(), .groups = "drop") %>% mutate(prop = 100 * N / sum(N))
fig2d_fmi_alts_tbl_pval$local_ILC_v_NST_pval <- fmi_pri_chi_pval$p.value
fig2d_fmi_alts_tbl_pval$distant_ILC_v_NST_pval <- fmi_met_chi_pval$p.value

assign("fig2d_fmi_alts", fig2d_fmi_alts, envir = .GlobalEnv)
assign("fig2d_fmi_alts_tbl", fig2d_fmi_alts_tbl, envir = .GlobalEnv)
assign("fig2d_fmi_alts_tbl_pval", fig2d_fmi_alts_tbl_pval, envir = .GlobalEnv)
message("  ✓ Fig 2D complete (fig2d_fmi_alts, fig2d_fmi_alts_tbl, fig2d_fmi_alts_tbl_pval assigned).\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/09_Fig2D_Tumor_local_distant_Alteration_barplots.R")
#
# # Script runs on source. Requires FILES$fmi_breast_diag.
#
# fig2d_fmi_alts, fig2d_fmi_alts_tbl, fig2d_fmi_alts_tbl_pval assigned; save inline in Rmd.







