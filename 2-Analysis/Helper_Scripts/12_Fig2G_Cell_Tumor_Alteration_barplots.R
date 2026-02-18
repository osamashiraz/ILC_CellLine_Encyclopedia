# ==============================================================================
# Script 12: Figure 2G - Cell vs Tumor CDH1 Alteration Barplots
# ==============================================================================
# Description: Generates CL/TCGA/FMI CDH1 alteration barplots and Fig 2F-style
#              CDH1 allele frequency jitter plots. Compares alteration frequencies
#              across cell lines and tumor cohorts.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - BRCA_CL_GAM: Cell line genomic alteration matrix
#   - TCGA_BRCA_GAM: TCGA genomic alteration matrix
#   - Foundation Medicine data files (via config paths)
#
# Output:
#   - Figure 2G barplot visualizations
#   - Figure 2F allele frequency plots
#
# Author: Osama Shiraz Shah
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(reshape2); library(readxl); library(ggplot2); library(ggthemes); library(gt) })

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}
if (!exists("BRCA_CL_GAM", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = TRUE)
}
# %notin% operator is defined in Helper_Functions.R
if (!exists("%notin%", envir = .GlobalEnv)) {
  `%notin%` <- Negate(`%in%`)
}
ensure_dir(DIRS$results_sub$cdh1)
out_dir <- DIRS$results_sub$cdh1

message("\n========================================")
message("Figure 2G: Cell vs Tumor CDH1 alteration barplots")
message("========================================")


# Cell Lines
cellids <- subset(CL_Annots, overlapWCCLE != "Y" & Name %in% colnames(BRCA_CL_GAM))$Name

CL_mat = data.frame(ID = cellids, Alt = BRCA_CL_GAM['CDH1', cellids],
                    Histology = CL_Annots[cellids, "Histology"])
CL_mat$Alt = gsub(x = CL_mat$Alt, pattern = ";", replacement = "+")
CL_mat$Alt[CL_mat$Alt %in% ""] = "WT"
CL_mat$Alt[CL_mat$Alt %in% "MUT+DEL"] = "DEL"
CL_mat$Alt[CL_mat$Alt %in% "MUT+AMP"] = "AMP"
CL_mat$Alt[CL_mat$Alt %in% "MUT+GAIN"] = "MUT"
CL_mat$Alt = gsub(x = CL_mat$Alt, pattern = "AMP", replacement = "GAIN")

# CL_mat$Alt_simple = ifelse(CL_mat$Alt %in% c("DEL", "MUT+DEL", "MUT+LOH"), "Dual Loss", "WT")

fig2g_cl_cdh1_alts_tbl <-  suppressMessages(CL_mat %>% group_by(Histology, Alt) %>% summarise(N=n()))

fig2g_cl_cdh1_alts_tbl$Histology = factor(fig2g_cl_cdh1_alts_tbl$Histology,levels = c("ILC", "ILC-like", "NST"))
fig2g_cl_cdh1_alts_tbl$Alt = factor(fig2g_cl_cdh1_alts_tbl$Alt, c("MUT+LOH", "DEL", "MUT", "LOH", "WT", "GAIN"))

fig2g_cl_cdh1_alts <- ggplot(fig2g_cl_cdh1_alts_tbl, aes(y = N, x = Histology, fill = Alt)) +
  geom_bar(stat = "identity", position = "fill", show.legend = FALSE) + theme_clean(20) +
  guides(fill = guide_legend(title = "CDH1 Alterations")) +
  scale_fill_manual("Alterations", values = annot_cols$Alt) + ylab("Frequency") + xlab("") +
  theme(panel.border = element_rect(color = "black", linewidth = 1))
assign("fig2g_cl_cdh1_alts", fig2g_cl_cdh1_alts, envir = .GlobalEnv)
assign("fig2g_cl_cdh1_alts_tbl", fig2g_cl_cdh1_alts_tbl, envir = .GlobalEnv)


## pvalue
CL_mat_sub = subset(CL_mat, Alt %notin% c("MUT","GAIN", "LOH"))
CL_mat_sub$Histology <- as.character(CL_mat_sub$Histology)
# CL_mat_sub$Histology[CL_mat_sub$Histology %in% c("ILC-like", "ILC")] = "ICLE"
CL_mat_sub$Event = ifelse(CL_mat_sub$Alt %notin% "WT", "Bi-allelic Inactivation", "WT")

cl_ilc_chi_pval <- chisq.test(table(subset(CL_mat_sub, Histology %in% c("ILC", "NST"))$Event, subset(CL_mat_sub, Histology %in% c("ILC", "NST"))$Histology))
cl_ilc_like_chi_pval <- chisq.test(table(subset(CL_mat_sub, Histology %in% c("ILC-like", "NST"))$Event, subset(CL_mat_sub, Histology %in% c("ILC-like", "NST"))$Histology))

fig2g_cl_cdh1_alts_tbl_pval <- suppressMessages(CL_mat_sub %>% group_by(Histology, Event) %>% 
  summarise(N = n()) %>% mutate(prop = 100*N/sum(N)))
fig2g_cl_cdh1_alts_tbl_pval[, "ILC_v_NST_pval"] <- cl_ilc_chi_pval$p.value
fig2g_cl_cdh1_alts_tbl_pval[, "ILClike_v_NST_pval"] <- cl_ilc_like_chi_pval$p.value

assign("fig2g_cl_cdh1_alts_tbl_pval", fig2g_cl_cdh1_alts_tbl_pval, envir = .GlobalEnv)






##### TCGA
TCGA_mat = data.frame(ID = colnames(TCGA_BRCA_GAM), Alt = TCGA_BRCA_GAM['CDH1',], 
                             Histology = TCGA_Annots[colnames(TCGA_BRCA_GAM), "Final Pathology"])
TCGA_mat$Alt[TCGA_mat$Alt %in% ""] = "WT"
TCGA_mat$Alt[TCGA_mat$Alt %in% "MUT+DEL"] = "DEL"
TCGA_mat$Alt[TCGA_mat$Alt %in% "MUT+AMP"] = "AMP"
TCGA_mat$Alt[TCGA_mat$Alt %in% "MUT+GAIN"] = "MUT"
TCGA_mat$Alt = gsub(x = TCGA_mat$Alt, pattern = ";", replacement = "+")

TCGA_mat = subset(TCGA_mat, Histology %in% c("ILC", "NST")) 

table(TCGA_mat$Alt, TCGA_mat$Histology)
TCGA_mat$Alt = gsub(x = TCGA_mat$Alt, pattern = "AMP", replacement = "GAIN")
TCGA_mat$Alt[TCGA_mat$Alt %in% "MUT+DEL"] = "DEL"

TCGA_mat$Alt_simple = ifelse(TCGA_mat$Alt %in% c("DEL", "MUT+DEL", "MUT+LOH"), "Dual Loss", "WT")

  
fig2g_tcga_cdh1_alts_tbl = suppressMessages(TCGA_mat %>% group_by(Histology, Alt) %>% summarise(N=n()))

fig2g_tcga_cdh1_alts_tbl$Histology = factor(fig2g_tcga_cdh1_alts_tbl$Histology,levels = c("ILC", "ILC-like", "NST"), labels = c("ILC", "ILC-like","NST"))
fig2g_tcga_cdh1_alts_tbl$Alt = factor(fig2g_tcga_cdh1_alts_tbl$Alt, levels = c("MUT+LOH", "DEL", "MUT", "LOH", "WT", "GAIN"))

fig2g_tcga_cdh1_alts <- ggplot(fig2g_tcga_cdh1_alts_tbl, aes(y = N, x = Histology, fill = Alt)) +
  geom_bar(stat = "identity", position = "fill", show.legend = FALSE) + theme_clean(20) +
  guides(fill = guide_legend(title = "CDH1 Alterations")) +
  scale_fill_manual("Alterations", values = annot_cols$Alt) + ylab("Frequency") + xlab("") +
  theme(panel.border = element_rect(color = "black", linewidth = 1))
assign("fig2g_tcga_cdh1_alts", fig2g_tcga_cdh1_alts, envir = .GlobalEnv)
assign("fig2g_tcga_cdh1_alts_tbl", fig2g_tcga_cdh1_alts_tbl, envir = .GlobalEnv)


## pvalue
TCGA_mat_sub = subset(TCGA_mat, Alt %notin% c("MUT","GAIN", "LOH"))
TCGA_mat_sub$Event = ifelse(TCGA_mat_sub$Alt %notin% c("WT"), "Bi-allelic Inactivation", "WT")


tcga_chi_pval <- chisq.test(table(TCGA_mat_sub$Event, TCGA_mat_sub$Histology))


fig2g_tcga_cdh1_alts_tbl_pval = suppressMessages(TCGA_mat_sub %>% group_by(Histology, Event) %>% 
  summarise(N = n()) %>% mutate(prop = 100*N/sum(N)))
fig2g_tcga_cdh1_alts_tbl_pval[, "ILC_v_NST_pval"] <- tcga_chi_pval$p.value

assign("fig2g_tcga_cdh1_alts_tbl_pval", fig2g_tcga_cdh1_alts_tbl_pval, envir = .GlobalEnv)


message("  ✓ Fig 2G complete (fig2g_cl_cdh1_alts, fig2g_tcga_cdh1_alts, fig2g_cl_cdh1_alts_tbl_pval, fig2g_tcga_cdh1_alts_tbl_pval assigned).\n")


# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/11_Fig2G_Cell_Tumor_Alteration_barplots.R")
#
# load_all_icle_data(load_external = TRUE)
# # Script runs on source. Assigns cdh1_af, CL_cdh1_alt, TCGA_cdh1_alt, FMI_*.
# # Save figures inline in Rmd (Fig 2G chunk).













