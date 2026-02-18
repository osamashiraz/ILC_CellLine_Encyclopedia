# ==============================================================================
# Script 06: Supplementary Figure 5 - Key ILC vs NST Alterations (Patient Tumors)
# ==============================================================================
# Description: Builds frequency table and barplot of LumA ILC vs NST alterations
#              from TCGA GAM. Compares alteration frequencies between ILC and NST
#              tumor subtypes.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - TCGA_BRCA_GAM: TCGA genomic alteration matrix
#   - TCGA_Annots: TCGA tumor annotations
#
# Output:
#   - supfigs5_tumor_alterations: Barplot (assigned to .GlobalEnv)
#   - freq_tbl: Frequency table (assigned to .GlobalEnv)
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
  library(reshape2)
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
})

# ------------------------------------------------------------------------------
# Note: build_freq_tbl_from_GAM is now defined in Helper_Functions.R
#       (loaded via Helper_Functions.R at the top of this script)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SupFig 5 analysis
# ------------------------------------------------------------------------------
message("\n========================================")
message("SupFig 5: ILC vs NST alteration frequencies (LumA tumors)")
message("========================================")

# Check if required data objects exist, load if needed
required_objects <- c("BRCA_CL_GAM", "TCGA_BRCA_GAM", "TCGA_Annots")
missing_objects <- required_objects[!sapply(required_objects, exists, envir = .GlobalEnv)]

if (length(missing_objects) > 0) {
  message("  Loading missing data objects: ", paste(missing_objects, collapse = ", "))
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = TRUE)
}

BRCA_CL_GAM_Simple <- BRCA_CL_GAM
for (alt in c("LOH", "GAIN")) BRCA_CL_GAM_Simple[BRCA_CL_GAM_Simple == alt] <- ""

TCGA_BRCA_GAM_Simple <- TCGA_BRCA_GAM
for (alt in c("LOH", "GAIN")) TCGA_BRCA_GAM_Simple[TCGA_BRCA_GAM_Simple == alt] <- ""

ILC_genes <- c("CDH1", "PTEN", "RUNX1", "FOXA1", "TBX3")
NST_genes <- c("TP53", "MYC", "GRIN2A", "CSMD1", "DDX5", "PPM1D", "MAP2K4", "MAP3K1", "GATA3")

tcga_luminal_subset <- subset(TCGA_Annots, Case.ID %in% colnames(TCGA_BRCA_GAM_Simple) &
  PAM50 %in% "LumA" & `Final Pathology` %in% c("ILC", "NST"))

freq_tbl <- build_freq_tbl_from_GAM(
  tumor_gam    = TCGA_BRCA_GAM_Simple[, tcga_luminal_subset$Case.ID],
  genes       = c(ILC_genes, NST_genes),
  tumor_labels = setNames(tcga_luminal_subset$`Final Pathology`, tcga_luminal_subset$Case.ID),
  ilc_label   = "ILC",
  nst_label   = "NST",
  min_count   = 2,
  fisher_alternative = "two.sided",
  debug       = TRUE,
  allowed_classes   = c("MUT", "MUT;LOH", "MUT;AMP", "MUT;DEL", "MUT;GAIN", "AMP", "DEL")
)
rownames(freq_tbl) <- freq_tbl$gene

ILC_LumA_N <- nrow(subset(tcga_luminal_subset, `Final Pathology` == "ILC"))
NST_LumA_N <- nrow(subset(tcga_luminal_subset, `Final Pathology` == "NST"))

ilc_nst_alt_df <- reshape2::melt(TCGA_BRCA_GAM[c(ILC_genes, NST_genes), tcga_luminal_subset$Case.ID])
ilc_nst_alt_df$Histology <- TCGA_Annots[as.character(ilc_nst_alt_df$Var2), "Final Pathology"]
ilc_nst_alt_df <- subset(ilc_nst_alt_df, value %notin% c("", "LOH", "GAIN"))
names(ilc_nst_alt_df) <- c("Gene", "ID", "Alt", "Histology")

ilc_nst_alt_df <- ilc_nst_alt_df[, c(1, 3, 4)] %>%
  group_by(Histology, Gene, Alt) %>%
  mutate(count = n()) %>%
  as.data.frame()
ilc_nst_alt_df <- ilc_nst_alt_df[!duplicated(paste0(ilc_nst_alt_df$Gene, ilc_nst_alt_df$Alt, ilc_nst_alt_df$Histology)), ]

ilc_nst_alt_df$Freq <- ifelse(ilc_nst_alt_df$Histology == "ILC",
  100 * ilc_nst_alt_df$count / ILC_LumA_N, 100 * ilc_nst_alt_df$count / NST_LumA_N)
ilc_nst_alt_df$Freq_rel <- ifelse(ilc_nst_alt_df$Histology == "ILC",
  100 * ilc_nst_alt_df$count / ILC_LumA_N, -100 * ilc_nst_alt_df$count / NST_LumA_N)

geneAltCount <- ilc_nst_alt_df %>% group_by(Gene) %>% summarise_at(.vars = "count", .funs = sum) %>% as.data.frame()
geneAltCount$Gene <- as.character(geneAltCount$Gene)
rownames(geneAltCount) <- geneAltCount$Gene
geneAltCount$count <- ifelse(geneAltCount$Gene %in% ILC_genes, geneAltCount$count, -geneAltCount$count)

events_ordered <- geneAltCount[order(geneAltCount$count, decreasing = FALSE), "Gene"]
ilc_nst_alt_df$Gene <- factor(ilc_nst_alt_df$Gene, levels = events_ordered, labels = events_ordered)

ilc_nst_alt_df <- ilc_nst_alt_df %>%
  mutate(Alt = gsub(";", "+", Alt), Alt = ifelse(Alt %in% c("MUT+AMP", "MUT+DEL", "MUT+GAIN", "MUT+LOH"), "MUT", Alt)) %>%
  mutate(Histology = factor(Histology, levels = c("NST", "ILC")))

supfigs5_tumor_alterations <- ggplot(ilc_nst_alt_df, aes(x = Gene, y = Freq_rel, fill = Alt)) +
  geom_bar(stat = "identity", width = 0.6) + coord_flip() +
  ylab("Alteration Frequency") + xlab("") +
  scale_fill_manual("Alteration", values = annot_cols$Alt_Simple) +
  theme_clean(base_size = 20) +
  theme(panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.5, linetype = "dashed")) +
  facet_grid(. ~ Histology, scales = "free")

assign("supfigs5_tumor_alterations", supfigs5_tumor_alterations, envir = .GlobalEnv)
assign("freq_tbl", freq_tbl, envir = .GlobalEnv)
message("  ✓ SupFig 5 complete (supfigs5_tumor_alterations, freq_tbl assigned).\n")

# ==============================================================================
# Usage Example 
# ==============================================================================
# source("config.R")
# source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
# source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
# load_all_icle_data(load_external = TRUE)
# source(file.path(DIRS$scripts$helpers, "06_SupFig5_ILC_NST_Alterations.R"))
#
# ggsave(file.path(DIRS$results, "SupFig5_BRCA_Tumor_Top_Alterations.pdf"), SupFigS5, width = 8, height = 5)
# ensure_dir(DIRS$results_sub$molecular_resemblance)
# write.table(freq_tbl, file.path(DIRS$results_sub$molecular_resemblance, "SupTable_ILC_NST_alteration_frequencies.tsv"),
#             sep = "\t", quote = FALSE, row.names = FALSE)
# Outputs: supfigs5_tumor_alterations, freq_tbl in .GlobalEnv; save PDF and TSV as above.
