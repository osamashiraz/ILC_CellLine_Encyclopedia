# ==============================================================================
# Script 14: Supplementary Figure 8 - TMB vs SV Preparation
# ==============================================================================
# Description: Prepares chrom_sizes, sv_tmb_summary, alt_count_chr, tmb_per_chr,
#              and SupFig8A-D plot objects. Analyzes relationship between tumor
#              mutational burden (TMB) and structural variations (SV).
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - ICLE_SV: Structural variation data
#   - TMB data (via config paths)
#   - BRCA_CL_MAF: Mutation data
#
# Output:
#   - chrom_sizes: Chromosome size data (assigned to .GlobalEnv)
#   - sv_tmb_summary: SV-TMB summary statistics (assigned to .GlobalEnv)
#   - alt_count_chr: Alteration counts per chromosome (assigned to .GlobalEnv)
#   - tmb_per_chr: TMB per chromosome (assigned to .GlobalEnv)
#   - SupFig8A-D: Plot objects (assigned to .GlobalEnv)
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
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  library(ggpubr)
  library(maftools)
  library(matrixStats)
  library(tidyr)
})

# %notin% operator is defined in Helper_Functions.R
if (!exists("%notin%", envir = .GlobalEnv)) {
  `%notin%` <- Negate(`%in%`)
}


message("\n========================================")
message("SupFig 8: TMB vs SV Comparison")
message("========================================")

# Check if required data objects exist, load if needed
if (!exists("BRCA_CL_MAF", envir = .GlobalEnv) || !exists("ICLE_SV", envir = .GlobalEnv) || !exists("CL_Annots", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_snv_data()
  load_sv_data()
}

# Load chromosome sizes if not already loaded
if (!exists("chrom_sizes", envir = .GlobalEnv)) {
  chrom_sizes <- read.delim(FILES$hg19_chrom_sizes, header = FALSE)
} else {
  chrom_sizes <- get("chrom_sizes", envir = .GlobalEnv)
}
names(chrom_sizes) <- c("Chromosome", "size_bp")
chrom_sizes$size_Mb <- chrom_sizes$size_bp / 1e6
rownames(chrom_sizes) <- chrom_sizes$Chromosome
chrom_sizes <- chrom_sizes[grep("_", chrom_sizes$Chromosome, value = TRUE, invert = TRUE), ]

hg19_wgs_size_mb <- sum(chrom_sizes$size_Mb)

sv_tmb_summary <- ICLE_SV %>% group_by(Sample) %>% summarise(SV = n(), .groups = "drop") %>%
  mutate(SV_Burden = SV / hg19_wgs_size_mb) %>% as.data.frame()
rownames(sv_tmb_summary) <- sv_tmb_summary$Sample

CL_Annots$SV <- sv_tmb_summary[CL_Annots$Sample, "SV"]
CL_Annots$SV_Burden <- sv_tmb_summary[CL_Annots$Sample, "SV_Burden"]

TMB_df <- read.delim(file.path(DIRS$icle$wes, "5_TMB", "TMB.tsv"))
rownames(TMB_df) <- TMB_df$Tumor_Sample_Barcode
CL_Annots$MUT <- TMB_df[CL_Annots$Name, "total"]
CL_Annots$TMB <- TMB_df[CL_Annots$Name, "total_perMB"]

alt_count_chr <- as.data.frame(matrixStats::rowMedians(as.matrix(table(subset(BRCA_CL_MAF, Center == "ICLE")[, c("Chromosome", "Tumor_Sample_Barcode")]))))
colnames(alt_count_chr) <- c("MUT")
alt_count_chr$Chr <- rownames(alt_count_chr)
alt_count_chr$Chr <- ifelse(!grepl("chr", alt_count_chr$Chr), paste0("chr", alt_count_chr$Chr), alt_count_chr$Chr)


alt_count_chr$SV <- as.data.frame(matrixStats::rowMedians(as.matrix(table(ICLE_SV$RefcontigID1, ICLE_SV$Sample))))[rownames(alt_count_chr), ]
alt_count_chr$SV_transloc <- as.data.frame(matrixStats::rowMedians(as.matrix(table(subset(ICLE_SV, Type %in% c("Tranloc-inter", "Tranloc-intra"))$RefcontigID1,
                                                                                   subset(ICLE_SV, Type %in% c("Tranloc-inter", "Tranloc-intra"))$Sample))))[rownames(alt_count_chr), ]
alt_count_chr$SV_notransloc <- as.data.frame(matrixStats::rowMedians(as.matrix(table(subset(ICLE_SV, Type %notin% c("Tranloc-inter", "Tranloc-intra"))$RefcontigID1,
                                                                                    subset(ICLE_SV, Type %notin% c("Tranloc-inter", "Tranloc-intra"))$Sample))))[rownames(alt_count_chr), ]

alt_count_chr <- alt_count_chr[!is.na(alt_count_chr$SV), ]

alt_count_chr$`Chr Size (MB)` <- round(chrom_sizes[rownames(alt_count_chr), "size_Mb"])
tmb_per_chr <- data.frame(TMB = alt_count_chr$MUT / alt_count_chr$`Chr Size (MB)`, row.names = alt_count_chr$Chr)
alt_count_chr$TMB <- tmb_per_chr[alt_count_chr$Chr, "TMB"]
alt_count_chr$SV_Burden <- alt_count_chr$SV / alt_count_chr$`Chr Size (MB)`
alt_count_chr$SV_noTrans_Burden <- alt_count_chr$SV_notransloc / alt_count_chr$`Chr Size (MB)`

SupFig8A <- ggplot(CL_Annots[grep("-I", CL_Annots$Name), ], aes(x = MUT, y = SV)) +
  geom_point(size = 6, alpha = 0.8, aes(color = `mRNA Subtypes`), show.legend = FALSE) +
  scale_color_manual(values = annot_cols$Subtypes) +
  ggrepel::geom_text_repel(aes(label = Sample), box.padding = 0, point.padding = 0, segment.color = "grey50", size = 4) +
  stat_cor(method = "pearson", label.x.npc = "center", label.y.npc = "bottom", size = 4, p.digits = 3) +
  labs(x = "Total Mutations", y = "Total Structural Variations") +
  theme_clean(20) + theme(panel.grid.major.x = element_line(color = "gray", linewidth = 0.1, size = 0.3, linetype = "dashed"),
                          panel.grid.major.y = element_line(color = "gray", linewidth = 0.1, size = 0.3, linetype = "dashed"))

SupFig8B <- ggplot(alt_count_chr, aes(x = SV_Burden, y = TMB)) +
  geom_point(size = 3, alpha = 0.8, show.legend = FALSE, color = "gray") +
  ggrepel::geom_text_repel(aes(label = gsub("chr", "", Chr)), box.padding = 0, point.padding = 0, segment.color = "grey50", size = 4) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "center", size = 4, p.digits = 3) +
  labs(x = "SV Burden", y = "Tumor Mutation Burder") +
  theme_clean(20) + theme(panel.grid.major.x = element_line(color = "gray", linewidth = 0.1, size = 0.3, linetype = "dashed"),
                          panel.grid.major.y = element_line(color = "gray", linewidth = 0.1, size = 0.3, linetype = "dashed"))

SupFig8B_no_outliers <- ggplot(subset(alt_count_chr, Chr %notin% c("chr19", "chr8")), aes(x = SV_Burden, y = TMB)) +
  geom_point(size = 3, alpha = 0.8, show.legend = FALSE, color = "gray") +
  ggrepel::geom_text_repel(aes(label = gsub("chr", "", Chr)), box.padding = 0, point.padding = 0, segment.color = "grey50", size = 4) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4, p.digits = 3) +
  labs(x = "SV Burden", y = "Tumor Mutation Burder") +
  theme_clean(20) + theme(panel.grid.major.x = element_line(color = "gray", linewidth = 0.1, size = 0.3, linetype = "dashed"),
                          panel.grid.major.y = element_line(color = "gray", linewidth = 0.1, size = 0.3, linetype = "dashed"))

plot_data <- alt_count_chr %>%
  select(Chr, `Chr Size (MB)`, SV_transloc, SV_notransloc) %>%
  tidyr::pivot_longer(cols = c(SV_transloc, SV_notransloc), names_to = "Variant_Type", values_to = "Count")

SupFig8C <- ggplot(plot_data, aes(y = reorder(Chr, Count), x = Count, fill = Variant_Type)) +
  geom_col() +
  labs(y = " ", x = "SV Count") + theme_clean(20) + theme(panel.grid.major.y = element_blank()) +
  scale_fill_manual(values = c("gray", "#4DB6AC"))

SupFig8D <- ggplot(alt_count_chr, aes(y = reorder(Chr, MUT), x = MUT)) +
  geom_col() +
  labs(y = " ", x = "Mutation Count") + theme_clean(20) + theme(panel.grid.major.y = element_blank())

assign("chrom_sizes", chrom_sizes, envir = .GlobalEnv)
assign("alt_count_chr", alt_count_chr, envir = .GlobalEnv)
assign("sv_tmb_summary", sv_tmb_summary, envir = .GlobalEnv)
assign("tmb_per_chr", tmb_per_chr, envir = .GlobalEnv)
assign("SupFig8A", SupFig8A, envir = .GlobalEnv)
assign("SupFig8B", SupFig8B, envir = .GlobalEnv)
assign("SupFig8B_no_outliers", SupFig8B_no_outliers, envir = .GlobalEnv)
assign("SupFig8C", SupFig8C, envir = .GlobalEnv)
assign("SupFig8D", SupFig8D, envir = .GlobalEnv)
message("  ✓ SupFig 8 prep complete (chrom_sizes, alt_count_chr, SupFig8A–D assigned).\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# load_all_icle_data(load_external = TRUE)
# source("Helper_Scripts/13_SupFig8_TMB_SV_Preparation.R")
#
# ensure_dir(DIRS$results_sub$ogm)
# ggsave(file.path(DIRS$results_sub$ogm, "SupFig8A_TMB_SV_Burden_Correlation.pdf"), SupFig8A, width = 5.5, height = 5)
# write.table(sv_tmb_summary, ...)
# etc.
