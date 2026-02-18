# ==============================================================================
# Script 13: Figure 2H - CDH1 Alteration Landscape (TCGA + Cell Lines)
# ==============================================================================
# Description: Generates multi-omic CDH1 heatmaps showing mutation, copy number,
#              RNA expression, protein expression, and DNA methylation patterns
#              across TCGA tumors and cell lines.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - TCGA_Annots: TCGA tumor annotations
#   - TCGA_BRCA_MAF: TCGA mutation data
#   - TCGA_BRCA_GAM: TCGA genomic alteration matrix
#   - BRCA_CL_GAM: Cell line genomic alteration matrix
#   - BRCA_CL_EXP: RNA expression data
#   - BRCA_CL_RPPA: Protein expression data
#   - BRCA_CL_DNAm: DNA methylation data
#   - CL_Annots: Cell line annotations
#
# Output:
#   - Figure 2H multi-omic heatmap
#
# Author: Osama Shiraz Shah
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ComplexHeatmap) })
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
message("Figure 2H: CDH1 alteration landscape")
message("========================================")

# Check if required data objects exist, load if needed
required_objects <- c("TCGA_Annots", "TCGA_BRCA_MAF", "TCGA_BRCA_GAM", "BRCA_CL_GAM", "BRCA_CL_EXP", 
                      "BRCA_CL_RPPA", "BRCA_CL_DNAm", "CL_Annots")
missing_objects <- required_objects[!sapply(required_objects, exists, envir = .GlobalEnv)]

if (length(missing_objects) > 0) {
  message("  Loading missing data objects: ", paste(missing_objects, collapse = ", "))
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = TRUE)
}

ht_opt(
  heatmap_row_names_gp = gpar(fontfamily = "Helvetica"),
  heatmap_column_names_gp = gpar(fontfamily = "Helvetica"),
  heatmap_row_title_gp = gpar(fontfamily = "Helvetica"),
  heatmap_column_title_gp = gpar(fontfamily = "Helvetica"),
  legend_title_gp = gpar(fontfamily = "Helvetica"),
  legend_labels_gp = gpar(fontfamily = "Helvetica")
)

TCGA_ids <- subset(TCGA_Annots, PAM50 %notin% "Basal" & `Final Pathology` %in% c("NST", "ILC"))$Case.ID

AF_mat <- subset(TCGA_BRCA_MAF, Hugo_Symbol == "CDH1" & Tumor_Sample_Barcode %in% TCGA_ids)[,c("Tumor_Sample_Barcode", "Broad_Variant_Classification", "Variant_Classification", "t_AF")] %>% 
  filter(Broad_Variant_Classification != "Other") %>%
  mutate(severity_score = case_when(Broad_Variant_Classification == "Truncating" ~ 1,Broad_Variant_Classification == "Missense"  ~ 2,TRUE ~ 3)) %>% 
  group_by(Tumor_Sample_Barcode) %>% arrange(severity_score, .by_group = TRUE) %>% dplyr::slice(1) %>% ungroup() %>% select(-severity_score) %>% as.data.frame()
rownames(AF_mat) <- AF_mat$Tumor_Sample_Barcode
AF_mat <- AF_mat[TCGA_ids, ]
rownames(AF_mat) <- TCGA_ids
AF_mat$Tumor_Sample_Barcode <- TCGA_ids
AF_mat[is.na(AF_mat$Broad_Variant_Classification), "Broad_Variant_Classification"] <- "WT"
AF_mat[is.na(AF_mat$Variant_Classification), "Variant_Classification"] <- "WT"


TCGA_CDH1_df <- data.frame(ID = TCGA_ids, Histology = TCGA_Annots[TCGA_ids,"Final Pathology"], PAM50 = TCGA_Annots[TCGA_ids, "PAM50"],
                           CN = TCGA_BRCA_GAM["CDH1", match(TCGA_ids, colnames(TCGA_BRCA_GAM))],
                           CN_logRR = as.numeric(TCGA_BRCA_CN["CDH1", match(TCGA_ids, colnames(TCGA_BRCA_CN))]),
                           Mutation = TCGA_BRCA_GAM["CDH1", match(TCGA_ids, colnames(TCGA_BRCA_GAM))],
                           Mutation_AF = AF_mat[TCGA_ids,"t_AF"], 
                           Mutation_type_broad = AF_mat[TCGA_ids,"Broad_Variant_Classification"], 
                           Mutation_type = AF_mat[TCGA_ids,"Variant_Classification"], 
                           RNA = scale(TCGA_BRCA_Log2CPM["CDH1",])[match(TCGA_ids, colnames(TCGA_BRCA_Log2CPM)),], 
                           Protein = scale(as.matrix(TCGA_BRCA_RPPA)["E-CADHERIN",])[match(TCGA_ids, colnames(TCGA_BRCA_RPPA)),],
                           DNAm = scale(TCGA_BRCA_DNAm[c("cg01251360"), ])[match(TCGA_ids, colnames(TCGA_BRCA_DNAm)), ]) %>% 
                           mutate(CN = case_when(grepl("AMP", CN) ~ 2, grepl("GAIN", CN) ~ 1, grepl("LOH", CN) ~ -1, grepl("DEL", CN) ~ -2, is.na(CN) ~ NA, grepl("MUT|", CN) ~ 0), 
                                  Mutation = ifelse(grepl("MUT", Mutation), 1, 0)) %>% 
                           filter(Histology %in% c("ILC", "NST"))
col_order <- order(TCGA_CDH1_df$RNA)

TCGA_MUT_ht <- Heatmap(t(TCGA_CDH1_df$Mutation), name = "Mutations", col = c("1"="black", "0"="white"), column_split = TCGA_CDH1_df$Histology, column_order = col_order,
                       border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)
TCGA_CN_ht <- Heatmap(t(TCGA_CDH1_df$CN), name = "Copy Number", col = annot_cols$CN, column_split = TCGA_CDH1_df$Histology, column_order = col_order, 
                      border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)
TCGA_RNA_ht <- Heatmap(t(TCGA_CDH1_df$RNA), name = "RNA", col = annot_cols$RNA_zscore, column_split = TCGA_CDH1_df$Histology, column_order = col_order, 
                       border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)
TCGA_RPPA_ht <- Heatmap(t(TCGA_CDH1_df$Protein), name = "Protein", col = annot_cols$RPPA, column_split = TCGA_CDH1_df$Histology, column_order = col_order, 
                        border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)
# TCGA_DNAm_ht <- Heatmap(t(TCGA_CDH1_df$DNAm), name = "DNAm", col = annot_cols$DNAm_zscore, column_split = TCGA_CDH1_df$Histology, column_order = col_order,
#                         border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)

fig2h_tcga <- TCGA_MUT_ht %v% TCGA_CN_ht %v% TCGA_RNA_ht %v% TCGA_RPPA_ht
assign("fig2h_tcga", fig2h_tcga, envir = .GlobalEnv)
assign("TCGA_CDH1_df", TCGA_CDH1_df, envir = .GlobalEnv)





CL_ids <- unique(c(subset(CL_Annots, overlapWCCLE != "Y" & `mRNA Subtypes` %notin% c("Basal"))$Name, grep("-I", CL_Annots$Name, value = T)))

BRCA_CL_CN[, "SKBR3-I"] = BRCA_CL_CN[, "SKBR3-C"] # borrow information from marcotte for missing cell line data


CL_AF_mat <- subset(BRCA_CL_MAF, Hugo_Symbol == "CDH1" & Tumor_Sample_Barcode %in% CL_ids)[,c("Tumor_Sample_Barcode", "Broad_Variant_Classification", "Variant_Classification", "t_AF")] %>%
  filter(Broad_Variant_Classification != "Other") %>%
  mutate(severity_score = case_when(Broad_Variant_Classification == "Truncating" ~ 1,Broad_Variant_Classification == "Missense"  ~ 2,TRUE ~ 3)) %>% 
  group_by(Tumor_Sample_Barcode) %>% arrange(severity_score, .by_group = TRUE) %>% slice(1) %>% ungroup() %>% select(-severity_score) %>% as.data.frame()
rownames(CL_AF_mat) <- CL_AF_mat$Tumor_Sample_Barcode
CL_AF_mat <- CL_AF_mat[CL_ids, ]
rownames(CL_AF_mat) <- CL_ids
CL_AF_mat$Tumor_Sample_Barcode <- CL_ids
CL_AF_mat[is.na(CL_AF_mat$Broad_Variant_Classification), "Broad_Variant_Classification"] <- "WT"
CL_AF_mat[is.na(CL_AF_mat$Variant_Classification), "Variant_Classification"] <- "WT"


CL_CDH1_df <- data.frame(ID = CL_ids, Histology = CL_Annots[CL_ids, "Histology"], Subtype = CL_Annots[CL_ids, "mRNA Subtypes"],
                         CN = BRCA_CL_GAM["CDH1", match(CL_ids, colnames(BRCA_CL_GAM))],
                         CN_logRR = as.numeric(as.matrix(BRCA_CL_CN)["CDH1", match(CL_ids, colnames(BRCA_CL_CN))]),
                         Mutation = BRCA_CL_GAM["CDH1", match(CL_ids, colnames(BRCA_CL_GAM))],
                         Mutation_AF = CL_AF_mat[CL_ids, "t_AF"],
                         Mutation_type_broad = CL_AF_mat[CL_ids,"Broad_Variant_Classification"], 
                         Mutation_type = CL_AF_mat[CL_ids,"Variant_Classification"], 
                         RNA = scale(BRCA_CL_EXP["CDH1",])[match(CL_ids, colnames(BRCA_CL_EXP)),], 
                         Protein = scale(as.matrix(BRCA_CL_RPPA)["E-CADHERIN",])[match(CL_ids, colnames(BRCA_CL_RPPA)),],
                         DNAm = scale(BRCA_CL_DNAm[c("cg01251360"), ])[match(CL_ids, colnames(BRCA_CL_DNAm)), ]) %>% 
  mutate(CN = case_when(grepl("AMP", CN) ~ 2, grepl("GAIN", CN) ~ 1, grepl("LOH", CN) ~ -1, grepl("DEL", CN) ~ -2, is.na(CN) ~ NA, grepl("MUT|", CN) ~ 0), 
         Mutation = ifelse(grepl("MUT", Mutation), 1, 0))

col_order <- order(CL_CDH1_df$RNA)

CL_MUT_ht <- Heatmap(t(CL_CDH1_df$Mutation), name = "Mutations", col = c("1"="black", "0"="white"), column_split = CL_CDH1_df$Histology, column_order = col_order, 
                       border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)
CL_CN_ht <- Heatmap(t(CL_CDH1_df$CN), name = "Copy Number", col = annot_cols$CN, column_split = CL_CDH1_df$Histology, column_order = col_order,
                      border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)
CL_RNA_ht <- Heatmap(t(CL_CDH1_df$RNA), name = "RNA", col = annot_cols$RNA_zscore, column_split = CL_CDH1_df$Histology, column_order = col_order,
                       border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)
CL_RPPA_ht <- Heatmap(t(CL_CDH1_df$Protein), name = "Protein", col = annot_cols$RPPA, column_split = CL_CDH1_df$Histology, column_order = col_order,
                        border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)
# CL_DNAm_ht <- Heatmap(t(CL_CDH1_df$DNAm), name = "DNAm", col = annot_cols$DNAm_zscore, column_split = CL_CDH1_df$Histology, column_order = col_order,
#                         border = T, height = unit(1, "cm"), width = unit(10, "cm"), heatmap_legend_param = heatmap_legend_param)

fig2h_cl <- CL_MUT_ht %v% CL_CN_ht %v% CL_RNA_ht %v% CL_RPPA_ht
assign("fig2h_cl", fig2h_cl, envir = .GlobalEnv)
assign("CL_CDH1_df", CL_CDH1_df, envir = .GlobalEnv)
message("  ✓ Fig 2H complete (fig2h_tcga, fig2h_cl, TCGA_CDH1_df, CL_CDH1_df assigned).\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/12_Fig2H_CDH1_Alteration_Landscape.R")
#
# load_all_icle_data(load_external = TRUE)
# # Script runs on source. Requires TCGA/CL GAM, MAF, RNA, RPPA, DNAm, annots.
#
# Outputs: Fig2H_CellLine_CDH1_Alteration_Heatmap.pdf in DIRS$results_sub$cdh1
