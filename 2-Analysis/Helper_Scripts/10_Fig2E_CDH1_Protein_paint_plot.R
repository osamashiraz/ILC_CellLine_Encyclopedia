# ==============================================================================
# Script 10: Figure 2E - CDH1 Protein Paint Query
# ==============================================================================
# Description: Prepares CDH1 mutation tables for Protein Paint visualization
#              (cell lines + TCGA). Formats mutation data for external
#              visualization tool.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - CL_Annots: Cell line annotations
#   - BRCA_CL_MAF: Cell line mutation data
#   - TCGA_Annots: TCGA tumor annotations
#   - TCGA_BRCA_MAF: TCGA mutation data
#
# Output:
#   - CDH1 mutation tables formatted for Protein Paint
#
# Author: Osama Shiraz Shah
# ==============================================================================

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
message("Figure 2E: CDH1 Protein Paint query")
message("========================================")

# Check if required data objects exist, load if needed
required_objects <- c("CL_Annots", "BRCA_CL_MAF", "TCGA_Annots", "TCGA_BRCA_MAF")
missing_objects <- required_objects[!sapply(required_objects, exists, envir = .GlobalEnv)]

if (length(missing_objects) > 0) {
  message("  Loading missing data objects: ", paste(missing_objects, collapse = ", "))
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = TRUE)
}

# protein paint mutation class dictionary
# mutation_class_dict <- list(
#   MISSENSE = "M",
#   EXON = "E",
#   FRAMESHIFT = "F",
#   NONSENSE = "N",
#   SILENT = "S",
#   PROTEINDEL = "D",
#   PROTEININS = "I",
#   PROTEINALTERING = "ProteinAltering",
#   SPLICE_REGION = "P",
#   SPLICE = "L",
#   INTRON = "Intron",
#   `Not tested` = "Blank",
#   Wildtype = "WT",
#   UTR_3 = "Utr3",
#   UTR_5 = "Utr5",
#   NONSTANDARD = "X",
#   NONCODING = "non-coding",
#   SNV = "snv",
#   MNV = "mnv",
#   `Sequence insertion` = "insertion",
#   `Sequence deletion` = "deletion"
# )

mutation_class_dict <- list(
  Truncating_ILC = "M",
  Missense_ILC = "E",
  Truncating_ILC_like = "F",
  Missense_ILC_like = "N",
  Truncating_NST = "D",
  Missense_NST = "I",
  Silent_NST = "X"
)



CL_ids <- unique(c(subset(CL_Annots, overlapWCCLE != "Y")$Name, grep("-I", CL_Annots$Name, value = T))) #  & `mRNA Subtypes` %notin% c("Basal")
CL_CDH1_DF <- subset(BRCA_CL_MAF, Tumor_Sample_Barcode %in% CL_ids & Hugo_Symbol == "CDH1")[, c("Tumor_Sample_Barcode", "Histology", "Chromosome", "Start_Position", "End_Position", "Broad_Variant_Classification","t_AF", "HGVSp_Short")]
CL_CDH1_DF <- subset(CL_CDH1_DF, Broad_Variant_Classification %notin% c("Other")) # remove non-truncating mutations

CL_CDH1_DF <- CL_CDH1_DF %>% group_by(Broad_Variant_Classification) %>% 
  mutate(Protein_Paint_Query = paste0(HGVSp_Short, ", ", Chromosome, ":", Start_Position, ", ", mutation_class_dict[paste0(Broad_Variant_Classification, "_", gsub("-","_", Histology))], ", ", Tumor_Sample_Barcode))

write.table(CL_CDH1_DF, file = file.path(out_dir, "Fig2E_CellLines_CDH1_Mutations_Protein_Paint_Query.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
message("  ✓ Saved Fig2E_CellLines_CDH1_Mutations_Protein_Paint_Query.tsv")
#


TCGA_ids <- subset(TCGA_Annots,  `Final Pathology` %in% c("NST", "ILC"))$Case.ID
TCGA_BRCA_MAF$Histology = TCGA_Annots[TCGA_BRCA_MAF$Tumor_Sample_Barcode, "Final Pathology"]
TCGA_BRCA_MAF$t_AF = TCGA_BRCA_MAF$t_alt_count/(TCGA_BRCA_MAF$t_alt_count+TCGA_BRCA_MAF$t_ref_count)

TCGA_CDH1_DF <- subset(TCGA_BRCA_MAF, Tumor_Sample_Barcode %in% TCGA_ids & Hugo_Symbol == "CDH1")[, c("Tumor_Sample_Barcode", "Histology", "Chromosome", "Start_Position", "End_Position", "Broad_Variant_Classification","t_AF", "HGVSp_Short")]
TCGA_CDH1_DF <- subset(TCGA_CDH1_DF, Broad_Variant_Classification %notin% c("Other", "Inframe")) # remove non-truncating mutations in MDAMB453

TCGA_CDH1_DF <- TCGA_CDH1_DF %>% group_by(Broad_Variant_Classification) %>% 
  mutate(Protein_Paint_Query = paste0(HGVSp_Short, ", chr", Chromosome, ":", Start_Position, ", ", mutation_class_dict[paste0(Broad_Variant_Classification, "_", gsub("-","_", Histology))], ", ", Tumor_Sample_Barcode))


write.table(TCGA_CDH1_DF, file = file.path(out_dir, "Fig2E_TCGA_CDH1_Mutations_Protein_Paint_Query.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
message("  ✓ Saved Fig2E_TCGA_CDH1_Mutations_Protein_Paint_Query.tsv\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/10_Fig2E_CDH1_Protein_paint_plot.R")
#
# load_all_icle_data(load_external = TRUE)
# # Script runs on source. Requires CL_Annots, BRCA_CL_MAF, TCGA_Annots, TCGA_BRCA_MAF.
#
# Outputs: CellLines_CDH1_Mutations_Protein_Paint_Query.tsv, TCGA_CDH1_Mutations_Protein_Paint_Query.tsv in DIRS$results_sub$cdh1