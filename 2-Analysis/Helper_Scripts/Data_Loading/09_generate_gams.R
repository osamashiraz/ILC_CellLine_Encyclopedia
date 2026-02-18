# ==============================================================================
# Generate Genomic Alteration Matrices (GAMs)
# ==============================================================================
# Builds TCGA, MSK, and ICLE GAMs from MAF + GISTIC. Used by 01_load_all_data.R.
#
# Author: Osama Shiraz Shah
# ==============================================================================

suppressPackageStartupMessages({
  library(maftools)
  library(dplyr)
  library(readxl)
  library(data.table)
  library(vroom)
})

if (!exists("DIRS", envir = .GlobalEnv)) {
  stop("config.R must be loaded before 09_generate_gams. Run source('config.R') first.")
}
if (!exists("broad_variant_classes", envir = .GlobalEnv)) {
  hf <- file.path(DIRS$scripts$helpers, "Helper_Functions.R")
  if (file.exists(hf)) source(hf) else stop("Helper_Functions.R required for broad_variant_classes")
}
ensure_dir(DIRS$results_sub$gam)

# Read and clean TCGA clinical annotations
load_TCGA_annotations <- function(filepath) {
  annot <- read_xlsx(filepath) %>% as.data.frame()
  rownames(annot) <- annot$Case.ID
  annot$Tumor_Sample_Barcode <- annot$Case.ID
  
  # Clean ER/HER2 values
  annot$`ER IHC`[!annot$`ER IHC` %in% c("Positive", "Negative")] <- NA
  annot$`HER2 IHC`[!annot$`HER2 IHC` %in% c("Positive", "Negative", "Equivocal")] <- NA
  
  # Clean Histology values
  annot$`Final Pathology`[annot$`Final Pathology` %in% "Mixed.IDC.ILC"] <- "mDLC"
  
  # Add simplified column names
  annot$ER <- annot$`ER IHC`
  annot$HER2 <- annot$`HER2 IHC`
  annot$Histology <- annot$`Final Pathology`
  return(annot)
}

# Generate mutation-based GAM from MAF
generate_mutation_GAM <- function(maf, isTCGA = FALSE) {
  mut_mat <- t(mutCountMatrix(read.maf(maf, isTCGA = isTCGA, rmFlags = TRUE), includeSyn = FALSE, removeNonMutated = FALSE))
  mut_mat[mut_mat > 1] <- 1
  
  mut_mat <- apply(mut_mat, c(1, 2), function(x) ifelse(x == 1, "MUT", ""))
  mut_df <- as.data.frame(mut_mat)
  
  return(mut_df)
}

# Generate CNV-based GAM from GISTIC file
generate_cnv_GAM <- function(cnv_df) {
  
  cnv_df <- cnv_df[!duplicated(cnv_df$Hugo_Symbol), ]
  
  rownames(cnv_df) <- cnv_df$Hugo_Symbol
  cnv_df <- cnv_df[, -(1:2)]
    
  cnv_mat <- t(as.matrix(cnv_df))
  
  cnv_mat <- apply(cnv_mat, c(1,2), function(x) {
    switch(as.character(x),
           "2" = "AMP", "1" = "GAIN",
           "0" = "", "-1" = "LOH",
           "-2" = "DEL", "")
  })
  
  return(cnv_mat)
}

# Merge CNV + SNV GAMs
merge_GAMs <- function(mutation_gam, cnv_gam) {
  cnv_gam = as.data.frame(cnv_gam)
  mutation_gam = as.data.frame(mutation_gam)
  
  common_samples <- intersect(rownames(mutation_gam), rownames(cnv_gam))
  features <- unique(c(colnames(mutation_gam), colnames(cnv_gam)))
  
  reindex_columns <- function(df, required_cols) {
    missing <- setdiff(required_cols, colnames(df))
    df[missing] <- ""
    df[, required_cols, drop = FALSE]
  }

  # Fill missing features in either GAM with empty strings
  temp1 = reindex_columns(cnv_gam[common_samples, , drop = FALSE], features)
  temp2 <- reindex_columns(mutation_gam[common_samples, , drop = FALSE], features)
  temp2[is.na(temp2)] <- ""
  temp1[is.na(temp1)] <- ""
  
  temp1$ID <- rownames(temp1)
  temp2$ID <- rownames(temp2)
  temp2 <- temp2[rownames(temp1), ]
  
  merged <- bind_rows(temp2, temp1) %>%
    group_by(ID) %>%
    summarise(across(everything(), ~ paste0(.x, collapse = ";")), .groups = "drop") %>%
    as.data.frame()
  
  rownames(merged) <- merged$ID
  merged <- merged[,-1] %>% t() %>% as.matrix()
  
  # Clean-up artifacts like "MUT;", ";DEL", etc.
  merged[merged == "MUT;"] <- "MUT"
  merged[merged == ";DEL"] <- "DEL"
  merged[merged == ";AMP"] <- "AMP"
  merged[merged == ";GAIN"] <- "GAIN"
  merged[merged == ";LOH"] <- "LOH"
  merged[merged == ";"] <- ""
  
  return(merged)
}


# ----------------------------
# Run Pipeline TCGA
# ----------------------------
if (exists("TCGA_BRCA_GAM", envir = .GlobalEnv) && exists("TCGA_BRCA_MAF", envir = .GlobalEnv)) {
  message("TCGA GAM data already loaded (skipping)")
} else if (!file.exists(FILES$tcga_gam)) {
  message("Processing TCGA GAMs...")
  TCGA_BRCA_MAF <- read.delim(FILES$tcga_maf)
  TCGA_BRCA_MAF$Tumor_Sample_Barcode <- substr(TCGA_BRCA_MAF$Tumor_Sample_Barcode, 1, 15)
  TCGA_BRCA_MAF$Broad_Variant_Classification <- broad_variant_classes[TCGA_BRCA_MAF$Variant_Classification]
  TCGA_BRCA_MAF$t_AF <- TCGA_BRCA_MAF$t_alt_count / (TCGA_BRCA_MAF$t_alt_count + TCGA_BRCA_MAF$t_ref_count)
  save(TCGA_BRCA_MAF, file = FILES$tcga_maf_rdata)
  TCGA_BRCA_GISTIC <- vroom(FILES$tcga_gistic, show_col_types = FALSE) %>% as.data.frame()
  colnames(TCGA_BRCA_GISTIC) <- substr(colnames(TCGA_BRCA_GISTIC), 1, 15)
  gam_mut <- generate_mutation_GAM(subset(TCGA_BRCA_MAF, Broad_Variant_Classification %notin% c("Other", "Silent")))
  gam_cnv <- generate_cnv_GAM(TCGA_BRCA_GISTIC)
  TCGA_BRCA_GAM <- merge_GAMs(gam_mut, gam_cnv)
  save(TCGA_BRCA_GAM, file = FILES$tcga_gam)
  assign("TCGA_BRCA_GAM", TCGA_BRCA_GAM, envir = .GlobalEnv)
  assign("TCGA_BRCA_MAF", TCGA_BRCA_MAF, envir = .GlobalEnv)
  message("  ✓ TCGA GAM saved: ", FILES$tcga_gam)
} else {
  message("Loading TCGA GAM...")
  load(FILES$tcga_gam, envir = .GlobalEnv)
  load(FILES$tcga_maf_rdata, envir = .GlobalEnv)
}




# ----------------------------
# Run Pipeline MSK
# ----------------------------
if (exists("MSK_BRCA_GAM", envir = .GlobalEnv) && exists("MSK_BRCA_MAF", envir = .GlobalEnv)) {
  message("MSK GAM data already loaded (skipping)")
} else if (!file.exists(FILES$msk_gam)) {
  message("Processing MSK GAMs...")
  MSK_BRCA_MAF <- read.delim(FILES$msk_maf)
  MSK_BRCA_MAF$Broad_Variant_Classification <- broad_variant_classes[MSK_BRCA_MAF$Variant_Classification]
  MSK_BRCA_MAF$t_AF <- MSK_BRCA_MAF$t_alt_count / (MSK_BRCA_MAF$t_alt_count + MSK_BRCA_MAF$t_ref_count)
  save(MSK_BRCA_MAF, file = FILES$msk_maf_rdata)
  MSK_BRCA_GISTIC <- vroom(FILES$msk_gistic, show_col_types = FALSE) %>% as.data.frame()
  gam_mut <- generate_mutation_GAM(subset(MSK_BRCA_MAF, Broad_Variant_Classification %notin% c("Other", "Silent")))
  gam_cnv <- generate_cnv_GAM(MSK_BRCA_GISTIC)
  MSK_BRCA_GAM <- merge_GAMs(gam_mut, gam_cnv)
  save(MSK_BRCA_GAM, file = FILES$msk_gam)
  assign("MSK_BRCA_GAM", MSK_BRCA_GAM, envir = .GlobalEnv)
  assign("MSK_BRCA_MAF", MSK_BRCA_MAF, envir = .GlobalEnv)
  message("  ✓ MSK GAM saved: ", FILES$msk_gam)
} else {
  message("Loading MSK GAM...")
  load(FILES$msk_gam, envir = .GlobalEnv)
  load(FILES$msk_maf_rdata, envir = .GlobalEnv)
}




# ----------------------------
# Run Pipeline BRCA Cell
# ----------------------------
if (exists("BRCA_CL_GAM", envir = .GlobalEnv) && exists("BRCA_CL_MAF", envir = .GlobalEnv)) {
  message("Cell line GAM data already loaded (skipping)")
} else if (!file.exists(FILES$cl_gam)) {
  message("Processing ICLE+CCLE+MARCOTTE GAMs...")
  BRCA_CL_MAF <- read.delim(FILES$cl_maf)
  save(BRCA_CL_MAF, file = FILES$cl_maf_rdata)
  cnv <- vroom(FILES$cl_gistic, show_col_types = FALSE) %>% as.data.frame()
  colnames(cnv)[c(1, 2)] <- c("Hugo_Symbol", "Entrez_Gene_Id")
  cnv <- cnv[, -3]
  colnames(cnv) <- gsub("-M", "-C", colnames(cnv))
  colnames(cnv) <- gsub("MM134-I", "MDAMB134VI-I", colnames(cnv))
  colnames(cnv) <- gsub("^MM", "MDAMB", colnames(cnv))
  colnames(cnv) <- gsub("MPE600", "600MPE", colnames(cnv))
  colnames(cnv) <- gsub("SUM44-C", "SUM44PE-C", colnames(cnv))
  skbr3_i <- data.frame("SKBR3-I" = cnv$`SKBR3-C`, check.names = FALSE)
  colnames(skbr3_i) <- gsub("[.]", "-", colnames(skbr3_i))
  cnv <- cbind(cnv, skbr3_i)
  gam_mut <- generate_mutation_GAM(subset(BRCA_CL_MAF, Broad_Variant_Classification %notin% c("Other", "Silent")))
  gam_cnv <- generate_cnv_GAM(cnv)
  BRCA_CL_GAM <- merge_GAMs(gam_mut, gam_cnv)
  save(BRCA_CL_GAM, file = FILES$cl_gam)
  assign("BRCA_CL_GAM", BRCA_CL_GAM, envir = .GlobalEnv)
  assign("BRCA_CL_MAF", BRCA_CL_MAF, envir = .GlobalEnv)
  message("  ✓ Cell line GAM saved: ", FILES$cl_gam)
} else {
  message("Loading cell line GAM...")
  load(FILES$cl_gam, envir = .GlobalEnv)
  load(FILES$cl_maf_rdata, envir = .GlobalEnv)
}


# Add manual SV/Exonic deletions if objects exist
if (exists("BRCA_CL_GAM", envir = .GlobalEnv)) {
  BRCA_CL_GAM <- get("BRCA_CL_GAM", envir = .GlobalEnv)
  BRCA_CL_GAM["CDH1", c("MDAMB134VI-I", "OCUBM-I")] <- "DEL"  # Exonic deletion
  BRCA_CL_GAM["CDH1", c("HCC2185-I", "HCC2218-I", "SKBR3-I")] <- "DEL"  # Large deletion
  BRCA_CL_GAM["CTNNA1", "MDAMB468-I"] <- "DEL"  # Exonic deletion
}

message("  ✓ Manual GAM annotations applied")


### Usuage
# genes = c("CDH1","TP53","PTEN","ERBB2","PIK3CA","GATA3","CTNNA1","TBX3", "FOXA1")
# setdiff(genes, rownames(BRCA_CL_GAM))
# 
# mat = BRCA_CL_GAM[genes,]
# 
# for(rm_alt in c("LOH","GAIN")){
#   mat[mat %in% rm_alt] = ""
# }
# 
# 
# library(ComplexHeatmap)
# col = c("GAIN"= "#FFCCBC", "LOH" = "#9FA8DA",  "AMP" = "#d50000", "DEL"= "#2962FF", "MUT" = "#008000")
# 
# alter_fun = list(
#   background = alter_graphic("rect", width = 1, fill = "#E0E0E0"),   ##E0E0E0
#   LOH = alter_graphic("rect", width = 1.01, height = 1, fill = col["LOH"]),
#   AMP = alter_graphic("rect",  width = 1.01, height = 1, fill = col["AMP"]),
#   GAIN = alter_graphic("rect",  width = 1.01, height = 1, fill = col["GAIN"]),
#   DEL = alter_graphic("rect",  width = 1.01, height = 1, fill = col["DEL"]),
#   MUT = alter_graphic("rect",  width = 1.01, height = 0.33, fill = col["MUT"])
# )
# 
# 
# 
# oncoPrint(mat, alter_fun = alter_fun, col = col, alter_fun_is_vectorized = FALSE,
#           heatmap_legend_param = list(column_gap = unit('0.1','mm')), show_column_names = T,
#           column_split = CL_Annots[match(gsub("-I|-C", "", colnames(mat)), CL_Annots$Sample), "Histology"],
#           bottom_annotation = HeatmapAnnotation(Histology = CL_Annots[match(gsub("-I|-C", "", colnames(mat)), CL_Annots$Sample), "Histology"],
#                                                 PAM50 = CL_Annots[colnames(mat), "PAM50"],
#                                                 col = annot_cols))















