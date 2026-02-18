# ==============================================================================
# Script 04: Figure 1D - Multiomics Overview Heatmap
# ==============================================================================
# Description: Creates a comprehensive heatmap showing multi-omic profiling of
#              ICLE cell lines including:
#              - Key genetic alterations (SNV, CNV)
#              - Genomic instability metrics (FGA, TMB, SV count, Fusions, DMI)
#              - Receptor protein expression (ER, PR, HER2)
#              - Top variable proteins (RPPA)
#              - Top variable gene expression (RNA-seq)
#              - Top variable DNA methylation patterns
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - CL_Annots: Cell line annotations with subtype information
#   - BRCA_CL_GAM: Gene alteration matrix
#   - BRCA_CL_RPPA: Protein expression data
#   - BRCA_CL_EXP: RNA expression data
#   - BRCA_CL_DNAm: DNA methylation data
#
# Output:
#   - fig1d_multiomics_overview: Multiomics overview heatmap (assigned to .GlobalEnv when script is run)
#   - Fig1D_Multiome_Overview.pdf written to DIRS$results
#
# Author: Osama Shiraz Shah
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(reshape2)
  library(ComplexHeatmap)
  library(CancerSubtypes)
  library(circlize)
  library(grid)
})

# Load configuration and helper functions
if (!exists("DIRS")) source("config.R")
if (!exists("annot_cols")) source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))

# Load main data objects (if not already loaded)
if (!exists("BRCA_CL_GAM")) {
  message("  Loading data objects...")
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data()
}

message("\n========================================")
message("Figure 1D: Multiomics Overview")
message("========================================\n")

# ==============================================================================
# 1. Load Additional Genomic Metrics
# ==============================================================================

load_additional_metrics <- function() {
  
  # DNA Methylation Instability (DMI)
  DMI <- subset(
    read.delim(file.path(DIRS$icle$dnam, "3_Integration", "DMI_Summary.tsv")),
    TissueSource %in% "Cell Lines"
  )
  rownames(DMI) <- DMI$SampleID
  
  # Tumor Mutational Burden (TMB)
  TMB <- read.delim(file.path(DIRS$icle$wes, "5_TMB", "TMB.tsv"))
  rownames(TMB) <- TMB$Tumor_Sample_Barcode
  
  # Fraction Genome Altered (FGA)
  FGA <- read.delim(file.path(DIRS$icle$cytosnp, "6_CIN_Analysis", "CINMetric_FGA.tsv"))
  FGA$sample_id <- gsub("MM", "MDAMB", FGA$sample_id)
  FGA$sample_id <- gsub("MDAMB134-I", "MDAMB134VI-I", FGA$sample_id)
  FGA$sample_id <- gsub("MPE600", "600MPE", FGA$sample_id)
  FGA$sample_id <- gsub("-M", "-C", FGA$sample_id)
  rownames(FGA) <- FGA$sample_id
  FGA["SKBR3-I", "fga"] <- FGA["SKBR3-C", "fga"] # borrow FGA for SKBR3-I from SKBR3-C

  # Structural Variations
  SV_data <- read.csv(file.path(DIRS$icle$bionano, "2_Structural_Variations", "ICLE_SV_Filtered.csv"))
  SV_count <- table(SV_data$Sample) %>% as.data.frame()
  colnames(SV_count) <- c("Sample", "Count")
  rownames(SV_count) <- paste0(SV_count$Sample, "-I")
  
  # Fusion events
  Fusions <- read.delim(file.path(DIRS$icle$bionano, "3_Fusions", "fusion_df.tsv"))
  Fusions_count <- table(Fusions$Sample) %>% as.data.frame()
  colnames(Fusions_count) <- c("Sample", "Count")
  rownames(Fusions_count) <- paste0(Fusions_count$Sample, "-I")

  metrics <- list(
    DMI = DMI,
    TMB = TMB,
    FGA = FGA,
    SV_count = SV_count,
    Fusions_count = Fusions_count
  )

  # Add metrics to annotations
  CL_Annots$FGA <- metrics$FGA[rownames(CL_Annots), "fga"]
  CL_Annots$SV <- metrics$SV_count[rownames(CL_Annots), "Count"]
  CL_Annots$Fusions <- metrics$Fusions_count[rownames(CL_Annots), "Count"]
  CL_Annots$DMI <- metrics$DMI[rownames(CL_Annots), "DMI_zscore"]
  CL_Annots$TMB <- metrics$TMB[rownames(CL_Annots), "total_perMB"]
  assign("CL_Annots", CL_Annots, envir = .GlobalEnv)

  write.table(CL_Annots[,c("Name", "Sample", "Study", "overlapWCCLE", "FGA", "SV", "Fusions", "DMI", "TMB")], 
              file = file.path(DIRS$results, "Genomic_Instability_Metrics.tsv"), sep = "\t", row.names = F)
  message("  ✓ Saved Genomic Instability Metrics to File", file.path(DIRS$results, "Genomic_Instability_Metrics.tsv"))
  
  # Add metrics to annotations
  return(metrics)
}

# ==============================================================================
# 2. Prepare Alteration Matrix
# ==============================================================================

prepare_alteration_matrix <- function() {
  # ICLE cell lines (suffix -I)
  ICLE_cells <- grep('-I', CL_Annots$Name, value = TRUE)

  # Key genes for display
  genes <- c("CDH1", "TP53", "ERBB2", "PIK3CA", "CTNNA1",
             "PTEN", "TBX3", "FOXA1", "BRCA2")

  # Extract alteration matrix
  alt_mat <- BRCA_CL_GAM[
    intersect(genes, rownames(BRCA_CL_GAM)),
    ICLE_cells
  ]
  
  # Simplify alteration types
  alt_mat[alt_mat == "LOH"] <- ""
  alt_mat[alt_mat == "GAIN"] <- ""
  alt_mat[alt_mat == ""] <- "WT"
  alt_mat[alt_mat == "AMP"] <- "AMP"
  alt_mat[alt_mat == "DEL"] <- "DEL"
  alt_mat[grepl("MUT", alt_mat)] <- "MUT"
  
  return(list(alt_mat = alt_mat, ICLE_cells = ICLE_cells))
}

# ==============================================================================
# Helper Function: Load DNAm ProbeSet
# ==============================================================================
# Load and simplify HM450K ProbeSet annotation
load_DNAm_probeset <- function(probe_file) {
  # Check if HM450K_ProbeSet already exists in global environment
  if (!exists("HM450K_ProbeSet", envir = .GlobalEnv)) {
    load(probe_file)  # loads `HM450K_ProbeSet`
  } else {
    HM450K_ProbeSet <- get("HM450K_ProbeSet", envir = .GlobalEnv)
  }
  
  # Copy the original region column
  HM450K_ProbeSet$region_simple <- HM450K_ProbeSet$region
  
  # Simplify multi-region annotations
  HM450K_ProbeSet$region_simple[
    HM450K_ProbeSet$region %in% c(
      "exon-fiveUTRs-threeUTRs-TSS_1k-TES_1k", "exon-fiveUTRs-threeUTRs-TSS_1k", 
      "exon-fiveUTRs-threeUTRs-TES_1k", "exon-fiveUTRs-threeUTRs", 
      "exon-threeUTRs-TSS_1k-TES_1k", "exon-threeUTRs-TES_1k", 
      "exon-fiveUTRs-TSS_1k-TES_1k", "exon-TSS_1k-TES_1k", 
      "exon-threeUTRs-TSS_1k", "TSS_1k-TES_1k", 
      "exon-fiveUTRs-TES_1k", "exon-fiveUTRs-TSS_1k"
    )
  ] <- "multiple"
  
  # Map specific region combinations to individual labels
  HM450K_ProbeSet$region_simple[HM450K_ProbeSet$region == "exon-TES_1k"] <- "TES_1k"
  HM450K_ProbeSet$region_simple[HM450K_ProbeSet$region == "exon-threeUTRs"] <- "threeUTRs"
  HM450K_ProbeSet$region_simple[HM450K_ProbeSet$region == "exon-TSS_1k"] <- "TSS_1k"
  HM450K_ProbeSet$region_simple[HM450K_ProbeSet$region == "exon-fiveUTRs"] <- "fiveUTRs"
  
  return(HM450K_ProbeSet)
}

# ==============================================================================
# 3. Select Top Variable Features
# ==============================================================================

select_top_features <- function(ICLE_cells) {
  
  # RNA: Top 6000 by MAD (with conditional plot suppression)
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    select_features_rna <- rownames(
      quiet_run(CancerSubtypes::FSbyMAD(
        as.matrix(BRCA_CL_EXP[, ICLE_cells]),
        cut.type = "topk",
        value = 6000
      ))
    )
  } else {
    select_features_rna <- rownames(
      CancerSubtypes::FSbyMAD(
        as.matrix(BRCA_CL_EXP[, ICLE_cells]),
        cut.type = "topk",
        value = 6000
      )
    )
  }
  mrna_mat <- t(scale(t(BRCA_CL_EXP[select_features_rna, ICLE_cells])))
  message("  ✓ RNA: ", nrow(mrna_mat), " genes selected")
  
  # RPPA: Top 50 by variance (with conditional plot suppression)
  icle_rppa_samples <- intersect(ICLE_cells, colnames(BRCA_CL_RPPA))
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    select_features_RPPA <- rownames(
      quiet_run(CancerSubtypes::FSbyVar(
        as.matrix(BRCA_CL_RPPA[, icle_rppa_samples]),
        cut.type = "topk",
        value = 50
      ))
    )
  } else {
    select_features_RPPA <- rownames(
      CancerSubtypes::FSbyVar(
        as.matrix(BRCA_CL_RPPA[, icle_rppa_samples]),
        cut.type = "topk",
        value = 50
      )
    )
  }
  rppa_mat <- t(t(BRCA_CL_RPPA[
    select_features_RPPA,
    match(ICLE_cells, colnames(BRCA_CL_RPPA))
  ]))
  colnames(rppa_mat)[is.na(colnames(rppa_mat))] <- c("WCRC25-I", "ZR7530-I")
  message("  ✓ RPPA: ", nrow(rppa_mat), " proteins selected")
  
  # # Receptor RPPA
  # receptor_mat <- t(t(BRCA_CL_RPPA[
  #   c("ER-A", "PR", "HER2", "HER2-PY1248"),
  #   match(ICLE_cells, colnames(BRCA_CL_RPPA))
  # ]))
  # colnames(receptor_mat)[is.na(colnames(receptor_mat))] <- c("WCRC25-I", "ZR7530-I")
  # message("  ✓ Receptor proteins: ", nrow(receptor_mat), " markers")
  
  # DNAm: Top 6000 promoter probes by variance
  HM450K_ProbeSet <- load_DNAm_probeset(
    file.path(DIRS$external$misc, "HM450K_ProbeSet_Unique_Gene_Promoter_Mapping.Rdata")
  )
  upstream_probes <- subset(
    HM450K_ProbeSet,
    promoter_refTSS_canonical == TRUE &
    region_simple %in% c("TSS_1k", 'fiveUTRs', "exon")
  )$probeID
  
  icle_dnam_samples <- intersect(ICLE_cells, colnames(BRCA_CL_DNAm))
  # DNAm: Top 6000 by variance (with conditional plot suppression)
  if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
      get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
    select_features_dnam <- rownames(
      quiet_run(CancerSubtypes::FSbyVar(
        as.matrix(BRCA_CL_DNAm[
          intersect(rownames(BRCA_CL_DNAm), upstream_probes),
          icle_dnam_samples
        ]),
        cut.type = "topk",
        value = 6000
      ))
    )
  } else {
    select_features_dnam <- rownames(
      CancerSubtypes::FSbyVar(
        as.matrix(BRCA_CL_DNAm[
          intersect(rownames(BRCA_CL_DNAm), upstream_probes),
          icle_dnam_samples
        ]),
        cut.type = "topk",
        value = 6000
      )
    )
  }
  
  dna_mat <- t(scale(t(BRCA_CL_DNAm[
    select_features_dnam,
    match(ICLE_cells, colnames(BRCA_CL_DNAm))
  ])))
  colnames(dna_mat)[is.na(colnames(dna_mat))] <- "SKBR3-I"
  dna_mat <- dna_mat[rowSums(is.na(dna_mat)) < 5, ]
  
  row_labs <- paste0(
    HM450K_ProbeSet[select_features_dnam, "transcript_geneName_noENS_NA"],
    " | ",
    HM450K_ProbeSet[select_features_dnam, "region_simple"]
  )
  message("  ✓ DNAm: ", nrow(dna_mat), " probes selected")
  
  return(list(
    mrna_mat = mrna_mat,
    rppa_mat = rppa_mat,
    # receptor_mat = receptor_mat,
    dna_mat = dna_mat,
    dna_row_labs = row_labs
  ))
}

# ==============================================================================
# 4. Create Heatmap Annotations
# ==============================================================================

myAnnoBarplot <- function(x, col = "#B0BEC5") {
  anno_barplot(
    x,
    gp = gpar(fill = col),
    axis = FALSE,
    border_gp = gpar(col = "gray"),
    height = unit(4, "mm"),
    bar_width = 0.4
  )
}

create_top_annotation <- function(alt_mat, metrics) {
  # message("4. Creating heatmap annotations...")
  
  # Create top annotation
  top_annotation <- HeatmapAnnotation(
    Subtype = CL_Annots[colnames(alt_mat), "mRNA Subtypes"],
    FGA = myAnnoBarplot(
      CL_Annots[gsub("SKBR3-I", "SKBR3-C", colnames(alt_mat)), "FGA"],
      "#F48FB1"
    ),
    SV = myAnnoBarplot(metrics$SV_count[colnames(alt_mat), "Count"], "#B3E5FC"),
    Fusions = myAnnoBarplot(metrics$Fusions_count[colnames(alt_mat), "Count"], "#C8E6C9"),
    TMB = myAnnoBarplot(metrics$TMB[colnames(alt_mat), "total_perMB"], "#D1C4E9"),
    DMI = myAnnoBarplot(metrics$DMI[colnames(alt_mat), "DMI_zscore"], "#FFCC80"),
    annotation_name_gp = gpar(
      fontsize = 12,
      fontface = "bold",
      col = "black",
      border = TRUE,
      fontfamily = helv
    ),
    gap = unit(0.5, "mm"),
    gp = gpar(col = NA, lwd = 0.6),
    annotation_legend_param = heatmap_legend_param,
    simple_anno_size = unit(4, "mm"),
    annotation_name_side = "right",
    col = list(Subtype = annot_cols$Subtypes)
  )
  
  message("  ✓ Top annotation created")
  
  return(top_annotation)
}

# ==============================================================================
# 5. Create Individual Heatmaps
# ==============================================================================

# Helper function to create common gpar objects for heatmaps
get_heatmap_gpars <- function() {
  list(
    column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
    column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
    row_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
    row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv)
  )
}

# Helper function to get common heatmap base parameters
# Returns a list that can be used with do.call() or passed individually
get_heatmap_base_params <- function(mat) {
  list(
    border = TRUE,
    border_gp = gpar(col = "black"),
    column_labels = gsub("-I", "", colnames(mat)),
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_column_names = TRUE,
    column_title = " ",
    na_col = "gray",
    heatmap_legend_param = heatmap_legend_param,
    column_split = CL_Annots[colnames(mat), "mRNA Subtypes"],
    width = unit(7, "cm")
  )
}

create_alteration_heatmap <- function(alt_mat, top_annotation) {
  # message("5. Creating alteration heatmap...")
  
  set.seed(123)
  gpars <- get_heatmap_gpars()
  base_params <- get_heatmap_base_params(alt_mat)
  
  ht <- Heatmap(
    as.matrix(alt_mat),
    col = annot_cols$Alt_Simple,
    row_title = " ",
    name = "Alterations",
    top_annotation = top_annotation,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    height = unit(3.3, "cm"),
    column_names_gp = gpars$column_names_gp,
    column_title_gp = gpars$column_title_gp,
    row_names_gp = gpars$row_names_gp,
    row_title_gp = gpars$row_title_gp,
    border = base_params$border,
    border_gp = base_params$border_gp,
    column_labels = base_params$column_labels,
    show_column_dend = base_params$show_column_dend,
    show_row_dend = base_params$show_row_dend,
    show_column_names = base_params$show_column_names,
    column_title = base_params$column_title,
    na_col = base_params$na_col,
    heatmap_legend_param = base_params$heatmap_legend_param,
    column_split = base_params$column_split,
    width = base_params$width
  )
  
  message("  ✓ Alteration heatmap created")
  return(ht)
}

create_receptor_heatmap <- function(receptor_mat) {
  # message("6. Creating receptor RPPA heatmap...")
  
  set.seed(123)
  gpars <- get_heatmap_gpars()
  base_params <- get_heatmap_base_params(receptor_mat)
  
  ht <- Heatmap(
    receptor_mat,
    col = annot_cols$RPPA,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = FALSE,
    name = "Protein",
    row_title = " ",
    show_row_names = TRUE,
    height = unit(1.5, "cm"),
    column_names_gp = gpars$column_names_gp,
    column_title_gp = gpars$column_title_gp,
    row_names_gp = gpars$row_names_gp,
    row_title_gp = gpars$row_title_gp,
    border = base_params$border,
    border_gp = base_params$border_gp,
    column_labels = base_params$column_labels,
    show_column_dend = base_params$show_column_dend,
    show_row_dend = base_params$show_row_dend,
    show_column_names = base_params$show_column_names,
    column_title = base_params$column_title,
    na_col = base_params$na_col,
    heatmap_legend_param = base_params$heatmap_legend_param,
    column_split = base_params$column_split,
    width = base_params$width
  )
  
  message("  ✓ Receptor heatmap created")
  return(ht)
}

create_rppa_heatmap <- function(rppa_mat) {
  # message("7. Creating RPPA expression heatmap...")
  
  set.seed(123)
  gpars <- get_heatmap_gpars()
  base_params <- get_heatmap_base_params(rppa_mat)
  
  ht <- Heatmap(
    rppa_mat,
    col = annot_cols$RPPA,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = TRUE,
    name = "Protein",
    row_title = "Top 50",
    show_row_names = FALSE,
    height = unit(2.5, "cm"),
    column_names_gp = gpars$column_names_gp,
    column_title_gp = gpars$column_title_gp,
    row_names_gp = gpars$row_names_gp,
    row_title_gp = gpars$row_title_gp,
    border = base_params$border,
    border_gp = base_params$border_gp,
    column_labels = base_params$column_labels,
    show_column_dend = base_params$show_column_dend,
    show_row_dend = base_params$show_row_dend,
    show_column_names = base_params$show_column_names,
    column_title = base_params$column_title,
    na_col = base_params$na_col,
    heatmap_legend_param = base_params$heatmap_legend_param,
    column_split = base_params$column_split,
    width = base_params$width
  )
  
  message("  ✓ RPPA heatmap created")
  return(ht)
}

create_rna_heatmap <- function(mrna_mat) {
  # message("8. Creating RNA expression heatmap...")
  
  set.seed(123)
  gpars <- get_heatmap_gpars()
  base_params <- get_heatmap_base_params(mrna_mat)
  
  ht <- Heatmap(
    mrna_mat,
    col = annot_cols$RNA_zscore,
    name = "RNA",
    row_title = "Top 6000",
    show_row_names = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    height = unit(2.5, "cm"),
    use_raster = TRUE,
    raster_quality = 5,
    column_names_gp = gpars$column_names_gp,
    column_title_gp = gpars$column_title_gp,
    row_names_gp = gpars$row_names_gp,
    row_title_gp = gpars$row_title_gp,
    border = base_params$border,
    border_gp = base_params$border_gp,
    column_labels = base_params$column_labels,
    show_column_dend = base_params$show_column_dend,
    show_row_dend = base_params$show_row_dend,
    show_column_names = base_params$show_column_names,
    column_title = base_params$column_title,
    na_col = base_params$na_col,
    heatmap_legend_param = base_params$heatmap_legend_param,
    column_split = base_params$column_split,
    width = base_params$width
  )
  
  message("  ✓ RNA heatmap created")
  return(ht)
}

create_dnam_heatmap <- function(dna_mat, row_labs) {
  # message("9. Creating DNA methylation heatmap...")
  
  set.seed(123)
  gpars <- get_heatmap_gpars()
  base_params <- get_heatmap_base_params(dna_mat)
  
  ht <- Heatmap(
    dna_mat,
    col = annot_cols$DNAm_zscore,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    cluster_rows = TRUE,
    name = "DNAm",
    row_title = "Top 6000",
    show_row_names = FALSE,
    row_labels = row_labs,
    height = unit(2.5, "cm"),
    use_raster = TRUE,
    raster_quality = 5,
    column_names_gp = gpars$column_names_gp,
    column_title_gp = gpars$column_title_gp,
    row_names_gp = gpars$row_names_gp,
    row_title_gp = gpars$row_title_gp,
    border = base_params$border,
    border_gp = base_params$border_gp,
    column_labels = base_params$column_labels,
    show_column_dend = base_params$show_column_dend,
    show_row_dend = base_params$show_row_dend,
    show_column_names = base_params$show_column_names,
    column_title = base_params$column_title,
    na_col = base_params$na_col,
    heatmap_legend_param = base_params$heatmap_legend_param,
    column_split = base_params$column_split,
    width = base_params$width
  )
  
  message("  ✓ DNAm heatmap created")
  return(ht)
}

# ==============================================================================
# 6. Main Execution
# ==============================================================================

main <- function() {
  # Load data
  message("  Step 1/6: Loading additional genomic metrics...")
  metrics <- load_additional_metrics()
  message("  ✓ Metrics loaded\n")
  
  # Prepare matrices
  message("  Step 2/6: Preparing genetic alteration matrix...")
  alt_results <- prepare_alteration_matrix()
  alt_mat <- alt_results$alt_mat
  ICLE_cells <- alt_results$ICLE_cells
  message("  ✓ Alteration matrix prepared: ", nrow(alt_mat), " genes x ", 
          ncol(alt_mat), " samples\n")
  
  # Select features
  message("  Step 3/6: Selecting top variable features...")
  features <- select_top_features(ICLE_cells)
  
  # Create annotations
  message("\n  Step 4/6: Creating annotations...")
  top_annotation <- create_top_annotation(alt_mat, metrics)
  
  # Create individual heatmaps
  message("\n  Step 5/6: Creating individual heatmaps...")
  ALT_ht <- create_alteration_heatmap(alt_mat, top_annotation)
  # receptor_ht <- create_receptor_heatmap(features$receptor_mat)
  RPPA_ht <- create_rppa_heatmap(features$rppa_mat)
  mRNA_ht <- create_rna_heatmap(features$mrna_mat)
  DNAm_ht <- create_dnam_heatmap(features$dna_mat, features$dna_row_labs)
  
  
  message("\n========================================")
  message("Figure 1D Completed")
  message("========================================\n")
  # receptor_ht
  return(ALT_ht %v% RPPA_ht %v% mRNA_ht %v% DNAm_ht)
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/04_Fig1D_Multiomics_Overview.R")
#
# load_all_icle_data(load_external = FALSE)
#
# fig1d_multiomics_overview <- main()
# draw(fig1d_multiomics_overview, merge_legend = TRUE,
#      annotation_legend_side = "right", heatmap_legend_side = "right")
#
# main() saves Fig1D to DIRS$results/Fig1D_Multiome_Overview.pdf.

# Run main function and assign to .GlobalEnv for downstream use
fig1d_multiomics_overview <- suppressWarnings(main())
