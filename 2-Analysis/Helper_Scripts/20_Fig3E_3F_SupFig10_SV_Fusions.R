# ==============================================================================
# Script 21: Fig 3E, 3F & SupFig 10 - SV Fusions Analysis
# ==============================================================================
# Description: Identifies and visualizes gene fusions from Bionano structural
#              variations. Generates functional fusion heatmaps and circos plots
#              for ICLE cell lines and TCGA tumors.
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - ICLE_SV: Structural variation data
#   - BRCA_CL_EXP: RNA expression data
#   - CL_Annots: Cell line annotations
#   - oncoVar: OncoVar gene annotations (loaded via Helper_Functions.R)
#
# Output Objects (assigned to .GlobalEnv):
#   - fusions_df: Complete fusion data frame with annotations
#   - fig3e_1: ICLE fusion distribution barplot
#   - fig3e_2: Tumor fusion distribution barplot
#   - supfig10a_fusion_breakpoints_ht: Fusion breakpoint heatmap
#   - supfig10b_recurring_fusions_ht: Recurring fusions expression heatmap
#   - fig3f_left_goe_fusions_circos: GOE fusions circos plot
#   - fig3f_right_loe_fusions_circos: LOE fusions circos plot
#   - recurring_gene_circos: List of circos plots for top recurring genes
#
# Author: Osama Shiraz Shah
# ==============================================================================

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.1 Resolve Complex Fusions
# ------------------------------------------------------------------------------

#' Resolve complex fusion strings with multiple candidates
#' 
#' @param fusion_string Character string with fusion candidates separated by ";"
#' @param gene_driver_scores Named vector of driver scores (default: global)
#' @return Single resolved fusion string
#' @details
#'   When a fusion has multiple candidates (e.g., "GENE1-GENE2;GENE3-GENE4"),
#'   this function selects the best candidate based on OncoVar driver scores.
#'   For genes with multiple candidates on one side, prioritizes driver genes.
resolve_complex_fusion <- function(fusion_string, gene_driver_scores = NULL) {
  if (is.null(gene_driver_scores)) {
    if (!exists("gene_driver_scores", envir = .GlobalEnv)) {
      stop("gene_driver_scores not found. Load oncoVar first.")
    }
    # Use global variable directly via .GlobalEnv$gene_driver_scores
  }
  
  # Split and get unique fusions
  candidates <- unique(unlist(strsplit(fusion_string, split = ";")))
  if (length(candidates) == 1) return(candidates)
  
  # Parse genes for all candidates
  temp_df <- data.frame(
    fusion = candidates,
    geneA = sapply(strsplit(candidates, "-"), `[`, 1),
    geneB = sapply(strsplit(candidates, "-"), `[`, 2),
    stringsAsFactors = FALSE
  )
  
  # Determine which side has multiple genes
  unique_A <- unique(temp_df$geneA)
  unique_B <- unique(temp_df$geneB)
  multiple_A <- length(unique_A) > 1
  multiple_B <- length(unique_B) > 1
  
  # Calculate scores based on "Multiple" rule
  # For the geneA or geneB which are multiple, only take the one that is oncovar driver
  # Use global variable if parameter is NULL
  driver_scores <- if (is.null(gene_driver_scores)) .GlobalEnv$gene_driver_scores else gene_driver_scores
  temp_df$score <- mapply(function(gA, gB) {
    s <- 0
    if (multiple_A) s <- s + get_oncovar_score(gA, driver_scores)
    if (multiple_B) s <- s + get_oncovar_score(gB, driver_scores)
    return(s)
  }, temp_df$geneA, temp_df$geneB)
  
  # Pick winner (highest score)
  winner_idx <- order(temp_df$score, decreasing = TRUE)[1]
  return(temp_df$fusion[winner_idx])
}

#' Get OncoVar driver score for a gene
#' 
#' @param gene Gene symbol
#' @param db Named vector of driver scores (default: global gene_driver_scores)
#' @return Driver score (0 if not found)
get_oncovar_score <- function(gene, db = NULL) {
  if (is.null(db)) {
    if (!exists("gene_driver_scores", envir = .GlobalEnv)) {
      return(0)
    }
    # Use global variable directly
    db <- .GlobalEnv$gene_driver_scores
  }
  score <- db[gene]
  if (is.na(score)) return(0)
  return(score)
}

# ------------------------------------------------------------------------------
# 1.2 Prepare Fusion Data Frame
# ------------------------------------------------------------------------------

#' Prepare fusion data frame from ICLE_SV data
#' 
#' @param ICLE_SV Structural variation data frame
#' @return Data frame with fusion information
prepare_fusion_dataframe <- function(ICLE_SV) {
  # Extract all putative fusions
  fusions <- subset(ICLE_SV, PutativeGeneFusion != "-")$PutativeGeneFusion
  
  # Create initial fusion data frame
  fusions_df <- data.frame(
    fusions = fusions,
    Sample = ICLE_SV[match(fusions, ICLE_SV$PutativeGeneFusion), "Sample"],
    type = ICLE_SV[match(fusions, ICLE_SV$PutativeGeneFusion), "Type"],
    geneA = sapply(strsplit(fusions, split = "-"), FUN = function(x) x[1]),
    chrA = ICLE_SV[match(fusions, ICLE_SV$PutativeGeneFusion), "RefcontigID1"],
    StartA = ICLE_SV[match(fusions, ICLE_SV$PutativeGeneFusion), "RefStartPos"],
    EndA = ICLE_SV[match(fusions, ICLE_SV$PutativeGeneFusion), "RefStartPos"],
    geneB = sapply(strsplit(fusions, split = "-"), FUN = function(x) x[2]),
    chrB = ICLE_SV[match(fusions, ICLE_SV$PutativeGeneFusion), "RefcontigID2"],
    StartB = ICLE_SV[match(fusions, ICLE_SV$PutativeGeneFusion), "RefEndPos"],
    EndB = ICLE_SV[match(fusions, ICLE_SV$PutativeGeneFusion), "RefEndPos"],
    stringsAsFactors = FALSE
  )
  
  # Remove duplicates based on unique coordinates
  fusions_df$unique_id <- paste0(fusions_df$StartA, fusions_df$EndA, 
                                   fusions_df$StartB, fusions_df$EndB, 
                                   fusions_df$Sample)
  fusions_df <- fusions_df[!duplicated(fusions_df$unique_id), ]
  rownames(fusions_df) <- fusions_df$fusions
  
  return(fusions_df)
}

# ------------------------------------------------------------------------------
# 1.3 Add Cytoband Information
# ------------------------------------------------------------------------------

#' Add cytoband information to fusion data frame
#' 
#' @param fusions_df Fusion data frame
#' @return Fusion data frame with cytoBandA and cytoBandB columns
add_cytoband_info <- function(fusions_df) {
  library(data.table)
  
  # Load cytoband data if not available
  if (!exists("cytoBand")) {
    data(cytoBand, package = "DNAcopy")
  }
  
  # Prepare cytoband data table
  cytoBand_dt <- as.data.table(DNAcopy::cytoBand[, c("chromNum", "chromStart", "chromEnd", "bandname")])
  cytoBand_dt[, chromNum := gsub("chr0?(\\d+|[XY])", "chr\\1", chromNum)]
  fusions_dt <- as.data.table(fusions_df)
  
  # Add cytoband for geneA
  fusions_dt[cytoBand_dt,
             on = .(chrA = chromNum, StartA >= chromStart, StartA <= chromEnd),
             cytoBandA := i.bandname]
  
  # Add cytoband for geneB
  fusions_dt[cytoBand_dt,
             on = .(chrB = chromNum, StartB >= chromStart, StartB <= chromEnd),
             cytoBandB := i.bandname]
  
  fusions_df <- as.data.frame(fusions_dt)
  return(fusions_df)
}

# ------------------------------------------------------------------------------
# 1.4 Calculate Fusion Gene Expression
# ------------------------------------------------------------------------------

#' Calculate expression status for fusion genes
#' 
#' @param fusions_df Fusion data frame
#' @param BRCA_CL_EXP RNA expression matrix
#' @param CL_Annots Cell line annotations
#' @return Fusion data frame with geneA_exp, geneB_exp, and exp columns
calculate_fusion_expression <- function(fusions_df, BRCA_CL_EXP, CL_Annots) {
  # Calculate z-scores for ICLE samples only
  icle_samples <- CL_Annots[CL_Annots$overlapWCCLE != "Y", "Name"]
  RNA_zscore <- t(scale(t(BRCA_CL_EXP[, icle_samples])))
  
  # Initialize expression columns
  fusions_df$geneA_exp <- NA
  fusions_df$geneB_exp <- NA
  
  # Calculate expression for each fusion
  for (fusion in rownames(fusions_df)) {
    # Clean gene names (remove antisense suffixes)
    geneA <- gsub(x = fusions_df[fusion, "geneA"], pattern = "-AS|-AS1|-AS2|-IT1", replacement = "")
    geneB <- gsub(x = fusions_df[fusion, "geneB"], pattern = "-AS|-AS1|-AS2|-IT1", replacement = "")
    Sample <- paste0(fusions_df[fusion, "Sample"], "-I")
    
    # Gene A expression
    if (geneA %in% rownames(BRCA_CL_EXP) && Sample %in% colnames(RNA_zscore)) {
      zscore_A <- RNA_zscore[geneA, Sample]
      fusions_df[fusion, "geneA_exp"] <- ifelse(zscore_A > 1, "up",
                                                 ifelse(zscore_A < -1, "dn", NA))
    }
    
    # Gene B expression
    if (geneB %in% rownames(BRCA_CL_EXP) && Sample %in% colnames(RNA_zscore)) {
      zscore_B <- RNA_zscore[geneB, Sample]
      fusions_df[fusion, "geneB_exp"] <- ifelse(zscore_B > 1, "up",
                                                ifelse(zscore_B < -1, "dn", NA))
    }
  }
  
  # Determine overall fusion expression status
  fusions_df$geneA_exp[is.na(fusions_df$geneA_exp)] <- "-"
  fusions_df$geneB_exp[is.na(fusions_df$geneB_exp)] <- "-"
  
  fusions_df$exp <- ifelse(fusions_df$geneA_exp == "up" | fusions_df$geneB_exp == "up", 
                          "up-regulated",
                          ifelse(fusions_df$geneA_exp == "dn" | fusions_df$geneB_exp == "dn", 
                                "down-regulated", "-"))
  
  # Exclude conflicting expression patterns
  fusions_df$exp[fusions_df$geneA_exp == "dn" & fusions_df$geneB_exp == "up"] <- "-"
  fusions_df$exp[fusions_df$geneA_exp == "up" & fusions_df$geneB_exp == "dn"] <- "-"
  
  # Add expression color for visualization
  fusions_df$exp_col <- ifelse(fusions_df$exp == "up-regulated", "red",
                               ifelse(fusions_df$exp == "down-regulated", "blue", "gray"))
  
  # Mark functional fusions
  fusions_df$functional_fusions <- ifelse(fusions_df$exp != "-", "Y", "N")
  
  return(fusions_df)
}

# ------------------------------------------------------------------------------
# 1.5 Add OncoVar Annotations
# ------------------------------------------------------------------------------

#' Add OncoVar annotations to fusion data frame
#' 
#' @param fusions_df Fusion data frame
#' @param oncoVar OncoVar data frame
#' @return Fusion data frame with OncoKB annotations
add_oncovar_annotations <- function(fusions_df, oncoVar) {
  # Add OncoKB status for geneA
  fusions_df$geneA_OncoTS <- ifelse(
    grepl(x = oncoVar[fusions_df$geneA, "OncoKB"], pattern = "Y"), NA,
    ifelse(grepl(x = oncoVar[fusions_df$geneA, "OncoKB"], pattern = "TSG"), "TSG",
           ifelse(grepl(x = oncoVar[fusions_df$geneA, "OncoKB"], pattern = "Oncogene"), 
                 "Oncogene", NA))
  )
  
  # Add OncoKB status for geneB
  fusions_df$geneB_OncoTS <- ifelse(
    grepl(x = oncoVar[fusions_df$geneB, "OncoKB"], pattern = "Y"), NA,
    ifelse(grepl(x = oncoVar[fusions_df$geneB, "OncoKB"], pattern = "TSG"), "TSG",
           ifelse(grepl(x = oncoVar[fusions_df$geneB, "OncoKB"], pattern = "Oncogene"), 
                 "Oncogene", NA))
  )
  
  return(fusions_df)
}

# ==============================================================================
# SECTION 2: DATA PREPARATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 Setup and Load Dependencies
# ------------------------------------------------------------------------------

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(ggthemes)
  library(data.table)
})

# %notin% operator is defined in Helper_Functions.R
if (!exists("%notin%", envir = .GlobalEnv)) {
  `%notin%` <- Negate(`%in%`)
}

# Check if required data objects exist, load if needed
required_objects <- c("ICLE_SV", "BRCA_CL_EXP", "CL_Annots")
missing_objects <- required_objects[!sapply(required_objects, exists, envir = .GlobalEnv)]

if (length(missing_objects) > 0) {
  message("  Loading missing data objects: ", paste(missing_objects, collapse = ", "))
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = TRUE)
}

# Load oncoVar
oncoVar <- read.delim(FILES$oncovar)
rownames(oncoVar) <- oncoVar$Gene_symbol

# Create gene driver scores lookup
gene_driver_scores <- unlist(split(oncoVar$Driver.Level, oncoVar$Gene_symbol))
assign("gene_driver_scores", gene_driver_scores, envir = .GlobalEnv)

# Load reference gene database
refGene_db <- read.delim(FILES$refgene, header = TRUE) %>% distinct(genename, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 2.2 Prepare ICLE Fusion Data
# ------------------------------------------------------------------------------

message("  Preparing ICLE fusion data...")

# Create initial fusion data frame
fusions_df <- prepare_fusion_dataframe(ICLE_SV)
rownames(fusions_df) <- fusions_df$fusions

# Add cytoband information
fusions_df <- add_cytoband_info(fusions_df)

# Resolve complex fusions (multiple candidates)
fusions_df$fusions <- sapply(fusions_df$fusions, resolve_complex_fusion)
fusions_df$geneA <- sapply(strsplit(fusions_df$fusions, split = "-"), FUN = function(x) x[1])
fusions_df$geneB <- sapply(strsplit(fusions_df$fusions, split = "-"), FUN = function(x) x[2])
fusions_df$fusions_inv <- paste0(fusions_df$geneB, "-", fusions_df$geneA)

# ------------------------------------------------------------------------------
# 2.3 Load and Prepare Tumor Fusion Data
# ------------------------------------------------------------------------------

message("  Loading tumor fusion data...")

TumorFusions <- readxl::read_excel(FILES$tumor_fusions, skip = 2) %>% as.data.frame()
TumorFusions$FusionPair <- gsub(x = TumorFusions$FusionPair, pattern = "__", replacement = "-")
TumorFusions$FusionPair <- gsub(x = TumorFusions$FusionPair, pattern = "__", replacement = "-")

# Clean and standardize tumor fusion data
TumorFusions <- TumorFusions[!is.na(TumorFusions$chrA) & !is.na(TumorFusions$chrB), ]
TumorFusions$exp_col <- "gray"
TumorFusions$StartA <- TumorFusions$posA
TumorFusions$EndA <- TumorFusions$posA
TumorFusions$StartB <- TumorFusions$posB
TumorFusions$EndB <- TumorFusions$posB
TumorFusions$chrA <- paste0("chr", TumorFusions$chrA)
TumorFusions$chrB <- paste0("chr", TumorFusions$chrB)

# Assign to global environment for use in visualization
assign("TumorFusions", TumorFusions, envir = .GlobalEnv)

# Match ICLE fusions with tumor fusions
fusion_match <- c(
  intersect(TumorFusions$FusionPair, fusions_df$fusions),
  intersect(TumorFusions$FusionPair, fusions_df$fusions_inv)
)
fusions_df$`Tumor Fusion Match` <- ifelse(fusions_df$fusions %in% fusion_match, "Y", "N")

# Match fusion partners (genes involved in fusions)
fusion_partner_match <- subset(fusions_df, 
                               (geneA %in% TumorFusions$`gene 1` | geneB %in% TumorFusions$`gene 1`) | 
                               (geneA %in% TumorFusions$`RNAseq gene 2` | geneB %in% TumorFusions$`RNAseq gene 2`)
)$fusions
fusions_df$`Tumor Fusion Partner Match` <- ifelse(fusions_df$fusions %in% fusion_partner_match, "Y", "N")

# ------------------------------------------------------------------------------
# 2.4 Calculate Fusion Expression
# ------------------------------------------------------------------------------

message("  Calculating fusion gene expression...")
fusions_df <- calculate_fusion_expression(fusions_df, BRCA_CL_EXP, CL_Annots)

# ------------------------------------------------------------------------------
# 2.5 Add OncoVar Annotations
# ------------------------------------------------------------------------------

message("  Adding OncoVar annotations...")
fusions_df <- add_oncovar_annotations(fusions_df, oncoVar)

# ------------------------------------------------------------------------------
# 2.6 Identify Recurring Fusions
# ------------------------------------------------------------------------------

message("  Identifying recurring fusions...")

functional_fusions <- subset(fusions_df, functional_fusions == "Y")

# Count how many samples each gene appears in
gene_sample_table <- table(
  c(functional_fusions$geneB, functional_fusions$geneA),
  c(functional_fusions$Sample, functional_fusions$Sample)
)
sample_counts <- rowSums(gene_sample_table > 0)
recurGenes <- head(sort(sample_counts, decreasing = TRUE), 15)

# Mark recurring fusions
recurring_fusions <- subset(functional_fusions, 
                           geneA %in% names(recurGenes) | geneB %in% names(recurGenes))
fusions_df$recurGene <- ifelse(
  fusions_df$geneA %in% names(recurGenes), fusions_df$geneA,
  ifelse(fusions_df$geneB %in% names(recurGenes), fusions_df$geneB, "")
)
fusions_df$recurring <- ifelse(rownames(fusions_df) %in% rownames(recurring_fusions), "Y", "N")

# Assign to global environment
assign("fusions_df", fusions_df, envir = .GlobalEnv)
assign("recurGenes", recurGenes, envir = .GlobalEnv)
assign("functional_fusions", functional_fusions, envir = .GlobalEnv)

message("  ✓ Fusion data preparation complete")

# ==============================================================================
# SECTION 3: VISUALIZATION FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Create Fusion Distribution Plot
# ------------------------------------------------------------------------------

#' Create fusion distribution barplot by chromosome
#' 
#' @param fusions_df Fusion data frame
#' @param top_n Number of top chromosomes to display (default: 10)
#' @param include_expression Logical, whether to include expression grouping
#' @return ggplot object
create_fusion_distribution_plot <- function(fusions_df, top_n = 10, include_expression = TRUE) {
  if (include_expression) {
    # Prepare data with expression grouping
    df_plot <- fusions_df %>% 
      dplyr::select(geneA, chrA, geneB, chrB, exp) %>% 
      {bind_rows(
        dplyr::select(., gene = geneA, chrom = chrA, exp),
        dplyr::select(., gene = geneB, chrom = chrB, exp)
      )} %>%
      dplyr::mutate(group = factor(dplyr::case_when(
        exp == "up-regulated" ~ "Up-Regulated",
        exp == "down-regulated" ~ "Down-Regulated",
        TRUE ~ "No-Change"
      ), levels = rev(c("Up-Regulated", "Down-Regulated", "No-Change"))))
    
    # Get top chromosomes
    top_chroms <- df_plot %>% 
      dplyr::count(chrom, sort = TRUE) %>% 
      dplyr::slice_head(n = top_n) %>% 
      dplyr::pull(chrom)
    
    # Create plot
    plt <- df_plot %>% 
      dplyr::filter(chrom %in% top_chroms) %>%
      dplyr::count(chrom, group) %>% 
      ggplot(aes(x = n, y = reorder(chrom, n), fill = group)) +
      geom_bar(position = "stack", stat = "identity", show.legend = FALSE) +
      scale_fill_manual(values = c("Down-Regulated" = "blue", 
                                   "Up-Regulated" = "red", 
                                   "No-Change" = "grey70")) +
      labs(y = " ", x = "Fusion Count") + 
      theme_clean(base_size = 20) + 
      theme(panel.grid.major.y = element_blank())
  } else {
    # Prepare data without expression grouping (for TumorFusions)
    # Check if this is TumorFusions (has 'gene 1' column) or regular fusions_df
    if ("gene 1" %in% colnames(fusions_df)) {
      df_plot <- fusions_df %>% 
        dplyr::select(`gene 1`, chrA, `gene 2`, chrB) %>% 
        {bind_rows(
          dplyr::select(., gene = `gene 1`, chrom = chrA),
          dplyr::select(., gene = `gene 2`, chrom = chrB)
        )}
    } else {
      # For regular fusions_df without expression, use geneA/geneB
      df_plot <- fusions_df %>% 
        dplyr::select(geneA, chrA, geneB, chrB) %>% 
        {bind_rows(
          dplyr::select(., gene = geneA, chrom = chrA),
          dplyr::select(., gene = geneB, chrom = chrB)
        )}
    }
    
    # Get top chromosomes
    top_chroms <- df_plot %>% 
      dplyr::count(chrom, sort = TRUE) %>% 
      dplyr::slice_head(n = top_n) %>% 
      dplyr::pull(chrom)
    
    # Create plot
    plt <- df_plot %>% 
      dplyr::filter(chrom %in% top_chroms) %>%
      dplyr::count(chrom) %>% 
      ggplot(aes(x = n, y = reorder(chrom, n))) +
      geom_bar(position = "stack", stat = "identity", show.legend = FALSE) +
      labs(y = " ", x = "Fusion Count") + 
      theme_clean(base_size = 20) + 
      theme(panel.grid.major.y = element_blank())
  }
  
  return(plt)
}

# ------------------------------------------------------------------------------
# 3.2 Create Fusion Breakpoint Heatmap
# ------------------------------------------------------------------------------

#' Create heatmap of fusion breakpoints by chromosome
#' 
#' @param fusions_df Fusion data frame
#' @param heatmap_legend_param Legend parameters
#' @return ComplexHeatmap object
create_fusion_breakpoint_heatmap <- function(fusions_df, heatmap_legend_param) {
  # Define chromosome order
  chroms <- paste0("chr", c(1:22, "X", "Y"))
  fusions_df$chrA <- factor(fusions_df$chrA, levels = chroms)
  fusions_df$chrB <- factor(fusions_df$chrB, levels = chroms)
  
  # Create count matrix (only functional fusions)
  mat <- as.matrix(table(
    subset(fusions_df, exp != "-")$chrA,
    subset(fusions_df, exp != "-")$chrB
  ))
  
  # Color function
  col_fun <- colorRamp2(c(0, max(mat)), c("white", "#F57F17"))
  
  # Row annotation with counts
  row_ha <- rowAnnotation(
    "Count" = anno_barplot(rowSums(mat), border = FALSE, 
                          gp = gpar(fill = "gray", col = "black"),
                          axis_param = list(side = "bottom")),
    width = unit(1, "cm")
  )
  
  # Clean chromosome labels
  rownames(mat) <- gsub("chr", "", rownames(mat))
  colnames(mat) <- gsub("chr", "", colnames(mat))
  
  # Create heatmap
  ht <- Heatmap(
    mat, 
    name = "Count", 
    col = col_fun, 
    heatmap_legend_param = heatmap_legend_param,
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    rect_gp = gpar(col = "gray", lwd = 0.5), 
    right_annotation = row_ha,
    row_names_side = "left", 
    column_names_side = "bottom",
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (mat[i, j] > 0) {
        grid.text(mat[i, j], x, y, gp = gpar(fontsize = 15))
      }
    }
  )
  
  return(ht)
}

# ------------------------------------------------------------------------------
# 3.3 Create Recurring Fusions Expression Heatmap
# ------------------------------------------------------------------------------

#' Create heatmap of recurring fusion gene expression
#' 
#' @param fusions_df Fusion data frame
#' @param heatmap_legend_param Legend parameters
#' @param annot_cols Annotation colors
#' @return ComplexHeatmap object
create_recurring_fusions_heatmap <- function(fusions_df, heatmap_legend_param, annot_cols) {
  rownames(fusions_df) <- fusions_df$fusions
  recurring_fusions <- subset(fusions_df, recurring == "Y")
  
  # Prepare expression matrix
  htMat <- t(recurring_fusions[, c("geneA_exp", "geneB_exp")])
  htMat <- htMat[, rownames(recurring_fusions)]
  rownames(htMat) <- c("Gene A Expression", "Gene B Expression")
  
  # Create heatmap
  set.seed(123)
  ht <- Heatmap(
    htMat,
    col = c("up" = "#ff760d", "dn" = "#99cc19", "-" = "gray"),
    name = "Expression",
    column_order = order(recurring_fusions$exp),
    heatmap_legend_param = heatmap_legend_param,
    column_split = factor(
      recurring_fusions$Sample,
      levels = rev(names(table(recurring_fusions$Sample))[
        order(table(recurring_fusions$Sample))
      ])
    ),
    column_title_rot = 90,
    show_column_names = TRUE,
    height = unit(2, "cm"),
    width = unit(20, "cm"),
    border = TRUE,
    border_gp = gpar(col = "black"),
    top_annotation = HeatmapAnnotation(
      simple_anno_size = unit(4, "mm"),
      border = TRUE,
      annotation_legend_param = heatmap_legend_param,
      SV = recurring_fusions$type,
      `Fusion Type` = recurring_fusions$exp,
      col = list(
        SV = annot_cols$SV,
        `Fusion Type` = c("down-regulated" = "blue", "up-regulated" = "red")
      )
    )
  )
  
  return(ht)
}

# ------------------------------------------------------------------------------
# 3.4 Create Circos Plot for Fusions
# ------------------------------------------------------------------------------

#' Create circos plot for fusion events
#' 
#' @param sv_list Data frame with fusion events
#' @param refGene_db Reference gene database (optional, uses global if not provided)
#' @param oncoVar OncoVar annotations (optional, uses global if not provided)
#' @param showLabs Logical, whether to show chromosome labels
#' @param driverLevels Vector of driver levels to include (default: c(2,3,4))
#' @param allChr Logical, whether to show all chromosomes (default: FALSE)
#' @return Recorded plot object
makeCircos <- function(sv_list, refGene_db = NULL, oncoVar = NULL, showLabs = TRUE, 
                      driverLevels = c(2, 3, 4), allChr = FALSE) {
  # Use global variables if not provided - reference them directly without assignment
  if (is.null(refGene_db)) {
    if (!exists("refGene_db", envir = .GlobalEnv)) {
      stop("refGene_db not found. Load reference gene database first.")
    }
  }
  
  if (is.null(oncoVar)) {
    if (!exists("oncoVar", envir = .GlobalEnv)) {
      stop("oncoVar not found. Load OncoVar annotations first.")
    }
  }
  
  # Use parameter if provided, otherwise use global directly
  refGene_db_use <- if (is.null(refGene_db)) .GlobalEnv$refGene_db else refGene_db
  oncoVar_use <- if (is.null(oncoVar)) .GlobalEnv$oncoVar else oncoVar
  
  chroms <- paste0("chr", c(1:22, "X", "Y"))
  sv_list$chrA <- factor(sv_list$chrA, levels = chroms)
  sv_list$chrB <- factor(sv_list$chrB, levels = chroms)
  
  # Get genes involved in fusions
  bed <- subset(refGene_db_use, genename %in% unique(c(sv_list$geneA, sv_list$geneB)))
  bed <- bed[!duplicated(bed$genename), ]
  rownames(bed) <- bed$genename
  
  # Filter to driver genes
  bed <- subset(bed, genename %in% unique(
    subset(oncoVar_use, Driver.Level %in% driverLevels)$Gene_symbol
  ))
  
  # Add OncoKB status
  bed$OncoTS <- ifelse(
    grepl(x = oncoVar_use[bed$genename, "OncoKB"], pattern = "Y"), "darkgray",
    ifelse(grepl(x = oncoVar_use[bed$genename, "OncoKB"], pattern = "TSG"), "blue",
           ifelse(grepl(x = oncoVar_use[bed$genename, "OncoKB"], pattern = "Oncogene"), 
                 "red", "darkgray"))
  )
  bed$chrom <- paste0("chr", bed$chrom)
  
  # Prepare coordinates
  coordA <- sv_list[, c("chrA", "StartA", "EndA")]
  coordB <- sv_list[, c("chrB", "StartB", "EndB")]
  
  library(circlize)
  
  # 1. Open a null device to suppress output
  pdf(NULL) 
  dev.control(displaylist = "enable")
  
  circos.clear()
  
  # Determine chromosomes to display
  chrs <- unique(c(sv_list$chrA, sv_list$chrB))
  if (allChr) {
    chrs <- paste0("chr", unique(refGene_db$chrom))
    chrs <- grep(x = chrs, pattern = "g|U|ss|hap|random|g1|M|Y", value = TRUE, invert = TRUE)
    chrs <- factor(chrs, levels = chroms)
  }
  
  # Initialize circos
  if (showLabs) {
    circos.initializeWithIdeogram(
      chromosome.index = chrs,
      species = "hg19",
      plotType = "labels",
      major.by = 0.2,
      sort.chr = FALSE
    )
  } else {
    circos.initializeWithIdeogram(
      plotType = NULL,
      chromosome.index = chrs,
      major.by = 0.1,
      species = "hg19"
    )
  }
  
  # Add gene labels
  circos.genomicLabels(bed, labels.column = 4, cex = 1.5, col = bed$OncoTS,
                       side = "outside", niceFacing = TRUE)
  
  # Add fusion links
  if (is.null(unique(sv_list$exp_col))) {
    circos.genomicLink(coordA, coordB, border = NA, lty = 1, lwd = 0.5)
  } else {
    circos.genomicLink(coordA, coordB, col = sv_list$exp_col, 
                      border = NA, lty = 1, lwd = 0.5)
  }
  
  circos.genomicIdeogram()
  circosP <- recordPlot()
  
  circos.clear()
  dev.off()
  
  return(circosP)
}

# ------------------------------------------------------------------------------
# 3.5 Create Circos Plot for Specific Gene
# ------------------------------------------------------------------------------

#' Create circos plot for fusions involving a specific gene
#' 
#' @param gene_name Gene symbol
#' @param fusions_df Fusion data frame
#' @param refGene_db Reference gene database (optional, uses global if not provided)
#' @param oncoVar OncoVar annotations (optional, uses global if not provided)
#' @param showLabs Logical, whether to show labels
#' @param driverLevels Vector of driver levels
#' @return Recorded plot object
create_gene_circos <- function(gene_name, fusions_df, refGene_db = NULL, oncoVar = NULL,
                               showLabs = FALSE, driverLevels = c(0, 1, 2, 3, 4)) {
  gene_fusions <- subset(fusions_df, 
                         geneA %in% gene_name | geneB %in% gene_name)
  circos_plot <- makeCircos(sv_list = gene_fusions, refGene_db = refGene_db,
                            oncoVar = oncoVar, showLabs = showLabs,
                            driverLevels = driverLevels, allChr = FALSE)
  return(circos_plot)
}

# ==============================================================================
# SECTION 4: GENERATE FIGURES
# ==============================================================================

message("  Generating fusion visualizations...")

# ------------------------------------------------------------------------------
# 4.1 Figure 3E: ICLE Fusion Distribution (Left)
# ------------------------------------------------------------------------------
message("\n========================================")
message("Figure 3E: ICLE Fusion Distribution (Left)")
message("========================================")

fig3e_1 <- create_fusion_distribution_plot(fusions_df, top_n = 10, include_expression = TRUE)
assign("fig3e_1", fig3e_1, envir = .GlobalEnv)
message("  ✓ Fig 3E (ICLE) complete")

# ------------------------------------------------------------------------------
# 4.2 Figure 3E: Tumor Fusion Distribution (Right)
# ------------------------------------------------------------------------------
message("\n========================================")
message("Figure 3E: Tumor Fusion Distribution (Right)")
message("========================================")

fig3e_2 <- create_fusion_distribution_plot(TumorFusions, top_n = 10, include_expression = FALSE)
assign("fig3e_2", fig3e_2, envir = .GlobalEnv)
message("  ✓ Fig 3E (Tumor) complete")

# ------------------------------------------------------------------------------
# 4.3 SupFig 10A: Fusion Breakpoint Heatmap
# ------------------------------------------------------------------------------
message("\n========================================")
message("SupFig 10A: Fusion Breakpoint Heatmap")
message("========================================")

supfig10a_fusion_breakpoints_ht <- create_fusion_breakpoint_heatmap(fusions_df, heatmap_legend_param)
assign("supfig10a_fusion_breakpoints_ht", supfig10a_fusion_breakpoints_ht, envir = .GlobalEnv)
message("  ✓ SupFig 10A complete")

# ------------------------------------------------------------------------------
# 4.4 SupFig 10B: Recurring Fusions Expression Heatmap
# ------------------------------------------------------------------------------
message("\n========================================")
message("SupFig 10B: Recurring Fusions Expression Heatmap")
message("========================================")

supfig10b_recurring_fusions_ht <- create_recurring_fusions_heatmap(
  fusions_df, heatmap_legend_param, annot_cols
)
assign("supfig10b_recurring_fusions_ht", supfig10b_recurring_fusions_ht, envir = .GlobalEnv)
message("  ✓ SupFig 10B complete")

# ------------------------------------------------------------------------------
# 4.5 Figure 3F: GOE and LOE Fusions Circos
# ------------------------------------------------------------------------------
message("\n========================================")
message("Figure 3F: GOE and LOE Fusions Circos")
message("========================================")

# GOE (Gain of Expression) fusions
fig3f_left_goe_fusions <- subset(fusions_df, exp == "up-regulated")
# Filter to top 3 chromosomes by count
top_chrA_goe <- names(tail(sort(table(fig3f_left_goe_fusions$chrA)), 3))
top_chrB_goe <- names(tail(sort(table(fig3f_left_goe_fusions$chrB)), 3))
fig3f_left_goe_fusions <- subset(fig3f_left_goe_fusions,
                                 chrA %in% top_chrA_goe & chrB %in% top_chrB_goe)

fig3f_left_goe_fusions_circos <- makeCircos(
  sv_list = fig3f_left_goe_fusions,
  showLabs = FALSE,
  driverLevels = c(1, 2, 3, 4)
)
assign("fig3f_left_goe_fusions_circos", fig3f_left_goe_fusions_circos, envir = .GlobalEnv)
message("  ✓ Fig 3F (GOE) complete")

# LOE (Loss of Expression) fusions
fig3f_right_loe_fusions <- subset(fusions_df, exp == "down-regulated")
# Filter to top 3 chromosomes by count
top_chrA_loe <- names(tail(sort(table(fig3f_right_loe_fusions$chrA)), 3))
top_chrB_loe <- names(tail(sort(table(fig3f_right_loe_fusions$chrB)), 3))
fig3f_right_loe_fusions <- subset(fig3f_right_loe_fusions,
                                  chrA %in% top_chrA_loe & chrB %in% top_chrB_loe)

fig3f_right_loe_fusions_circos <- makeCircos(
  sv_list = fig3f_right_loe_fusions,
  showLabs = FALSE,
  driverLevels = c(1, 2, 3, 4)
)
assign("fig3f_right_loe_fusions_circos", fig3f_right_loe_fusions_circos, envir = .GlobalEnv)
message("  ✓ Fig 3F (LOE) complete")

# ------------------------------------------------------------------------------
# 4.6 SupFig 10C: Recurring Gene Circos Plots
# ------------------------------------------------------------------------------
message("\n========================================")
message("SupFig 10C: Recurring Gene Circos Plots")
message("========================================")

# Create circos plots for top 5 recurring genes
recurring_gene_circos <- list()
top_genes <- names(recurGenes)[1:min(5, length(recurGenes))]

for (gene in top_genes) {
  gene_circos <- create_gene_circos(
    gene_name = gene,
    fusions_df = fusions_df,
    showLabs = FALSE,
    driverLevels = c(0, 1, 2, 3, 4)
  )
  recurring_gene_circos[[gene]] <- gene_circos
}

assign("recurring_gene_circos", recurring_gene_circos, envir = .GlobalEnv)
message("  ✓ SupFig 10C (recurring genes) complete")

# message("═══════════════════════════════════════════════════════")
# message("  Fusion analysis complete")
# message("═══════════════════════════════════════════════════════")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/21_Fig3E_3F_SupFig10_SV_Fusions.R")
#
# load_all_icle_data(load_external = TRUE)
# # All figures are generated and assigned to global environment
# # Save figures in Main_Data_Analysis.Rmd
