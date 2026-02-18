# ==============================================================================
# Script 22: Fig 4 & SupFig 11 - DNA Methylation Alterations (DMI, ILC vs NST)
# ==============================================================================
# Description: Analyzes DNA methylation instability (DMI) and generates
#              SupFig 11 and Fig 4A-D showing DNAm alterations across tissue types
#              (cell lines, primary tumors, metastatic tumors).
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
# Differential analysis functions (run_deseq2, run_diff_methylation, etc.) are inlined below.
#
# Required Data Objects:
#   - BRCA_CL_DNAm: Cell line DNA methylation data
#   - TCGA_BRCA_DNAm: TCGA DNA methylation data
#   - TCGA_BRCA_FC: TCGA RNA-seq count matrix
#   - TCGA_BRCA_Log2CPM: TCGA RNA-seq log2CPM matrix
#   - BRCA_CL_CTS: Cell line RNA-seq count matrix
#   - BRCA_CL_EXP: Cell line RNA-seq expression matrix
#   - CL_Annots: Cell line annotations
#   - TCGA_Annots: TCGA tumor annotations
#
# Output:
#   - dnam_dir: output path (assigned to .GlobalEnv)
#   - DMI summary statistics (DMI_Summary.tsv)
#   - SupFig11 DNAm visualizations (tumor_pam50, cl_pam50, tumor_cl_histology_LumA)
#   - Fig4A-D DNAm alteration plots (tissue_dmi, dnma_lfc_plt, fig4c_ht, fig4e)
#   - TCGA_DNAm_Regulated_Genes.csv
#
# Author: Osama Shiraz Shah
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: Load Configuration and Dependencies
# ------------------------------------------------------------------------------

# Check for config and load if needed
if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) {
    source("config.R")
  } else {
    source("../config.R")
  }
}

# Load helper functions if not already loaded
if (!exists("annot_cols", envir = .GlobalEnv)) {
  source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))
}

# ------------------------------------------------------------------------------
# Differential analysis functions (inlined from former 22_Identify_Differential_RNA_DNAm.R)
# ------------------------------------------------------------------------------
# ==============================================================================
# Utility Functions: HM450K ProbeSet Loading
# ==============================================================================

# ==============================================================================
# FUNCTION: Load and Simplify HM450K ProbeSet
# ==============================================================================
#' Load and simplify HM450K probe annotation data frame
#' 
#' Loads the HM450K probe annotation file and simplifies complex region
#' annotations for easier downstream analysis.
#' 
#' @param probe_file Path to Rdata file containing HM450K_ProbeSet object
#' @return HM450K_ProbeSet data frame with simplified region annotations
#' 
#' @details
#'   - Checks if HM450K_ProbeSet already exists in global environment
#'   - Simplifies multi-region annotations to "multiple"
#'   - Maps specific region combinations to individual labels
load_DNAm_probeset <- function(probe_file) {
  # Load HM450K_ProbeSet from file
  load(probe_file)  # loads `HM450K_ProbeSet`
  
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
  HM450K_ProbeSet$region_simple[
    HM450K_ProbeSet$region == "exon-TES_1k"
  ] <- "TES_1k"
  
  HM450K_ProbeSet$region_simple[
    HM450K_ProbeSet$region == "exon-threeUTRs"
  ] <- "threeUTRs"
  
  HM450K_ProbeSet$region_simple[
    HM450K_ProbeSet$region == "exon-TSS_1k"
  ] <- "TSS_1k"
  
  HM450K_ProbeSet$region_simple[
    HM450K_ProbeSet$region == "exon-fiveUTRs"
  ] <- "fiveUTRs"
  
  return(HM450K_ProbeSet)
}

# ------------------------------------------------------------------------------
# Utility Functions
# ------------------------------------------------------------------------------

# ==============================================================================
# FUNCTION: Beta to M-value Conversion
# ==============================================================================
#' Convert DNA methylation beta values to M-values
#' 
#' @param beta Vector or matrix of beta values (0-1 range)
#' @param offset Small offset to prevent log(0) (default: 1e-6)
#' @return M-values (log2 transformed)
#' @details
#'   M-values are logit-transformed beta values, providing better statistical
#'   properties for differential analysis. Formula: log2(beta / (1 - beta))
beta_to_M <- function(beta, offset = 1e-6) {
  # Clip beta values to valid range [offset, 1-offset]
  beta <- pmin(pmax(beta, offset), 1 - offset)
  # Convert to M-values: log2(beta / (1 - beta))
  log2(beta / (1 - beta))
}

# ==============================================================================
# FUNCTION: Run DESeq2 Differential Expression Analysis
# ==============================================================================
#' Run DESeq2 differential expression analysis
#' 
#' Performs differential gene expression analysis using DESeq2 with optional
#' covariate adjustment (e.g., tumor purity).
#' 
#' @param count_matrix Raw counts matrix (rows = genes, columns = samples)
#' @param group_labels Vector of group labels (length = ncol(count_matrix))
#' @param covariate Optional numeric vector for covariate adjustment (e.g., purity)
#' @param alpha Adjusted p-value cutoff (default: 0.05)
#' @param lfc_threshold Log2 fold-change threshold for DEGs (default: 1)
#' @return List with:
#'   - all_results: Full DESeq2 results data frame
#'   - degs: Significant differentially expressed genes
#'   - dds: DESeqDataSet object
#'   - coef_used: Coefficient name used for contrast
#' 
#' @details
#'   - Filters low-count genes (rowSum >= 10)
#'   - Uses apeglm shrinkage for log2FC when nlevels == 2
#'   - Returns results ordered by adjusted p-value
run_deseq2 <- function(count_matrix, group_labels, 
                       covariate = NULL, 
                       alpha = 0.05, lfc_threshold = 1) {
  
  # Validate DESeq2 package availability
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Please install DESeq2: BiocManager::install('DESeq2')")
  }
  
  # Ensure group_labels is a factor
  group_labels <- factor(group_labels)
  
  # ----------------------------------------------------------------------------
  # Step 1: Handle covariate (e.g., tumor purity) - drop samples with NA
  # ----------------------------------------------------------------------------
  if (!is.null(covariate)) {
    if (length(covariate) != ncol(count_matrix)) {
      stop("Length of 'covariate' must match number of samples (columns in count_matrix).")
    }
    
    na_samples <- is.na(covariate)
    if (any(na_samples)) {
      message("Dropping samples with NA in covariate: ", 
              paste(colnames(count_matrix)[na_samples], collapse = ", "))
      count_matrix <- count_matrix[, !na_samples, drop = FALSE]
      group_labels <- group_labels[!na_samples]
      covariate <- covariate[!na_samples]
    }
  }
  
  # ----------------------------------------------------------------------------
  # Step 2: Build sample metadata and design formula
  # ----------------------------------------------------------------------------
  if (is.null(covariate)) {
    col_data <- data.frame(condition = group_labels, row.names = colnames(count_matrix))
    design_formula <- ~ condition
  } else {
    col_data <- data.frame(condition = group_labels, covariate = covariate,
                           row.names = colnames(count_matrix))
    design_formula <- ~ covariate + condition
  }
  
  message("Design formula: ", deparse(design_formula))
  
  # ----------------------------------------------------------------------------
  # Step 3: Create DESeqDataSet and filter low-count genes
  # ----------------------------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = col_data,
                                design = design_formula)
  
  # Filter out genes with very low counts (rowSum >= 10)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  message("Genes kept after filtering: ", nrow(dds))
  
  # ----------------------------------------------------------------------------
  # Step 4: Run DESeq2 pipeline
  # ----------------------------------------------------------------------------
  dds <- DESeq(dds)
  
  # Get available coefficients for contrast
  coef_names <- resultsNames(dds)
  message("Available coefficients: ", paste(coef_names, collapse = ", "))
  
  # Extract condition coefficient
  condition_coef <- grep("condition", coef_names, value = TRUE)
  res <- results(dds, name = condition_coef, alpha = alpha)
  
  # Apply log2FC shrinkage for better estimates (when 2 groups)
  if (nlevels(group_labels) == 2) {
    res <- lfcShrink(dds, coef = condition_coef, res = res, type = "apeglm")
  }
  
  coef_used <- condition_coef
  
  # ----------------------------------------------------------------------------
  # Step 5: Convert to data frame and filter significant DEGs
  # ----------------------------------------------------------------------------
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Filter significant DEGs (padj < alpha & |log2FC| >= threshold)
  degs <- subset(res_df, padj < alpha & abs(log2FoldChange) >= lfc_threshold)
  degs <- degs[order(degs$padj), ]
  
  return(list(
    all_results = res_df,
    degs = degs,
    dds = dds,
    coef_used = coef_used
  ))
}

# ==============================================================================
# FUNCTION: Correlate Probe Methylation and Gene Expression
# ==============================================================================
#' Correlate DNA methylation probe values with gene expression
#' 
#' Computes correlation between DNAm beta values (probes) and gene expression
#' levels for probe-gene pairs. Used to identify probes whose methylation
#' correlates with expression changes.
#' 
#' @param probe_info Data frame with columns "Probe" and "gene"
#' @param beta_matrix Matrix of probe beta values (rows = probes, columns = samples)
#' @param expression_matrix Matrix of gene expression (rows = genes, columns = samples)
#' @param method Correlation method: "spearman" (default) or "pearson"
#' @return Data frame with correlation coefficient and p-value per probe-gene pair
#' 
#' @details
#'   - Matches probes to genes based on probe_info
#'   - Requires at least 3 samples with non-NA values
#'   - Suppresses warnings from cor.test for cleaner output
correlate_probe_expression <- function(probe_info, beta_matrix, expression_matrix, method = "spearman") {
  
  # Match probes to genes and check availability
  common_probes <- intersect(probe_info$Probe, rownames(beta_matrix))
  probe_info <- subset(probe_info, Probe %in% common_probes)
  
  # Filter expression genes
  common_genes <- intersect(probe_info$gene, rownames(expression_matrix))
  probe_info <- subset(probe_info, gene %in% common_genes)
  
  # Initialize result holder
  cor_results <- lapply(seq_len(nrow(probe_info)), function(i) {
    probe <- probe_info$Probe[i]
    gene <- probe_info$gene[i]
    
    # Extract vectors
    meth_vals <- as.numeric(beta_matrix[probe, ])
    expr_vals <- as.numeric(expression_matrix[gene, ])
    
    # Ensure matching samples and no NA
    if (length(meth_vals) != length(expr_vals)) return(NULL)
    df <- na.omit(data.frame(meth = meth_vals, expr = expr_vals))
    if (nrow(df) < 3) return(NULL)
    
    # Correlation with suppressed warnings
    suppressWarnings({
      cor_test <- cor.test(df$meth, df$expr, method = method, use = "complete.obs")
    })
    
    return(data.frame(
      Probe = probe,
      gene = gene,
      correlation = cor_test$estimate,
      pvalue = cor_test$p.value,
      method = method
    ))
  })
  
  # Combine and return
  do.call(rbind, cor_results)
}

# ==============================================================================
# FUNCTION: Differential Probe Methylation Analysis (Limma)
# ==============================================================================
#' Run Limma differential methylation analysis
#' 
#' Performs differential DNA methylation analysis using limma with optional
#' covariate adjustment (e.g., tumor purity). Converts beta to M-values for
#' better statistical properties.
#' 
#' @param beta_matrix Matrix of beta values (rows = probes, cols = samples)
#' @param group_labels Vector of group labels (length = ncol(beta_matrix))
#' @param purity Optional numeric vector for purity covariate
#' @param p_cutoff Adjusted p-value cutoff (default: 0.05)
#' @param logFC_cut Log2 fold-change threshold (default: 0.5)
#' @param use_mvalues If TRUE, convert beta values to M-values for analysis (default: TRUE)
#' @return List with:
#'   - all_results: Full limma results data frame
#'   - significant_probes: Significant differentially methylated probes
#'   - design: Design matrix used
#'   - fit: Fitted limma model object
#' 
#' @details
#'   - Converts beta to M-values for better statistical properties
#'   - Uses empirical Bayes moderation (eBayes)
#'   - Calculates delta-beta (difference in mean beta) for interpretation
run_diff_methylation <- function(beta_matrix, group_labels, 
                                 purity = NULL,
                                 p_cutoff = 0.05, logFC_cut = 0.5,
                                 use_mvalues = TRUE) {
  
  # Validate limma package availability
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Please install the limma package: BiocManager::install('limma')")
  }
  
  # Ensure group labels are a factor
  group_labels <- factor(group_labels)
  
  # Save original beta matrix for delta-beta calculation
  beta_orig <- beta_matrix
  
  # ----------------------------------------------------------------------------
  # Step 1: Convert Beta values to M-values (if needed)
  # ----------------------------------------------------------------------------
  message("Converting to M-values...")
  m_matrix <- beta_matrix
  if (use_mvalues) {
    m_matrix[m_matrix <= 0] <- 1e-6
    m_matrix[m_matrix >= 1] <- 1 - 1e-6
    m_matrix <- log2(m_matrix / (1 - m_matrix))
  }
  
  # ----------------------------------------------------------------------------
  # Step 2: Build design matrix (with optional purity covariate)
  # ----------------------------------------------------------------------------
  if (is.null(purity)) {
    design <- model.matrix(~ 0 + group_labels)
    colnames(design) <- levels(group_labels)
  } else {
    if (length(purity) != ncol(beta_matrix)) {
      stop("Length of 'purity' must match number of samples (columns in beta_matrix).")
    }
    design <- model.matrix(~ 0 + group_labels + purity)
    colnames(design) <- make.names(colnames(design))
  }
  
  # ----------------------------------------------------------------------------
  # Step 3: Fit linear model
  # ----------------------------------------------------------------------------
  message("Performing linear model fit...")
  fit <- lmFit(m_matrix, design)
  
  # ----------------------------------------------------------------------------
  # Step 4: Set up contrasts (assumes exactly 2 levels in group_labels)
  # ----------------------------------------------------------------------------
  if (nlevels(group_labels) != 2) {
    stop("This function currently supports only two groups.")
  }
  
  contrast_name <- paste0("group_labels", levels(group_labels)[2], "-", "group_labels", levels(group_labels)[1])
  
  if (is.null(purity)) {
    contrast_matrix <- makeContrasts(
      contrasts = paste0(levels(group_labels)[2], "-", levels(group_labels)[1]),
      levels = design
    )
  } else {
    # Match the encoded group variable names from design
    contrast_matrix <- makeContrasts(
      contrasts = paste0("group_labels", levels(group_labels)[2], "-", "group_labels", levels(group_labels)[1]),
      levels = design
    )
  }
  
  # Apply contrasts and empirical Bayes moderation
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # ----------------------------------------------------------------------------
  # Step 5: Extract results
  # ----------------------------------------------------------------------------
  res <- topTable(fit2, adjust.method = "BH", number = Inf)
  res$Probe <- rownames(res)
  
  # Calculate delta-beta (using original beta values) for biological interpretation
  beta_diff <- rowMeans(beta_orig[, group_labels == levels(group_labels)[2], drop = FALSE], na.rm = TRUE) -
    rowMeans(beta_orig[, group_labels == levels(group_labels)[1], drop = FALSE], na.rm = TRUE)
  res$delta_beta <- beta_diff[res$Probe]
  
  # ----------------------------------------------------------------------------
  # Step 6: Filter significant probes
  # ----------------------------------------------------------------------------
  sig_res <- subset(res, adj.P.Val < p_cutoff & abs(logFC) >= logFC_cut)
  sig_res <- sig_res[order(sig_res$adj.P.Val), ]
  
  return(list(
    all_results = res,
    significant_probes = sig_res,
    design = design,
    fit = fit2
  ))
}

# ==============================================================================
# Wrapper Functions for Specific Analyses
# ==============================================================================

# ==============================================================================
# FUNCTION: Run TCGA DESeq2 Analysis
# ==============================================================================
#' Run DESeq2 differential expression analysis on TCGA LumA ILC vs NST tumors
#' 
#' Wrapper function that sets up TCGA-specific sample selection and runs DESeq2
#' with purity correction. Results are automatically saved to files.
#' 
#' @param TCGA_BRCA_FC TCGA RNA-seq count matrix
#' @param TCGA_Annots TCGA tumor annotations
#' @param output_file Path to save results (.Rdata file)
#' @param csv_file Path to save DEGs CSV file
#' @param alpha Adjusted p-value threshold (default: 0.05)
#' @param lfc_threshold Log2 fold-change threshold (default: 1)
#' @return List with DESeq2 results (same structure as run_deseq2 output)
#' 
#' @details
#'   - Filters to LumA PAM50 subtype
#'   - Requires non-NA tumor purity (TumorPurity_CPE)
#'   - Uses purity as covariate in DESeq2 model
run_tcga_deseq_analysis <- function(TCGA_BRCA_FC, TCGA_Annots, 
                                     output_file, csv_file,
                                     alpha = 0.05, lfc_threshold = 1) {
  # ----------------------------------------------------------------------------
  # Step 1: Define sample groups (LumA ILC vs NST with purity data)
  # ----------------------------------------------------------------------------
  ILC_ids <- intersect(
    subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("ILC"))$Case.ID,
    colnames(TCGA_BRCA_FC)
  )
  NST_ids <- intersect(
    subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("NST"))$Case.ID,
    colnames(TCGA_BRCA_FC)
  )
  
  if (length(ILC_ids) == 0 || length(NST_ids) == 0) {
    stop("Insufficient samples for TCGA DESeq2 analysis. ILC: ", length(ILC_ids), ", NST: ", length(NST_ids))
  }
  
  message("  ILC samples: ", length(ILC_ids), ", NST samples: ", length(NST_ids))
  
  # ----------------------------------------------------------------------------
  # Step 2: Prepare data matrices and labels
  # ----------------------------------------------------------------------------
  TCGA_de_mat <- TCGA_BRCA_FC[, c(ILC_ids, NST_ids)]
  group_labels <- factor(c(rep("ILC", length(ILC_ids)), rep("NST", length(NST_ids))),
                         levels = c("NST", "ILC"))
  purity_covar <- TCGA_Annots[TCGA_Annots$Case.ID %in% c(ILC_ids, NST_ids), "TumorPurity_CPE"]
  
  # ----------------------------------------------------------------------------
  # Step 3: Run DESeq2 with purity correction
  # ----------------------------------------------------------------------------
  TCGA_DESeq_Results <- run_deseq2(TCGA_de_mat, group_labels, purity_covar, 
                                    alpha = alpha, lfc_threshold = lfc_threshold)
  TCGA_DESeq_Results$Run_Notes <- "TCGA BRCA LumA ILC vs NST Purity Corrected"
  
  # ----------------------------------------------------------------------------
  # Step 4: Save results
  # ----------------------------------------------------------------------------
  ensure_dir(dirname(output_file))
  save(TCGA_DESeq_Results, file = output_file)
  message("  ✓ Saved TCGA DESeq2 results to: ", output_file)
  
  if (!is.null(csv_file)) {
    ensure_dir(dirname(csv_file))
    write.csv(TCGA_DESeq_Results$degs, csv_file, row.names = FALSE)
    message("  ✓ Saved TCGA DEGs to: ", csv_file)
  }
  
  return(TCGA_DESeq_Results)
}

# ==============================================================================
# FUNCTION: Run TCGA Limma Analysis
# ==============================================================================
#' Run Limma differential methylation analysis on TCGA LumA ILC vs NST tumors
#' 
#' Wrapper function that sets up TCGA-specific sample selection and runs Limma
#' with purity correction. Adds gene annotations and RNA-DNAm correlations.
#' 
#' @param TCGA_BRCA_DNAm TCGA DNA methylation beta matrix
#' @param TCGA_Annots TCGA tumor annotations
#' @param HM450K_ProbeSet Probe annotation data frame
#' @param allProbes Vector of probe IDs to analyze
#' @param output_file Path to save results (.Rdata file)
#' @param csv_file Path to save significant probes CSV file
#' @param p_cutoff Adjusted p-value threshold (default: 0.05)
#' @param logFC_cut Log2 fold-change threshold (default: 0.5)
#' @param use_mvalues Whether to use M-values for analysis (default: TRUE)
#' @return List with Limma results (same structure as run_diff_methylation output)
#' 
#' @details
#'   - Filters to LumA PAM50 subtype with purity data
#'   - Adds gene names and probe regions from HM450K_ProbeSet
#'   - Computes RNA-DNAm correlations if expression data available
run_tcga_limma_analysis <- function(TCGA_BRCA_DNAm, TCGA_Annots, HM450K_ProbeSet,
                                    allProbes, output_file, csv_file,
                                    p_cutoff = 0.05, logFC_cut = 0.5, use_mvalues = TRUE) {
  # ----------------------------------------------------------------------------
  # Step 1: Define sample groups (LumA ILC vs NST with purity data)
  # ----------------------------------------------------------------------------
  ILC_ids <- intersect(
    subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("ILC"))$Case.ID,
    colnames(TCGA_BRCA_DNAm)
  )
  NST_ids <- intersect(
    subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("NST"))$Case.ID,
    colnames(TCGA_BRCA_DNAm)
  )
  
  if (length(ILC_ids) == 0 || length(NST_ids) == 0) {
    stop("Insufficient samples for TCGA Limma analysis. ILC: ", length(ILC_ids), ", NST: ", length(NST_ids))
  }
  
  message("  ILC samples: ", length(ILC_ids), ", NST samples: ", length(NST_ids))
  
  # ----------------------------------------------------------------------------
  # Step 2: Prepare data matrices and labels
  # ----------------------------------------------------------------------------
  TCGA_limma_mat <- TCGA_BRCA_DNAm[allProbes, c(ILC_ids, NST_ids)]
  group_labels <- factor(c(rep("ILC", length(ILC_ids)), rep("NST", length(NST_ids))),
                         levels = c("NST", "ILC"))
  purity_covar <- TCGA_Annots[TCGA_Annots$Case.ID %in% c(ILC_ids, NST_ids), "TumorPurity_CPE"]
  
  # ----------------------------------------------------------------------------
  # Step 3: Run differential methylation analysis
  # ----------------------------------------------------------------------------
  TCGA_DPM_Results <- run_diff_methylation(TCGA_limma_mat, group_labels, purity_covar,
                                           p_cutoff = p_cutoff, use_mvalues = use_mvalues, logFC_cut = logFC_cut)
  
  # ----------------------------------------------------------------------------
  # Step 4: Add gene and region annotations from HM450K_ProbeSet
  # ----------------------------------------------------------------------------
  TCGA_DPM_Results$all_results$gene <- HM450K_ProbeSet[TCGA_DPM_Results$all_results$Probe, "transcript_geneName_noENS_NA"]
  TCGA_DPM_Results$significant_probes$gene <- HM450K_ProbeSet[TCGA_DPM_Results$significant_probes$Probe, "transcript_geneName_noENS_NA"]
  TCGA_DPM_Results$all_results$Probe_region <- HM450K_ProbeSet[TCGA_DPM_Results$all_results$Probe, "region_simple"]
  TCGA_DPM_Results$significant_probes$Probe_region <- HM450K_ProbeSet[TCGA_DPM_Results$significant_probes$Probe, "region_simple"]
  TCGA_DPM_Results$Run_Notes <- "TCGA BRCA LumA ILC vs NST"
  
  # ----------------------------------------------------------------------------
  # Step 5: Add RNA-DNAm correlation results (if expression data available)
  # ----------------------------------------------------------------------------
  if (exists("TCGA_BRCA_Log2CPM", envir = .GlobalEnv) && exists("correlate_probe_expression", envir = .GlobalEnv)) {
    commonIDs <- intersect(colnames(TCGA_limma_mat), colnames(TCGA_BRCA_Log2CPM))
    if (length(commonIDs) > 0 && exists("beta_to_M", envir = .GlobalEnv)) {
      message("  Computing RNA-DNAm correlations...")
      TCGA_RNA_DNAm_Correlation_Res <- correlate_probe_expression(
        probe_info = TCGA_DPM_Results$significant_probes,
        beta_matrix = beta_to_M(TCGA_BRCA_DNAm[TCGA_DPM_Results$significant_probes$Probe, commonIDs]),
        expression_matrix = TCGA_BRCA_Log2CPM[, commonIDs]
      )
      rownames(TCGA_RNA_DNAm_Correlation_Res) <- TCGA_RNA_DNAm_Correlation_Res$Probe
      TCGA_DPM_Results$significant_probes$RNA_cor <- TCGA_RNA_DNAm_Correlation_Res[TCGA_DPM_Results$significant_probes$Probe, "correlation"]
      TCGA_DPM_Results$significant_probes$RNA_cor_pvalue <- TCGA_RNA_DNAm_Correlation_Res[TCGA_DPM_Results$significant_probes$Probe, "pvalue"]
      
      # Rank by lowest p-value, then highest absolute correlation
      TCGA_RNA_DNAm_Correlation_Res <- TCGA_RNA_DNAm_Correlation_Res[order(
        TCGA_RNA_DNAm_Correlation_Res$gene,
        TCGA_RNA_DNAm_Correlation_Res$pvalue,
        -abs(TCGA_RNA_DNAm_Correlation_Res$correlation)
      ), ]
      
      # Select top entry per gene (best correlated probe per gene)
      TCGA_RNA_DNAm_Correlation_Res_unique <- TCGA_RNA_DNAm_Correlation_Res[!duplicated(TCGA_RNA_DNAm_Correlation_Res$gene), ]
      TCGA_DPM_Results$significant_probes[, "top_cor_probe"] <- ifelse(
        TCGA_DPM_Results$significant_probes$Probe %in% TCGA_RNA_DNAm_Correlation_Res_unique$Probe, TRUE, FALSE
      )
      message("  ✓ Added RNA-DNAm correlations")
    }
  }
  
  # ----------------------------------------------------------------------------
  # Step 6: Save results
  # ----------------------------------------------------------------------------
  ensure_dir(dirname(output_file))
  save(TCGA_DPM_Results, file = output_file)
  message("  ✓ Saved TCGA Limma results to: ", output_file)
  
  if (!is.null(csv_file)) {
    ensure_dir(dirname(csv_file))
    write.csv(TCGA_DPM_Results$significant_probes, csv_file, row.names = FALSE)
    message("  ✓ Saved TCGA significant probes to: ", csv_file)
  }
  
  return(TCGA_DPM_Results)
}

# ==============================================================================
# FUNCTION: Run Cell Line DESeq2 Analysis
# ==============================================================================
#' Run DESeq2 differential expression analysis on ICLE non-basal cell lines
#' 
#' Wrapper function that sets up cell line-specific sample selection (non-basal,
#' ILC/ILC-like vs NST) and runs DESeq2. Results are automatically saved.
#' 
#' @param BRCA_CL_CTS Cell line RNA-seq count matrix
#' @param CL_Annots Cell line annotations
#' @param cellsID Vector of cell line IDs to include
#' @param output_file Path to save results (.Rdata file)
#' @param csv_file Path to save DEGs CSV file
#' @param alpha Adjusted p-value threshold (default: 0.3, more lenient for cell lines)
#' @param lfc_threshold Log2 fold-change threshold (default: 1)
#' @return List with DESeq2 results (same structure as run_deseq2 output)
#' 
#' @details
#'   - Filters to non-basal cell lines (excludes Basal subtype)
#'   - Combines ILC and ILC-like into single group
#'   - Excludes DCIS cell lines (HCC202, HCC1500)
run_cl_deseq_analysis <- function(BRCA_CL_CTS, CL_Annots, cellsID,
                                  output_file, csv_file,
                                  alpha = 0.3, lfc_threshold = 1) {
  # ----------------------------------------------------------------------------
  # Step 1: Define sample groups (non-basal ILC/ILC-like vs NST)
  # ----------------------------------------------------------------------------
  nonBasal_cells <- subset(CL_Annots[cellsID, ], TopCall %notin% c("Basal") & overlapWCCLE != "Y")$Name
  ILC_cl_ids <- intersect(
    subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("ILC"))$Name,
    colnames(BRCA_CL_CTS)
  )
  ILC_like_cl_ids <- intersect(
    subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("ILC-like"))$Name,
    colnames(BRCA_CL_CTS)
  )
  NST_cl_ids <- intersect(
    subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("NST") & Sample != c("HCC202", "HCC1500"))$Name,
    colnames(BRCA_CL_CTS)
  )
  
  if (length(c(ILC_cl_ids, ILC_like_cl_ids)) == 0 || length(NST_cl_ids) == 0) {
    stop("Insufficient samples for Cell Line DESeq2 analysis. ILC: ", 
         length(c(ILC_cl_ids, ILC_like_cl_ids)), ", NST: ", length(NST_cl_ids))
  }
  
  message("  ILC samples: ", length(c(ILC_cl_ids, ILC_like_cl_ids)), ", NST samples: ", length(NST_cl_ids))
  
  # ----------------------------------------------------------------------------
  # Step 2: Prepare data matrices and labels
  # ----------------------------------------------------------------------------
  de_mat <- as.matrix(BRCA_CL_CTS[, c(c(ILC_cl_ids, ILC_like_cl_ids), NST_cl_ids)])
  mode(de_mat) <- "integer"  # Ensure integer counts for DESeq2
  group_labels <- factor(c(rep("ILC", length(c(ILC_cl_ids, ILC_like_cl_ids))), rep("NST", length(NST_cl_ids))),
                         levels = c("NST", "ILC"))
  
  # ----------------------------------------------------------------------------
  # Step 3: Run DESeq2 (no covariate for cell lines)
  # ----------------------------------------------------------------------------
  BRCAL_CL_DESeq_Results <- run_deseq2(de_mat, group_labels, alpha = alpha, lfc_threshold = lfc_threshold)
  BRCAL_CL_DESeq_Results$Run_Notes <- "BRCA Non-Basal ICLE vs NST Cell lines"
  
  # ----------------------------------------------------------------------------
  # Step 4: Save results
  # ----------------------------------------------------------------------------
  ensure_dir(dirname(output_file))
  save(BRCAL_CL_DESeq_Results, file = output_file)
  message("  ✓ Saved Cell Line DESeq2 results to: ", output_file)
  
  if (!is.null(csv_file)) {
    ensure_dir(dirname(csv_file))
    write.csv(BRCAL_CL_DESeq_Results$degs, csv_file, row.names = FALSE)
    message("  ✓ Saved Cell Line DEGs to: ", csv_file)
  }
  
  return(BRCAL_CL_DESeq_Results)
}

# ==============================================================================
# FUNCTION: Run Cell Line Limma Analysis
# ==============================================================================
#' Run Limma differential methylation analysis on ICLE non-basal cell lines
#' 
#' Wrapper function that sets up cell line-specific sample selection and runs
#' Limma. Adds gene annotations and RNA-DNAm correlations.
#' 
#' @param BRCA_CL_DNAm Cell line DNA methylation beta matrix
#' @param CL_Annots Cell line annotations
#' @param HM450K_ProbeSet Probe annotation data frame
#' @param allProbes Vector of probe IDs to analyze
#' @param cellsID Vector of cell line IDs to include
#' @param output_file Path to save results (.Rdata file)
#' @param csv_file Path to save significant probes CSV file
#' @param p_cutoff Adjusted p-value threshold (default: 0.3, more lenient for cell lines)
#' @param logFC_cut Log2 fold-change threshold (default: 0.5)
#' @param use_mvalues Whether to use M-values for analysis (default: TRUE)
#' @return List with Limma results (same structure as run_diff_methylation output)
#' 
#' @details
#'   - Filters to non-basal cell lines
#'   - Removes probes with >5 missing values across samples
#'   - Computes RNA-DNAm correlations if expression data available
run_cl_limma_analysis <- function(BRCA_CL_DNAm, CL_Annots, HM450K_ProbeSet,
                                   allProbes, cellsID, output_file, csv_file,
                                   p_cutoff = 0.3, logFC_cut = 0.5, use_mvalues = TRUE) {
  # ----------------------------------------------------------------------------
  # Step 1: Define sample groups (non-basal ILC/ILC-like vs NST)
  # ----------------------------------------------------------------------------
  nonBasal_cells <- subset(CL_Annots[cellsID, ], TopCall %notin% c("Basal") & overlapWCCLE != "Y")$Name
  ILC_cl_ids <- intersect(
    subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("ILC"))$Name,
    colnames(BRCA_CL_DNAm)
  )
  ILC_like_cl_ids <- intersect(
    subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("ILC-like"))$Name,
    colnames(BRCA_CL_DNAm)
  )
  NST_cl_ids <- intersect(
    subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("NST") & Sample != c("HCC202", "HCC1500"))$Name,
    colnames(BRCA_CL_DNAm)
  )
  
  if (length(c(ILC_cl_ids, ILC_like_cl_ids)) == 0 || length(NST_cl_ids) == 0) {
    stop("Insufficient samples for Cell Line Limma analysis. ILC: ", 
         length(c(ILC_cl_ids, ILC_like_cl_ids)), ", NST: ", length(NST_cl_ids))
  }
  
  message("  ILC samples: ", length(c(ILC_cl_ids, ILC_like_cl_ids)), ", NST samples: ", length(NST_cl_ids))
  
  # ----------------------------------------------------------------------------
  # Step 2: Prepare data matrices and filter probes with too many NAs
  # ----------------------------------------------------------------------------
  limma_mat <- BRCA_CL_DNAm[allProbes, c(ILC_cl_ids, ILC_like_cl_ids, NST_cl_ids)]
  limma_mat <- limma_mat[rowSums(is.na(limma_mat)) <= 5, ]  # Keep probes with <=5 NAs
  group_labels <- factor(c(rep("ILC", length(c(ILC_cl_ids, ILC_like_cl_ids))), rep("NST", length(NST_cl_ids))),
                         levels = c("NST", "ILC"))
  
  # ----------------------------------------------------------------------------
  # Step 3: Run differential methylation analysis (no covariate for cell lines)
  # ----------------------------------------------------------------------------
  BRCA_CL_DPM_Results <- run_diff_methylation(limma_mat, group_labels, use_mvalues = use_mvalues, 
                                              logFC_cut = logFC_cut, p_cutoff = p_cutoff)
  
  # ----------------------------------------------------------------------------
  # Step 4: Add gene and region annotations from HM450K_ProbeSet
  # ----------------------------------------------------------------------------
  BRCA_CL_DPM_Results$all_results$gene <- HM450K_ProbeSet[BRCA_CL_DPM_Results$all_results$Probe, "transcript_geneName_noENS_NA"]
  BRCA_CL_DPM_Results$significant_probes$gene <- HM450K_ProbeSet[BRCA_CL_DPM_Results$significant_probes$Probe, "transcript_geneName_noENS_NA"]
  BRCA_CL_DPM_Results$all_results$Probe_region <- HM450K_ProbeSet[BRCA_CL_DPM_Results$all_results$Probe, "region_simple"]
  BRCA_CL_DPM_Results$significant_probes$Probe_region <- HM450K_ProbeSet[BRCA_CL_DPM_Results$significant_probes$Probe, "region_simple"]
  BRCA_CL_DPM_Results$Run_Notes <- "BRCA Non-Basal ICLE vs NST Cell lines"
  
  # ----------------------------------------------------------------------------
  # Step 5: Add RNA-DNAm correlation results (if expression data available)
  # ----------------------------------------------------------------------------
  if (exists("BRCA_CL_EXP", envir = .GlobalEnv) && exists("correlate_probe_expression", envir = .GlobalEnv)) {
    commonIDs <- intersect(colnames(limma_mat), colnames(BRCA_CL_EXP))
    if (length(commonIDs) > 0 && exists("beta_to_M", envir = .GlobalEnv)) {
      message("  Computing RNA-DNAm correlations...")
      BRCA_CL_RNA_DNAm_Correlation_Res <- correlate_probe_expression(
        probe_info = BRCA_CL_DPM_Results$significant_probes,
        beta_matrix = beta_to_M(BRCA_CL_DNAm[BRCA_CL_DPM_Results$significant_probes$Probe, commonIDs]),
        expression_matrix = BRCA_CL_EXP[, commonIDs]
      )
      rownames(BRCA_CL_RNA_DNAm_Correlation_Res) <- BRCA_CL_RNA_DNAm_Correlation_Res$Probe
      BRCA_CL_DPM_Results$significant_probes$RNA_cor <- BRCA_CL_RNA_DNAm_Correlation_Res[BRCA_CL_DPM_Results$significant_probes$Probe, "correlation"]
      BRCA_CL_DPM_Results$significant_probes$RNA_cor_pvalue <- BRCA_CL_RNA_DNAm_Correlation_Res[BRCA_CL_DPM_Results$significant_probes$Probe, "pvalue"]
      
      # Rank by lowest p-value, then highest absolute correlation
      BRCA_CL_RNA_DNAm_Correlation_Res <- BRCA_CL_RNA_DNAm_Correlation_Res[order(
        BRCA_CL_RNA_DNAm_Correlation_Res$gene,
        BRCA_CL_RNA_DNAm_Correlation_Res$pvalue,
        -abs(BRCA_CL_RNA_DNAm_Correlation_Res$correlation)
      ), ]
      
      # Select top entry per gene (best correlated probe per gene)
      BRCA_CL_RNA_DNAm_Correlation_Res_unique <- BRCA_CL_RNA_DNAm_Correlation_Res[!duplicated(BRCA_CL_RNA_DNAm_Correlation_Res$gene), ]
      BRCA_CL_DPM_Results$significant_probes[, "top_cor_probe"] <- ifelse(
        BRCA_CL_DPM_Results$significant_probes$Probe %in% BRCA_CL_RNA_DNAm_Correlation_Res_unique$Probe, TRUE, FALSE
      )
      message("  ✓ Added RNA-DNAm correlations")
    }
  }
  
  # ----------------------------------------------------------------------------
  # Step 6: Save results
  # ----------------------------------------------------------------------------
  ensure_dir(dirname(output_file))
  save(BRCA_CL_DPM_Results, file = output_file)
  message("  ✓ Saved Cell Line Limma results to: ", output_file)
  
  if (!is.null(csv_file)) {
    ensure_dir(dirname(csv_file))
    write.csv(BRCA_CL_DPM_Results$significant_probes, csv_file, row.names = FALSE)
    message("  ✓ Saved Cell Line significant probes to: ", csv_file)
  }
  
  return(BRCA_CL_DPM_Results)
}

ensure_dir(DIRS$results_sub$dna_methylation)
dnam_dir <- DIRS$results_sub$dna_methylation
assign("dnam_dir", dnam_dir, envir = .GlobalEnv)

# Load required packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(readxl)
  library(dplyr)
  library(ComplexHeatmap)
  library(sesame)
  library(patchwork)
  library(tidytext)
  library(forcats)
})

# ------------------------------------------------------------------------------
# SECTION 1: Load Data
# ------------------------------------------------------------------------------

message("=", strrep("=", 78))
message("Section 1: Loading Data")
message("=", strrep("=", 78))

# Load Oncovar annotations
message("Loading Oncovar annotations...")
oncoVar <- read.delim(FILES$oncovar)

# Load and simplify HM450K probe set
message("Loading HM450K probe set...")
HM450K_ProbeSet <- load_DNAm_probeset(FILES$hm450k_probeset)

# Load RNA-seq data
message("Loading RNA-seq data...")
load(FILES$rnaseq_cts)
load(FILES$rnaseq_data)

# Load cell line DNA methylation data
message("Loading cell line DNA methylation data...")
load(FILES$icle_dnam)
colnames(BRCA_CL_DNAm) <- gsub(x = colnames(BRCA_CL_DNAm), pattern = "-S", replacement = "-C")

# Filter to cell lines with annotations
cellsID <- intersect(colnames(BRCA_CL_DNAm), CL_Annots$Name)
BRCA_CL_DNAm <- BRCA_CL_DNAm[, cellsID]

# Check if TCGA DNA methylation data exists, load if needed
if (!exists("TCGA_BRCA_DNAm", envir = .GlobalEnv)) {
  message("Loading TCGA DNA methylation data...")
  # Try to load from external data
  if (exists("load_external_data", mode = "function")) {
    source(file.path(DIRS$scripts$helpers, "Data_Loading", "08_load_external_data.R"))
    load_tcga_data()
  }
  
  # If still not found, try loading directly from file
  if (!exists("TCGA_BRCA_DNAm", envir = .GlobalEnv)) {
    tcga_dnam_file <- file.path(DIRS$external$tcga, "DNAm", "TCGA_BRCA_DNAm.Rdata")
    if (file.exists(tcga_dnam_file)) {
      load(tcga_dnam_file, envir = .GlobalEnv)
      message("  ✓ Loaded TCGA_BRCA_DNAm from file")
    } else {
      stop("TCGA_BRCA_DNAm not found. Please ensure external data is loaded.")
    }
  }
} else {
  message("TCGA DNA methylation data already loaded")
}

# Identify common probes across datasets
allProbes <- intersect(intersect(HM450K_ProbeSet$probeID, rownames(BRCA_CL_DNAm)), rownames(TCGA_BRCA_DNAm))
message("  Common probes across datasets: ", length(allProbes))

# ------------------------------------------------------------------------------
# SECTION 2: Calculate DMI (DNA Methylation Index)
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("Section 2: Calculating DMI (DNA Methylation Index)")
message("=", strrep("=", 78))

# ------------------------------------------------------------------------------
# Step 1: Calculate median DNA methylation in normal samples (identified by "-11")
# ------------------------------------------------------------------------------
message("  Step 1: Calculating median DNA methylation in normal samples...")

# Extract normal tissue samples
beta_normal_matrix <- TCGA_BRCA_DNAm[allProbes, grepl("-11", colnames(TCGA_BRCA_DNAm))]

# Calculate probe-wise median beta values for normals
median_normal_beta <- apply(beta_normal_matrix[allProbes, ], 1, function(x) median(x, na.rm = TRUE))
message("    Calculated median beta for ", length(median_normal_beta), " probes")

# ------------------------------------------------------------------------------
# Step 2: Calculate Global Methylation Aberration (GMA) scores
# ------------------------------------------------------------------------------
message("  Step 2: Calculating Global Methylation Aberration (GMA) scores...")

# GMA scores for patient tumors (sum of absolute deviation from normal medians)
patient_GMA_scores <- colSums(abs(TCGA_BRCA_DNAm[allProbes, ] - as.numeric(median_normal_beta)), na.rm = TRUE)

# GMA scores for cell lines
cell_line_GMA_scores <- colSums(abs(BRCA_CL_DNAm[allProbes, ] - as.numeric(median_normal_beta)), na.rm = TRUE)
message("  Patient tumors: ", length(patient_GMA_scores), " samples")
message("  Cell lines: ", length(cell_line_GMA_scores), " samples")

# ------------------------------------------------------------------------------
# Step 3: Combine patient and cell line data into a single DataFrame
# ------------------------------------------------------------------------------
message("  Step 3: Combining patient and cell line data...")

DMI_df <- data.frame(
  DMI = c(as.numeric(patient_GMA_scores), as.numeric(cell_line_GMA_scores)),
  SampleID = c(names(patient_GMA_scores), cellsID),
  TissueSource = c(rep("Patient Tumors", length(names(patient_GMA_scores))), 
                   rep("Cell Lines", length(names(cell_line_GMA_scores)))),
  Histology = c(TCGA_Annots[names(patient_GMA_scores), "Final Pathology"], 
                as.character(CL_Annots[cellsID, "Histology"])), 
  PAM50 = c(TCGA_Annots[names(patient_GMA_scores), "PAM50"], 
            as.character(CL_Annots[cellsID, "TopCall"]))
)

# Remove duplicates and set row names
DMI_df <- DMI_df[!duplicated(DMI_df$SampleID), ]
rownames(DMI_df) <- DMI_df$SampleID

# Fix PAM50 labels for cell lines with multiple subtypes
idx <- grep(pattern = ";", DMI_df$PAM50)
DMI_df[idx, "PAM50"] <- CL_Annots[DMI_df[idx, "SampleID"], "PCAPAM50"]

# Update TissueSource to "Normal" for normal samples
DMI_df$TissueSource[grepl("-11", DMI_df$SampleID)] <- "Normal"
DMI_df$TissueSource[grepl("-01", DMI_df$SampleID)] <- "Patient Tumors"

# ------------------------------------------------------------------------------
# Step 4: Define group variables for plotting and subsetting
# ------------------------------------------------------------------------------
message("  Step 4: Defining group variables...")

# Group labels by histology and tissue source
DMI_df$HistologyGroup <- paste0(DMI_df$Histology, "-", DMI_df$TissueSource)

# Group labels by histology and PAM50 subtype
DMI_df$HistologySubtype <- paste0(DMI_df$Histology, "-", DMI_df$PAM50)

# Filter for relevant histology types
# DMI_df <- subset(DMI_df, Histology %notin% c("Other", "mDLC"))
DMI_df <- DMI_df[grep("-20|-06", DMI_df$SampleID, invert = T), ] # remove -20 and -06 - not primary tumors
# DMI_df <- subset(DMI_df, (TissueSource == "Patient Tumors" & !is.na(Histology)) | TissueSource != "Patient Tumors")
message("  Filtered to ", nrow(DMI_df), " samples (", paste0(names(table(DMI_df$TissueSource)), ": ", table(DMI_df$TissueSource), ". "), ")")

# ------------------------------------------------------------------------------
# Step 5: Compute DMI z-scores and normalize to [0,1]
# ------------------------------------------------------------------------------
message("  Step 5: Computing DMI z-scores...")

# Extract DMI scores for normal samples
normal_DMI_scores <- subset(DMI_df, TissueSource == "Normal")$DMI

# Z-score calculation: (DMI - median_normal) / standard deviation_normal
DMI_df$DMI_zscore <- (DMI_df$DMI - median(normal_DMI_scores)) / sqrt(var(normal_DMI_scores))

# ------------------------------------------------------------------------------
# Step 6: Remove outlier normals (z-score > ±1.5)
# ------------------------------------------------------------------------------
message("  Step 6: Removing outlier normals...")

outlier_normals <- subset(DMI_df, TissueSource == "Normal" & abs(DMI_zscore) > 1.5)$SampleID
DMI_df <- subset(DMI_df, !(SampleID %in% outlier_normals))
if (length(outlier_normals) > 0) {
  message("  Removed ", length(outlier_normals), " outlier normal samples")
}

# ------------------------------------------------------------------------------
# Step 7: Generate combined histology-type labels
# ------------------------------------------------------------------------------
message("  Step 7: Generating combined histology-type labels...")

DMI_df$Histology_TissueSource <- paste0(DMI_df$Histology, " | ", DMI_df$TissueSource)
DMI_df$Histology_TissueSource[DMI_df$Histology_TissueSource == "Normal | Normal"] <- "Normal"

rownames(DMI_df) <- DMI_df$SampleID

# Save DMI summary
write.table(DMI_df, file = file.path(dnam_dir, "DMI_Summary.tsv"), 
            quote = FALSE, row.names = FALSE, sep = "\t")
message("  ✓ Saved DMI summary to: ", file.path(dnam_dir, "DMI_Summary.tsv"))

# ------------------------------------------------------------------------------
# SECTION 3: DMI Comparisons and Visualizations (SupFig11, Fig4A)
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("Section 3: DMI Comparisons and Visualizations")
message("=", strrep("=", 78))

# ------------------------------------------------------------------------------
# SupFig 11A - DMI Patient Tumors By PAM50
# ------------------------------------------------------------------------------
message("\n========================================")
message("SupFig 11A: DMI patient tumors by PAM50")
message("========================================")

DMI_df_t <- subset(DMI_df, TissueSource == "Patient Tumors" & PAM50 != "Normal")

tumor_pam50 <- ggplot(DMI_df_t, aes(x = reorder(PAM50, DMI_zscore), y = DMI_zscore, color = PAM50)) +
  ggbeeswarm::geom_quasirandom(size = 1, alpha = 1, width = 0.4) + 
  geom_boxplot(alpha = 0.6, linewidth = 0.4, colour = "black", outlier.color = NA, width = 0.2) +
  scale_color_manual(values = c(annot_cols$PAM50), guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  xlab(" ") + ylab("DNAm Index") + theme_bw(20) + theme(axis.text.x = element_blank())

# Suppress warnings from pairwise.t.test
suppressWarnings({
  SupFig11A_pval_table <- pairwise.t.test(x = DMI_df_t$DMI, g = DMI_df_t$PAM50)$p.value
})

# ------------------------------------------------------------------------------
# SupFig 11B - DMI Cell Lines By PAM50
# ------------------------------------------------------------------------------
message("\n========================================")
message("SupFig 11B: DMI cell lines by PAM50")
message("========================================")

DMI_df_cl <- subset(DMI_df, TissueSource == "Cell Lines")
DMI_df_cl <- DMI_df_cl[!is.na(DMI_df_cl$PAM50), ]

cl_pam50 <- ggplot(DMI_df_cl, aes(x = reorder(PAM50, DMI_zscore), y = DMI_zscore, color = PAM50)) +
  ggbeeswarm::geom_quasirandom(size = 1, alpha = 1, width = 0.4) +
  geom_boxplot(alpha = 0.6, linewidth = 0.4, colour = "black", outlier.color = NA, width = 0.2) +
  scale_color_manual(values = c(annot_cols$PAM50), guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  xlab(" ") + ylab("DNAm Index") + theme_bw(20) + theme(axis.text.x = element_blank())

suppressWarnings({
  SupFig11B_pval_table <- pairwise.t.test(x = DMI_df_cl$DMI, g = DMI_df_cl$PAM50)$p.value
})

# ------------------------------------------------------------------------------
# SupFig 11C - DMI by Histology (LumA tumors and non-basal cell lines)
# ------------------------------------------------------------------------------
message("\n========================================")
message("SupFig 11C: DMI by histology")
message("========================================")

# Add RNA subtypes for cell lines
DMI_df$RNA_Subtypes <- CL_Annots[DMI_df$SampleID, "mRNA Subtypes"]

# Filter to LumA tumors and non-basal cell lines
DMI_df_luminal <- DMI_df %>% filter((TissueSource == "Patient Tumors" & PAM50 == "LumA" & Histology %in% c("NST", "ILC")) | (TissueSource == "Normal") | 
                                    (TissueSource == "Cell Lines" & RNA_Subtypes %notin% c("Basal")))

tumor_cl_histology_LumA <- ggplot(DMI_df_luminal, aes(x = reorder(Histology_TissueSource, DMI_zscore), y = DMI_zscore)) +
  ggbeeswarm::geom_quasirandom(size = 1, alpha = 1, width = 0.4, color = "black") + 
  geom_boxplot(alpha = 0.6, linewidth = 0.4, colour = "black", outlier.color = NA, width = 0.2) +
  xlab(" ") + ylab("DNAm Index") + theme_bw(20) + theme(axis.text.x = element_blank())

suppressWarnings({
  SupFig11C_pval_table <- pairwise.t.test(x = DMI_df_luminal$DMI, g = DMI_df_luminal$Histology_TissueSource)$p.value
})

# ------------------------------------------------------------------------------
# Fig 4A - DMI by Tissue Source (Normal, Patient Tumors, Cell Lines)
# ------------------------------------------------------------------------------
message("\n========================================")
message("Figure 4A: DMI by tissue source")
message("========================================")

# Filter to Luminal A tumors and non-basal cell lines for final DMI_df
normal_subset <- DMI_df$TissueSource == "Normal"
LumA_tumors_subset <- (DMI_df$TissueSource == "Patient Tumors") & (DMI_df$PAM50 == "LumA")
Non_basal_cl_subset <- (DMI_df$TissueSource == "Cell Lines") & (DMI_df$PAM50 != "Basal")
DMI_df <- subset(DMI_df, normal_subset | LumA_tumors_subset | Non_basal_cl_subset)

tissue_dmi <- ggplot(
  subset(DMI_df, TissueSource %in% c("Patient Tumors", "Cell Lines", "Normal")) %>% 
    mutate(TissueSource = forcats::fct_relevel(TissueSource, "Normal", "Patient Tumors", "Cell Lines")),
  aes(x = reorder(TissueSource, DMI_zscore), y = DMI_zscore)
) +
  ggbeeswarm::geom_quasirandom(aes(color = TissueSource), size = 1, alpha = 1, width = 0.4) +
  geom_boxplot(alpha = 0.6, linewidth = 0.4, colour = "black", outlier.color = NA, width = 0.2) +
  scale_color_manual("Tissue Source", values = annot_cols$Type, guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
  theme_bw(20) + xlab(" ") + ylab("DNAm Index") + theme(axis.text.x = element_blank())

# Suppress warnings from t.test
suppressWarnings({
  t.test(subset(DMI_df, TissueSource == "Normal")$DMI_zscore, 
         subset(DMI_df, TissueSource == "Patient Tumors")$DMI_zscore)$p.value
  t.test(subset(DMI_df, TissueSource == "Cell Lines")$DMI_zscore, 
         subset(DMI_df, TissueSource == "Patient Tumors")$DMI_zscore)$p.value
})

p1 <- t.test(subset(DMI_df, TissueSource == "Normal")$DMI_zscore, subset(DMI_df, TissueSource == "Patient Tumors")$DMI_zscore)$p.value
p2 <- t.test(subset(DMI_df, TissueSource == "Cell Lines")$DMI_zscore, subset(DMI_df, TissueSource == "Patient Tumors")$DMI_zscore)$p.value

library(gt)
tissue_dmi_pval <- tibble(`Normal vs Tumor` = p1, `Cell Line vs Tumor` = p2) %>%
  gt() %>%
  fmt_scientific(columns = everything(), decimals = 3) %>%
  tab_header(title = "DMI Z-Score Comparison (p-values)")

message("  ✓ SupFig 11 and Fig 4A DMI figures assigned to global environment")

# ------------------------------------------------------------------------------
# SECTION 4: Differential Analysis (TCGA and Cell Lines)
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("Section 4: Differential Analysis (TCGA and Cell Lines)")
message("=", strrep("=", 78))

# ------------------------------------------------------------------------------
# TCGA DESeq2 Analysis (LumA ILC vs NST)
# ------------------------------------------------------------------------------
message("TCGA DESeq2 Analysis (LumA ILC vs NST)...")

# Define sample groups
ILC_ids <- intersect(
  subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("ILC"))$Case.ID, 
  colnames(TCGA_BRCA_FC)
)
NST_ids <- intersect(
  subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("NST"))$Case.ID, 
  colnames(TCGA_BRCA_FC)
)

TCGA_de_mat <- TCGA_BRCA_FC[, c(ILC_ids, NST_ids)]

# Check if results exist, otherwise run analysis
if (!file.exists(FILES$tcga_deseq)) {
  message("  TCGA DESeq2 results not found. Running analysis...")
  TCGA_DESeq_Results <- run_tcga_deseq_analysis(
    TCGA_BRCA_FC = TCGA_BRCA_FC,
    TCGA_Annots = TCGA_Annots,
    output_file = FILES$tcga_deseq,
    csv_file = file.path(dirname(FILES$tcga_deseq), "TCGA_BRCA_DESeq_LumA_ILC_vs_NST_p0.05_lfc1.csv"),
    alpha = 0.05,
    lfc_threshold = 1
  )
} else {
  message("  Loading existing TCGA DESeq2 results...")
  load(FILES$tcga_deseq)
}

# Add group labels
TCGA_DESeq_Results$degs$group <- ifelse(TCGA_DESeq_Results$degs$log2FoldChange > 0, "ILC", "NST")
message("  ✓ TCGA DESeq2: ", nrow(TCGA_DESeq_Results$degs), " significant DEGs")

# ------------------------------------------------------------------------------
# TCGA Limma Analysis (LumA ILC vs NST)
# ------------------------------------------------------------------------------
message("TCGA Limma Analysis (LumA ILC vs NST)...")

# Define sample groups for DNAm
ILC_ids_dnam <- intersect(
  subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("ILC"))$Case.ID, 
  colnames(TCGA_BRCA_DNAm)
)
NST_ids_dnam <- intersect(
  subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("NST"))$Case.ID, 
  colnames(TCGA_BRCA_DNAm)
)

TCGA_limma_mat <- TCGA_BRCA_DNAm[allProbes, c(ILC_ids_dnam, NST_ids_dnam)]

# Check if results exist, otherwise run analysis
if (!file.exists(FILES$tcga_limma)) {
  message("  TCGA Limma results not found. Running analysis...")
  TCGA_DPM_Results <- run_tcga_limma_analysis(
    TCGA_BRCA_DNAm = TCGA_BRCA_DNAm,
    TCGA_Annots = TCGA_Annots,
    HM450K_ProbeSet = HM450K_ProbeSet,
    allProbes = allProbes,
    output_file = FILES$tcga_limma,
    csv_file = file.path(dirname(FILES$tcga_limma), "TCGA_BRCA_Limma_LumA_ILC_vs_NST_p0.05_betacut0.1.csv"),
    p_cutoff = 0.05,
    logFC_cut = 0.5,
    use_mvalues = TRUE
  )
} else {
  message("  Loading existing TCGA Limma results...")
  load(FILES$tcga_limma)
  
  # Add correlation results if not already present
  if (!"RNA_cor" %in% colnames(TCGA_DPM_Results$significant_probes)) {
    message("  Adding RNA-DNAm correlation results...")
    commonIDs <- intersect(colnames(TCGA_limma_mat), colnames(TCGA_de_mat))
    TCGA_RNA_DNAm_Correlation_Res <- correlate_probe_expression(
      probe_info = TCGA_DPM_Results$significant_probes,
      beta_matrix = beta_to_M(TCGA_BRCA_DNAm[TCGA_DPM_Results$significant_probes$Probe, commonIDs]),
      expression_matrix = TCGA_BRCA_Log2CPM[, commonIDs]
    )
    rownames(TCGA_RNA_DNAm_Correlation_Res) <- TCGA_RNA_DNAm_Correlation_Res$Probe
    TCGA_DPM_Results$significant_probes$RNA_cor <- TCGA_RNA_DNAm_Correlation_Res[TCGA_DPM_Results$significant_probes$Probe, "correlation"]
    TCGA_DPM_Results$significant_probes$RNA_cor_pvalue <- TCGA_RNA_DNAm_Correlation_Res[TCGA_DPM_Results$significant_probes$Probe, "pvalue"]
    
    # Rank by lowest p-value, then highest absolute correlation
    TCGA_RNA_DNAm_Correlation_Res <- TCGA_RNA_DNAm_Correlation_Res[order(
      TCGA_RNA_DNAm_Correlation_Res$gene,
      TCGA_RNA_DNAm_Correlation_Res$pvalue,
      -abs(TCGA_RNA_DNAm_Correlation_Res$correlation)
    ), ]
    
    # Select top entry per gene
    TCGA_RNA_DNAm_Correlation_Res_unique <- TCGA_RNA_DNAm_Correlation_Res[!duplicated(TCGA_RNA_DNAm_Correlation_Res$gene), ]
    TCGA_DPM_Results$significant_probes[, "top_cor_probe"] <- ifelse(
      TCGA_DPM_Results$significant_probes$Probe %in% TCGA_RNA_DNAm_Correlation_Res_unique$Probe, TRUE, FALSE
    )
  }
}

# Clean up results
TCGA_DPM_Results$all_results <- TCGA_DPM_Results$all_results[!is.na(TCGA_DPM_Results$all_results$adj.P.Val), ]
TCGA_DPM_Results$significant_probes$group <- ifelse(TCGA_DPM_Results$significant_probes$logFC > 0, "ILC", "NST")
message("  ✓ TCGA Limma: ", nrow(TCGA_DPM_Results$significant_probes), " significant probes")

# ------------------------------------------------------------------------------
# Cell Line DESeq2 Analysis (Non-Basal ILC vs NST)
# ------------------------------------------------------------------------------
message("Cell Line DESeq2 Analysis (Non-Basal ILC vs NST)...")

# Define sample groups
nonBasal_cells <- subset(CL_Annots[cellsID, ], TopCall %notin% c("Basal") & overlapWCCLE != "Y")$Name
ILC_cl_ids <- intersect(
  subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("ILC"))$Name, 
  colnames(BRCA_CL_CTS)
)
ILC_like_cl_ids <- intersect(
  subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("ILC-like"))$Name, 
  colnames(BRCA_CL_CTS)
)
NST_cl_ids <- intersect(
  subset(CL_Annots[nonBasal_cells, ], TopCall != "Basal" & Histology %in% c("NST") & Sample != c("HCC202", "HCC1500"))$Name, 
  colnames(BRCA_CL_CTS)
)

de_mat <- as.matrix(BRCA_CL_CTS[, c(c(ILC_cl_ids, ILC_like_cl_ids), NST_cl_ids)])

# Check if results exist, otherwise run analysis
if (!file.exists(FILES$cl_deseq)) {
  message("  Cell Line DESeq2 results not found. Running analysis...")
  BRCAL_CL_DESeq_Results <- run_cl_deseq_analysis(
    BRCA_CL_CTS = BRCA_CL_CTS,
    CL_Annots = CL_Annots,
    cellsID = cellsID,
    output_file = FILES$cl_deseq,
    csv_file = file.path(dirname(FILES$cl_deseq), "BRCA_CL_DEseq2_NonBasal_ILC_vs_NST_p0.3_lfc1.csv"),
    alpha = 0.3,
    lfc_threshold = 1
  )
} else {
  message("  Loading existing Cell Line DESeq2 results...")
  load(FILES$cl_deseq)
}

# Add group labels
BRCAL_CL_DESeq_Results$degs$group <- ifelse(BRCAL_CL_DESeq_Results$degs$log2FoldChange > 0, "ILC", "NST")
message("  ✓ Cell Line DESeq2: ", nrow(BRCAL_CL_DESeq_Results$degs), " significant DEGs")

# ------------------------------------------------------------------------------
# Cell Line Limma Analysis (Non-Basal ILC vs NST)
# ------------------------------------------------------------------------------
message("Cell Line Limma Analysis (Non-Basal ILC vs NST)...")

# Define sample groups for DNAm
limma_mat <- BRCA_CL_DNAm[allProbes, c(ILC_cl_ids, ILC_like_cl_ids, NST_cl_ids)]
limma_mat <- limma_mat[rowSums(is.na(limma_mat)) <= 5, ]

# Check if results exist, otherwise run analysis
if (!file.exists(FILES$cl_limma)) {
  message("  Cell Line Limma results not found. Running analysis...")
  BRCA_CL_DPM_Results <- run_cl_limma_analysis(
    BRCA_CL_DNAm = BRCA_CL_DNAm,
    CL_Annots = CL_Annots,
    HM450K_ProbeSet = HM450K_ProbeSet,
    allProbes = allProbes,
    cellsID = cellsID,
    output_file = FILES$cl_limma,
    csv_file = file.path(dirname(FILES$cl_limma), "BRCA_CL_Limma_NonBasal_ILC_vs_NST_p0.3_betacut0.1.csv"),
    p_cutoff = 0.3,
    logFC_cut = 0.5,
    use_mvalues = TRUE
  )
} else {
  message("  Loading existing Cell Line Limma results...")
  load(FILES$cl_limma)
  
  # Add correlation results if not already present
  if (!"RNA_cor" %in% colnames(BRCA_CL_DPM_Results$significant_probes)) {
    message("  Adding RNA-DNAm correlation results...")
    commonIDs <- intersect(colnames(limma_mat), colnames(de_mat))
    BRCA_CL_RNA_DNAm_Correlation_Res <- correlate_probe_expression(
      probe_info = BRCA_CL_DPM_Results$significant_probes,
      beta_matrix = beta_to_M(BRCA_CL_DNAm[BRCA_CL_DPM_Results$significant_probes$Probe, commonIDs]),
      expression_matrix = BRCA_CL_EXP[, commonIDs]
    )
    rownames(BRCA_CL_RNA_DNAm_Correlation_Res) <- BRCA_CL_RNA_DNAm_Correlation_Res$Probe
    BRCA_CL_DPM_Results$significant_probes$RNA_cor <- BRCA_CL_RNA_DNAm_Correlation_Res[BRCA_CL_DPM_Results$significant_probes$Probe, "correlation"]
    BRCA_CL_DPM_Results$significant_probes$RNA_cor_pvalue <- BRCA_CL_RNA_DNAm_Correlation_Res[BRCA_CL_DPM_Results$significant_probes$Probe, "pvalue"]
    
    # Rank by lowest p-value, then highest absolute correlation
    BRCA_CL_RNA_DNAm_Correlation_Res <- BRCA_CL_RNA_DNAm_Correlation_Res[order(
      BRCA_CL_RNA_DNAm_Correlation_Res$gene,
      BRCA_CL_RNA_DNAm_Correlation_Res$pvalue,
      -abs(BRCA_CL_RNA_DNAm_Correlation_Res$correlation)
    ), ]
    
    # Select top entry per gene
    BRCA_CL_RNA_DNAm_Correlation_Res_unique <- BRCA_CL_RNA_DNAm_Correlation_Res[!duplicated(BRCA_CL_RNA_DNAm_Correlation_Res$gene), ]
    BRCA_CL_DPM_Results$significant_probes[, "top_cor_probe"] <- ifelse(
      BRCA_CL_DPM_Results$significant_probes$Probe %in% BRCA_CL_RNA_DNAm_Correlation_Res_unique$Probe, TRUE, FALSE
    )
  }
}

# Clean up results
BRCA_CL_DPM_Results$all_results <- BRCA_CL_DPM_Results$all_results[!is.na(BRCA_CL_DPM_Results$all_results$adj.P.Val), ]
BRCA_CL_DPM_Results$significant_probes$group <- ifelse(BRCA_CL_DPM_Results$significant_probes$logFC > 0, "ILC", "NST")
message("  ✓ Cell Line Limma: ", nrow(BRCA_CL_DPM_Results$significant_probes), " significant probes")

# ------------------------------------------------------------------------------
# SECTION 5: Identify Consensus Events Between Cell Lines and Tumors
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("Section 5: Identifying Consensus Events Between Cell Lines and Tumors")
message("=", strrep("=", 78))

# ------------------------------------------------------------------------------
# Consensus for Differential Gene Expression (DGE)
# ------------------------------------------------------------------------------
message("Identifying consensus DGE events...")

commonGenes <- intersect(TCGA_DESeq_Results$degs$gene, BRCAL_CL_DESeq_Results$degs$gene)
DGE_matches <- intersect(BRCAL_CL_DESeq_Results$degs$gene, TCGA_DESeq_Results$degs$gene)

BRCAL_CL_DESeq_Results$degs[DGE_matches, "consensus_DGE"] <- 
  BRCAL_CL_DESeq_Results$degs[DGE_matches, "group"] == TCGA_DESeq_Results$degs[DGE_matches, "group"]

TCGA_DESeq_Results$all_results[, "consensus_DGE"] <- FALSE
TCGA_DESeq_Results$all_results[commonGenes, "consensus_DGE"] <- 
  BRCAL_CL_DESeq_Results$degs[commonGenes, "consensus_DGE"]

TCGA_DESeq_Results$degs[, "consensus_DGE"] <- FALSE
TCGA_DESeq_Results$degs[commonGenes, "consensus_DGE"] <- 
  BRCAL_CL_DESeq_Results$degs[commonGenes, "consensus_DGE"]

message("  Consensus DGE events: ", sum(BRCAL_CL_DESeq_Results$degs$consensus_DGE, na.rm = TRUE))

# ------------------------------------------------------------------------------
# Consensus for Differential Probe Methylation (DPM)
# ------------------------------------------------------------------------------
message("Identifying consensus DPM events...")

DPM_matches <- intersect(BRCA_CL_DPM_Results$significant_probes$Probe, TCGA_DPM_Results$significant_probes$Probe)

BRCA_CL_DPM_Results$significant_probes[DPM_matches, "consensus_DPM"] <- 
  BRCA_CL_DPM_Results$significant_probes[DPM_matches, "group"] == 
  TCGA_DPM_Results$significant_probes[DPM_matches, "group"]

TCGA_DPM_Results$all_results[, "consensus_DPM"] <- FALSE
TCGA_DPM_Results$all_results[DPM_matches, "consensus_DPM"] <- 
  BRCA_CL_DPM_Results$significant_probes[DPM_matches, "consensus_DPM"]

TCGA_DPM_Results$significant_probes[, "consensus_DPM"] <- FALSE
TCGA_DPM_Results$significant_probes[DPM_matches, "consensus_DPM"] <- 
  BRCA_CL_DPM_Results$significant_probes[DPM_matches, "consensus_DPM"]

message("  Consensus DPM events: ", sum(BRCA_CL_DPM_Results$significant_probes$consensus_DPM, na.rm = TRUE))

# Add cell line correlation to TCGA results
TCGA_DPM_Results$significant_probes$RNA_cor_CL <- 
  BRCA_CL_DPM_Results$significant_probes[TCGA_DPM_Results$significant_probes$Probe, "RNA_cor"]

# ------------------------------------------------------------------------------
# Merge RNA and DNAm Results
# ------------------------------------------------------------------------------
message("Merging RNA and DNAm results...")

# Background merge (all genes with DNAm data)
TCGA_merged_results_bg <- TCGA_DESeq_Results$all_results %>%
  inner_join(
    TCGA_DPM_Results$all_results %>% group_by(gene) %>%
      summarise(
        n_probes = n(),
        avg_delta_beta = mean(delta_beta, na.rm = TRUE),
        top_DNAm_logFoldChange = max(logFC, na.rm = TRUE),
        Probe_Status = ifelse(all(logFC > 0), "hypermethylated",
                           ifelse(all(logFC < 0), "hypomethylated", "mixed"))
      ), 
    by = "gene"
  )

# Probe summary for significant probes with RNA correlation
TCGA_probe_summary <- subset(TCGA_DPM_Results$significant_probes, RNA_cor < -0.2) %>%
  group_by(gene) %>%
  summarise(
    n_probes = n(),
    probes = paste(Probe, collapse = ","),
    probe_regions = paste(Probe_region, collapse = ","),
    top_probe = Probe[which.max(logFC)],
    top_probe_region = Probe_region[which.max(logFC)],
    top_probe_logFoldChange = logFC[which.max(logFC)],
    top_cor_probe = Probe[top_cor_probe],
    top_cor_probe_region = Probe_region[top_cor_probe],
    top_cor_probe_RNAcor = RNA_cor[which.max(RNA_cor)],
    top_cor_probe_RNAcor_CL = RNA_cor_CL[which.max(RNA_cor)],
    avg_delta_beta = mean(delta_beta, na.rm = TRUE),
    consensus_DPM = {
      vals <- consensus_DPM[!is.na(consensus_DPM) & consensus_DPM != FALSE]
      if (length(vals) > 0) vals[1] else NA
    },
    Probe_Status = ifelse(all(logFC > 0), "hypermethylated",
                       ifelse(all(logFC < 0), "hypomethylated", "mixed"))
  )

TCGA_probe_summary <- as.data.frame(TCGA_probe_summary)
rownames(TCGA_probe_summary) <- TCGA_probe_summary$gene

# Merge DEG results with probe summary
TCGA_deg_results <- TCGA_DESeq_Results$degs
TCGA_merged_results <- TCGA_deg_results %>%
  inner_join(TCGA_probe_summary, by = "gene")

# Filter to canonical DNAm-RNA relationships
TCGA_DNAm_Regulated_Genes <- TCGA_merged_results %>%
  filter(
    (avg_delta_beta > 0 & log2FoldChange < 0) |   # hypermethylated & downregulated
    (avg_delta_beta < 0 & log2FoldChange > 0)      # hypomethylated & upregulated
  ) %>% 
  mutate(Regulation = ifelse(
    (avg_delta_beta > 0 & log2FoldChange < 0), "HyperDNAm -> Down-regulated",
    ifelse((avg_delta_beta < 0 & log2FoldChange > 0), "HypoDNAm -> Up-regulated", "Non-canonical")
  ))

message("  ✓ Identified ", nrow(TCGA_DNAm_Regulated_Genes), " DNAm-regulated genes")

# ------------------------------------------------------------------------------
# Add Annotations and Expression/Methylation Values
# ------------------------------------------------------------------------------
message("Adding annotations and expression/methylation values...")

# Add OncoKB and Driver Level annotations
rownames(oncoVar) <- oncoVar$Gene_symbol
TCGA_DNAm_Regulated_Genes[, "OncoKB"] <- gsub("[|]Pan", "", oncoVar[TCGA_DNAm_Regulated_Genes$gene, "OncoKB"])
TCGA_DNAm_Regulated_Genes[, "Driver Level"] <- oncoVar[TCGA_DNAm_Regulated_Genes$gene, "Driver.Level"]

# Add expression fold changes
TCGA_DNAm_Regulated_Genes$EXP_LFC <- TCGA_DESeq_Results$degs[TCGA_DNAm_Regulated_Genes$gene, "log2FoldChange"]
TCGA_DNAm_Regulated_Genes$EXP_LFC_CL <- BRCAL_CL_DESeq_Results$all_results[TCGA_DNAm_Regulated_Genes$gene, "log2FoldChange"]

# Add DNAm fold changes
TCGA_DNAm_Regulated_Genes$DNAm_LFC <- TCGA_DPM_Results$significant_probes[TCGA_DNAm_Regulated_Genes$top_cor_probe, "logFC"]
TCGA_DNAm_Regulated_Genes$DNAm_LFC_CL <- BRCA_CL_DPM_Results$significant_probes[TCGA_DNAm_Regulated_Genes$top_cor_probe, "logFC"]
TCGA_DNAm_Regulated_Genes$DNAm_LFC_CL[is.na(TCGA_DNAm_Regulated_Genes$DNAm_LFC_CL)] <- 
  BRCA_CL_DPM_Results$significant_probes[TCGA_DNAm_Regulated_Genes$top_probe[is.na(TCGA_DNAm_Regulated_Genes$DNAm_LFC_CL)], "logFC"]

# Add consensus flags
TCGA_DESeq_Results$all_results[, "consensus_DGE"] <- FALSE
TCGA_DESeq_Results$all_results[commonGenes, "consensus_DGE"] <- 
  BRCAL_CL_DESeq_Results$degs[commonGenes, "consensus_DGE"]

TCGA_DESeq_Results$degs[, "consensus_DGE"] <- FALSE
TCGA_DESeq_Results$degs[commonGenes, "consensus_DGE"] <- 
  BRCAL_CL_DESeq_Results$degs[commonGenes, "consensus_DGE"]

common_probes <- intersect(TCGA_DPM_Results$significant_probes$Probe, BRCA_CL_DPM_Results$significant_probes$Probe)
TCGA_DPM_Results$all_results[, "consensus_DPM"] <- FALSE
TCGA_DPM_Results$all_results[common_probes, "consensus_DPM"] <- 
  BRCA_CL_DPM_Results$significant_probes[common_probes, "consensus_DPM"]

TCGA_DPM_Results$significant_probes[, "consensus_DPM"] <- FALSE
TCGA_DPM_Results$significant_probes[common_probes, "consensus_DPM"] <- 
  BRCA_CL_DPM_Results$significant_probes[common_probes, "consensus_DPM"]

# Add expression values for different groups
normal_ids <- setdiff(grep("-11", value = TRUE, colnames(TCGA_BRCA_Log2CPM)), outlier_normals)
ILC_ids_exp <- intersect(
  subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("ILC"))$Case.ID, 
  colnames(TCGA_BRCA_FC)
)
NST_ids_exp <- intersect(
  subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("NST"))$Case.ID, 
  colnames(TCGA_BRCA_FC)
)

TCGA_DNAm_Regulated_Genes$Exp_ILC_Tumors <- log2(rowMeans(2^TCGA_BRCA_Log2CPM[TCGA_DNAm_Regulated_Genes$gene, ILC_ids_exp], na.rm = TRUE))
TCGA_DNAm_Regulated_Genes$Exp_NST_Tumors <- log2(rowMeans(2^TCGA_BRCA_Log2CPM[TCGA_DNAm_Regulated_Genes$gene, NST_ids_exp], na.rm = TRUE))
TCGA_DNAm_Regulated_Genes$Exp_Normal_Adjacent <- log2(rowMeans(2^TCGA_BRCA_Log2CPM[TCGA_DNAm_Regulated_Genes$gene, normal_ids], na.rm = TRUE))
TCGA_DNAm_Regulated_Genes$Exp_ILC_CL <- log2(rowMeans(2^BRCA_CL_EXP[match(TCGA_DNAm_Regulated_Genes$gene, rownames(BRCA_CL_EXP)), ILC_cl_ids], na.rm = TRUE))
TCGA_DNAm_Regulated_Genes$Exp_NST_CL <- log2(rowMeans(2^BRCA_CL_EXP[match(TCGA_DNAm_Regulated_Genes$gene, rownames(BRCA_CL_EXP)), NST_cl_ids], na.rm = TRUE))
TCGA_DNAm_Regulated_Genes$Exp_ILC_like_CL <- log2(rowMeans(2^BRCA_CL_EXP[match(TCGA_DNAm_Regulated_Genes$gene, rownames(BRCA_CL_EXP)), ILC_like_cl_ids], na.rm = TRUE))

# Add DNAm values for different groups
normal_ids_dnam <- setdiff(grep("-11", value = TRUE, colnames(TCGA_BRCA_DNAm)), outlier_normals)
ILC_ids_dnam_final <- intersect(
  subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("ILC"))$Case.ID, 
  colnames(TCGA_BRCA_DNAm)
)
NST_ids_dnam_final <- intersect(
  subset(TCGA_Annots, !is.na(TumorPurity_CPE) & PAM50 == "LumA" & `Final Pathology` %in% c("NST"))$Case.ID, 
  colnames(TCGA_BRCA_DNAm)
)

TCGA_DNAm_Regulated_Genes$DNAm_ILC_Tumors <- rowMeans(TCGA_BRCA_DNAm[TCGA_DNAm_Regulated_Genes$top_cor_probe, ILC_ids_dnam_final], na.rm = TRUE)
TCGA_DNAm_Regulated_Genes$DNAm_NST_Tumors <- rowMeans(TCGA_BRCA_DNAm[TCGA_DNAm_Regulated_Genes$top_cor_probe, NST_ids_dnam_final], na.rm = TRUE)
TCGA_DNAm_Regulated_Genes$DNAm_Normal_Adjacent <- rowMeans(TCGA_BRCA_DNAm[TCGA_DNAm_Regulated_Genes$top_cor_probe, normal_ids_dnam], na.rm = TRUE)
TCGA_DNAm_Regulated_Genes$DNAm_ILC_CL <- rowMeans(BRCA_CL_DNAm[TCGA_DNAm_Regulated_Genes$top_cor_probe, ILC_cl_ids], na.rm = TRUE)
TCGA_DNAm_Regulated_Genes$DNAm_NST_CL <- rowMeans(BRCA_CL_DNAm[TCGA_DNAm_Regulated_Genes$top_cor_probe, NST_cl_ids], na.rm = TRUE)
TCGA_DNAm_Regulated_Genes$DNAm_ILC_like_CL <- rowMeans(BRCA_CL_DNAm[TCGA_DNAm_Regulated_Genes$top_cor_probe, ILC_like_cl_ids], na.rm = TRUE)

# ------------------------------------------------------------------------------
# Identify Consensus Genes
# ------------------------------------------------------------------------------
message("Identifying consensus genes...")

rownames(TCGA_DNAm_Regulated_Genes) <- TCGA_DNAm_Regulated_Genes$gene

# Add consensus flags from merged results
TCGA_DNAm_Regulated_Genes$consensus_DGE <- TCGA_DESeq_Results$degs[TCGA_DNAm_Regulated_Genes$gene, "consensus_DGE"]
TCGA_DNAm_Regulated_Genes$consensus_DPM <- TCGA_DPM_Results$significant_probes[TCGA_DNAm_Regulated_Genes$top_cor_probe, "consensus_DPM"]

# Define consensus based on multiple criteria
TCGA_DNAm_Regulated_Genes$Consensus <- ifelse(
  (TCGA_DNAm_Regulated_Genes$consensus_DGE == TRUE | TCGA_DNAm_Regulated_Genes$consensus_DPM == TRUE) & 
  !is.na(TCGA_DNAm_Regulated_Genes$top_cor_probe_RNAcor_CL) & 
  !is.na(TCGA_DNAm_Regulated_Genes$EXP_LFC_CL) &
  !is.na(TCGA_DNAm_Regulated_Genes$EXP_LFC) &
  (sign(TCGA_DNAm_Regulated_Genes$top_cor_probe_RNAcor) == sign(TCGA_DNAm_Regulated_Genes$top_cor_probe_RNAcor_CL)) & 
  (sign(TCGA_DNAm_Regulated_Genes$EXP_LFC) == sign(TCGA_DNAm_Regulated_Genes$EXP_LFC_CL)) &
  (sign(TCGA_DNAm_Regulated_Genes$DNAm_LFC) == sign(TCGA_DNAm_Regulated_Genes$DNAm_LFC_CL)) &
  TCGA_DNAm_Regulated_Genes$top_cor_probe_RNAcor < -0.1 & 
  TCGA_DNAm_Regulated_Genes$top_cor_probe_RNAcor_CL < -0.1, 
  "Y", "N"
)

TCGA_DNAm_Regulated_Genes$Consensus[is.na(TCGA_DNAm_Regulated_Genes$Consensus)] <- "N"

TCGA_DNAm_Regulated_Genes$Category <- ifelse(
  TCGA_DNAm_Regulated_Genes$Consensus == "Y",
  paste0(TCGA_DNAm_Regulated_Genes$Regulation, " (consensus)"),
  TCGA_DNAm_Regulated_Genes$Regulation
)

# Sort by absolute log2FoldChange
TCGA_DNAm_Regulated_Genes <- TCGA_DNAm_Regulated_Genes[order(abs(TCGA_DNAm_Regulated_Genes$log2FoldChange), decreasing = TRUE), ]

consensus_genes <- subset(TCGA_DNAm_Regulated_Genes, Consensus == "Y")
message("  ✓ Identified ", nrow(consensus_genes), " consensus genes")

# ------------------------------------------------------------------------------
# SECTION 6: Generate Fig 4B - DNAm-mRNA Alterations Scatter Plot
# ------------------------------------------------------------------------------

message("\n========================================")
message("Figure 4B: DNAm-mRNA alterations scatter plot")
message("========================================")

# Select top genes for labeling
top_genes <- rbind(
  head(subset(TCGA_DNAm_Regulated_Genes, Category == "HypoDNAm -> Up-regulated" & Consensus == "N" & 
              abs(top_probe_logFoldChange) >= 0.5 & padj < 0.05), 8),
  head(subset(TCGA_DNAm_Regulated_Genes, Category == "HyperDNAm -> Down-regulated" & Consensus == "N" & 
              abs(top_probe_logFoldChange) >= 0.5 & padj < 0.05), 8)
)

# Create scatter plot
dnma_lfc_plt <- ggplot(TCGA_merged_results_bg, aes(x = avg_delta_beta, y = log2FoldChange)) +
  geom_point(color = "gray90", size = 1, shape = 20, show.legend = FALSE) + 
  geom_point(aes(color = Category, size = -log10(padj)), shape = 20, stroke = 0.1,
             data = subset(TCGA_DNAm_Regulated_Genes, gene %notin% c(consensus_genes$gene)), show.legend = FALSE) +
  geom_point(aes(color = Category, size = -log10(padj)), shape = 17, alpha = 0.6, 
             data = consensus_genes, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = rbind(top_genes, consensus_genes),
                           aes(label = gene), segment.color = "gray50",
                           box.padding = 0, point.padding = 0, min.segment.length = 0,
                           max.overlaps = Inf, force = 60,
                           size = 5, show.legend = FALSE) +
  scale_color_manual(values = c(
    "HyperDNAm -> Down-regulated" = "#FFCCBC", 
    "HyperDNAm -> Down-regulated (consensus)" = "#99000D",
    "HypoDNAm -> Up-regulated (consensus)" = "#08306B", 
    "HypoDNAm -> Up-regulated" = "#B3E5FC"
  )) + 
  theme_bw(20) +
  xlab("Average Delta Beta (DNAm)") + ylab("log2 Fold Change (RNA)") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.2, color = "#263238", alpha = 1) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2, color = "#263238", alpha = 1) + 
  scale_size_continuous(range = c(2, 15))

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("dnma_lfc_plt", dnma_lfc_plt, envir = .GlobalEnv)
message("  ✓ Fig 4B complete (dnma_lfc_plt assigned).")

# ------------------------------------------------------------------------------
# SECTION 7: Generate Fig 4C - Consensus DNAm-mRNA Alterations Heatmap
# ------------------------------------------------------------------------------

message("\n========================================")
message("Figure 4C: Consensus DNAm-mRNA alterations heatmap")
message("========================================")

# Extract consensus genes
TCGA_DNAm_Regulated_Genes_Consensus <- subset(TCGA_DNAm_Regulated_Genes, Consensus == "Y")
rownames(TCGA_DNAm_Regulated_Genes_Consensus) <- TCGA_DNAm_Regulated_Genes_Consensus$gene

# Prepare z-score matrices
colNames <- c("ILC (Patient Tumors)", "NST (Patient Tumors)", "ILC (Cell Lines)", "ILC-like (Cell Lines)", "NST (Cell Lines)")

RNA_zcore_mat <- cbind(
  t(scale(t(TCGA_DNAm_Regulated_Genes_Consensus[, c("Exp_ILC_Tumors", "Exp_NST_Tumors")]))),
  t(scale(t(TCGA_DNAm_Regulated_Genes_Consensus[, c("Exp_ILC_CL", "Exp_ILC_like_CL", "Exp_NST_CL")])))
)
colnames(RNA_zcore_mat) <- colNames

DNAm_zscore_mat <- cbind(
  t(scale(t(TCGA_DNAm_Regulated_Genes_Consensus[, c("DNAm_ILC_Tumors", "DNAm_NST_Tumors")]))),
  t(scale(t(TCGA_DNAm_Regulated_Genes_Consensus[, c("DNAm_ILC_CL", "DNAm_ILC_like_CL", "DNAm_NST_CL")])))
)
colnames(DNAm_zscore_mat) <- colNames

heatmap_legend_param_h <- modifyList(heatmap_legend_param, list(direction = "horizontal"))

# Create DNAm heatmap
DNAm_zscore_ht <- Heatmap(
  t((t(DNAm_zscore_mat))), 
  show_row_names = TRUE, 
  heatmap_legend_param = heatmap_legend_param_h, 
  height = unit(4.5, "cm"), 
  width = unit(6, "cm"),
  row_labels = TCGA_DNAm_Regulated_Genes_Consensus[, "top_cor_probe"], 
  row_names_side = "left",
  left_annotation = HeatmapAnnotation(
    which = "row", 
    col = list(Region = annot_cols$DNAm_regions),
    annotation_legend_param = heatmap_legend_param_h, 
    simple_anno_size = unit(4, "mm"),
    Region = HM450K_ProbeSet[TCGA_DNAm_Regulated_Genes_Consensus$top_cor_probe, "region_simple"]
  ),
  name = "DNAm z-scores", 
  border = TRUE, 
  column_split = factor(colNames, colNames), 
  row_title = " ", 
  column_title = " ",
  cluster_column_slices = FALSE, 
  cluster_row_slices = FALSE, 
  column_title_gp = gpar(fontface = "bold", fontsize = 14),
  column_names_gp = gpar(fontsize = 15, fontface = "bold"), 
  row_names_gp = gpar(fontsize = 15),
  cluster_columns = FALSE, 
  cluster_rows = FALSE, 
  show_row_dend = FALSE, 
  show_column_dend = FALSE,
  row_split = TCGA_DNAm_Regulated_Genes_Consensus$Regulation, 
  col = annot_cols$DNAm_zscore
)

# Create RNA heatmap
RNA_zscore_ht <- Heatmap(
  t((t(RNA_zcore_mat))), 
  heatmap_legend_param = heatmap_legend_param_h,
  row_names_side = "right", 
  height = unit(4.5, "cm"), 
  width = unit(6.2, "cm"),
  name = "RNA z-scores", 
  border = TRUE, 
  column_split = factor(colNames, colNames), 
  row_title = " ", 
  column_title = " ",
  cluster_column_slices = FALSE, 
  cluster_row_slices = FALSE, 
  column_title_gp = gpar(fontface = "bold", fontsize = 14),
  column_names_gp = gpar(fontsize = 15, fontface = "bold"), 
  row_names_gp = gpar(fontsize = 15),
  cluster_columns = FALSE, 
  cluster_rows = FALSE, 
  show_row_dend = FALSE, 
  show_column_dend = FALSE,
  row_split = TCGA_DNAm_Regulated_Genes_Consensus$Regulation, 
  col = annot_cols$RNA_zscore
)

# Combine heatmaps
fig4c_ht <- DNAm_zscore_ht + RNA_zscore_ht

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("fig4c_ht", fig4c_ht, envir = .GlobalEnv)
message("  ✓ Fig 4C complete (fig4c_ht assigned).")

# ------------------------------------------------------------------------------
# SECTION 8: Generate Fig 4D - MSI1 & TFAP2B Expression and DNAm Barplots
# ------------------------------------------------------------------------------

message("\n========================================")
message("Figure 4D: MSI1 and TFAP2B expression and DNAm barplots")
message("========================================")

# Prepare data for cell lines
nonBasal_cells_final <- subset(
  CL_Annots[cellsID, ], 
  TopCall %notin% c("Basal") & overlapWCCLE != "Y" & Sample %notin% c("HCC202", "HCC1500")
)$Name

DNAm_zscore_mat_cl <- BRCA_CL_DNAm[TCGA_DNAm_Regulated_Genes_Consensus$top_cor_probe, nonBasal_cells_final]
RNA_zcore_mat_cl <- BRCA_CL_EXP[TCGA_DNAm_Regulated_Genes_Consensus$gene, nonBasal_cells_final]

# Melt both matrices into long format
expr_long <- reshape2::melt(t(scale(t(RNA_zcore_mat_cl))), varnames = c("Gene", "CellLine"), value.name = "Expression")
meth_long <- reshape2::melt(t(scale(t(DNAm_zscore_mat_cl))), varnames = c("Gene", "CellLine"), value.name = "Methylation")
meth_long$Gene <- HM450K_ProbeSet[as.character(meth_long$Gene), "transcript_geneName_noENS_NA"]

# Merge expression and methylation data
merged_df <- merge(expr_long, meth_long, by = c("Gene", "CellLine"))
merged_df <- merged_df[!is.na(merged_df$Methylation), ]
merged_df$Histology <- CL_Annots[as.character(merged_df$CellLine), "Histology"]

# Melt into long format for combined plot
plot_df <- merged_df %>%
  select(Gene, CellLine, Histology, Expression, Methylation) %>%
  tidyr::pivot_longer(cols = c("Expression", "Methylation"), names_to = "Measure", values_to = "Value")

plot_df$Measure <- factor(plot_df$Measure, levels = c("Methylation", "Expression"))

# Calculate order from methylation
ordering_info <- plot_df %>%
  filter(Measure == "Methylation") %>%
  group_by(Gene) %>%
  arrange(desc(Value)) %>%
  mutate(order = row_number(), Gene_CellLine = reorder_within(CellLine, -order, Gene)) %>%
  ungroup() %>%
  select(Gene, CellLine, order, Gene_CellLine)

# Filter to MSI1 and TFAP2B
gene1 <- "MSI1"
gene2 <- "TFAP2B"

plot_df_sub <- plot_df %>%
  left_join(ordering_info, by = c("Gene", "CellLine")) %>%
  filter(Gene %in% c(gene1, gene2)) %>%
  mutate(
    Gene_CellLine = reorder_within(CellLine, -order, Gene),
    facet_group = factor(
      paste0(Gene, " ", Measure), 
      levels = c(
        paste0(gene1, " Methylation"), paste0(gene2, " Methylation"), 
        paste0(gene1, " Expression"), paste0(gene2, " Expression")
      )
    )
  )

# Create barplot
fig4e <- ggplot(plot_df_sub, aes(x = Gene_CellLine, y = Value, fill = Histology)) +
  geom_col(color = NA, width = 0.8) +
  facet_wrap(~ facet_group, scales = "free", ncol = length(unique(plot_df_sub$Gene))) +
  scale_x_reordered() +
  scale_fill_manual(values = annot_cols$Histology) +
  theme_minimal(base_size = 12) + 
  ylim(c(-2, 2)) +
  labs(y = "Z-score", x = "") + 
  theme_bw(20, base_family = "Helvetica") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top", 
    strip.background = element_blank()
  )

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("fig4e", fig4e, envir = .GlobalEnv)
message("  ✓ Fig 4D complete (fig4e assigned).")

# ------------------------------------------------------------------------------
# SECTION 9: Save Results
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("Section 9: Saving Results")
message("=", strrep("=", 78))

# Save TCGA_DNAm_Regulated_Genes to CSV
write.csv(TCGA_DNAm_Regulated_Genes, 
          file = file.path(dnam_dir, "TCGA_DNAm_Regulated_Genes.csv"), 
          row.names = FALSE)
message("  ✓ Saved TCGA_DNAm_Regulated_Genes.csv to: ", file.path(dnam_dir, "TCGA_DNAm_Regulated_Genes.csv"))

message("\n", "=", strrep("=", 78))
message("Script 22 completed successfully!")
message("=", strrep("=", 78))
message("  All figures assigned to global environment for saving in Main_Data_Analysis.Rmd")
message("  - SupFig11: tumor_pam50, cl_pam50, tumor_cl_histology_LumA")
message("  - Fig4A: tissue_dmi")
message("  - Fig4B: dnma_lfc_plt")
message("  - Fig4C: fig4c_ht")
message("  - Fig4D: fig4e")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/21_Fig4_SupFig11_DNAm_Alterations.R")
#
# load_all_icle_data(load_external = TRUE)
# # Script runs on source. Uses FILES$hm450k_probeset, rnaseq_cts, icle_dnam, etc.
#
# Outputs: 
#   - DMI_Summary.tsv in DIRS$results_sub$dna_methylation
#   - TCGA_DNAm_Regulated_Genes.csv in DIRS$results_sub$dna_methylation
#   - All figures assigned to global environment for saving in Main_Data_Analysis.Rmd
