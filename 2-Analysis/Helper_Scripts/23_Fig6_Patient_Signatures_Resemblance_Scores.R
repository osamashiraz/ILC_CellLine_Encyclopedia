# ==============================================================================
# Script 24: Figure 6 - Patient Signatures & Resemblance Scores
# ==============================================================================
# Description: Generates ILC vs NST multi-omic signatures from TCGA patient
#              tumors and calculates cell-line resemblance scores. Compares ICLE
#              cell lines to patient tumor signatures across multiple data
#              modalities (CNV, Alterations, RNA, DNAm, RPPA).
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - TCGA_Annots: TCGA tumor annotations
#   - TCGA_BRCA_GAM: TCGA genomic alteration matrix
#   - TCGA_BRCA_Log2CPM: TCGA RNA expression data
#   - TCGA_BRCA_DNAm: TCGA DNA methylation data
#   - TCGA_BRCA_RPPA: TCGA RPPA protein data
#   - TCGA_BRCA_CN: TCGA copy number data (genes x samples)
#   - BRCA_CL_GAM: Cell line genomic alteration matrix
#   - BRCA_CL_EXP: Cell line RNA expression data
#   - BRCA_CL_DNAm: Cell line DNA methylation data
#   - BRCA_CL_RPPA: Cell line RPPA protein data
#   - BRCA_CL_CN: Cell line copy number data (genes x samples)
#   - CL_Annots: Cell line annotations
#
# Output:
#   - Part 1: Patient-derived signatures saved to files:
#     * ILC_CNV_Genes.tsv
#     * ILC_alterations.tsv
#     * ILC_DEGs.tsv
#     * ILC_DNAm_Probes.tsv
#     * ILC_DEPs.tsv
#   - Part 2: Resemblance scores and consensus:
#     * Resemblance_Scores.tsv
#   - Assigned to .GlobalEnv: cn_sig_ht, mut_sig_ht, dnam_sign_ht, rna_sign_ht, rppa_sig_ht, fig6_resemblance_ht
#
# Author: Osama Shiraz Shah
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: Load Configuration and Dependencies
# ------------------------------------------------------------------------------

DIRS$results_sub$molecular_resemblance <- DIRS$results_sub$molecular_resemblance
signatures_dir <- file.path(DIRS$results_sub$molecular_resemblance, "ILC_Signatures")
dir.create(signatures_dir, showWarnings = FALSE, recursive = TRUE)

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(limma)
  library(vroom)
})


# Check for CN data (may need to be loaded separately)
if (!exists("TCGA_BRCA_CN", envir = .GlobalEnv) || !exists("BRCA_CL_CN", envir = .GlobalEnv)) {
  message("  Loading CN data...")
  # Try to load from CNTools files if they exist
  tcga_cn_file <- file.path(DIRS$external$tcga, "CN", "TCGA_BRCA_LogRR_CNTools.tsv")
  cl_cn_file <- file.path(DIRS$icle$cytosnp, "7_CNTools_logRR", "BRCA_CL_LogRR_CNTools.tsv")
  
  if (file.exists(tcga_cn_file) && !exists("TCGA_BRCA_CN", envir = .GlobalEnv)) {
    TCGA_BRCA_CN <- vroom::vroom(tcga_cn_file, delim = "\t", show_col_types = FALSE) %>% as.data.frame()
    rownames(TCGA_BRCA_CN) <- TCGA_BRCA_CN$genename
    TCGA_BRCA_CN <- TCGA_BRCA_CN[, -c(1:5)]
    assign("TCGA_BRCA_CN", TCGA_BRCA_CN, envir = .GlobalEnv)
  }
  
  if (file.exists(cl_cn_file) && !exists("BRCA_CL_CN", envir = .GlobalEnv)) {
    BRCA_CL_CN <- vroom::vroom(cl_cn_file, delim = "\t", show_col_types = FALSE) %>% as.data.frame()
    rownames(BRCA_CL_CN) <- BRCA_CL_CN$genename
    BRCA_CL_CN <- BRCA_CL_CN[, -c(1:5)]
    colnames(BRCA_CL_CN) <- gsub("-M", "-C", colnames(BRCA_CL_CN))
    assign("BRCA_CL_CN", BRCA_CL_CN, envir = .GlobalEnv)
  }
}

# ------------------------------------------------------------------------------
# SECTION 1: Helper Functions
# ------------------------------------------------------------------------------

# ==============================================================================
# Utility Functions
# ==============================================================================

#' Rank normalization to [0,1]
rank01 <- function(x) (rank(x, na.last = "keep", ties.method = "average") - 1) / (sum(!is.na(x)) - 1)

#' Scale values to [0,1]
scale01 <- function(v) (v + 1) / 2

# ==============================================================================
# FUNCTION: Find Differential Proteins (RPPA)
# ==============================================================================
#' Find differentially expressed proteins using limma
#' 
#' @param rppa_mat Protein expression matrix (proteins x samples)
#' @param groups Named vector of group labels (ILC/NST)
#' @param contrast Contrast string (default: "ILC - NST")
#' @return Data frame with differential protein results
find_diff_proteins <- function(rppa_mat, groups, contrast = "ILC - NST") {
  stopifnot(all(colnames(rppa_mat) %in% names(groups)))
  
  # Subset to matched samples
  common_samples <- intersect(colnames(rppa_mat), names(groups))
  expr <- rppa_mat[, common_samples, drop = FALSE]
  g <- factor(groups[common_samples])
  
  # Design and fit
  design <- model.matrix(~ 0 + g)
  colnames(design) <- levels(g)
  fit <- lmFit(expr, design)
  
  # Contrast
  cont <- makeContrasts(contrasts = contrast, levels = design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)
  
  # Results
  res <- topTable(fit2, number = Inf, sort.by = "P")
  res$Protein <- rownames(res)
  rownames(res) <- NULL
  return(res)
}

# ==============================================================================
# FUNCTION: Compute Resemblance (RNA, DNAm, RPPA)
# ==============================================================================
#' Compute resemblance scores between tumor centroids and cell lines
#' 
#' @param tumor_mat Tumor data matrix (features x samples)
#' @param tumor_samples Vector of tumor sample IDs
#' @param cell_mat Cell line data matrix (features x samples)
#' @param cell_samples Vector of cell line IDs
#' @param features Vector of feature names to use
#' @param n_iter Number of bootstrap iterations (default: 1000)
#' @param n_feature Number of features to sample per iteration (default: NULL = all)
#' @param cor_method Correlation method (default: "spearman")
#' @return Data frame with BootMedian scores per cell line
compute_resemblance <- function(tumor_mat, tumor_samples, cell_mat, cell_samples, 
                                features, n_iter = 1000, n_feature = NULL, 
                                cor_method = "spearman") {
  cell_samples <- intersect(colnames(cell_mat), cell_samples)
  tumor_samples <- intersect(colnames(tumor_mat), tumor_samples)
  
  features <- intersect(features, intersect(rownames(tumor_mat), rownames(cell_mat)))
  centroid <- rowMeans(tumor_mat[features, tumor_samples])
  
  # Bootstrap across features
  res <- lapply(cell_samples, function(cl) {
    cors <- replicate(n_iter, {
      if (!is.null(n_feature)) {
        feats <- sample(features, n_feature)
      } else {
        feats <- features
      }
      stats::cor(centroid[feats], cell_mat[feats, cl], use = "complete.obs", method = cor_method)
    })
    data.frame(
      CellLine = cl,
      BootMedian = median(cors, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, res)
  rownames(out) <- out$CellLine
  return(out)
}

# ==============================================================================
# FUNCTION: Compute Resemblance CNV
# ==============================================================================
#' Compute resemblance scores for copy number data (per-chromosome bootstrap)
#' 
#' @param tumor_cn Tumor CN matrix (genes x samples)
#' @param cell_cn Cell line CN matrix (genes x samples)
#' @param ga Gene annotation data frame (genename, chrom)
#' @param tumor_samples Vector of tumor sample IDs
#' @param cell_samples Vector of cell line IDs
#' @param B Number of bootstrap iterations (default: 1000)
#' @param method Correlation method (default: "spearman")
#' @param seed Random seed (default: NULL)
#' @param return_distributions Return full distributions (default: FALSE)
#' @return Data frame with BootMean and BootMedian scores per cell line
compute_resemblance_cnv <- function(tumor_cn, cell_cn, ga, tumor_samples, 
                                    cell_samples, B = 1000, method = "spearman",
                                    seed = NULL, return_distributions = FALSE) {
  method <- match.arg(method, c("pearson", "spearman"))
  
  tumor_cn <- as.matrix(tumor_cn)
  cell_cn <- as.matrix(cell_cn)
  
  cell_samples <- intersect(colnames(cell_cn), cell_samples)
  tumor_samples <- intersect(colnames(tumor_cn), tumor_samples)
  
  common_genes <- intersect(rownames(tumor_cn), intersect(rownames(cell_cn), ga$genename))
  ga <- ga[ga$genename %in% common_genes, , drop = FALSE]
  
  genes_by_chr <- split(ga$genename, ga$chrom)
  chromosomes <- as.character(ga$chrom)
  
  pools <- lapply(chromosomes, function(ch) {
    gpool <- genes_by_chr[[ch]]
    if (is.null(gpool)) character(0) else intersect(gpool, common_genes)
  })
  names(pools) <- chromosomes
  
  n_chr <- length(chromosomes)
  n_cells <- length(cell_samples)
  
  # Storage
  cor_mat <- matrix(NA_real_, nrow = B, ncol = n_cells,
                    dimnames = list(paste0("iter", seq_len(B)), cell_samples))
  selected_genes <- matrix(NA_character_, nrow = n_chr, ncol = B,
                           dimnames = list(chromosomes, paste0("iter", seq_len(B))))
  
  if (!is.null(seed)) set.seed(seed)
  
  # Main bootstrap loop
  for (b in seq_len(B)) {
    # For each chromosome choose one gene
    chosen <- vapply(pools, function(pool) {
      if (length(pool) == 0L) NA_character_ else sample(pool, 1)
    }, FUN.VALUE = character(1))
    selected_genes[, b] <- chosen
    
    # Build tumor centroid vector (per-chromosome mean of chosen gene)
    tumor_centroid <- vapply(chosen, function(g) {
      if (is.na(g)) return(NA_real_)
      mean(as.numeric(tumor_cn[g, tumor_samples]), na.rm = TRUE)
    }, FUN.VALUE = numeric(1))
    
    # For each cell line extract the chosen-gene values and correlate
    for (ci in seq_along(cell_samples)) {
      cl <- cell_samples[ci]
      cell_vals <- vapply(chosen, function(g) {
        if (is.na(g)) return(NA_real_)
        as.numeric(cell_cn[g, cl])
      }, FUN.VALUE = numeric(1))
      
      ok <- complete.cases(tumor_centroid, cell_vals)
      if (sum(ok) >= 2) {
        cor_mat[b, ci] <- cor(tumor_centroid[ok], cell_vals[ok], use = "complete.obs", method = method)
      } else {
        cor_mat[b, ci] <- NA_real_
      }
    }
  }
  
  # Summaries per cell line
  summary_df <- data.frame(
    CellLine = cell_samples,
    BootMean = apply(cor_mat, 2, function(x) mean(x, na.rm = TRUE)),
    BootMedian = apply(cor_mat, 2, function(x) median(x, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  if (return_distributions) {
    list(summary = summary_df, correlations = cor_mat, selected_genes = selected_genes)
  } else {
    return(summary_df)
  }
}

# ==============================================================================
# FUNCTION: Compute Alteration Resemblance
# ==============================================================================
#' Compute alteration resemblance scores based on ILC/NST signatures
#' 
#' @param mut_mat Alteration matrix (genes x cell lines)
#' @param ILC_specs ILC-associated alteration specs (e.g., "CDH1_MUT_MUT;LOH")
#' @param NST_specs NST-associated alteration specs
#' @param freq_tbl Frequency table with p_ILC and p_NST
#' @param scale01 Scale combined score to [0,1] (default: TRUE)
#' @param min_genes_per_side Minimum genes per side (default: 1)
#' @return Data frame with resemblance scores per cell line
compute_alteration_resemblance <- function(mut_mat, ILC_specs = NULL, NST_specs = NULL,
                                           freq_tbl = NULL, scale01 = TRUE,
                                           min_genes_per_side = 1) {
  mut_mat <- as.matrix(mut_mat)
  cell_lines <- colnames(mut_mat)
  
  if (!is.null(freq_tbl)) {
    stopifnot(all(c("gene", "class", "p_ILC", "p_NST") %in% names(freq_tbl)))
  }
  
  # Parse specs: "GENE_CLASS1_CLASS2" -> gene, classes
  parse_specs <- function(specs) {
    if (is.null(specs) || length(specs) == 0) {
      return(data.frame(gene = character(0), alts = I(list()), stringsAsFactors = FALSE))
    }
    out <- lapply(specs, function(s) {
      parts <- strsplit(as.character(s), "_", fixed = TRUE)[[1]]
      gene <- parts[1]
      alts <- if (length(parts) > 1) parts[-1] else NA_character_
      list(gene = gene, alts = unname(alts))
    })
    data.frame(gene = vapply(out, `[[`, character(1), "gene"),
               alts = I(lapply(out, `[[`, "alts")),
               stringsAsFactors = FALSE)
  }
  
  ILC_df <- parse_specs(ILC_specs)
  NST_df <- parse_specs(NST_specs)
  
  # Build lookup maps from freq_tbl
  key <- function(g, cl) paste0(g, "||", cl)
  freq_p_ILC <- freq_p_NST <- list()
  if (!is.null(freq_tbl) && nrow(freq_tbl) > 0) {
    kk <- key(freq_tbl$gene, freq_tbl$class)
    freq_p_ILC <- setNames(freq_tbl$p_ILC, kk)
    freq_p_NST <- setNames(freq_tbl$p_NST, kk)
  }
  
  # Helper: contribution for one spec row in one cell line
  contribution_one <- function(gene, allowed_alts, cl, side) {
    if (!gene %in% rownames(mut_mat)) return(0)
    val <- as.character(mut_mat[gene, cl])
    if (is.na(val) || identical(val, "")) return(0)
    
    if (all(is.na(allowed_alts))) {
      kk <- key(gene, val)
      p_base <- if (side == "ILC") {
        if (kk %in% names(freq_p_NST)) freq_p_NST[[kk]] else 0
      } else {
        if (kk %in% names(freq_p_ILC)) freq_p_ILC[[kk]] else 0
      }
      return(1 - p_base)
    } else {
      if (val %in% allowed_alts) {
        kk <- key(gene, val)
        p_base <- if (side == "ILC") {
          if (kk %in% names(freq_p_NST)) freq_p_NST[[kk]] else 0
        } else {
          if (kk %in% names(freq_p_ILC)) freq_p_ILC[[kk]] else 0
        }
        return(1 - p_base)
      } else {
        return(0)
      }
    }
  }
  
  # Compute per-cell-line components
  out_list <- lapply(cell_lines, function(cl) {
    # ILC side
    if (nrow(ILC_df) > 0) {
      ilc_vals <- mapply(function(g, a) contribution_one(g, a, cl, "ILC"),
                         ILC_df$gene, ILC_df$alts, SIMPLIFY = TRUE)
      part_ilc <- mean(ilc_vals)
      ilc_hits <- sum(ilc_vals > 0)
    } else {
      part_ilc <- NA_real_
      ilc_hits <- NA_integer_
    }
    
    # NST side
    if (nrow(NST_df) > 0) {
      nst_vals <- mapply(function(g, a) contribution_one(g, a, cl, "NST"),
                         NST_df$gene, NST_df$alts, SIMPLIFY = TRUE)
      part_nst <- mean(nst_vals)
      nst_hits <- sum(nst_vals > 0)
    } else {
      part_nst <- NA_real_
      nst_hits <- NA_integer_
    }
    
    # Enforce minimum genes per side
    if (nrow(ILC_df) < min_genes_per_side) { part_ilc <- NA_real_; ilc_hits <- NA_integer_ }
    if (nrow(NST_df) < min_genes_per_side) { part_nst <- NA_real_; nst_hits <- NA_integer_ }
    
    combined <- (ifelse(is.na(part_ilc), 0, part_ilc)) - (ifelse(is.na(part_nst), 0, part_nst))
    
    data.frame(CellLine = cl,
               ILC_tests = nrow(ILC_df), ILC_hits = ilc_hits, ILC_component = round(as.numeric(part_ilc), 3),
               NST_tests = nrow(NST_df), NST_hits = nst_hits, NST_component = round(as.numeric(part_nst), 3),
               Combined_score = round(as.numeric(combined), 3),
               stringsAsFactors = FALSE)
  })
  
  out <- do.call(rbind, out_list)
  rownames(out) <- out$CellLine
  
  # Optional 0..1 scaling
  if (scale01) {
    rng <- range(out$Combined_score, na.rm = TRUE)
    if (is.finite(rng[1]) && diff(rng) > 0) {
      out$Combined01 <- round((out$Combined_score - rng[1]) / diff(rng), 3)
    } else {
      out$Combined01 <- NA_real_
    }
  }
  return(out)
}

# ==============================================================================
# FUNCTION: Compute Consensus Score
# ==============================================================================
#' Compute consensus/combined score from multiple resemblance scores
#' 
#' @param scores Data frame with resemblance scores (rows = cell lines, cols = modalities)
#' @param cols Column names to use (default: all columns)
#' @return Vector of consensus scores
consensus_equal_weight <- function(scores, cols = NULL) {
  if (is.null(cols)) cols <- colnames(scores)
  mat <- as.matrix(scores[, cols, drop = FALSE])
  apply(mat, 1, function(x) {
    mask <- !is.na(x)
    k <- sum(mask)
    if (k == 0) return(NA_real_)
    sum(x[mask] * (1 / k))
  })
}

# ------------------------------------------------------------------------------
# SECTION 2: Select Patient and Cell Line Subsets
# ------------------------------------------------------------------------------

message("=", strrep("=", 78))
message("Selecting Patient and Cell Line Subsets")
message("=", strrep("=", 78))

# Select LumA ILC and NST tumors
LumA_tumors_subset_df <- subset(TCGA_Annots, PAM50 %in% "LumA" & `Final Pathology` %in% c("ILC", "NST"))
ILC_tumors <- subset(LumA_tumors_subset_df, `Final Pathology` == "ILC")$Case.ID
NST_tumors <- subset(LumA_tumors_subset_df, `Final Pathology` == "NST")$Case.ID

# Select non-basal cell lines (with exceptions)

if (is.null(CL_Annots$'mRNA Subtype')){
  ms <- read.delim("../3-Results/Molecular_Subtyping/Molecular_Subtyping_Summary.tsv")[,c("Name","mRNA.Subtypes")]
  rownames(ms) <- ms$Name
  CL_Annots[ms$Name, 'mRNA Subtypes'] <- ms$mRNA.Subtypes  
}

CL_Annots$`mRNA Subtypes` <- factor(CL_Annots$`mRNA Subtypes`, levels = SAMPLE_GROUPS$subtype_order)
cell_ids <- subset(CL_Annots, `mRNA Subtypes` != "Basal" | Sample %in% c("MDAMB468", "HCC1187"))$Name

message("  ILC tumors: ", length(ILC_tumors))
message("  NST tumors: ", length(NST_tumors))
message("  Cell lines: ", length(cell_ids))

# ------------------------------------------------------------------------------
# PART 1: Compute Patient-Derived Signatures and Save to Files
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("PART 1: Computing Patient-Derived Signatures from TCGA")
message("=", strrep("=", 78))

# Load OncoVar for driver gene annotations
oncoVar <- read.delim(FILES$oncovar)
rownames(oncoVar) <- oncoVar$Gene_symbol

# ------------------------------------------------------------------------------
# 1.1 CNV Signature
# ------------------------------------------------------------------------------

message("\n1.1 Computing CNV signature...")

# Load reference gene database
refGene_db <- read.delim(FILES$refgene)

# Define ILC patient CN signature (driver genes only)
chr_gene_annot <- subset(refGene_db, genename %in% subset(oncoVar, Driver.Level %in% c(2, 3, 4))$Gene_symbol)[, c("genename", "chrom")]
chr_gene_annot$chrom <- as.character(chr_gene_annot$chrom)

# Compute per-chromosome mean CN for ILC tumors
tumor_ids_cnv <- intersect(colnames(TCGA_BRCA_CN), ILC_tumors)
BRCA_CN_Mean <- t(sapply(unique(chr_gene_annot$chrom), FUN = function(x) {
  colMeans(TCGA_BRCA_CN[subset(chr_gene_annot, chrom == x)$genename, tumor_ids_cnv], na.rm = TRUE)
}))

# CNV Signature heatmap
cn_sig_ht = Heatmap(t((t(BRCA_CN_Mean))), col = annot_cols$CN_logRR, border = T, row_order = order(as.numeric(unique(chr_gene_annot$chrom))),
                    width = unit(4, "cm"), height = unit(5, "cm"), show_column_dend = F, show_row_dend = F, name = "CN",cluster_column_slices = F,
                    column_split = LumA_tumors_subset_df[tumor_ids_cnv,]$`Final Pathology`, heatmap_legend_param = heatmap_legend_param,
                    show_column_names = F, show_row_names = T)

# Save CNV gene annotation
write.table(chr_gene_annot, 
            file = file.path(signatures_dir, "ILC_CNV_Genes.tsv"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
message("  ✓ Saved CNV gene annotation to: ", file.path(signatures_dir, "ILC_CNV_Genes.tsv"))



# ------------------------------------------------------------------------------
# 1.2 Alterations (MUT+GISTIC) Signature
# ------------------------------------------------------------------------------

message("\n1.2 Computing Alterations signature...")

# Simplify GAMs (remove simple LOH/GAIN)
BRCA_CL_GAM_Simple <- BRCA_CL_GAM
for (alt in c("LOH", "GAIN")) {
  BRCA_CL_GAM_Simple[BRCA_CL_GAM_Simple == alt] <- ""
}

TCGA_BRCA_GAM_Simple <- TCGA_BRCA_GAM
for (alt in c("LOH", "GAIN")) {
  TCGA_BRCA_GAM_Simple[TCGA_BRCA_GAM_Simple == alt] <- ""
}


# Define ILC and NST enriched genes
ilc_alterations <- c("CDH1_MUT_MUT;LOH_MUT;DEL_DEL", "PTEN_MUT_MUT;LOH_DEL", 
               "FOXA1_MUT_MUT;GAIN_MUT;AMP", "TBX3_MUT_MUT;DEL_MUT;LOH") # TCGA BRCA 2015 Study, 
ilc_alterations <- c("CTNNA1_MUT_MUT;LOH_MUT;DEL", ilc_alterations) # CTNNA1 is bona fida driver of lobular phenotype based on in vivo studies
nst_alterations <- c("TP53_MUT_MUT;LOH_MUT;DEL_DEL", "MYC_AMP_MUT;AMP",  
               "GRIN2A_AMP_MUT;AMP", "EFCAB1_AMP_MUT;AMP", 
               "CSMD1_MUT_MUT;LOH_MUT;DEL", "HLA-DRB1_MUT_MUT;LOH_MUT;DEL", 
               "GATA3_MUT", "MAP3K1_MUT_MUT;GAIN_MUT;AMP", 
               "MAP2K4_MUT_MUT;GAIN_MUT;AMP") # TCGA BRCA 2015 Study


# ilc_alterations <- c("CDH1", "PTEN", "FOXA1", "TBX3") # TCGA BRCA 2015 Study
# nst_alterations <- c("TP53", "MYC",  "GRIN2A", "EFCAB1", "CSMD1", "HLA-DRB1", "GATA3", "MAP3K1", "MAP2K4") # TCGA BRCA 2015 Study

tumor_ids_alt <- intersect(colnames(TCGA_BRCA_GAM_Simple), c(ILC_tumors, NST_tumors))

# Build frequency table from GAM
freq_tbl <- build_freq_tbl_from_GAM(
  tumor_gam = TCGA_BRCA_GAM_Simple[, tumor_ids_alt], 
  genes = c(ilc_alterations, nst_alterations),
  tumor_labels = setNames(LumA_tumors_subset_df$`Final Pathology`, LumA_tumors_subset_df$Case.ID),
  ilc_label = "ILC",
  nst_label = "NST",
  min_count = 3,
  fisher_alternative = "two.sided",
  debug = FALSE,
  allowed_classes = c("MUT", "MUT;LOH", "MUT;AMP", "MUT;DEL", "MUT;GAIN", "AMP", "DEL")
)


# Alteration Signature Heatmap
genes_sorted <- unique(freq_tbl$gene[order(freq_tbl$p_ILC, decreasing = TRUE)])
t_alt_mat <- TCGA_BRCA_GAM_Simple[genes_sorted, intersect(tumor_ids_alt, ILC_tumors)]

t_alt_mat[t_alt_mat == ""] <- "WT"
t_alt_mat[grepl("^MUT", t_alt_mat)] <- "MUT"

map_values <- c("WT" = 0, "DEL" = 1, "AMP" = 2, "MUT" = 3)
t_alt_mat_numeric <- matrix(map_values[t_alt_mat], 
                            nrow = nrow(t_alt_mat), 
                            dimnames = dimnames(t_alt_mat))

mut_sig_ht = Heatmap(t_alt_mat_numeric, border = T, col = c("3"="black", "2"="#B71C1C", "1"="#0D47A1",   "0"="white" ),
                   width = unit(4, "cm"), height = unit(5, "cm"), show_column_dend = F, show_row_dend = F, name = "Alterations", cluster_column_slices = F,
                   column_split = LumA_tumors_subset_df[intersect(tumor_ids_alt, ILC_tumors),]$`Final Pathology`, heatmap_legend_param = heatmap_legend_param,
                   show_column_names = F, show_row_names = T)

# Save alterations signature
write.table(freq_tbl[, c("gene", "class", "k_ILC", "n_ILC", "p_ILC", "k_NST", "n_NST", "p_NST")], 
            file = file.path(signatures_dir, "ILC_alterations.tsv"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
message("  ✓ Saved alterations signature to: ", file.path(signatures_dir, "ILC_Alterations.tsv"))

# ------------------------------------------------------------------------------
# 1.3 RNA Signature
# ------------------------------------------------------------------------------

message("\n1.3 Computing RNA signature...")

# Load TCGA DESeq results
load(FILES$tcga_deseq)

# Prepare z-scored matrices
BRCA_CL_EXP_Z <- as.matrix(t(scale(t(BRCA_CL_EXP))))
TCGA_BRCA_EXP_Z <- as.matrix(t(scale(t(TCGA_BRCA_Log2CPM))))

# Extract top ILC-enriched DEGs
TCGA_DESeq_Results$degs$ILC_mean <- rowMeans(TCGA_BRCA_EXP_Z[TCGA_DESeq_Results$degs$gene, ILC_tumors], na.rm = TRUE)
degs <- head(TCGA_DESeq_Results$degs[order(abs(TCGA_DESeq_Results$degs$ILC_mean), decreasing = TRUE), ], 
             nrow(TCGA_DESeq_Results$degs) * 5 / 100)

# RNA Signature Heatmap
tumor_ids_rna <- intersect(colnames(TCGA_BRCA_EXP_Z), ILC_tumors)
rowlabs <- degs$gene
rowlabs[rowlabs %notin% head(degs[order(degs$padj),],10)$gene] = ""
rna_sign_ht = Heatmap(t((t(TCGA_BRCA_EXP_Z[degs$gene, tumor_ids_rna]))), col = annot_cols$RNA_zscore, border = T, row_labels = rowlabs,
                   row_split = factor(ifelse(degs$log2FoldChange > 0, "up", "dn"), levels = c("up", "dn")), cluster_row_slices = F,
                   width = unit(4, "cm"), height = unit(5, "cm"), show_column_dend = F, show_row_dend = F, name = "RNA",cluster_column_slices = F,
                   column_split = LumA_tumors_subset_df[tumor_ids_rna,]$`Final Pathology`, heatmap_legend_param = heatmap_legend_param,
                   show_column_names = F, show_row_names = T)


# Save RNA signature
write.table(degs[, c("log2FoldChange", "pvalue", "padj", "ILC_mean", "gene")], 
            file = file.path(signatures_dir, "ILC_DEGs.tsv"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
message("  ✓ Saved RNA signature to: ", file.path(signatures_dir, "ILC_DEGs.tsv"))


# ------------------------------------------------------------------------------
# 1.4 DNA Methylation Signature
# ------------------------------------------------------------------------------

message("\n1.4 Computing DNAm signature...")

# Load TCGA Limma results and DNAm regulated genes
TCGA_DNAm_Regulated_Genes <- read.csv(FILES$tcga_dnam_regulated)

# Extract probes with strong RNA-DNAm correlation
probes_df <- TCGA_DNAm_Regulated_Genes %>% 
  filter(abs(top_cor_probe_RNAcor) > 0.3) %>%
  transmute(probe = as.character(top_cor_probe), 
            gene = as.character(gene), 
            direction = ifelse(top_probe_logFoldChange > 0, "up", "dn"), 
            region = top_probe_region) %>%
  filter(!is.na(probe), probe != "", !is.na(gene), gene != "") %>%
  distinct(gene, .keep_all = TRUE)


# DNAm Signature Heatmap
tumor_ids_dnam <- intersect(colnames(TCGA_BRCA_DNAm), ILC_tumors)
rowlabs <- paste0(probes_df$gene, "(", probes_df$region, ")")
rowlabs[probes_df$gene %notin% degs$gene] = ""

dnam_sign_ht = Heatmap(t(scale(t(TCGA_BRCA_DNAm[probes_df$probe, ])))[, tumor_ids_dnam], col = annot_cols$DNAm_zscore, border = T, row_labels = rowlabs,
                    row_split = factor(probes_df$direction,  levels = c("up", "dn")), cluster_row_slices = F,
                    width = unit(4, "cm"), height = unit(5, "cm"), show_column_dend = F, show_row_dend = F, cluster_column_slices = F,
                    column_split = LumA_tumors_subset_df[tumor_ids_dnam,]$`Final Pathology`, name = "DNAm", heatmap_legend_param = heatmap_legend_param,
                    show_column_names = F, show_row_names = T)



# Save DNAm signature
write.table(probes_df, 
            file = file.path(signatures_dir, "ILC_DNAm_Probes.tsv"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
message("  ✓ Saved DNAm signature to: ", file.path(signatures_dir, "ILC_DNAm_Probes.tsv"))

# ------------------------------------------------------------------------------
# 1.5 RPPA Signature
# ------------------------------------------------------------------------------

message("\n1.5 Computing RPPA signature...")

# Normalize protein names
normalize <- function(x) toupper(gsub("[^A-Z0-9]", "", x))

TCGA_BRCA_RPPA_ <- TCGA_BRCA_RPPA
norm1 <- normalize(rownames(TCGA_BRCA_RPPA_))
rownames(TCGA_BRCA_RPPA_) <- norm1

BRCA_CL_RPPA_ <- BRCA_CL_RPPA
norm2 <- normalize(rownames(BRCA_CL_RPPA_))
rownames(BRCA_CL_RPPA_) <- norm2

commonProteins <- intersect(norm1, norm2)

# Find differential proteins
groups_rppa <- setNames(LumA_tumors_subset_df$`Final Pathology`, LumA_tumors_subset_df$Case.ID)
match_idx <- match(LumA_tumors_subset_df$Case.ID, colnames(TCGA_BRCA_RPPA_))
groups_rppa <- groups_rppa[!is.na(match_idx)]

diff_prots <- find_diff_proteins(TCGA_BRCA_RPPA_[commonProteins, match_idx[!is.na(match_idx)]], groups_rppa)
rownames(diff_prots) <- diff_prots$Protein
diff_prots$ILC_mean <- rowMeans(TCGA_BRCA_RPPA_[diff_prots$Protein, intersect(colnames(TCGA_BRCA_RPPA_), ILC_tumors)], na.rm = TRUE)

# Filter significant differential proteins
diff_prots_sig <- subset(diff_prots, abs(logFC) > 0.4 & adj.P.Val < 0.05)

# RPPA Signature Heatmap
tumor_ids_rppa <- intersect(colnames(TCGA_BRCA_RPPA_), ILC_tumors)
rppa_sig_ht = Heatmap(t(scale(t(TCGA_BRCA_RPPA_)))[diff_prots_sig$Protein, tumor_ids_rppa], col = annot_cols$RPPA, border = T,
                    row_split = factor(ifelse(diff_prots_sig$logFC >0, "up", "dn") ,  levels = c("up", "dn")), cluster_row_slices = F,
                    width = unit(4, "cm"), height = unit(5, "cm"), show_column_dend = F, show_row_dend = F, name = "Protein",cluster_column_slices = F,
                    column_split = LumA_tumors_subset_df[tumor_ids_rppa,]$`Final Pathology`, heatmap_legend_param = heatmap_legend_param,
                    show_column_names = F, show_row_names = T)


# Save RPPA signature
write.table(diff_prots_sig[, c("logFC", "P.Value", "adj.P.Val", "ILC_mean", "Protein")], 
            file = file.path(signatures_dir, "ILC_DEPs.tsv"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
message("  ✓ Saved RPPA signature to: ", file.path(signatures_dir, "ILC_DEPs.tsv"))

message("\n", "=", strrep("=", 78))
message("PART 1 Complete: All patient-derived signatures saved to files")
message("=", strrep("=", 78))

# ------------------------------------------------------------------------------
# PART 2: Score Cell Lines and Compute Consensus
# ------------------------------------------------------------------------------

message("\n", "=", strrep("=", 78))
message("PART 2: Scoring Cell Lines and Computing Consensus")
message("=", strrep("=", 78))

# ------------------------------------------------------------------------------
# 2.1 CNV Resemblance Scores
# ------------------------------------------------------------------------------

message("\n2.1 Computing CNV resemblance scores...")

cnv_res <- compute_resemblance_cnv(
  tumor_cn = t(scale(t(TCGA_BRCA_CN[chr_gene_annot$genename, ]))),
  cell_cn = t(scale(t(BRCA_CL_CN[chr_gene_annot$genename, ]))),
  ga = chr_gene_annot,
  tumor_samples = ILC_tumors,
  cell_samples = unique(c(cell_ids, "SKBR3-C")),
  B = 1000,
  method = "spearman",
  seed = 123,
  return_distributions = FALSE
)

# cnv_res$CellLine <- gsub('SKBR3-C', "SKBR3-I", cnv_res$CellLine)
rownames(cnv_res) <- cnv_res$CellLine
message("  ✓ CNV scores computed for ", nrow(cnv_res), " cell lines")

# ------------------------------------------------------------------------------
# 2.2 Alterations Resemblance Scores
# ------------------------------------------------------------------------------

message("\n2.2 Computing Alterations resemblance scores...")

mut_res <- compute_alteration_resemblance(
  mut_mat = BRCA_CL_GAM_Simple[, intersect(colnames(BRCA_CL_GAM_Simple), cell_ids)],
  ILC_specs = ilc_alterations,
  NST_specs = nst_alterations,
  freq_tbl = freq_tbl,
  scale01 = TRUE,
  min_genes_per_side = 1
)

mut_res$CombinedRank <- rank01(mut_res$Combined_score)
message("  ✓ Alterations scores computed for ", nrow(mut_res), " cell lines")

# ------------------------------------------------------------------------------
# 2.3 RNA Resemblance Scores
# ------------------------------------------------------------------------------

message("\n2.3 Computing RNA resemblance scores...")

rna_res <- compute_resemblance(
  TCGA_BRCA_EXP_Z, ILC_tumors, BRCA_CL_EXP_Z, cell_ids, degs$gene, 
  n_iter = 1000, n_feature = 0.8 * length(degs$gene), "spearman"
)
message("  ✓ RNA scores computed for ", nrow(rna_res), " cell lines")

# ------------------------------------------------------------------------------
# 2.4 DNA Methylation Resemblance Scores
# ------------------------------------------------------------------------------

message("\n2.4 Computing DNAm resemblance scores...")

dnam_res <- compute_resemblance(
  t(scale(t(TCGA_BRCA_DNAm[probes_df$probe, ]))), ILC_tumors,
  t(scale(t(BRCA_CL_DNAm[probes_df$probe, ]))), cell_ids, probes_df$probe,
  n_iter = 1000, n_feature = 0.8 * length(probes_df$probe), "spearman"
)
message("  ✓ DNAm scores computed for ", nrow(dnam_res), " cell lines")

# ------------------------------------------------------------------------------
# 2.5 RPPA Resemblance Scores
# ------------------------------------------------------------------------------

message("\n2.5 Computing RPPA resemblance scores...")

rppa_res <- compute_resemblance(
  t(scale(t(TCGA_BRCA_RPPA_))), ILC_tumors,
  t(scale(t(BRCA_CL_RPPA_))), cell_ids, diff_prots_sig$Protein,
  n_iter = 1000, n_feature = 0.8 * length(diff_prots_sig$Protein), "spearman"
)
message("  ✓ RPPA scores computed for ", nrow(rppa_res), " cell lines")

# ------------------------------------------------------------------------------
# 2.6 Combine Scores and Compute Consensus
# ------------------------------------------------------------------------------

message("\n2.6 Combining scores and computing consensus...")

# Combine all scores
scores <- data.frame(
  RNA = round(rna_res[cell_ids, "BootMedian"], 1), 
  DNAm = round(dnam_res[cell_ids, "BootMedian"], 1), 
  CNV = round(cnv_res[cell_ids, "BootMedian"], 1), 
  RPPA = round(rppa_res[cell_ids, "BootMedian"], 1), 
  MUT = round(mut_res[cell_ids, "Combined_score"], 1)
)

# Rank normalize each modality
molecular_resemblance_scores <- as.data.frame(apply(scores, 2, rank01))
rownames(molecular_resemblance_scores) <- cell_ids

# Compute consensus score (equal weight average)
molecular_resemblance_scores$CONSENSUS <- rank01(round(consensus_equal_weight(
  molecular_resemblance_scores, 
  cols = colnames(molecular_resemblance_scores)
), 1))

# Save resemblance scores
write.table(molecular_resemblance_scores, quote = FALSE, col.names = NA,
            file = file.path(DIRS$results_sub$molecular_resemblance, "Resemblance_Scores.tsv"), 
            row.names = TRUE, sep = "\t")
message("  ✓ Saved resemblance scores to: ", file.path(DIRS$results_sub$molecular_resemblance, "Resemblance_Scores.tsv"))

# ------------------------------------------------------------------------------
# 2.7 Generate Fig 6 - Resemblance Scores Heatmap
# ------------------------------------------------------------------------------

message("\n========================================")
message("Figure 6: Resemblance scores heatmap")
message("========================================")

# Prepare matrix for heatmap (filter to ICLE cell lines with "-I" suffix)
mat <- as.matrix(molecular_resemblance_scores[, c("CNV", "MUT", "DNAm", "RPPA", "RNA", "CONSENSUS")])
mat <- mat[grep("-I", rownames(mat)), ]

# Order cell lines by consensus score
row_order <- order(mat[, 6], decreasing = TRUE)

# Clean rownames (remove "-I" suffix)
rownames(mat) <- gsub("-I", "", rownames(mat))

# Main heatmap: scores 1-5
resem_ht <- Heatmap(
  mat[, -6, drop = FALSE],
  name = "Score",
  col = annot_cols$Score,
  na_col = "gray",
  height = unit(10, "cm"),
  width = unit(6, "cm"),
  row_order = row_order,
  row_split = CL_Annots[paste0(rownames(mat), "-I"), "mRNA Subtypes"],
  show_row_dend = TRUE,
  column_title = "Resemblance scores",
  row_title = "Cell lines",
  border = TRUE,
  border_gp = gpar(col = "black"),
  heatmap_legend_param = heatmap_legend_param,
  row_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
  row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_names_side = "right",
  column_names_side = "bottom",
  left_annotation = rowAnnotation(
    Subtype = CL_Annots[paste0(rownames(mat), "-I"), "mRNA Subtypes"],
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold", col = "black", fontfamily = helv),
    annotation_legend_param = heatmap_legend_param,
    gp = gpar(col = NA, lwd = 0.6),
    simple_anno_size = unit(4, "mm"),
    col = list(Subtype = annot_cols$Subtypes)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- mat[i, j]
    if (!is.na(val)) {
      lab <- sprintf("%.1f", val)
      grid.text(lab, x, y, gp = gpar(fontsize = 10, col = "black"))
    }
  }
)

# Side heatmap: CONSENSUS column
final <- mat[, 6, drop = FALSE]
colnames(final) <- "CONSENSUS"

combined_resem_ht <- Heatmap(
  final,
  name = "Score",
  col = annot_cols$Score,
  na_col = "gray",
  height = unit(10, "cm"),
  width = unit(1.5, "cm"),
  row_order = row_order,
  row_split = CL_Annots[paste0(rownames(final), "-I"), "mRNA Subtypes"],
  show_row_dend = FALSE,
  column_title = "Combined",
  row_title = " ",
  border = TRUE,
  border_gp = gpar(col = "black"),
  heatmap_legend_param = heatmap_legend_param,
  row_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black", fontfamily = helv),
  column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
  row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = helv),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_names_side = "right",
  column_names_side = "bottom",
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- final[i, j]
    if (!is.na(val)) {
      lab <- sprintf("%.1f", val)
      grid.text(lab, x, y, gp = gpar(fontsize = 10, col = "black"))
    }
  }
)

# Combine heatmaps
fig6_resemblance_ht <- resem_ht + combined_resem_ht

# Assign to global environment (saved in Main_Data_Analysis.Rmd)
assign("cn_sig_ht", cn_sig_ht, envir = .GlobalEnv)
assign("mut_sig_ht", mut_sig_ht, envir = .GlobalEnv)
assign("dnam_sign_ht", dnam_sign_ht, envir = .GlobalEnv)
assign("rna_sign_ht", rna_sign_ht, envir = .GlobalEnv)
assign("rppa_sig_ht", rppa_sig_ht, envir = .GlobalEnv)
assign("fig6_resemblance_ht", fig6_resemblance_ht, envir = .GlobalEnv)
message("  ✓ Fig 6 complete (fig6_resemblance_ht and signature heatmaps assigned).\n")

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/23_Fig6_Patient_Signatures_Resemblance_Scores.R")
#
# load_all_icle_data(load_external = TRUE)
# # Script runs on source. Uses FILES$oncovar, refgene, tcga_deseq, tcga_limma, etc.
#
# Outputs: 
#   - Patient signatures in DIRS$results_sub$molecular_resemblance/ILC_Signatures/
#   - Resemblance_Scores.tsv in DIRS$results_sub$molecular_resemblance/
#   - Fig6 heatmap assigned to global environment for saving in Main_Data_Analysis.Rmd
