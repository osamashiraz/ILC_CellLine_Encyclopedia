# ==============================================================================
# Helper Functions for ICLE Analysis Pipeline
# ==============================================================================
# Description: Centralized utility functions, color schemes, themes, and helper
#              functions used across the ICLE analysis pipeline.
#
# Dependencies:
#   - config.R (must be sourced first for FILES and DIRS)
#   - Various R packages (loaded as needed)
#
# Usage:
#   source("config.R")
#   source("Helper_Scripts/Helper_Functions.R")
#
# Author: Osama Shiraz Shah
# ==============================================================================

# ==============================================================================
# SECTION 1: PACKAGE LOADING & INITIALIZATION
# ==============================================================================

library(ggpubr)
library(ggthemes)
library(ComplexHeatmap)
library(dplyr)

# Load reference datasets from config paths
# Note: Requires config.R to be sourced first
if (exists("FILES")) {
  load(FILES$biomart_genemap)
  oncoVar <- read.delim(FILES$oncovar)
  rownames(oncoVar) <- oncoVar$Gene_symbol
} else {
  warning("config.R not loaded. Please source config.R before using Helper_Functions.R")
}

# ==============================================================================
# SECTION 2: GLOBAL UTILITIES & OPERATORS
# ==============================================================================
# Used in: All scripts (general utility).
# Sankey: Fig 1B uses inline code in 03_Fig1B_Multiomics_Sankey.R; additional
#         Sankey helpers are in SECTION 12 below.
# ------------------------------------------------------------------------------

# Negation operator for %in%
`%notin%` <- Negate(`%in%`)

# Font definition
helv <- "Helvetica"

# ==============================================================================
# SECTION 3: COMPLEXHEATMAP CONFIGURATION
# ==============================================================================
# Used in: All scripts generating ComplexHeatmap visualizations
#          - 02_Fig1_SupFig2_3_4_Molecular_Subtyping.R
#          - 04_Fig1D_Multiomics_Overview.R
#          - 21_Fig4_SupFig11_DNAm_Alterations.R
#          - 21_Fig3E_3F_SupFig10_SV_Fusions.R
#          - And many others
# ------------------------------------------------------------------------------

heatmap_legend_param <- list(
  title_gp   = gpar(fontsize = 15, fontface = "bold", fontfamily = helv),
  labels_gp  = gpar(fontsize = 12, fontfamily = helv), 
  grid_width = unit(4, "mm"), 
  grid_height = unit(4, "mm"),
  border     = "black", 
  legend_gp  = gpar(col = "black", lwd = 0.6)
)

# ==============================================================================
# SECTION 4: COLOR SCHEMES & ANNOTATION COLORS
# ==============================================================================
# Used in: All scripts generating heatmaps and visualizations
#          - annot_cols: Used in ComplexHeatmap annotations across all figures
#          - Color functions: Used for data type-specific color mappings
# ------------------------------------------------------------------------------

annot_cols <- list(
  # Study/cohort annotations
  Study = c("ICLE" = "black", "Marcotte" = "gray", "CCLE" = "lightgray", "Other" = "gray"),
  
  # Sample type annotations
  Type = c("Patient Tumors" = "#548134", "Cell Lines" = "#a8d08e", "Normal" = "gray"),
  
  # Histology annotations
  Histology = c("ILC" = "#D32F2F", "ILC-like" = "#AA00FF", "NST" = "#0066ff", 
                "Normal" = "gray", "Other" = "#616161", "Immortalized normal" = "#9c7a3b"),
  
  # PAM50 subtyping annotations
  PAM50 = c("LumA" = "#0D47A1", "LumB" = "#00BCD4", "LumA;LumB" = "#2196F3", 
            "LumA;Normal" = "#CDDC39", "LumA;Her2" = "#B39DDB", "Her2" = "#ffcdd2", 
            "Normal" = "#C8E6C9", "Basal" = "#FF5722"),
  
  # mRNA subtype annotations
  Subtypes = c("Lum" = "#039BE5", "HER2" = "#ffcdd2", "Lum/HER2" = "#009688", 
               "Basal" = "#FF5722"),
  
  # Receptor status annotations
  ER = c("ER-" = "white", "ER+" = "gray30"),
  HER2 = c("HER2-" = "white", "HER2+" = "gray30"),
  
  # DNA methylation region annotations
  DNAm_regions = c("intron" = "#BDBDBD", "exon" = "#AED581", "multiple" = "#006064", 
                   "fiveUTRs" = "#1976D2", "TSS_1k" = "#0D47A1", "threeUTRs" = "#E91E63", 
                   "TES_1k" = "#AD1457"),
  
  # Structural variant type annotations
  SV = c("Del" = "#E64B35FF", "Dup" = "#4DBBD5FF", "Ins" = "#00A087FF",
         "Inv" = "#3C5488FF", "Tranloc-inter" = "#F39B7FFF", 
         "Tranloc-intra" = "#8491B4FF"),
  
  # Mutation type annotations
  Mutation = c("Truncating" = "black", "Missense" = "#43A047", "Inframe" = "#377EB8", 
               "Silent" = "gray", "Other" = "#9575CD"),
  
  # Alteration type annotations (simplified)
  Alt_Simple = c("MUT" = "black", "AMP" = "#B71C1C", "DEL" = "#0D47A1", "WT" = "white"),
  
  # Alteration type annotations (detailed)
  Alt = c("LOH" = "#64b5f6", "GAIN" = "#f06292", "AMP" = "#B71C1C", 
          "DEL" = "#0D47A1", "MUT+LOH" = "black", "MUT+GAIN" = "#FB8C00", 
          "MUT" = "#bf7537", "WT" = "gray"),
  
  # Alteration type annotations (FMI)
  Alt_FMI = c("DEL" = "#0D47A1", "GAIN" = "#f06292", "RE" = "#4DB6AC", 
              "MUT" = "#bf7537", "WT" = "gray"),
  
  # Continuous color mappings using circlize colorRamp2
  CN = circlize::colorRamp2(c(-2, -1, 0, 1, 2), 
                           colors = c("#0D47A1", "#64B5F6", "#FAFAFA", "#F06292", "red")),
  
  CN_logRR = circlize::colorRamp2(c(-1, 0, 1), 
                                  colors = c("#0D47A1", "#FAFAFA", "red")),
  
  RNA_zscore = circlize::colorRamp2(c(-2, -1, 0, 1, 2), 
                                    colors = c('#003d30', "#7dc8bf", "#EEEEEE", "#ca9446", '#543005'), 
                                    space = 'RGB'),
  
  RNA = circlize::colorRamp2(c(0, 3, 6, 9, 12), 
                            colors = c('#003d30', "#7dc8bf", "#EEEEEE", "#ca9446", '#543005'), 
                            space = 'RGB'),
  
  RPPA = circlize::colorRamp2(c(-2, -1, 0, 1, 2), 
                              colors = c('#1B5E20', "#43A047", "white", "#E53935", "#EF5350"), 
                              space = 'RGB'),
  
  DNAm = circlize::colorRamp2(c(0, 0.15, 0.3, 0.5, 0.7, 0.85, 1), 
                             colors = c("#08306B", "#4292C6", "#B3E5FC", "white", 
                                       "#FFCCBC", "#FB6A4A", "#99000D")),
  
  DNAm_zscore = circlize::colorRamp2(c(-2, -1, 0, 1, 2), 
                                     colors = c("#08306B", "#4292C6", "white", "#FB6A4A", "#99000D")),
  
  Score = circlize::colorRamp2(c(0, 0.5, 1), 
                              colors = c("#4292C6", "white", "#EF9A9A")),
  
  SET = circlize::colorRamp2(c(-2, 0, 2), c("black", "white", "gold"), space = 'RGB'),
  
  Correlation = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), 
                                     colors = c("blue", "#81D4FA", "white", "#EF9A9A", "red"))
)

# ==============================================================================
# SECTION 5: MUTATION CLASSIFICATION
# ==============================================================================
# Used in: Scripts processing mutation data
#          - 13_Fig2H_CDH1_Alteration_Landscape.R
#          - 08_Fig2C_CDH1_exonic_deletions.R
#          - And other mutation analysis scripts
# ------------------------------------------------------------------------------

# Mapping of TCGA mutation variant classifications to simplified categories
broad_variant_classes <- c(
  "In_Frame_Del" = "Inframe",
  "In_Frame_Ins" = "Inframe",
  "Silent" = "Silent",
  "Targeted_Region" = "Inframe",
  "Missense_Mutation" = "Missense",
  "Frame_Shift" = "Truncating",
  "Frame_Shift_Del" = "Truncating",
  "Frame_Shift_Ins" = "Truncating",
  "Nonsense_Mutation" = "Truncating",
  "Nonstop_Mutation" = "Truncating",
  "Splice_Region" = "Truncating",
  "Splice_Site" = "Truncating",
  "3'Flank" = "Other",
  "3'UTR" = "Other",
  "5'Flank" = "Other",
  "5'UTR" = "Other",
  "De_novo_Start_InFrame" = "Other",
  "De_novo_Start_OutOfFrame" = "Other",
  "Fusion" = "Other",
  "IGR" = "Other",
  "Intron" = "Other",
  "lincRNA" = "Other",
  "RNA" = "Other",
  "Start_Codon_Del" = "Other",
  "Start_Codon_Ins" = "Other",
  "Start_Codon_SNP" = "Other",
  "Stop_Codon_Del" = "Other",
  "Stop_Codon_Ins" = "Other",
  "Translation_Start_Site" = "Other",
  "Unknown" = "Other"
)

# ==============================================================================
# SECTION 6: PLOT SUPPRESSION FUNCTIONS
# ==============================================================================
# Used in: Scripts that generate intermediate plots during feature selection
#          and consensus clustering
#          - 02_Fig1_SupFig2_3_4_Molecular_Subtyping.R (FSbyMAD, FSbyVar, ConsensusClusterPlus)
#          - 04_Fig1D_Multiomics_Overview.R (FSbyMAD, FSbyVar)
# ------------------------------------------------------------------------------

# Optional: suppress package startup messages and warnings (if flags are set in Main_Data_Analysis.Rmd)
if (exists("SUPPRESS_WARNINGS", envir = .GlobalEnv) && get("SUPPRESS_WARNINGS", envir = .GlobalEnv)) {
  old_warn <- options(warn = -1)
  on.exit(options(old_warn), add = TRUE)
}
if (exists("SUPPRESS_PKG_MESSAGES", envir = .GlobalEnv) && get("SUPPRESS_PKG_MESSAGES", envir = .GlobalEnv)) {
  lib_orig <- base::library
  assign("library", function(..., warn.conflicts = TRUE) {
    suppressPackageStartupMessages(lib_orig(..., warn.conflicts = warn.conflicts))
  }, envir = .GlobalEnv)
}

#' Helper function to suppress plots by redirecting to null device
#' 
#' @param expr Expression to evaluate (lazy evaluation)
#' @return Result of expression evaluation
#' @details
#'   Opens a null PDF device, evaluates the expression, then closes the device.
#'   This suppresses any plots generated during evaluation.
#'   Uses lazy evaluation to capture the unevaluated expression.
#' 
#' @examples
#'   quiet_run(CancerSubtypes::FSbyMAD(mat, cut.type = "topk", value = 6000))
quiet_run <- function(expr) {
  pdf(NULL)                 # Open dummy device
  on.exit(dev.off())        # Ensure it closes when done
  # Use substitute to capture the call, then eval in parent frame
  result <- eval(substitute(expr), envir = parent.frame())
  return(result)
}

# Conditional wrappers for feature selection functions
if (SUPPRESS_FEATURE_SELECTION_PLOTS) {
  # Wrapper that conditionally suppresses plots
  FSbyMAD_quiet <- function(...) {
    if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
        get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
      return(quiet_run(CancerSubtypes::FSbyMAD(...)))
    } else {
      return(CancerSubtypes::FSbyMAD(...))
    }
  }
  
  FSbyVar_quiet <- function(...) {
    if (exists("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv) && 
        get("SUPPRESS_FEATURE_SELECTION_PLOTS", envir = .GlobalEnv)) {
      return(quiet_run(CancerSubtypes::FSbyVar(...)))
    } else {
      return(CancerSubtypes::FSbyVar(...))
    }
  }
  
  # Assign to global environment for use in downstream scripts
  assign("FSbyMAD_quiet", FSbyMAD_quiet, envir = .GlobalEnv)
  assign("FSbyVar_quiet", FSbyVar_quiet, envir = .GlobalEnv)
} else {
  # If suppression is disabled, create passthrough functions
  FSbyMAD_quiet <- function(...) CancerSubtypes::FSbyMAD(...)
  FSbyVar_quiet <- function(...) CancerSubtypes::FSbyVar(...)
  assign("FSbyMAD_quiet", FSbyMAD_quiet, envir = .GlobalEnv)
  assign("FSbyVar_quiet", FSbyVar_quiet, envir = .GlobalEnv)
}

# Conditional wrapper for ConsensusClusterPlus
if (SUPPRESS_CONSENSUS_CLUSTER_PLOTS) {
  ConsensusClusterPlus_quiet <- function(...) {
    if (exists("SUPPRESS_CONSENSUS_CLUSTER_PLOTS", envir = .GlobalEnv) && 
        get("SUPPRESS_CONSENSUS_CLUSTER_PLOTS", envir = .GlobalEnv)) {
      return(quiet_run(ConsensusClusterPlus::ConsensusClusterPlus(..., plot = NULL)))
    } else {
      return(ConsensusClusterPlus::ConsensusClusterPlus(...))
    }
  }
  assign("ConsensusClusterPlus_quiet", ConsensusClusterPlus_quiet, envir = .GlobalEnv)
} else {
  ConsensusClusterPlus_quiet <- function(...) ConsensusClusterPlus::ConsensusClusterPlus(...)
  assign("ConsensusClusterPlus_quiet", ConsensusClusterPlus_quiet, envir = .GlobalEnv)
}

# ==============================================================================
# SECTION 7: DNA METHYLATION FUNCTIONS
# ==============================================================================
# Used in: Scripts processing DNA methylation data
#          - 21_Fig4_SupFig11_DNAm_Alterations.R
#          - 04_Fig1D_Multiomics_Overview.R
# ------------------------------------------------------------------------------

#' Load and simplify HM450K probe set annotations
#' 
#' @param probe_file Path to probe set RData file (optional, uses FILES$hm450k_probeset if NULL)
#' @return Data frame with HM450K probe annotations including simplified region labels
#' @details
#'   Loads the HM450K probe set and creates a simplified region annotation
#'   by collapsing complex multi-region annotations into single categories.
#'   Used for DNA methylation analysis and visualization.
#' 
#' @examples
#'   HM450K_ProbeSet <- load_DNAm_probeset()
#'   HM450K_ProbeSet <- load_DNAm_probeset(FILES$hm450k_probeset)
load_DNAm_probeset <- function(probe_file = NULL) {
  # Use config path if not specified
  if (is.null(probe_file)) {
    if (!exists("FILES")) {
      stop("config.R must be loaded or probe_file path must be provided")
    }
    probe_file <- FILES$hm450k_probeset
  }
  
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

# ==============================================================================
# SECTION 8: CIRCOS PLOT FUNCTIONS
# ==============================================================================
# Used in: Scripts generating interactive circos plots
#          - 20_Fig3D_Circos_Selected_Samples.R
# ------------------------------------------------------------------------------

#' Save interactive circos plot to PDF
#' 
#' @param tracks Circos plot tracks object (from interacCircos)
#' @param sample_name Sample name for output file naming
#' @param output_dir Directory to save PDF file
#' @param width HTML widget width in pixels (default: 800)
#' @param height HTML widget height in pixels (default: 800)
#' @return HTML widget object (invisibly)
#' @details
#'   Converts interactive HTML circos plots (from interacCircos) to static PDF
#'   files. First saves as HTML using htmlwidgets, then converts to PDF using
#'   webshot2. The intermediate HTML file is automatically cleaned up.
#' 
#' @examples
#'   save_circos_to_pdf(trk_bck4, "BCK4", fig3d_dir)
save_circos_to_pdf <- function(tracks, sample_name, output_dir, width = 800, height = 800) {
  # Create the circos plot widget
  circos_widget <- draw_circos_plot(tracks)
  
  # Save as HTML first (intermediate step)
  html_file <- file.path(output_dir, paste0("Fig3D_Circos_", sample_name, ".html"))
  htmlwidgets::saveWidget(circos_widget, file = html_file, selfcontained = TRUE)
  
  # Convert HTML to PDF using webshot2
  pdf_file <- file.path(output_dir, paste0("Fig3D_Circos_", sample_name, ".pdf"))
  tryCatch({
    webshot2::webshot(
      url = html_file, 
      file = pdf_file, 
      vwidth = width,
      vheight = height,
      delay = 3,         # Wait for JS to render
      zoom = 2,          # 2x resolution for better quality
      selector = ".html-widget"
    )
    message("  ✓ Saved PDF: ", basename(pdf_file))
    # Clean up HTML file
    unlink(html_file)
  }, error = function(e) {
    warning("Failed to save PDF for ", sample_name, ": ", e$message)
  })
  
  return(circos_widget)
}

# ==============================================================================
# SECTION 9: GGPLOT2 THEMES
# ==============================================================================
# Used in: Scripts generating ggplot2 visualizations
#          - 21_Fig4_SupFig11_DNAm_Alterations.R (myTheme)
#          - 22_Fig5_RNAi_Differential_Dependencies.R (myTheme)
#          - 05_Fig1F_Alteration_barplots.R (myTheme)
#          - 22_Fig5_RNAi_Differential_Dependencies.R (myTheme_barplot via myEnrich_barplot)
# ------------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)

#' Custom ggplot2 theme for ICLE analysis figures
#' 
#' @param text_size Base text size (default: 20)
#' @param axis_text_size Axis text size (default: 15)
#' @param legend_text_size Legend text size (default: 15)
#' @param subtitle_size Subtitle text size (default: 8)
#' @return ggplot2 theme object
#' @details
#'   Provides consistent styling for ICLE analysis plots with Helvetica font,
#'   clean backgrounds, and customizable text sizes.
#' 
#' @examples
#'   ggplot(data, aes(x, y)) + geom_point() + myTheme(20)
#'   ggplot(data, aes(x, y)) + geom_point() + myTheme(15, 12, 12, 10)
myTheme <- function(text_size = 20, axis_text_size = 15, legend_text_size = 15, subtitle_size = 8) {
  theme_bw() +
    theme(
      text = element_text(size = text_size, family = "Helvetica"),
      axis.text = element_text(size = axis_text_size, face = "bold", color = "black", family = "Helvetica"),
      legend.text = element_text(size = legend_text_size, face = "bold", family = "Helvetica"),
      plot.subtitle = element_text(size = subtitle_size, face = "bold", family = "Helvetica"),
      legend.background = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_line(color = "gray60", linetype = "solid")
    )
}

#' Custom ggplot2 theme for barplot visualizations
#' 
#' @return ggplot2 theme object
#' @details
#'   Minimal theme for barplots with clean white background, no grid lines,
#'   and simplified axis styling. Used primarily in enrichment analysis plots.
#' 
#' @examples
#'   ggplot(data, aes(x, y)) + geom_bar() + myTheme_barplot
myTheme_barplot <- theme_bw() + theme(
  # TEXT
  text = element_text(size = 15, family = "Helvetica", face = "bold"), 
  legend.text = element_text(size = 10, face = "bold"),
  legend.title = element_text(size = 12, face = "bold"),
  rect = element_blank(),
  # AXES
  axis.line.y = element_line(colour = "#424242", size = 0.5),
  plot.title = element_text(hjust = 0.5),
  # PANEL
  panel.background = element_rect(fill = "white"),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank()
)

# ==============================================================================
# SECTION 10: GENE ENRICHMENT ANALYSIS FUNCTIONS
# ==============================================================================
# Used in: 22_Fig5_RNAi_Differential_Dependencies.R (doGSEA, myEnrich_barplot)
# ------------------------------------------------------------------------------

#' Convert gene names between different identifier types
#' 
#' @param geneNames Vector of gene identifiers
#' @param fromType Source identifier type (default: "SYMBOL")
#' @param toType Target identifier type (default: "ENTREZID")
#' @return Data frame with mapping between identifier types
#' @details
#'   Uses clusterProfiler and org.Hs.eg.db to convert between gene identifier
#'   types (e.g., gene symbols to Entrez IDs). Required for enrichment analysis.
#' 
#' @examples
#'   gene_map <- gene_name_conversion(c("CDH1", "TP53", "PIK3CA"), 
#'                                    fromType = "SYMBOL", toType = "ENTREZID")
gene_name_conversion <- function(geneNames, fromType = "SYMBOL", toType = "ENTREZID") {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  map <- bitr(geneID = geneNames, fromType = fromType, toType = toType, OrgDb = org.Hs.eg.db)
  map <- map[!duplicated(map[, 1]), ]
  rownames(map) <- map[, 1]
  geneNames <- intersect(map[, 1], geneNames)
  return(map)
}

#' Create enrichment analysis barplot
#' 
#' @param df Data frame with enrichment results (must have columns: label, pval/fdr, overlap)
#' @param order Logical, whether to order results (default: TRUE)
#' @param order_by Column to order by (default: "logpval")
#' @param plot_top_n Number of top results to plot (default: 10)
#' @param plot_fdr Logical, whether to plot FDR instead of p-value (default: TRUE)
#' @param overlap Minimum overlap threshold (default: 3)
#' @param abr Logical, whether to abbreviate labels (default: FALSE)
#' @param abr_len Minimum abbreviation length (default: 50)
#' @param title Plot title (default: "")
#' @param abline_val FDR/p-value cutoff for vertical line (default: 0.05)
#' @param title_color Title color (default: "#0097A7")
#' @param font_color Font color (default: "black")
#' @param bar_color Bar fill color (default: "#80DEEA")
#' @param add_cutoff Logical, whether to add vertical cutoff line (default: FALSE)
#' @return List with plot object (plt) and data frame (df)
#' @details
#'   Creates a horizontal barplot of enrichment analysis results. Supports both
#'   custom enrichment data frames and clusterProfiler output format.
#'   Used in conjunction with doGSEA for pathway enrichment visualization.
#' 
#' @examples
#'   enrich_plot <- myEnrich_barplot(enrich_df, plot_top_n = 10, 
#'                                   title = "KEGG Pathways", bar_color = "#29B6F6")
myEnrich_barplot <- function(df, order = TRUE, order_by = "logpval", plot_top_n = 10, 
                            plot_fdr = TRUE, overlap = 3, abr = FALSE, abr_len = 50, 
                            title = "", abline_val = 0.05, title_color = "#0097A7", 
                            font_color = "black", bar_color = "#80DEEA", add_cutoff = FALSE) {
  # Filter by overlap threshold
  df <- subset(df, overlap >= overlap)
  df <- df[!duplicated(df$label), ]
  
  # Handle clusterProfiler output format
  if (sum(names(df) %in% "label") == 0) {
    df$label <- df$Description
    df$pval <- df$pvalue
    df$fdr <- df$p.adjust
  }
  
  # Abbreviate labels if requested
  if (abr) {
    df$abr_label <- abbreviate(df[, "label"], minlength = abr_len)
  } else {
    df[, "abr_label"] <- df[, "label"]
  }
  
  # Calculate -log10 p-value or FDR
  if (plot_fdr) {
    df[, "logpval"] <- -log10(df[, "fdr"])
  } else {
    df[, "logpval"] <- -log10(df[, "pval"])
  }
  
  # Order results
  if (order) {
    df <- df[order(df[, order_by], decreasing = TRUE), ]
  }
  
  # Select top N results
  df_sub <- na.omit(df[1:plot_top_n, ])
  hide_NA <- NULL
  
  # Pad if fewer than plot_top_n results
  if (nrow(df_sub) < plot_top_n) {
    pad_length <- plot_top_n - nrow(df_sub)
    padding <- as.data.frame(matrix(data = NA, nrow = pad_length, ncol = ncol(df_sub)))
    colnames(padding) <- colnames(df_sub)
    padding$logpval <- 0
    padding$abr_label <- as.character(1:pad_length)
    hide_NA <- theme(axis.text.x = element_text(
      colour = c(rep("#424242", nrow(df_sub)), rep("#FAFAFA", pad_length))
    ))
    df_sub <- rbind(df_sub, padding)
  }
  
  # Add cutoff line if requested
  if (add_cutoff == TRUE) {
    geomVline <- geom_vline(xintercept = -log10(abline_val), linetype = 'dashed', 
                            color = "darkgray", alpha = 0.6, size = 1)
  } else {
    geomVline <- NULL
  }
  
  # Adjust font size based on label length
  fontsize <- ifelse(nchar(df_sub$label) > 20, 4, 5)
  
  # Create plot
  plt <- ggplot(df_sub, aes(y = reorder(abr_label, logpval), x = logpval, label = label)) + 
    geomVline +
    geom_bar(stat = "identity", fill = bar_color, alpha = 0.5, width = 0.75) + 
    myTheme_barplot + 
    ggtitle(title) +
    geom_text(aes(x = 0.05), size = fontsize, face = "bold", hjust = 0, vjust = 0.4, 
              alpha = 0.6, color = font_color) + 
    scale_x_continuous(expand = c(0, 0.2)) + 
    theme(title = element_text(size = 20, color = title_color)) +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
    hide_NA
  
  plt <- plt + theme(axis.title.x = element_text(size = 20)) + xlab("-log10 Pvalue")
  
  return(list(plt = plt, df = df))
}

#' Perform gene set enrichment analysis (GSEA) using hypeR
#' 
#' @param gset_db Database name (default: "KEGG")
#' @param gset_db_obj Gene set database object (from hypeR)
#' @param input_list Named list of gene vectors (e.g., list(NST = nst_genes, ILC = icle_genes))
#' @param plot_top_n Number of top pathways to plot (default: 10)
#' @param add_cutoff Logical, whether to add FDR cutoff line (default: FALSE)
#' @return List with enrichment plots for each input group (lobular, ductal)
#' @details
#'   Performs hypergeometric test-based enrichment analysis using hypeR and
#'   generates barplots for each group in the input list. Designed for comparing
#'   enrichment between ILC and NST samples.
#' 
#' @examples
#'   KEGG_ALL <- doGSEA(gset_db = "KEGG", KEGG, 
#'                      list(NST = NST_genes, ILC = ICLE_genes), 
#'                      plot_top_n = 7)
doGSEA <- function(gset_db = "KEGG", gset_db_obj, input_list = input_list, 
                  plot_top_n = 10, add_cutoff = FALSE) {
  hyp_obj <- hypeR(input_list, gset_db_obj, test = "hypergeometric",
                   pval = 0.05, plotting = FALSE)
  
  ductal_obj <- myEnrich_barplot(hyp_obj$data$`NST`$data, plot_top_n = plot_top_n, abr = TRUE, 
                                title = paste0("NST - ", gset_db), bar_color = "#29B6F6",
                                title_color = "#29B6F6", plot_fdr = FALSE, add_cutoff = add_cutoff)
  
  lobular_obj <- myEnrich_barplot(hyp_obj$data$ILC$data, plot_top_n = plot_top_n, abr = TRUE, 
                                 title = paste0("ICLE - ", gset_db), bar_color = "#EF5350",
                                 title_color = "#EF5350", plot_fdr = FALSE, add_cutoff = add_cutoff)
  
  return(list(lobular = lobular_obj, ductal = ductal_obj))
}

# ==============================================================================
# SECTION 11: GAM ANALYSIS FUNCTIONS
# ==============================================================================
# Used in: 06_SupFig5_ILC_NST_Alterations.R, 23_Fig6_Patient_Signatures_Resemblance_Scores.R
# ------------------------------------------------------------------------------

#' Build frequency table from Genomic Alteration Matrix (GAM)
#' 
#' Calculates alteration frequencies between ILC and NST groups from a GAM,
#' performs Fisher's exact tests, and returns a comprehensive results table.
#' 
#' @param tumor_gam Matrix of genomic alterations (genes x samples)
#' @param tumor_labels Named vector/factor with sample labels (ILC/NST)
#' @param genes Vector of gene names to analyze (optional, NULL = all genes)
#' @param ilc_label Label for ILC group (default: "ILC")
#' @param nst_label Label for NST group (default: "NST")
#' @param laplace Laplace smoothing parameter (default: 0.5)
#' @param fisher_alternative Alternative hypothesis for Fisher's test (default: "two.sided")
#' @param min_count Minimum count threshold for inclusion (default: 3)
#' @param allowed_classes Vector of allowed alteration classes (default: common classes)
#' @param debug Logical, whether to print debug messages (default: FALSE)
#' @return Data frame with frequency statistics, odds ratios, and p-values
#' 
#' @examples
#'   freq_tbl <- build_freq_tbl_from_GAM(
#'     tumor_gam = TCGA_BRCA_GAM_Simple,
#'     tumor_labels = setNames(tcga_subset$`Final Pathology`, tcga_subset$Case.ID),
#'     genes = c("CDH1", "TP53", "PTEN"),
#'     min_count = 2
#'   )
build_freq_tbl_from_GAM <- function(
    tumor_gam,
    tumor_labels,
    genes = NULL,
    ilc_label   = "ILC",
    nst_label   = "NST",
    laplace     = 0.5,
    fisher_alternative = c("two.sided", "greater", "less"),
    min_count   = 3,
    allowed_classes = c("MUT", "DEL", "AMP", "MUT;GAIN", "MUT;LOH", "MUT;AMP", "MUT;DEL", "LOH", "GAIN"),
    debug       = FALSE
) {
  fisher_alternative <- match.arg(fisher_alternative)
  tumor_gam <- as.matrix(tumor_gam)
  if (!is.numeric(tumor_gam)) tumor_gam[] <- trimws(tumor_gam)
  for (cl in setdiff(unique(as.character(tumor_gam)), allowed_classes)) tumor_gam[tumor_gam == cl] <- ""
  stopifnot(all(colnames(tumor_gam) %in% names(tumor_labels)))
  lbl <- factor(tumor_labels[colnames(tumor_gam)])
  ilc_idx <- which(lbl == ilc_label)
  nst_idx <- which(lbl == nst_label)
  n_ilc <- length(ilc_idx)
  n_nst <- length(nst_idx)
  if (n_ilc == 0 || n_nst == 0) stop("No samples for one or both groups in tumor_labels.")
  clean_vec <- function(x) {
    if (is.null(x)) return(character(0))
    x <- as.character(x)
    x <- trimws(x)
    x[!is.na(x) & nzchar(x)]
  }
  parse_specs_long <- function(spec_vec, tag) {
    spec_vec <- clean_vec(spec_vec)
    if (!length(spec_vec)) return(data.frame(gene = character(), class = character(), source = character()))
    out <- do.call(rbind, lapply(spec_vec, function(s) {
      p <- strsplit(s, "_", fixed = TRUE)[[1]]
      g <- p[1]
      if (length(p) == 1) data.frame(gene = g, class = "ANY", source = tag, stringsAsFactors = FALSE)
      else {
        cls <- unique(p[-1])
        data.frame(gene = g, class = cls, source = tag, stringsAsFactors = FALSE)
      }
    }))
    out
  }
  ilc_long <- parse_specs_long(genes, "ILC")
  nst_long <- parse_specs_long(genes, "NST")
  specs_long <- rbind(ilc_long, nst_long)
  if (debug) message(sprintf("Specs: %d rows, %d genes", nrow(specs_long), length(unique(specs_long$gene))))
  if (!nrow(specs_long)) stop("No specs provided after cleaning.")
  specs_long <- specs_long[specs_long$gene %in% rownames(tumor_gam), , drop = FALSE]
  if (!nrow(specs_long)) stop("None of the specified genes are present in tumor_gam.")
  is_any <- specs_long$class == "ANY"
  bad_cls <- !is_any & !(specs_long$class %in% allowed_classes)
  specs_long <- specs_long[is_any | (!is_any & specs_long$class %in% allowed_classes), , drop = FALSE]
  if (!nrow(specs_long)) stop("After validation, no (gene,class) rows remain.")
  src_tab <- reshape(transform(specs_long, flag = 1), timevar = "source", idvar = c("gene", "class"), direction = "wide")
  names(src_tab) <- sub("^flag\\.", "", names(src_tab))
  if (!"ILC" %in% names(src_tab)) src_tab$ILC <- NA
  if (!"NST" %in% names(src_tab)) src_tab$NST <- NA
  src_tab$in_ILC_spec <- !is.na(src_tab$ILC)
  src_tab$in_NST_spec <- !is.na(src_tab$NST)
  src_tab <- src_tab[, c("gene", "class", "in_ILC_spec", "in_NST_spec")]
  pairs <- unique(src_tab[, c("gene", "class")])
  count_exact <- function(vec, cl) sum(vec == cl, na.rm = TRUE)
  count_any   <- function(vec) sum(vec %in% allowed_classes, na.rm = TRUE)
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    g <- pairs$gene[i]; cl <- pairs$class[i]
    if (cl == "ANY") {
      k_ilc <- count_any(tumor_gam[g, ilc_idx]); k_nst <- count_any(tumor_gam[g, nst_idx])
    } else {
      k_ilc <- count_exact(tumor_gam[g, ilc_idx], cl); k_nst <- count_exact(tumor_gam[g, nst_idx], cl)
    }
    f_ilc <- if (n_ilc > 0) k_ilc / n_ilc else NA_real_
    f_nst <- if (n_nst > 0) k_nst / n_nst else NA_real_
    data.frame(gene = g, class = cl, k_ILC = k_ilc, n_ILC = n_ilc, f_ILC_raw = f_ilc,
               k_NST = k_nst, n_NST = n_nst, f_NST_raw = f_nst, stringsAsFactors = FALSE)
  })
  counts <- do.call(rbind, rows)
  keep <- pmax(counts$k_ILC, counts$k_NST) >= min_count
  counts <- counts[keep, , drop = FALSE]
  if (!nrow(counts)) { warning("No (gene,class) rows met min_count."); return(data.frame()) }
  smooth_freq <- function(k, n, a = laplace) (k + a) / (n + 2 * a)
  res <- do.call(rbind, lapply(seq_len(nrow(counts)), function(i) {
    a <- counts$k_ILC[i]; b <- counts$n_ILC[i] - a; c <- counts$k_NST[i]; d <- counts$n_NST[i] - c
    ft_p <- ft_or <- ci_low <- ci_high <- NA_real_
    if ((a + b > 0) && (c + d > 0)) {
      ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = fisher_alternative), error = function(e) NULL)
      if (!is.null(ft)) { ft_p <- unname(ft$p.value); ft_or <- if (!is.null(ft$estimate)) unname(ft$estimate) else NA_real_; if (!is.null(ft$conf.int)) { ci_low <- ft$conf.int[1]; ci_high <- ft$conf.int[2] } }
    }
    if (!is.finite(ft_or)) ft_or <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    p_ILC <- smooth_freq(counts$k_ILC[i], counts$n_ILC[i]); p_NST <- smooth_freq(counts$k_NST[i], counts$n_NST[i])
    data.frame(gene = counts$gene[i], class = counts$class[i], k_ILC = counts$k_ILC[i], n_ILC = counts$n_ILC[i],
               f_ILC_raw = counts$f_ILC_raw[i], p_ILC = p_ILC, k_NST = counts$k_NST[i], n_NST = counts$n_NST[i],
               f_NST_raw = counts$f_NST_raw[i], p_NST = p_NST, delta_freq = p_ILC - p_NST,
               OR = ft_or, OR_low = ci_low, OR_high = ci_high, p_value = ft_p,
               log2OR = if (is.finite(ft_or) && ft_or > 0) log2(ft_or) else NA_real_, stringsAsFactors = FALSE)
  }))
  out <- merge(res, src_tab, by = c("gene", "class"), all.x = TRUE)
  out$in_ILC_spec[is.na(out$in_ILC_spec)] <- FALSE
  out$in_NST_spec[is.na(out$in_NST_spec)] <- FALSE
  out$FDR <- p.adjust(out$p_value, method = "BH")
  out <- out[order(out$FDR, -abs(out$log2OR), -abs(out$delta_freq)), ]
  rownames(out) <- NULL
  out
}

# ==============================================================================
# SECTION 12: SANKEY DIAGRAM VISUALIZATION FUNCTIONS
# ==============================================================================
# Moved from Visualization/plot_sankey.R. Scripts sourcing Helper_Functions.R
# get these functions; no need to source plot_sankey.R.
# ------------------------------------------------------------------------------

if (requireNamespace("ggsankey", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggsankey))
}

plot_subtype_sankey <- function(data, source_col, target_col, title = NULL, colors = NULL, levels = NULL, base_size = 15, width = 7, height = 5) {
  if (!requireNamespace("ggsankey", quietly = TRUE)) stop("Package 'ggsankey' is required.")
  if (is.null(colors)) {
    if (exists("VIZ_PARAMS", envir = .GlobalEnv)) colors <- c(VIZ_PARAMS$subtypes, c("LumA" = "#E1CCC6", "LumB" = "#D4A5A5", "Her2" = "#D12C73", "Normal" = "#CCCCCC"))
    else colors <- c("Lum" = "#E1CCC6", "Lum/HER2" = "#EDB185", "HER2" = "#D12C73", "Basal" = "#57429E", "LumA" = "#E1CCC6", "LumB" = "#D4A5A5", "Her2" = "#D12C73", "Normal" = "#CCCCCC")
  }
  if (is.null(levels)) {
    if (exists("SAMPLE_GROUPS", envir = .GlobalEnv)) levels <- SAMPLE_GROUPS$cl_levels
    else levels <- unique(c(data[[source_col]], data[[target_col]]))
  }
  df <- data %>% filter(!is.na(.data[[source_col]]) & !is.na(.data[[target_col]])) %>%
    select(all_of(c(source_col, target_col))) %>% make_long(.data[[source_col]], .data[[target_col]]) %>%
    drop_na(node) %>% mutate(node = factor(node, levels = levels), next_node = factor(next_node, levels = levels))
  p <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) +
    geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) +
    scale_fill_manual("Molecular\nSubtypes", values = colors) + theme_sankey(base_size = base_size) + xlab("")
  if (!is.null(title)) p <- p + ggtitle(title)
  return(p)
}

plot_multiplatform_sankey <- function(data, columns, colors = NULL, levels = NULL, base_size = 15) {
  if (!requireNamespace("ggsankey", quietly = TRUE)) stop("Package 'ggsankey' is required.")
  if (is.null(colors) && exists("VIZ_PARAMS", envir = .GlobalEnv)) colors <- VIZ_PARAMS$subtypes
  if (is.null(levels) && exists("SAMPLE_GROUPS", envir = .GlobalEnv)) levels <- SAMPLE_GROUPS$cl_levels
  df <- data %>% select(all_of(columns)) %>% make_long(!!!rlang::syms(columns)) %>% drop_na(node) %>%
    mutate(node = factor(node, levels = levels), next_node = factor(next_node, levels = levels))
  ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) +
    geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) +
    scale_fill_manual("Molecular\nSubtypes", values = colors) + theme_sankey(base_size = base_size) + xlab("")
}

generate_all_subtype_sankeyplots <- function(annotations, save_plots = TRUE, output_dir = NULL) {
  icle_data <- annotations %>% filter(Study == "ICLE")
  if (is.null(output_dir) && exists("DIRS")) output_dir <- DIRS$results_sub$molecular_subtyping
  plots <- list(
    mrna = plot_subtype_sankey(icle_data, "mRNA Subtypes", "TopCall", title = "mRNA-based Classification (PAM50)"),
    rppa = plot_subtype_sankey(icle_data, "RPPA Subtypes", "TopCall", title = "RPPA-based Classification"),
    dnam = plot_subtype_sankey(icle_data, "DNAm Subtypes", "TopCall", title = "DNA Methylation-based Classification"),
    integrated = plot_multiplatform_sankey(icle_data, c("Sample", "mRNA Subtypes", "RPPA Subtypes", "DNAm Subtypes"))
  )
  if (save_plots && !is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    ggsave(file.path(output_dir, "SupFig2_RNA_PAM50.pdf"), plot = plots$mrna, width = 7, height = 5)
    ggsave(file.path(output_dir, "SupFig3_RPPA_Subtypes.pdf"), plot = plots$rppa, width = 7, height = 5)
    ggsave(file.path(output_dir, "SupFig4_DNAm_Subtypes.pdf"), plot = plots$dnam, width = 7, height = 5)
    ggsave(file.path(output_dir, "Fig1B_MolecularSubtypes.pdf"), plot = plots$integrated, width = 7, height = 5)
    message("✓ All subtype sankey plots saved to: ", output_dir)
  }
  return(plots)
}

.sankey_cl_levels <- rev(c("MDAMB134VI", "MDAMB330", "SUM44PE", "CAMA1", "BCK4", "HCC2185", "IPH926", "UACC3133", "HCC2218", "MDAMB453", "ZR7530", "WCRC25", "OCUBM", "600MPE", "SKBR3", "MDAMB468", "HCC1187", "Lum", "LumA", "LumB", "LumA;LumB", "LumA;Normal", "LumA;Her2", "Lum/HER2", "Her2", "HER2", "Normal", "Basal"))

create_mrna_pam50_sankey <- function(CL_Annots, annot_cols, output_file = NULL) {
  df <- CL_Annots %>% filter(Study == "ICLE") %>% select(`mRNA Subtypes`, TopCall) %>% make_long(`mRNA Subtypes`, TopCall) %>% drop_na(node) %>% mutate(node = factor(node, levels = .sankey_cl_levels), next_node = factor(next_node, levels = .sankey_cl_levels))
  p <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) + geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) + scale_fill_manual("Molecular\nSubtypes", values = c(annot_cols$PAM50, annot_cols$Subtypes)) + theme_sankey(base_size = 15) + xlab("") + labs(title = "mRNA Subtypes vs PAM50")
  if (!is.null(output_file)) { ggsave(output_file, plot = p, width = 7, height = 5); message("  ✓ Saved: ", output_file) }
  return(p)
}

create_rppa_pam50_sankey <- function(CL_Annots, annot_cols, output_file = NULL) {
  df <- CL_Annots %>% filter(Study == "ICLE") %>% select(`RPPA Subtypes`, TopCall) %>% make_long(`RPPA Subtypes`, TopCall) %>% drop_na(node) %>% mutate(node = factor(node, levels = .sankey_cl_levels), next_node = factor(next_node, levels = .sankey_cl_levels))
  p <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) + geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) + scale_fill_manual("Molecular\nSubtypes", values = c(annot_cols$PAM50, annot_cols$Subtypes)) + theme_sankey(base_size = 15) + xlab("") + labs(title = "RPPA Subtypes vs PAM50")
  if (!is.null(output_file)) { ggsave(output_file, plot = p, width = 7, height = 5); message("  ✓ Saved: ", output_file) }
  return(p)
}

create_dnam_pam50_sankey <- function(CL_Annots, annot_cols, output_file = NULL) {
  df <- CL_Annots %>% filter(Study == "ICLE") %>% select(`DNAm Subtypes`, TopCall) %>% make_long(`DNAm Subtypes`, TopCall) %>% drop_na(node) %>% mutate(node = factor(node, levels = .sankey_cl_levels), next_node = factor(next_node, levels = .sankey_cl_levels))
  p <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) + geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) + scale_fill_manual("Molecular\nSubtypes", values = c(annot_cols$PAM50, annot_cols$Subtypes)) + theme_sankey(base_size = 15) + xlab("") + labs(title = "DNAm Subtypes vs PAM50")
  if (!is.null(output_file)) { ggsave(output_file, plot = p, width = 7, height = 5); message("  ✓ Saved: ", output_file) }
  return(p)
}

create_multiomic_sankey <- function(CL_Annots, annot_cols, output_file = NULL) {
  df <- CL_Annots %>% filter(Study == "ICLE") %>% select(Sample, `mRNA Subtypes`, `RPPA Subtypes`, `DNAm Subtypes`) %>% make_long(Sample, `mRNA Subtypes`, `RPPA Subtypes`, `DNAm Subtypes`) %>% drop_na(node) %>% mutate(node = factor(node, levels = .sankey_cl_levels), next_node = factor(next_node, levels = .sankey_cl_levels))
  p <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) + geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) + scale_fill_manual("Molecular\nSubtypes", values = annot_cols$Subtypes) + theme_sankey(base_size = 15) + xlab("") + labs(title = "Multi-Omic Molecular Subtypes")
  if (!is.null(output_file)) { ggsave(output_file, plot = p, width = 7, height = 5); message("  ✓ Saved: ", output_file) }
  return(p)
}

generate_all_sankey_plots <- function(CL_Annots, annot_cols, output_dir = NULL) {
  if (is.null(output_dir) && exists("DIRS")) output_dir <- DIRS$results_sub$molecular_subtyping
  list(mrna_pam50 = create_mrna_pam50_sankey(CL_Annots, annot_cols, file.path(output_dir, "SupFig2C_RNA_PAM50.pdf")), rppa_pam50 = create_rppa_pam50_sankey(CL_Annots, annot_cols, file.path(output_dir, "SupFig3C_RPPA_PAM50.pdf")), dnam_pam50 = create_dnam_pam50_sankey(CL_Annots, annot_cols, file.path(output_dir, "SupFig4C_DNAm_PAM50.pdf")), multiomic = create_multiomic_sankey(CL_Annots, annot_cols, file.path(output_dir, "Fig1B_MolecularSubtypes.pdf")))
}

plot_subtype_sankey_annot <- function(annots, from_col, to_col, color_palette, filter_study = "ICLE", title = "") {
  if (!requireNamespace("ggsankey", quietly = TRUE)) stop("Package 'ggsankey' is required.")
  df <- annots %>% filter(Study == filter_study) %>% select(all_of(c(from_col, to_col))) %>% make_long(all_of(from_col), all_of(to_col)) %>% drop_na(node) %>% mutate(node = factor(node, levels = .sankey_cl_levels), next_node = factor(next_node, levels = .sankey_cl_levels))
  ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) + geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) + scale_fill_manual("Molecular\nSubtypes", values = color_palette) + theme_sankey(base_size = 15) + xlab("") + ggtitle(title)
}

plot_multicolumn_sankey <- function(annots, columns, color_palette, filter_study = "ICLE", title = "") {
  if (!requireNamespace("ggsankey", quietly = TRUE)) stop("Package 'ggsankey' is required.")
  df <- annots %>% filter(Study == filter_study) %>% select(all_of(columns)) %>% make_long(all_of(columns)) %>% drop_na(node) %>% mutate(node = factor(node, levels = .sankey_cl_levels), next_node = factor(next_node, levels = .sankey_cl_levels))
  ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1, smooth = 3, space = 0.6, width = 0.1) + geom_sankey_label(size = 4, color = 1, fill = "white", space = 0.6) + scale_fill_manual("Molecular\nSubtypes", values = color_palette) + theme_sankey(base_size = 15) + xlab("") + ggtitle(title)
}
