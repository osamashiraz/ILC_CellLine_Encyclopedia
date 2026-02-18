# ==============================================================================
# Script 05: Figure 1F - Alteration Frequency Barplots
# ==============================================================================
# Description: Creates barplots comparing alteration frequencies of key genes
#              across ICLE cell lines, TCGA primary ILC tumors, and MSK 
#              metastatic ILC tumors.
#
# Dependencies:
#   - config.R (paths)
#   - Helper_Functions.R (themes and colors)
#   - 01_load_all_data.R (GAM and annotation data)
#
# Required Data Objects:
#   - CL_Annots: Cell line annotations
#   - BRCA_CL_GAM: Cell line genomic alteration matrix
#   - TCGA_Annots, TCGA_BRCA_GAM: TCGA tumor data
#   - MSK_Annots, MSK_BRCA_GAM: MSK metastatic tumor data
#
# Output:
#   - fig1f_alteration_barplot: Alteration barplot (assigned to .GlobalEnv when script is run)
#   - PDFs and p-value tables in DIRS$results/Molecular_Resemblance
#
# Author: Osama Shiraz Shah
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggplot2)
  library(ggpattern)
  library(gt)
})

# Load configuration and helper functions
if (!exists("DIRS")) source("config.R")
if (!exists("annot_cols")) source(file.path(DIRS$scripts$helpers, "Helper_Functions.R"))

# Load data if not already loaded
if (!exists("BRCA_CL_GAM")) {
  message("Loading data objects...")
  source(file.path(DIRS$scripts$helpers, "Data_Loading", "01_load_all_data.R"))
  load_all_icle_data(load_external = TRUE)
}

message("\n========================================")
message("Figure 1F: Alteration Barplots")
message("========================================\n")

# ==============================================================================
# 1. Define Gene Sets
# ==============================================================================

define_gene_sets <- function() {
  ILC_genes <- c("RUNX1", "FOXA1", "TBX3", "CDH1", "PTEN", "CTNNA1")
  NST_genes <- c("BRIP1", "PPM1D", "GATA3", "MAP3K1", "MAP2K4", "MYC", "TP53")
  
  escat_genes <- c(
    "ERBB2", "TP53", "PIK3CA", "CDH1", "BRCA1", "BRCA2", "ESR1", "AKT1", 
    "PTEN", "ERBB3", "NF1", "MDM2", "ARID1A", "ARID1B", "ATR", "PALB2", 
    "ATM", "IGF1R", "MYC", "MAP2K4", "MAP3K1", "CCND1", "FGFR1", "PGR", "FGFR4"
  )
  
  events <- unique(c(escat_genes, NST_genes, ILC_genes))
  
  message("  ✓ Gene sets defined:")
  message("    - ILC genes: ", length(ILC_genes))
  message("    - NST genes: ", length(NST_genes))
  message("    - ESCAT genes: ", length(escat_genes))
  message("    - Total unique genes: ", length(events), "\n")
  
  return(list(
    ILC_genes = ILC_genes,
    NST_genes = NST_genes,
    escat_genes = escat_genes,
    events = events
  ))
}

# ==============================================================================
# 2. Calculate Alteration Frequencies
# ==============================================================================

calculate_alteration_freq <- function(GAM, sample_ids, events, dataset_name) {
  message("    Processing ", dataset_name, "...")
  
  # Melt GAM and filter
  alt_df <- reshape2::melt(GAM[intersect(events, rownames(GAM)), sample_ids])
  alt_df <- subset(alt_df, value %notin% c("", "LOH", "GAIN"))
  names(alt_df) <- c("Gene", "ID", "Alt")
  
  # Count alterations per gene
  alt_df$Gene <- as.character(alt_df$Gene)
  alt_df <- alt_df %>%
    group_by(Gene) %>%
    mutate(Count = n()) %>%
    as.data.frame()
  
  alt_df <- alt_df[!duplicated(alt_df$Gene), ]
  
  # Calculate frequency
  alt_df$Freq <- 100 * alt_df$Count / length(sample_ids)
  alt_df$Gene <- as.character(alt_df$Gene)
  rownames(alt_df) <- alt_df$Gene
  alt_df$Type <- dataset_name
  
  message("      - Samples: ", length(sample_ids))
  message("      - Altered genes: ", nrow(alt_df))
  
  return(alt_df)
}

get_alteration_frequencies <- function(gene_sets) {
  events <- gene_sets$events
  
  # ICLE Cell Lines (non-HER2, non-Basal)
  ICLE_cells <- subset(
    CL_Annots,
    `mRNA Subtypes` %notin% c("HER2", "Basal") & grepl("-I", Name)
  )$Name
  
  ICLE_Alt_Df <- calculate_alteration_freq(
    BRCA_CL_GAM, ICLE_cells, events, "Cell Lines"
  )
  
  # TCGA Primary ILC Tumors
  ILC_Pri <- subset(
    TCGA_Annots,
    `Final Pathology` %in% "ILC" & Case.ID %in% colnames(TCGA_BRCA_GAM)
  )$Case.ID
  
  TCGA_Alt_Df <- calculate_alteration_freq(
    TCGA_BRCA_GAM, ILC_Pri, events, "Primary Tumors"
  )
  
  # MSK Metastatic ILC Tumors
  ILC_Met <- subset(
    MSK_Annots,
    ONCOTREE_CODE %in% "ILC" & SAMPLE_ID %in% colnames(MSK_BRCA_GAM)
  )$SAMPLE_ID
  
  MSK_Alt_Df <- calculate_alteration_freq(
    MSK_BRCA_GAM, ILC_Met, events, "Metastatic Tumors"
  )
  
  message("")
  
  return(list(
    ICLE = ICLE_Alt_Df,
    TCGA = TCGA_Alt_Df,
    MSK = MSK_Alt_Df,
    sample_counts = list(
      Cell_Lines = length(ICLE_cells),
      Primary_Tumors = length(ILC_Pri),
      Metastatic_Tumors = length(ILC_Met)
    )
  ))
}

# ==============================================================================
# 3. Combine and Order Data
# ==============================================================================

combine_and_order_data <- function(alt_freq_list, gene_sets) {
  ILC_genes <- gene_sets$ILC_genes
  NST_genes <- gene_sets$NST_genes
  
  # Find common events
  common_events <- unique(c(
    intersect(
      alt_freq_list$ICLE$Gene,
      c(alt_freq_list$MSK$Gene, alt_freq_list$TCGA$Gene)
    ),
    "TBX3", "GATA3", "MAP3K1", "PTEN", "FOXA1"
  ))
  
  # Combine dataframes
  Combined_Alt_Df <- rbind(
    alt_freq_list$ICLE[, c("Gene", "Type", "Count", "Freq")] %>%
      filter(Gene %in% common_events),
    alt_freq_list$MSK[, c("Gene", "Type", "Count", "Freq")] %>%
      filter(Gene %in% common_events),
    alt_freq_list$TCGA[, c("Gene", "Type", "Count", "Freq")] %>%
      filter(Gene %in% common_events)
  )
  
  # Set factor levels
  Combined_Alt_Df$Type <- factor(
    Combined_Alt_Df$Type,
    levels = rev(c("Primary Tumors", "Metastatic Tumors", "Cell Lines"))
  )
  
  # Fill missing combinations with zeros
  Combined_Alt_Df <- Combined_Alt_Df %>%
    tidyr::complete(Gene, Type, fill = list(Count = 0, Freq = 0))
  
  # Order genes by cell line frequency
  gene_order <- Combined_Alt_Df %>%
    filter(Type %in% c("Cell Lines")) %>%
    group_by(Gene) %>%
    summarise(mean_freq_tumor = mean(Freq, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_freq_tumor) %>%
    pull(Gene)
  
  Combined_Alt_Df$Gene <- factor(Combined_Alt_Df$Gene, levels = gene_order)
  
  # Select genes to plot
  genes_to_plot <- unique(c(
    tail(levels(Combined_Alt_Df$Gene), 15),
    ILC_genes,
    NST_genes
  ))
  
  message("  ✓ Combined data prepared:")
  message("    - Common events: ", length(common_events))
  message("    - Genes to plot: ", length(genes_to_plot), "\n")
  
  return(list(
    data = Combined_Alt_Df,
    genes_to_plot = genes_to_plot
  ))
}

# ==============================================================================
# 4. Create Barplot
# ==============================================================================

create_alteration_barplot <- function(combined_data) {
  Combined_Alt_Df <- combined_data$data
  genes <- combined_data$genes_to_plot
  
  # Define visual parameters
  type_cols <- c(
    "Cell Lines"        = "#EF9A9A",
    "Primary Tumors"    = "#B71C1C",
    "Metastatic Tumors" = "#880E4F"
  )
  
  type_patterns <- c(
    "Cell Lines"        = "none",
    "Primary Tumors"    = "stripe",
    "Metastatic Tumors" = "stripe"
  )
  
  type_pattern_angles <- c(
    "Cell Lines"        = 0,
    "Primary Tumors"    = 45,
    "Metastatic Tumors" = -45
  )
  
  # Create plot
  plt_combined <- subset(Combined_Alt_Df, Gene %in% genes) %>%
    ggplot(aes(
      x = Gene, y = Freq,
      fill = Type, pattern = Type,
      pattern_angle = Type, pattern_color = Type
    )) +
    geom_col_pattern(
      stat = "identity",
      position = position_dodge(width = 0.72),
      alpha = 1, width = 0.70,
      color = "black", lwd = 0.2,
      pattern_key_scale_factor = 0.8,
      pattern_frequency = 0.1,
      pattern_size = 0.2,
      pattern_density = 0.005,
      pattern_spacing = 0.02
    ) +
    coord_flip() +
    ylab("Alteration Frequency") +
    xlab("") +
    scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0, 100)) +
    scale_pattern_angle_manual(values = type_pattern_angles) +
    scale_pattern_color_manual(values = c("white", "white", "white")) +
    scale_fill_manual("Type", values = type_cols) +
    scale_pattern_manual("Type", values = type_patterns) +
    myTheme(14, 14, 14, 14) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.y = element_text(size = 14),
      panel.grid.major.x = element_line(
        color = "#E0E0E0",
        linewidth = 0.3,
        linetype = "dotted"
      ),
      panel.grid.minor.x = element_blank()
    )
  
  message("  ✓ Barplot created\n")
  
  return(plt_combined)
}

# ==============================================================================
# 5. Statistical Testing
# ==============================================================================

perform_statistical_tests <- function(combined_data, sample_counts) {
  Combined_Alt_Df <- combined_data$data
  genes <- combined_data$genes_to_plot
  
  # Add sample totals and calculate counts
  Combined_Alt_Df <- subset(Combined_Alt_Df, Gene %in% genes) %>%
    mutate(
      Total = case_when(
        Type == "Cell Lines"        ~ sample_counts$Cell_Lines,
        Type == "Primary Tumors"    ~ sample_counts$Primary_Tumors,
        Type == "Metastatic Tumors" ~ sample_counts$Metastatic_Tumors
      ),
      Alt = round(Freq / 100 * Total)
    )
  
  # Safe wrapper for prop.test
  prop_pval <- function(alt, total) {
    if (length(alt) != 2L || length(total) != 2L ||
        any(is.na(alt)) || any(is.na(total))) {
      return(NA_real_)
    }
    prop.test(alt, total)$p.value
  }
  
  # Cell Lines vs Primary Tumors
  p_CL_PT <- Combined_Alt_Df %>%
    filter(Type %in% c("Cell Lines", "Primary Tumors")) %>%
    group_by(Gene) %>%
    summarise(
      p_CL_vs_PT = prop_pval(Alt, Total),
      .groups = "drop"
    )
  
  # Cell Lines vs Metastatic Tumors
  p_CL_MT <- Combined_Alt_Df %>%
    filter(Type %in% c("Cell Lines", "Metastatic Tumors")) %>%
    group_by(Gene) %>%
    summarise(
      p_CL_vs_MT = prop_pval(Alt, Total),
      .groups = "drop"
    )
  
  # Combine p-values
  pvals_all <- p_CL_PT %>%
    full_join(p_CL_MT, by = "Gene") %>%
    as.data.frame()
  
  rownames(pvals_all) <- pvals_all$Gene
  pvals_all <- pvals_all[, -1]
  
  # Create notation version
  pvals_notation <- pvals_all
  pvals_notation[pvals_all > 0.05] <- "ns"
  pvals_notation[pvals_all <= 0.05] <- "*"
  pvals_notation[pvals_all <= 0.01] <- "**"
  pvals_notation[pvals_all <= 0.001] <- "***"
  pvals_notation[pvals_all <= 0.0001] <- "****"
  
  message("  ✓ Statistical tests completed\n")
  
  return(list(
    pvals_numeric = pvals_all,
    pvals_notation = pvals_notation
  ))
}

# ==============================================================================
# 6. Save Results
# ==============================================================================

save_results <- function(plot, pvals, genes) {
  # Save main plot
  # output_file <- file.path(DIRS$results, "Molecular_Resemblance", "Fig1F_Alteration_Barplot.pdf")
  # ggsave(filename = output_file, plot, width = 6, height = 5)
  # message("  ✓ Barplot saved: ", output_file)
  
  # Save numeric p-values
  fig1F_pval <- gt::gt(data.frame(
    genes = rev(rownames(pvals$pvals_numeric)),
    signif(pvals$pvals_numeric, 1)[rev(rownames(pvals$pvals_numeric)), ]
  ))
  
  pval_file <- file.path(DIRS$results, "Molecular_Resemblance", "Fig1F_Alteration_Barplot_pval.pdf")
  suppressMessages(gt::gtsave(fig1F_pval, filename = pval_file))
  message("  ✓ P-values (numeric) saved: ", pval_file)
  
  # Save notation p-values
  fig1F_pval_nt <- gt::gt(data.frame(
    genes = rev(rownames(pvals$pvals_notation)),
    pvals$pvals_notation[rev(rownames(pvals$pvals_notation)), ]
  ))
  
  pval_nt_file <- file.path(DIRS$results, "Molecular_Resemblance", "Fig1F_Alteration_Barplot_pval_nt.pdf")
  suppressMessages(gt::gtsave(fig1F_pval_nt, filename = pval_nt_file))
  message("  ✓ P-values (notation) saved: ", pval_nt_file)
  
  message("  ✓ Figure 1F complete\n")
  
  return(invisible(NULL))
}

# ==============================================================================
# 7. Main Execution
# ==============================================================================

main <- function() {

  # Define gene sets
  message("  Step 1/5: Defining gene sets...")
  gene_sets <- define_gene_sets()
  # message("  ✓ Gene sets defined: ", length(gene_sets$genes), " genes\n")
  
  # Calculate frequencies
  message("  Step 2/5: Calculating alteration frequencies...")
  alt_freq_list <- get_alteration_frequencies(gene_sets)
  # message("  ✓ Frequencies calculated across all datasets\n")
  
  # Combine and order
  message("  Step 3/5: Combining and ordering data...")
  combined_data <- combine_and_order_data(alt_freq_list, gene_sets)
  # message("  ✓ Data combined: ", nrow(combined_data$plot_data), " gene-dataset combinations\n")
  
  # Create plot
  message("  Step 4/5: Creating alteration barplot...")
  fig1f_alteration_barplot <- create_alteration_barplot(combined_data)
  # message("  ✓ Plot generated\n")
  
  # Statistical tests
  message("  Step 5/5: Performing statistical tests...")
  pvals <- perform_statistical_tests(combined_data, alt_freq_list$sample_counts)
  # message("  ✓ Statistical tests completed\n")
  
  # Save results
  message("  Saving results p-values plots...")
  save_results(fig1f_alteration_barplot, pvals, combined_data$genes_to_plot)
  
  return(fig1f_alteration_barplot)
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Helper_Functions.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/05_Fig1F_Alteration_barplots.R")
#
# load_all_icle_data(load_external = TRUE)
# main()
#
# Outputs: Fig1F_Alteration_Barplot.pdf, pval PDFs in DIRS$results/Molecular_Resemblance

# Run main function and assign to .GlobalEnv for downstream use
fig1f_alteration_barplot <- main()
assign("fig1f_alteration_barplot", fig1f_alteration_barplot, envir = .GlobalEnv)
