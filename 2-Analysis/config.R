# ==============================================================================
# ICLE Project Configuration
# ==============================================================================
# This file centralizes project paths and analysis parameters for reproducibility
# Author: Osama Shiraz Shah
# Last Updated: 2025
# ==============================================================================

# ------------------------------------------------------------------------------
# Project Root Detection
# ------------------------------------------------------------------------------
get_project_root <- function() {
  # Try multiple strategies to find project root
  wd <- getwd()
  
  # Strategy 1: If running from 2-Analysis directory
  if (basename(wd) == "2-Analysis" && file.exists(file.path(wd, "..", "1-Datasets"))) {
    return(normalizePath(file.path(wd, "..")))
  }
  
  # Strategy 2: Check if we're in project root
  if (file.exists(file.path(wd, "1-Datasets")) && file.exists(file.path(wd, "2-Analysis"))) {
    return(normalizePath(wd))
  }
  
  # Strategy 3: Go up from Helper_Scripts or subdirectories
  if (grepl("2-Analysis", wd)) {
    # Navigate up to find 2-Analysis directory
    temp_wd <- wd
    while (!basename(temp_wd) == "2-Analysis" && temp_wd != dirname(temp_wd)) {
      temp_wd <- dirname(temp_wd)
    }
    if (basename(temp_wd) == "2-Analysis") {
      parent <- dirname(temp_wd)
      if (file.exists(file.path(parent, "1-Datasets"))) {
        return(normalizePath(parent))
      }
    }
  }
  
  # Strategy 4: Go up from current directory to find project root
  temp_wd <- wd
  for (i in 1:5) {  # Try going up 5 levels max
    if (file.exists(file.path(temp_wd, "1-Datasets")) && 
        file.exists(file.path(temp_wd, "2-Analysis"))) {
      return(normalizePath(temp_wd))
    }
    temp_wd <- dirname(temp_wd)
    if (temp_wd == dirname(temp_wd)) break  # Reached root
  }
  
  # Default: return current directory (will trigger warning)
  return(normalizePath(wd))
}

# Set project root
PROJECT_ROOT <- get_project_root()

# Verify root
if (!dir.exists(file.path(PROJECT_ROOT, "1-Datasets")) || 
    !dir.exists(file.path(PROJECT_ROOT, "2-Analysis"))) {
  warning("Project root detection may have failed. Please verify PROJECT_ROOT: ", PROJECT_ROOT)
}

# ------------------------------------------------------------------------------
# Directory Paths
# ------------------------------------------------------------------------------
DIRS <- list(
  root = PROJECT_ROOT,
  data = file.path(PROJECT_ROOT, "1-Datasets"),
  analysis = file.path(PROJECT_ROOT, "2-Analysis"),
  results = file.path(PROJECT_ROOT, "3-Results"),
  figures = file.path(PROJECT_ROOT, "4-Figures"),
  sharing = file.path(PROJECT_ROOT, "5-DataSharing"),
  
  # ICLE subdirectories (use base paths; add subdirs via file.path e.g. bionano + "2_Structural_Variations")
  icle = list(
    base = file.path(PROJECT_ROOT, "1-Datasets", "ICLE"),
    rppa = file.path(PROJECT_ROOT, "1-Datasets", "ICLE", "RPPA"),
    rnaseq = file.path(PROJECT_ROOT, "1-Datasets", "ICLE", "RNAseq"),
    bionano = file.path(PROJECT_ROOT, "1-Datasets", "ICLE", "Bionano"),
    cytosnp = file.path(PROJECT_ROOT, "1-Datasets", "ICLE", "CytoSNP"),
    wes = file.path(PROJECT_ROOT, "1-Datasets", "ICLE", "WES"),
    dnam = file.path(PROJECT_ROOT, "1-Datasets", "ICLE", "DNAm")
  ),
  
  # External data subdirectories
  external = list(
    base = file.path(PROJECT_ROOT, "1-Datasets", "External"),
    tcga = file.path(PROJECT_ROOT, "1-Datasets", "External", "TCGA"),
    msk = file.path(PROJECT_ROOT, "1-Datasets", "External", "MSK"),
    ccle = file.path(PROJECT_ROOT, "1-Datasets", "External", "CCLE"),
    sanger = file.path(PROJECT_ROOT, "1-Datasets", "External", "Sanger"),
    fmi = file.path(PROJECT_ROOT, "1-Datasets", "External", "FMI"),
    marcotte = file.path(PROJECT_ROOT, "1-Datasets", "External", "Marcotte"),
    misc = file.path(PROJECT_ROOT, "1-Datasets", "External", "Misc")
  ),
  
  # Scripts subdirectories
  scripts = list(
    helpers = file.path(PROJECT_ROOT, "2-Analysis", "Helper_Scripts"),
    data_prep = file.path(PROJECT_ROOT, "2-Analysis", "Helper_Scripts"),
    analysis = file.path(PROJECT_ROOT, "2-Analysis", "Helper_Scripts"),
    visualization = file.path(PROJECT_ROOT, "2-Analysis", "Helper_Scripts")
  ),
  
  # Results subdirectories
  results_sub = list(
    molecular_subtyping = file.path(PROJECT_ROOT, "3-Results", "Molecular_Subtyping"),
    dna_methylation = file.path(PROJECT_ROOT, "3-Results", "DNA_Methylation"),
    molecular_resemblance = file.path(PROJECT_ROOT, "3-Results", "Molecular_Resemblance"),
    ogm = file.path(PROJECT_ROOT, "3-Results", "Optical_Genome_Mapping"),
    gam = file.path(PROJECT_ROOT, "3-Results", "Genomic_Alteration_Matrix"),
    cdh1 = file.path(PROJECT_ROOT, "3-Results", "CDH1_Alteration_Landscape"),
    dependencies = file.path(PROJECT_ROOT, "3-Results", "Gene_Dependencies")
  )
)

# ------------------------------------------------------------------------------
# File Paths (Commonly Used Files)
# ------------------------------------------------------------------------------
FILES <- list(
  # Annotation files
  cl_annots_simple = file.path(DIRS$icle$base, "BRCA_CL_Annots_Simple.tsv"),
  cl_annotations = file.path(DIRS$icle$base, "BRCA_CL_Annots.tsv"),
  
  # ICLE data files
  rppa_data = file.path(DIRS$icle$rppa, "BRCA_CL_RPPA_COMBAT.Rdata"),
  rnaseq_data = file.path(DIRS$icle$rnaseq, "3_Counts", "BRCA_CL_EXP_Log2CPM.Rdata"),
  icle_genotype = file.path(DIRS$icle$cytosnp, "4_Genotyping", "ICLE_Genotype_Mat.Rdata"),
  sv_data = file.path(DIRS$icle$bionano, "2_Structural_Variations", "ICLE_SV_Filtered.csv"),
  cnv_seg = file.path(DIRS$icle$cytosnp, "3_Segmentation", "BRCA_CL_LogRR_DNACopy.seg"),
  
  # External data files
  marcotte_genotype = file.path(DIRS$external$marcotte, "SNP", "Marcotte_Genotype_Mat.Rdata"),
  pathway_dir = file.path(DIRS$external$tcga, "RPPA", "Akhbani Pathways"),
  set_genes = file.path(DIRS$external$misc, "SET", "SET_GENES.xlsx"),
  oncovar = file.path(DIRS$external$misc, "OncoVar_2021", "TCGA.BRCA.all.genes.OncoVar.tsv"),
  refgene = file.path(DIRS$external$misc, "refGene.hg19.sorted_dedup_geneSumarized.bed"),
  hg19_chrom_sizes = file.path(DIRS$external$misc, "hg19.chrom.sizes"),
  tumor_fusions = file.path(DIRS$external$tcga, "Fusions", "Tumor_Fusions.xlsx"),
  biomart_genemap = file.path(DIRS$external$misc, "bioMart_geneMap_w_length.Rdata"),
  hm450k_probeset = file.path(DIRS$external$misc, "HM450K_ProbeSet_Unique_Gene_Promoter_Mapping.Rdata"),

  # External Misc (pathways, DDR, fusion results)
  ddr_genes = file.path(DIRS$external$misc, "DDR_genes.txt"),
  maftools_pathways = file.path(DIRS$external$misc, "maftools_oncogenic_sig_pathways.tsv"),
  fusion_df = file.path(DIRS$icle$bionano, "3_Fusions", "fusion_df.tsv"),
  refgene_unique = file.path(DIRS$external$misc, "refGene.hg19.sorted.unique.bed"),
  fmi_breast_diag = file.path(DIRS$external$fmi, "Breast_Diagnostic_REs.xlsx"),

  # CDH1 exonic deletions (Fig2C)
  cdh1_gtf = file.path(DIRS$icle$wes, "6_CDH1_Reads", "CDH1_gtf.bed"),
  cdh1_pileups_dir = file.path(DIRS$icle$wes, "6_CDH1_Reads", "pileups"),

  # GAM pipeline inputs/outputs
  tcga_maf = file.path(DIRS$external$tcga, "SNV", "TCGA_BRCA_MAF.tsv"),
  tcga_maf_rdata = file.path(DIRS$external$tcga, "SNV", "TCGA_BRCA_MAF.Rdata"),
  tcga_gistic = file.path(DIRS$external$tcga, "CN", "TCGA_BRCA_GISTIC.txt"),
  msk_maf = file.path(DIRS$external$msk, "SNV", "MSK_BRCA_MAF.tsv"),
  msk_maf_rdata = file.path(DIRS$external$msk, "SNV", "MSK_BRCA_MAF.Rdata"),
  msk_gistic = file.path(DIRS$external$msk, "CN", "MSK_BRCA_GISTIC.txt"),
  cl_maf = file.path(DIRS$icle$wes, "4_MAF", "BRCA_CL_MAF.tsv"),
  cl_maf_rdata = file.path(DIRS$icle$wes, "4_MAF", "BRCA_CL_MAF.Rdata"),
  cl_gistic = file.path(DIRS$icle$cytosnp, "5_GISTIC", "342744", "all_thresholded.by_genes.txt"),
  tcga_gam = file.path(DIRS$results_sub$gam, "TCGA_BRCA_GAM.Rdata"),
  msk_gam = file.path(DIRS$results_sub$gam, "MSK_BRCA_MET_GAM.Rdata"),
  cl_gam = file.path(DIRS$results_sub$gam, "BRCA_CL_GAM.Rdata"),

  # Genomic metrics (Fig1D, etc.)
  tmb_tsv = file.path(DIRS$icle$wes, "5_TMB", "TMB.tsv"),
  fga_tsv = file.path(DIRS$icle$cytosnp, "6_CIN_Analysis", "CINMetric_FGA.tsv"),
  dmi_tsv = file.path(DIRS$icle$dnam, "3_Integration", "DMI_Summary.tsv"),

  # Fig4 / Fig6 (DNAm, DE, resemblance)
  rnaseq_cts = file.path(DIRS$icle$rnaseq, "3_Counts", "BRCA_CL_EXP_CTS.Rdata"),
  icle_dnam = file.path(DIRS$icle$dnam, "3_Integration", "BRCA_CL_DNAm.Rdata"),
  tcga_deseq = file.path(DIRS$external$tcga, "RNA", "TCGA_BRCA_DESeq_LumA_ILC_vs_NST.Rdata"),
  tcga_limma = file.path(DIRS$external$tcga, "DNAm", "TCGA_BRCA_Limma_LumA_ILC_vs_NST.Rdata"),
  cl_deseq = file.path(DIRS$icle$rnaseq, "BRCA_CL_DEseq2_NonBasal_ILC_vs_NST.Rdata"),
  cl_limma = file.path(DIRS$icle$dnam, "BRCA_CL_Limma_NonBasal_ILC_vs_NST.Rdata"),
  tcga_dnam_regulated = file.path(DIRS$results_sub$dna_methylation, "TCGA_DNAm_Regulated_Genes.csv"),

  # Fig5 (RNAi, dependencies)
  rnai_d2 = file.path(DIRS$external$marcotte, "RNAi", "D2_combined_gene_dep_scores.csv"),
  rnai_breast_screens = file.path(DIRS$external$marcotte, "RNAi", "breast_screens_with_weights.eset"),
  rnai_hairpin = file.path(DIRS$external$marcotte, "RNAi", "hairpin_annotations.txt"),
  rnai_simem = file.path(DIRS$results_sub$dependencies, "simem", "2025_ILC_vs_NST_allgenes.tsv"),
  rnai_simem_rdata = file.path(DIRS$results_sub$dependencies, "simem", "2025_ILC_vs_NST_allgenes.Rdata"),
  rnai_annots = file.path(DIRS$external$marcotte, "RNAi", "RNAi_CL_Annots.tsv"),
  dgidb = file.path(DIRS$external$misc, "DGIdb", "2024-01-02", "interactions.tsv")
)

# ------------------------------------------------------------------------------
# Analysis Parameters
# ------------------------------------------------------------------------------
ANALYSIS_PARAMS <- list(
  # Random seed for reproducibility
  random_seed = 123,
  
  # Parallel processing
  n_cores = parallel::detectCores() - 1,
  
  # Statistical thresholds
  min_samples_per_group = 3,
  p_value_threshold = 0.05,
  fdr_threshold = 0.1,
  
  # Sample size for random sampling
  random_sample_size = 10000,
  
  # Figure dimensions
  fig_width = 10,
  fig_height = 8,
  dpi = 300
)

# ------------------------------------------------------------------------------
# Visualization Parameters
# ------------------------------------------------------------------------------
VIZ_PARAMS <- list(
  # Font family
  font_family = "Helvetica",
  
  # Base font sizes
  base_fontsize = 12,
  title_fontsize = 14,
  
  # Color palettes (from Helper_Functions.R)
  # subtypes = c("Lum" = "#E1CCC6", "Lum/HER2" = "#EDB185", "HER2" = "#D12C73", "Basal" = "#57429E"),
  # histology = c("ILC" = "#CD5038", "ILC-like" = "#FFC266", "NST" = "#8E9BAE"),
  
  # Heatmap colors
  rna_zscore = circlize::colorRamp2(c(-2, -1, 0, 1, 2), c("#003d30", "#7dc8bf","white", "#ca9446", "#543005")),
  rppa_zscore = circlize::colorRamp2(c(-2, -1, 0, 1, 2), c("#003d30", "#7dc8bf","white", "#ca9446", "#543005")),
  
  # Heatmap legend parameters
  heatmap_legend = list(
    legend_direction = "horizontal",
    legend_width = ggplot2::unit(4, "cm"),
    title_gp = grid::gpar(fontsize = 10, fontface = "bold")
  )
)

# ------------------------------------------------------------------------------
# Sample Group Definitions
# ------------------------------------------------------------------------------
SAMPLE_GROUPS <- list(
  # Cell line levels for plotting
  cl_levels = c(
    "MDAMB134VI", "MDAMB330", "SUM44PE", "CAMA1", "BCK4",
    "HCC2185", "IPH926", "UACC3133", "HCC2218", "MDAMB453",
    "ZR7530", "WCRC25", "OCUBM", "600MPE", "SKBR3",
    "MDAMB468", "HCC1187", "Lum", "LumA", "LumB", "LumA;LumB",
    "LumA;Normal", "LumA;Her2", "Lum/HER2", "Her2", "HER2", 
    "Normal", "Basal"
  ),
  
  # Molecular subtype order
  subtype_order = c("Lum", "Lum/HER2", "HER2", "Basal"),
  
  # Histology order
  histology_order = c("ILC", "ILC-like", "NST")
)

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

#' Create output directory if it doesn't exist
#' @param path Directory path
#' @return NULL (creates directory as side effect)
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created directory: ", path)
  }
}

#' Get file path with automatic directory creation
#' @param ... Path components
#' @param create_dir Whether to create parent directory
#' @return Full file path
get_file_path <- function(..., create_dir = TRUE) {
  path <- file.path(...)
  if (create_dir) {
    ensure_dir(dirname(path))
  }
  return(path)
}

#' Save plot with standardized parameters
#' @param plot ggplot or ComplexHeatmap object
#' @param filename Output filename
#' @param width Plot width (inches)
#' @param height Plot height (inches)
#' @param dpi Resolution
#' @return NULL (saves plot as side effect)
save_plot <- function(plot, filename, width = ANALYSIS_PARAMS$fig_width, 
                     height = ANALYSIS_PARAMS$fig_height, dpi = ANALYSIS_PARAMS$dpi) {
  ensure_dir(dirname(filename))
  
  if (inherits(plot, "Heatmap") || inherits(plot, "HeatmapList")) {
    # ComplexHeatmap object
    pdf(file = filename, width = width, height = height)
    draw(plot, merge_legend = TRUE)
    dev.off()
  } else {
    # ggplot object
    ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi)
  }
  message("Saved plot: ", filename)
}

# ------------------------------------------------------------------------------
# Verify Critical Directories
# ------------------------------------------------------------------------------
verify_directories <- function(create_missing = FALSE) {
  critical_dirs <- c(
    DIRS$data,
    DIRS$analysis,
    DIRS$results,
    DIRS$figures,
    DIRS$icle$base,
    DIRS$icle$rppa,
    DIRS$icle$rnaseq,
    DIRS$icle$bionano,
    DIRS$icle$cytosnp,
    DIRS$icle$wes,
    DIRS$icle$dnam,
    DIRS$external$base,
    DIRS$external$tcga,
    DIRS$external$msk,
    DIRS$scripts$helpers
  )
  
  missing <- critical_dirs[!dir.exists(critical_dirs)]
  if (length(missing) > 0) {
    if (create_missing) {
      lapply(missing, dir.create, recursive = TRUE)
      message("Created missing directories:\n", paste(missing, collapse = "\n"))
    } else {
      warning("The following critical directories are missing:\n",
              paste(missing, collapse = "\n"))
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Verify critical files exist
#' @return Logical indicating if all files exist
verify_files <- function() {
  critical_files <- unlist(FILES)
  missing <- critical_files[!file.exists(critical_files)]
  if (length(missing) > 0) {
    warning("The following critical files are missing:\n",
            paste(missing, collapse = "\n"))
    return(FALSE)
  }
  return(TRUE)
}

# ------------------------------------------------------------------------------
# Initialization
# ------------------------------------------------------------------------------

# Create results subdirectories if they don't exist
lapply(DIRS$results_sub, ensure_dir)

# Set random seed for reproducibility
set.seed(ANALYSIS_PARAMS$random_seed)

# Print configuration status
cat("===============================================\n")
cat("ICLE Project Configuration Loaded\n")
cat("===============================================\n")
cat("Project Root:", DIRS$root, "\n")
cat("Directories verified:", verify_directories(), "\n")
cat("Files verified:", verify_files(), "\n")
cat("Random seed set:", ANALYSIS_PARAMS$random_seed, "\n")
cat("===============================================\n")
