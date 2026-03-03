# ==============================================================================
# ICLE Project - Required R Packages
# ==============================================================================
# This script lists all required R packages and provides installation functions.
# Packages are grouped by source (CRAN, Bioconductor, GitHub) and category.
#
# Author: Osama Shiraz Shah
# Last Updated: 2026-03-02
# ==============================================================================

# ------------------------------------------------------------------------------
# Required R Version
# ------------------------------------------------------------------------------
required_r_version <- "4.0.0"
if (getRversion() < required_r_version) {
  stop("This analysis requires R version ", required_r_version, " or higher. ",
       "Current version: ", getRversion())
}

# ------------------------------------------------------------------------------
# CRAN Packages
# ------------------------------------------------------------------------------
cran_packages <- c(
  # Data manipulation
  "dplyr",
  "tidyr",
  "data.table",
  "reshape2",
  "plyr",
  "magrittr",
  "vroom",
  "tidytext",

  # Visualization - ggplot2 ecosystem
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "ggthemes",
  "ggbeeswarm",
  "ggpattern",
  "ggsci",
  "ggsankey",
  "cowplot",
  "patchwork",
  "gridExtra",
  "viridis",
  "scales",

  # Heatmaps and circular plots
  "circlize",

  # Interactive/HTML visualisation
  "interacCircos",
  "htmlwidgets",
  "webshot2",

  # Tables
  "gt",

  # Statistical utilities
  "matrixStats",
  "EnvStats",

  # File I/O
  "readxl",

  # Parallel processing
  "parallel",
  "doParallel",

  # String manipulation
  "stringr",
  "stringi",
  "forcats",

  # Other utilities
  "devtools",
  "msigdbr"
)

# ------------------------------------------------------------------------------
# Bioconductor Packages
# ------------------------------------------------------------------------------
bioc_packages <- c(
  # Genomics / copy number
  "GenomicRanges",
  "rtracklayer",
  "DNAcopy",
  "CNTools",

  # Mutation / TMB
  "maftools",

  # DNA methylation
  "sesame",

  # Expression analysis
  "limma",
  "edgeR",
  "DESeq2",
  "sva",

  # Subtyping and clustering
  "ConsensusClusterPlus",
  "CancerSubtypes",

  # Pathway / enrichment
  "clusterProfiler",
  "hypeR",
  "org.Hs.eg.db",

  # Infrastructure
  "SummarizedExperiment",
  "Biobase",
  "hiAnnotator",

  # Visualization
  "ComplexHeatmap"
)

# ------------------------------------------------------------------------------
# GitHub Packages
# ------------------------------------------------------------------------------
# Install with: devtools::install_github("user/repo")
github_packages <- list(
  "CGATOxford/ShatterSeek",   # Chromothripsis detection
  "bhklab/BreastSubtypeR"     # Intrinsic breast cancer subtyping
)

# ------------------------------------------------------------------------------
# Package Versions (for reproducibility)
# ------------------------------------------------------------------------------
# Critical packages with specific version requirements
package_versions <- list(
  "ComplexHeatmap" = "2.10.0",
  "ggplot2"        = "3.4.0",
  "dplyr"          = "1.1.0"
)

# ------------------------------------------------------------------------------
# Installation Functions
# ------------------------------------------------------------------------------

#' Install CRAN packages
install_cran_packages <- function(packages = cran_packages, update = FALSE) {
  installed <- rownames(installed.packages())
  missing <- setdiff(packages, installed)

  if (length(missing) > 0) {
    message("Installing ", length(missing), " CRAN packages...")
    install.packages(missing, dependencies = TRUE)
  } else {
    message("All CRAN packages already installed.")
  }

  if (update) {
    message("Updating CRAN packages...")
    update.packages(oldPkgs = packages, ask = FALSE)
  }
}

#' Install Bioconductor packages
install_bioc_packages <- function(packages = bioc_packages, update = FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  installed <- rownames(installed.packages())
  missing <- setdiff(packages, installed)

  if (length(missing) > 0) {
    message("Installing ", length(missing), " Bioconductor packages...")
    BiocManager::install(missing, update = FALSE, ask = FALSE)
  } else {
    message("All Bioconductor packages already installed.")
  }

  if (update) {
    message("Updating Bioconductor packages...")
    BiocManager::install(packages, update = TRUE, ask = FALSE)
  }
}

#' Install GitHub packages
install_github_packages <- function(packages = github_packages) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }

  for (pkg in packages) {
    pkg_name <- basename(pkg)
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      message("Installing ", pkg, " from GitHub...")
      devtools::install_github(pkg)
    }
  }
}

#' Install all required packages
install_all_packages <- function(update = FALSE) {
  message("==============================================")
  message("Installing ICLE Project Dependencies")
  message("==============================================\n")

  # Install CRAN packages
  tryCatch({
    install_cran_packages(update = update)
  }, error = function(e) {
    warning("Error installing CRAN packages: ", e$message)
  })

  # Install Bioconductor packages
  tryCatch({
    install_bioc_packages(update = update)
  }, error = function(e) {
    warning("Error installing Bioconductor packages: ", e$message)
  })

  # Install GitHub packages
  if (length(github_packages) > 0) {
    tryCatch({
      install_github_packages()
    }, error = function(e) {
      warning("Error installing GitHub packages: ", e$message)
    })
  }

  message("\n==============================================")
  message("Package installation complete!")
  message("==============================================\n")
}

#' Check package versions
check_package_versions <- function() {
  message("Checking critical package versions...")

  for (pkg in names(package_versions)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      installed_version <- as.character(packageVersion(pkg))
      required_version <- package_versions[[pkg]]

      if (package_version(installed_version) < package_version(required_version)) {
        warning(pkg, " version ", installed_version, " is older than required version ",
                required_version, ". Consider updating.")
      } else {
        message(pkg, " version ", installed_version, " \u2713")
      }
    } else {
      warning(pkg, " is not installed.")
    }
  }
}

#' Load all required packages
load_packages <- function() {
  all_packages <- unique(c(cran_packages, bioc_packages))

  message("Loading packages...")
  success <- sapply(all_packages, function(pkg) {
    suppressPackageStartupMessages(
      requireNamespace(pkg, quietly = TRUE)
    )
  })

  failed <- names(success)[!success]
  if (length(failed) > 0) {
    warning("Failed to load packages: ", paste(failed, collapse = ", "))
    message("Run install_all_packages() to install missing packages.")
    return(FALSE)
  } else {
    message("All packages loaded successfully \u2713")
    return(TRUE)
  }
}

#' Generate session info report
generate_session_info <- function(output_file = NULL) {
  info <- sessionInfo()

  if (!is.null(output_file)) {
    sink(output_file)
    print(info)
    sink()
    message("Session info saved to: ", output_file)
  }

  return(info)
}

# ------------------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------------------

# Only run if script is executed directly (not sourced)
if (sys.nframe() == 0) {
  cat("ICLE Project - Package Management\n")
  cat("==================================\n\n")
  cat("Available functions:\n")
  cat("  install_all_packages()     - Install all required packages\n")
  cat("  check_package_versions()   - Verify package versions\n")
  cat("  load_packages()            - Load all packages\n")
  cat("  generate_session_info()    - Generate session info report\n\n")
  cat("Example usage:\n")
  cat("  source('requirements.R')\n")
  cat("  install_all_packages()\n")
  cat("  load_packages()\n")
}

# ------------------------------------------------------------------------------
# Package Documentation
# ------------------------------------------------------------------------------

#' Package Documentation
#'
#' Data Manipulation:
#'   - dplyr / tidyr / data.table: Core tabular data manipulation
#'   - reshape2 / plyr: Legacy reshape and split-apply-combine helpers
#'   - magrittr: Pipe operator
#'   - vroom: Fast delimited file I/O
#'
#' Visualisation (ggplot2 ecosystem):
#'   - ggplot2: Grammar of graphics
#'   - ggpubr / cowplot / patchwork / gridExtra: Plot arrangement and annotation
#'   - ggrepel: Non-overlapping text labels
#'   - ggthemes / ggsci / viridis / scales: Themes and colour palettes
#'   - ggbeeswarm / ggpattern: Specialised geoms
#'   - ggsankey: Sankey / alluvial diagrams
#'   - circlize / ComplexHeatmap: Circular and heatmap visualisation
#'   - interacCircos / htmlwidgets / webshot2: Interactive circos → PDF export
#'
#' Genomics:
#'   - GenomicRanges / rtracklayer: Genomic interval handling
#'   - DNAcopy: Circular binary segmentation for CNV
#'   - CNTools: Gene-level CNV from segmentation
#'   - maftools: MAF file parsing and TMB calculation
#'   - ShatterSeek (GitHub): Chromothripsis detection
#'
#' DNA Methylation:
#'   - sesame: EPIC/HM450K array preprocessing (SeSAMe pipeline)
#'
#' Expression Analysis:
#'   - limma / edgeR / DESeq2: Differential expression
#'   - sva: Surrogate variable analysis / ComBat batch correction
#'
#' Subtyping and Clustering:
#'   - ConsensusClusterPlus: Consensus clustering
#'   - CancerSubtypes: Multi-omics cancer subtyping
#'   - BreastSubtypeR (GitHub): Intrinsic breast cancer subtyping
#'
#' Pathway and Enrichment:
#'   - clusterProfiler: GO / KEGG enrichment
#'   - hypeR: Gene set over-representation analysis
#'   - msigdbr: MSigDB gene sets in R
#'   - org.Hs.eg.db: Human gene annotation
#'
#' For detailed help on any package use: help(package = "packageName")
