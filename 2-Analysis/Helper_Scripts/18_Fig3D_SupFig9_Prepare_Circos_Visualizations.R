# ==============================================================================
# Script 18: Fig 3D & SupFig 9 - Circos Plots (OGM SV, CNV, Chromothripsis)
# ==============================================================================
# Description: Prepares data and generates interactive Circos plots per sample
#              showing OGM SV, CNV, and chromothripsis events. Requires hg19_ideogram
#              in environment (e.g. from interacCircos or precomputed).
#
# Dependencies:
#   - config.R (paths and configuration)
#   - Helper_Functions.R (plotting parameters and color schemes)
#   - 01_load_all_data.R (data loading functions)
#
# Required Data Objects:
#   - ICLE_SV: Structural variation data
#   - CNV data (via config paths)
#   - Chromothripsis data (from ShatterSeek analysis)
#   - hg19_ideogram: Genome ideogram data
#
# Outputs (assigned to .GlobalEnv / available after sourcing):
#   - genome_tracks: Shared genome annotation tracks for circos
#   - prepare_sample_tracks: Function to build tracks for a sample
#   - draw_circos_plot: Function to render a circos plot
#   - draw_circos_plot_for_all: Function to generate circos for all samples (used by Main for SupFig 9)
#
# Author: Osama Shiraz Shah
# ==============================================================================

if (!exists("DIRS", envir = .GlobalEnv)) {
  if (file.exists("config.R")) source("config.R") else source("../config.R")
}

suppressPackageStartupMessages({
  library(hiAnnotator)
  library(interacCircos)
  library(GenomicRanges)
  library(htmlwidgets)
  library(webshot2)
})


# -------------------------------
# 2. Prepare Genome and Annotation Tracks
# -------------------------------
prepare_genome_tracks <- function(oncoVar_path, gene_bed_path, hg19_ideogram) {
  
  ## Cytoband track
  ideogram_track = CircosArc("arc01", data = hg19_ideogram, innerRadius = 218, outerRadius = 235, animationDisplay = F, 
                             animationDelay = 10, animationTime = 100)
  
  
  ## Gene track
  # Load oncoVar if not already loaded
  if (!exists("oncoVar", envir = .GlobalEnv)) {
    oncoVar <- read.delim(oncoVar_path, sep = "\t")
    rownames(oncoVar) <- oncoVar$Gene_symbol
  } else {
    oncoVar <- get("oncoVar", envir = .GlobalEnv)
  }

  # Load gene_bed if not already loaded
  if (!exists("gene_bed", envir = .GlobalEnv)) {
    gene_bed <- read.delim(gene_bed_path, header = FALSE)
  } else {
    gene_bed <- get("gene_bed", envir = .GlobalEnv)
  }
  gene_bed$V1 <- gsub("chr", "", gene_bed$V1)

  BRCA_genes <- subset(gene_bed, V4 %in% subset(oncoVar, Driver.Level %in% 1:4)$Gene_symbol)
  BRCA_genes_track <- BRCA_genes
  names(BRCA_genes_track) <- c("chr", "start", "end", "des")
  BRCA_genes_track$type <- oncoVar[BRCA_genes_track$des, "OncoKB"]
  BRCA_genes_track$color <- ifelse(BRCA_genes_track$type %in% "TSG|Pan", "blue",
                                   ifelse(BRCA_genes_track$type %in% "Oncogene|Pan", "red", "white"))
  BRCA_genes_track <- subset(BRCA_genes_track, color != "white")
  BRCA_genes_track <- subset(BRCA_genes_track, chr %in% unique(hg19_ideogram$chr))

  BRCA_OG_TSG_track <- CircosArc("arc02", data = BRCA_genes_track, innerRadius = 110, outerRadius = 122,
                                 animationDisplay = FALSE, animationDelay = 10, animationTime = 100)
  
  
  return(list(BRCA_OG_TSG_track = BRCA_OG_TSG_track, ideogram_track = ideogram_track))
}



# -------------------------------
# 3. Prepare Sample-Specific Circos Tracks
# -------------------------------
prepare_sample_tracks <- function(cellLine, sv_df, cnv_df, chromothripsis_df) {
  nonTransSVs <- subset(sv_df, Sample == cellLine & RefcontigID1 == RefcontigID2 & Type != "Tranloc-intra")
  nonTransSVs$RefcontigID1 <- gsub("chr", "", nonTransSVs$RefcontigID1)
  nonTransSVs$RefcontigID2 <- gsub("chr", "", nonTransSVs$RefcontigID2)
  nonTransSVs$ID <- paste0(nonTransSVs$RefcontigID1, nonTransSVs$RefStartPos, nonTransSVs$RefEndPos)
  nonTransSVs <- nonTransSVs[, c("RefcontigID1", "RefStartPos", "RefEndPos", "Type", "OverlapGenes")]
  nonTransSVs$name <- paste0(nonTransSVs$Type, " --> ", nonTransSVs$OverlapGenes)
  names(nonTransSVs) <- c("chr", "start", "end", "type", "genes", "name")
  nonTransSVs_list <- split(nonTransSVs, nonTransSVs$type)

  TransSVs <- subset(sv_df, Sample == cellLine & Type %in% c("Tranloc-intra", "Tranloc-inter"))
  TransSVs$RefcontigID1 <- gsub("chr", "", TransSVs$RefcontigID1)
  TransSVs$RefcontigID2 <- gsub("chr", "", TransSVs$RefcontigID2)
  TransSVs <- TransSVs[, c("RefcontigID1", "RefStartPos", "RefStartPos", "RefcontigID2", "RefEndPos", "RefEndPos", "Type", "OverlapGenes")]
  names(TransSVs) <- c("g1chr", "g1start", "g1end", "g2chr", "g2start", "g2end", "type", "fusion")
  TransSVs_list <- split(TransSVs, TransSVs$type)

  cn <- subset(cnv_df, Sample == cellLine)[,c("chromosome", "start", "end", "total_cn")]
  colnames(cn) <- c("chr", "start", "end", "value")
  cn$color <- ifelse(cn$value < 2, "blue", ifelse(cn$value > 2, "red", "gray"))

  thripsisData <- subset(chromothripsis_df, Sample == cellLine)[, c("chrom", "start", "end", "confidence")]
  names(thripsisData) <- c("chr", "start", "end", "des")
  thripsisData$color <- ifelse(thripsisData$des == "High", "#eded2a",
                               ifelse(thripsisData$des == "Low", "#a1d99b", "#EFEBE9"))
  thripsisData <- subset(thripsisData, des %in% c("High", "Low"))

  
  ## Chromothripsis track
  thripsis_track = CircosArc(data = thripsisData, 'thripsis_1', innerRadius = 199, outerRadius = 215, animationDisplay = TRUE, opacity = 1)
  
  
  ## SV track
  SV_tracks <- list(
    Del = if (!is.null(nonTransSVs_list$Del)) CircosScatter("SV_1", data = nonTransSVs_list$Del,
                                                            radius = 190, width = 0.5, animationDisplay = F, 
                                                            innerCircleColor = "#E64B35FF", outerCircleColor = "#E64B35FF", outerCircleSize = 2, 
                                                            displayLinkLabel = F, animationTime = 5000,animationDelay = 100) else NULL,
    Dup = if (!is.null(nonTransSVs_list$Dup)) CircosScatter("SV_2", data = nonTransSVs_list$Dup, 
                                                           radius = 175, width = 0.5, animationDisplay = F, 
                                                           innerCircleColor = "#4DBBD5FF", outerCircleColor = "#4DBBD5FF", outerCircleSize = 2, 
                                                           displayLinkLabel = F, animationTime = 5000,animationDelay = 100) else NULL,
    Ins = if (!is.null(nonTransSVs_list$Ins)) CircosScatter("SV_3", data = nonTransSVs_list$Ins, 
                                                            radius = 180, width = 0.5, animationDisplay = F, 
                                                            innerCircleColor = "#00A087FF", outerCircleColor = "#00A087FF", outerCircleSize = 2, 
                                                            displayLinkLabel = F, animationTime = 5000,animationDelay = 100) else NULL,
    Inv = if (!is.null(nonTransSVs_list$Inv)) CircosScatter("SV_4", data = nonTransSVs_list$Inv, 
                                                            radius = 185, width = 0.5, animationDisplay = F, 
                                                            innerCircleColor = "#3C5488FF", outerCircleColor = "#3C5488FF", outerCircleSize = 2, 
                                                            displayLinkLabel = F, animationTime = 5000,animationDelay = 100) else NULL
  )
  
  SV_track_bg = CircosBackground("bg02",axisShow = F, minRadius = 170, maxRadius = 195, borderSize = 0.5, 
                                 fillColors = "#ECEFF1", borderColors = "#CFD8DC", animationDisplay = F)

  ## Copy Number track
  cn_track = CircosCnv("cn1", data = cn, minRadius = 130, maxRadius = 160, animationDisplay = F, 
                       animationDelay = 1, animationTime = 50, strokeWidth = 0, width = 3, animationType = "linear")
  cn_track_bg = CircosBackground("bg01",axisShow = F, minRadius = 125, maxRadius = 165, borderSize = 0, fillColors = "#eeeeff",
                                 animationDisplay = F)

  ## Translocations track
  link_track01 <- if (!is.null(TransSVs_list$`Tranloc-intra`)) link_track01 = CircosLink("link01", data = TransSVs_list$`Tranloc-intra`,
                                                                                         radius = 105, width = 1, animationDisplay = T, fillColor = "#8491B4FF", 
                                                                                         displayLinkLabel = F, animationTime = 5000,animationDelay = 10) else NULL
  link_track02 <- if (!is.null(TransSVs_list$`Tranloc-inter`)) CircosLink("link02",data = TransSVs_list$`Tranloc-inter`, radius = 105, width = 1, 
                                                                          fillColor = "#F39B7FFF", animationDisplay = T, displayLinkLabel = F, 
                                                                          animationTime = 5000, animationDelay = 10) else NULL 

  return(list(
    thripsis_track = thripsis_track,
    SV_track_bg = SV_track_bg,
    SV_tracks = SV_tracks,
    cn_track = cn_track,
    cn_track_bg = cn_track_bg,
    link_track01 = link_track01,
    link_track02 = link_track02
  ))
}




# -------------------------------
# 4. Draw Circos Plot
# -------------------------------

draw_circos_plot <- function(tracks, genome = "hg19", outerRadius = 230) {
  plt = Circos(
    genome = genome, outerRadius = outerRadius,
    genomeFillColor = rep("white", 24), genomeBorderColor = "white", genomeBorderSize = 0,
    moduleList = tracks$ideogram_track + tracks$thripsis_track + tracks$SV_track_bg +
      tracks$SV_tracks$Del + tracks$SV_tracks$Dup + tracks$SV_tracks$Ins + tracks$SV_tracks$Inv +
      tracks$cn_track + tracks$cn_track_bg + tracks$BRCA_OG_TSG_track +
      tracks$link_track01 + tracks$link_track02, 
    ARCMouseOverDisplay = F, SCATTERMouseOverDisplay = F, CNVMouseOverDisplay = F, LINKMouseOverDisplay = F
  )
  return(plt)
}





# -------------------------------
# Make Circos Plots for All ICLE Samples using SV, CN and Chromothripsis data
# -------------------------------

# Prepare Genome Track
if (!exists("hg19_ideogram", envir = .GlobalEnv)) {
  if (exists("getIdeogramData", envir = asNamespace("interacCircos"))) {
    assign("hg19_ideogram", interacCircos::getIdeogramData(genome = "hg19"), envir = .GlobalEnv)
  } else {
    stop("hg19_ideogram not found. Create it or load interacCircos and ensure getIdeogramData exists.")
  }
}
genome_tracks <- prepare_genome_tracks(FILES$oncovar, FILES$refgene_unique, get("hg19_ideogram", envir = .GlobalEnv))



draw_circos_plot_for_all <- function(sv_df, cnv_df, chromothripsis_df, out_dir){
  
  # Ensure output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  for (cellLine in unique(sv_df$Sample)){
    
    print(paste0("Making Circos Plot for: ", cellLine))
    
    sample_tracks <- prepare_sample_tracks(cellLine, sv_df, cnv_df, chromothripsis_df)
    circos_widget <- draw_circos_plot(c(genome_tracks, sample_tracks))
    
    html_file <- file.path(out_dir, paste0(cellLine, "_circos_plot.html"))
    pdf_file  <- file.path(out_dir, paste0(cellLine, "_circos_plot.pdf"))
    
    # 3. Save as HTML first (Intermediate step)
    # selfcontained = FALSE can sometimes be faster/lighter if you delete it immediately, 
    # but TRUE is safer for portability.
    htmlwidgets::saveWidget(circos_widget, file = html_file, selfcontained = TRUE)
    
    # 4. Convert HTML to PDF using webshot2
    # delay: gives the JS time to render/animate before taking the snapshot
    # zoom: increases resolution (2 = 2x resolution)
    tryCatch({
      webshot2::webshot(url = html_file, file = pdf_file, vwidth = 800, vheight = 20, delay = 3, zoom = 2, selector = ".html-widget")
      message("  Saved PDF for sample: ", cellLine)
    }, error = function(e) {
      message("  Failed to save PDF for: ", cellLine, " - ", e$message)
    })
  
  }
}

# ==============================================================================
# Usage Example
# ==============================================================================
# source("config.R")
# source("Helper_Scripts/Data_Loading/01_load_all_data.R")
# source("Helper_Scripts/13_Fig3_Identify_Chromothripsis.R")
# source("Helper_Scripts/14_Fig3D_Make_Circos_Plots.R")
#
# load_all_icle_data(load_external = TRUE)
# shatterseek_outs <- run_shatterseek_main(FILES$cnv_seg, ICLE_SV, DIRS$results_sub$ogm)
#
# trk <- c(prepare_sample_tracks("BCK4", ICLE_SV, shatterseek_outs$cnv_df, shatterseek_outs$chromothripsis_df), genome_tracks)
# draw_circos_plot(trk)
#
# draw_circos_plot_for_all(shatterseek_outs$sv_df, shatterseek_outs$cnv_df, shatterseek_outs$chromothripsis_df,
#                          out_dir = file.path(DIRS$results_sub$ogm, "SupFig9_OGM_Circos"))


