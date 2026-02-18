# ICLE Analysis Script Index

Author: Osama Shiraz Shah

This document lists the order and roles of scripts used by `Main_Data_Analysis.Rmd` to avoid dependency conflicts when running scripts manually.

## 1. Configuration and helpers (run first)

- **config.R** – Project root, `DIRS`, `FILES`. Source from `2-Analysis/` (e.g. when knitting Main_Data_Analysis.Rmd).
- **Helper_Scripts/Helper_Functions.R** – Shared utilities, colors, themes (`annot_cols`, `heatmap_legend_param`, `myTheme`, etc.).

## 2. Data loading order (01_load_all_data.R)

`Main_Data_Analysis.Rmd` calls `load_all_icle_data(load_external = TRUE)` via `Data_Loading/01_load_all_data.R`, which runs, in order:

| Step | Script | Role |
|------|--------|------|
| 0 | Data_Loading/00_load_annotations.R | Cell line annotations (CL_Annots, CL_Annots_simple) |
| 1 | Data_Loading/02_load_rppa_data.R | RPPA (sources 1-Datasets/ICLE/RPPA/Prepare_RPPA_Data.R) |
| 2 | Data_Loading/03_load_rna_data.R | RNA-seq (loads from FILES$rnaseq_data) |
| 3 | Data_Loading/04_load_cnv_data.R | CNV (sources 1-Datasets/ICLE/CytoSNP/2_GenomeStudio/Prepare_SNP_Data.R) |
| 4 | Data_Loading/05_load_snv_data.R | SNV (sources 1-Datasets/ICLE/WES/Prepare_SNV_Data.R) |
| 5 | Data_Loading/06_load_dnam_data.R | DNAm (sources 1-Datasets/ICLE/DNAm/Prepare_DNAm_Data.R) |
| 6 | Data_Loading/07_load_sv_data.R | SV (sources 1-Datasets/ICLE/Bionano/2_Structural_Variations/Prepare_SV_Data.R) |
| 7 | Data_Loading/08_load_external_data.R | TCGA + MSK (sources 1-Datasets/External/TCGA/Load_TCGA_Data.R, Load_MSK_Data.R) |
| 8 | Data_Loading/09_generate_gams.R | Genomic alteration matrices |

## 3. Figure scripts (order in Main_Data_Analysis.Rmd)

| # | Script | Notes |
|---|--------|------|
| 01 | 01_SupFig1_Genotype_Similarity.R | Genotype similarity |
| 02 | 02_Fig1_SupFig2_3_4_Molecular_Subtyping.R | Fig 1 (subtyping, SupFig 2–4, Fig 1C SET, Fig 1B Sankey); sources 03_Fig1C_SET_Signature.R; assigns fig1b_sankey, fig1c_setheatmap |
| 04 | 04_Fig1D_Multiomics_Overview.R | Fig 1D |
| 05 | 05_Fig1F_Alteration_barplots.R | Fig 1F |
| 06 | 06_SupFig5_ILC_NST_Alterations.R | SupFig 5 |
| 07 | 07_SupFig6_Pathway_Alterations.R | SupFig 6 |
| 08 | 08_Fig2_CDH1_Alteration_Landscape_All.R | **Orchestrator**: sources 08_Fig2C, 09_Fig2D, 10_Fig2E, 11_Fig2F, 12_Fig2G, 13_Fig2H (Fig 2 CDH1) |
| 09–13 | 08_Fig2C … 13_Fig2H | Fig 2 panels (sourced by 08_Fig2_CDH1_Alteration_Landscape_All.R only) |
| 14 | 14_Fig3_SV_All.R | **Orchestrator**: sources 14_SupFig8, 15_Fig3A, 16_Fig3B, 17_Fig3C, 18_Fig3D_SupFig9, 19_Fig3D, 20_Fig3E_3F_SupFig10 |
| 15–20 | 14_SupFig8 … 20_Fig3E_3F_SupFig10 | Fig 3 + SupFig 8–10 (sourced by 14_Fig3_SV_All.R only); 17–18 sourced by 19–20 |
| 21 | 21_Fig4_SupFig11_DNAm_Alterations.R | Fig 4, SupFig 11 (differential RNA/DNAm functions inlined) |
| 22 | 22_Fig5_RNAi_Differential_Dependencies.R | Fig 5 (uses simem/ libs) |
| 23 | 23_Fig6_Patient_Signatures_Resemblance_Scores.R | Fig 6 |

## 4. Outputs (3-Results)

Figures and tables are written under `DIRS$results` and `DIRS$results_sub$...` (see config.R for exact paths). Main outputs include PDFs in `3-Results/` and subfolders (e.g. Molecular_Subtyping, CDH1_Alteration_Landscape, Optical_Genome_Mapping, DNA_Methylation, dependenciesGene_Dependencies, Molecular_Resemblance).