Single-cell RNA-seq Analysis Pipeline for Ovarian Endometriosis
Overview
This repository contains the complete R analysis pipeline for single-cell RNA-sequencing data from ovarian endometriosis (OEM) and normal ovarian tissue samples.

Study Design
Samples: 11 samples (6 normal, 5 OEM) from cortex and medulla regions
Technology: Single-cell RNA-sequencing
Subjects: Young adult females (Normal: n=3, OEM: n=3)
Analysis Workflow
1. Quality Control & Preprocessing
Cell-level QC filtering (nFeature_RNA > 300, percent.mt < 25%)
Doublet detection using DoubletFinder
Normalization and feature selection
2. Batch Effect Correction
Harmony integration across samples
Multiple resolution clustering
3. Cell Type Annotation (Level 1)
Identified 10 major cell types: Endometrial stroma cells, Mesenchymal cells, Smooth muscle cells, T/NK cells, Myeloid cells, Endothelial cells, Granulosa cells, Mast cells, B cells, Epithelial cells.
4. Differential Expression Analysis
OEM vs Normal comparison for each cell type
Identification of upregulated and downregulated genes
5. Cell Proportion Analysis
Sample-level and individual-level cell composition
Visualization with stacked bar plots and bubble plots
6. Sample Correlation Analysis
Expression correlation between samples
Hierarchical clustering and dendrogram visualization
7. Perturbation Analysis
Augur: Cell type-specific prioritization by AUC scores
scDist: Statistical distance measurement between conditions
Combined visualization of both metrics
