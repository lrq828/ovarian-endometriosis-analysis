# ==============================================================================
# Section 1: Quality Control and Preprocessing
# ==============================================================================

## Load required libraries ####
library(Seurat)
library(data.table)
library(DoubletFinder)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(pheatmap)
set.seed(123)

## Set working directory
setwd("project_dir/")

## Define data paths and sample names ####
data_dir <- "project_dir"
sample_names <- c('Normal.Y01.Cortex','Normal.Y01.Medulla','Normal.Y02.Cortex',
                  'Normal.Y03.Cortex','Normal.Y03.Medulla',
                  'OEM.Y01.Cortex','OEM.Y01.Medulla','OEM.Y02.Cortex',
                  'OEM.Y02.Medulla','OEM.Y03.Cortex','OEM.Y03.Medulla')

## Get CSV files from directory
csv_files <- list.files(path = data_dir, pattern = "*.csv", full.names = TRUE)

## Initialize storage lists and QC dataframe ####
sce_list <- list()
sce_unfilter_list <- list()

sample_info <- data.frame(
  sample = character(),
  `number_of_cells_unfilter` = integer(),
  `UMIs_0.75` = numeric(),
  `UMIs_0.25` = numeric(),
  `UMIs_Median` = integer(),
  `genes_0.75` = numeric(),
  `genes_0.25` = numeric(),
  `genes_Median` = integer(),
  `percent.mt_0.75` = numeric(),
  `percent.mt_0.25` = numeric(),
  `percent.mt_median` = numeric(),
  `Doublet` = integer(),
  `Singlet` = integer(),
  `Doublet_ratio` = numeric(),
  `number_of_cells_filter` = integer(),
  stringsAsFactors = FALSE
)

## Define doublet detection function ####
Find_doublet <- function(data) {
  # Perform parameter sweep to find optimal pK value
  sweep.res.list <- paramSweep(data, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
  
  # Calculate expected doublet rate
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  Doubletrate <- 0.08 * (length(rownames(data@meta.data)) / 10000)
  nExp_poi <- round(Doubletrate * ncol(data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Identify doublets
  data <- doubletFinder(data, PCs = 1:20, pN = 0.25, pK = pK, 
                        nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
  return(data)
}

## Create output directories ####
if (!dir.exists("./Result/01.L1/01.QC/Doublet/")) {
  dir.create("./Result/01.L1/01.QC/Doublet/", recursive = TRUE)
}
if (!dir.exists("./Result/01.L1/01.QC/feature_count_mt/")) {
  dir.create("./Result/01.L1/01.QC/feature_count_mt/", recursive = TRUE)
}

## Process each sample ####
for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  data_path <- csv_files[i]
  
  # Read and format data
  temp <- fread(file = data_path, header = TRUE)
  Cell_id <- temp$Cell_Index
  temp <- temp[, -1]
  temp <- as.sparse(t(temp))
  colnames(temp) <- paste0(sample_name, Cell_id)
  
  # Create Seurat object
  sce <- CreateSeuratObject(counts = temp, project = sample_name, 
                            min.cells = 3, min.features = 200)
  
  # Add metadata
  sce$sample <- sample_name
  sce$class <- strsplit(sample_name, split = "\\.")[[1]][1]
  sce$people <- paste0(strsplit(sample_name, split = "\\.")[[1]][1], '.', 
                       strsplit(sample_name, split = "\\.")[[1]][2])
  sce$Age <- ifelse(substring(strsplit(sample_name, split = "\\.")[[1]][2],1,1)=='Y', 
                    'Young', 'Middle')
  sce$histologic_type <- strsplit(sample_name, split = "\\.")[[1]][3]
  
  # Calculate mitochondrial percentage
  sce[['percent.mt']] <- PercentageFeatureSet(sce, pattern = '^MT-')
  
  # Store unfiltered data
  sce_unfilter_list[[sample_name]] <- sce
  
  # Plot QC metrics before filtering
  p <- VlnPlot(sce, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), 
               pt.size = 0, group.by = 'sample')
  ggsave(paste0("./Result/01.L1/01.QC/feature_count_mt/", sample_name, "_unfilter.pdf"), 
         plot = p)
  
  # Calculate QC statistics
  nCount_info <- quantile(sce$nCount_RNA, probs = c(0.25, 0.5, 0.75))
  nFeature_info <- quantile(sce$nFeature_RNA, probs = c(0.25, 0.5, 0.75))
  mt_info <- quantile(sce$percent.mt, probs = c(0.25, 0.5, 0.75))
  
  # Update sample info
  sample_info <- rbind(sample_info, data.frame(
    sample = sample_name,
    `number_of_cells_unfilter` = ncol(sce),
    `UMIs_0.75` = nCount_info[3],
    `UMIs_0.25` = nCount_info[1],
    `UMIs_Median` = nCount_info[2],
    `genes_0.75` = nFeature_info[3],
    `genes_0.25` = nFeature_info[1],
    `genes_Median` = nFeature_info[2],
    `percent.mt_0.75` = mt_info[3],
    `percent.mt_0.25` = mt_info[1],
    `percent.mt_median` = mt_info[2],
    Doublet = NA,
    Singlet = NA,
    `Doublet_ratio` = NA,
    `number_of_cells_filter` = NA
  ))
  
  # Perform QC filtering, normalization, and dimensionality reduction
  sce <- sce %>%
    subset(subset = nFeature_RNA > 300 & percent.mt < 25) %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:20, reduction = "pca", n.neighbors = 30L, min.dist = 0.3) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters() %>%
    Find_doublet()  # Identify doublets (computationally intensive)
  
  # Update doublet statistics
  sample_info[nrow(sample_info), 'Doublet'] <- sum(sce$doublet_info == "Doublet")
  sample_info[nrow(sample_info), 'Singlet'] <- sum(sce$doublet_info == "Singlet")
  sample_info[nrow(sample_info), 'Doublet_ratio'] <- sum(sce$doublet_info == "Doublet") / ncol(sce)
  
  # Save doublet UMAP plot
  p1 <- DimPlot(sce, reduction = "umap", group.by = "doublet_info")
  ggsave(paste0("./Result/01.L1/01.QC/Doublet/", sample_name, "_doublet.pdf"), plot = p1)
  
  # Filter doublets
  sce <- subset(sce, subset = doublet_info == 'Singlet')
  sample_info[nrow(sample_info), 'number_of_cells_filter'] <- ncol(sce)
  
  # Save filtered Seurat object
  saveRDS(sce, paste0('./Data/01.All_QC_Rds/', sample_name, '_filter.rds'))
  sce_list[[sample_name]] <- sce
}

## Merge all samples ####
sce.all <- merge(sce_list[[1]], y=sce_list[2:length(sce_list)], project ="ovary_OEM")

## Save QC summary ####
write.csv(sample_info, file = "./Result/01.L1/01.QC/sample_QC_info.csv", row.names = TRUE)

## Plot merged QC metrics ####
p2 <- VlnPlot(sce.all, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), 
              pt.size = 0, group.by = 'Sample',
              cols = c('#E5D2DD', '#F1BB72', '#53A85F', '#F3B1A0', '#D6E7A3', '#57C3F3', 
                       '#E59CC4', '#6778AE', '#3A6963', '#BD956A', '#8C549C', '#C1E6F3'))
ggsave("./Result/01.L1/01.QC/AllSample_feature_count_mt.pdf", plot = p2, 
       width = 22, height = 7)

saveRDS(sce.all, './Data/01.All_QC_Rds/merge_filter.rds')


# ==============================================================================
# Section 2: Batch Effect Correction with Harmony
# ==============================================================================

## Load required libraries ####
library(Seurat)
library(ggplot2)
library(future)
library(tidyverse)
library(harmony)
library(pheatmap)
library(ggpubr)
library(rstatix)
library(scCustomize)
library(cowplot)

## Load filtered data ####
path <- "project_dir/Data/01.All_QC_Rds/"
file <- list.files(path, pattern = 'rds')

sce <- list()
for (i in 1:length(file)) {
  dat <- readRDS(paste0(path, file[i]))
  sce[[i]] <- dat
}

sce <- merge(sce[[1]], y=sce[2:length(sce)], project ="ovary_OEM")

## Set factor levels for metadata ####
Sample_order <- c('Normal.Y01.Cortex','Normal.Y01.Medulla','Normal.Y02.Cortex',
                  'Normal.Y03.Cortex','Normal.Y03.Medulla',
                  'OEM.Y01.Cortex','OEM.Y01.Medulla','OEM.Y02.Cortex',
                  'OEM.Y02.Medulla','OEM.Y03.Cortex','OEM.Y03.Medulla')
People_order <- c("Normal.Y01","Normal.Y02","Normal.Y03","OEM.Y01","OEM.Y02","OEM.Y03")
Class_order <- c("Normal","OEM")
Histologic_type_order <- c("Cortex", "Medulla")

sce@meta.data[["orig.ident"]] <- factor(sce@meta.data[["orig.ident"]], levels = Sample_order)
sce@meta.data[["Sample"]] <- factor(sce@meta.data[["Sample"]], levels = Sample_order)
sce@meta.data[["People"]] <- factor(sce@meta.data[["People"]], levels = People_order)
sce@meta.data[["Class"]] <- factor(sce@meta.data[["Class"]], levels = Class_order)
sce@meta.data[["Histologic_type"]] <- factor(sce@meta.data[["Histologic_type"]], 
                                             levels = Histologic_type_order)

## Calculate percentages for specific gene types ####
erccs <- grep('^ERCC-', x= rownames(sce), value = T)
RP <- grep("^RP[SL][[:digit:]]", x= rownames(sce), value = T)
mt <- grep('^MT-', x= rownames(sce), value = T)
ncRNA <- grep("^[A-Z][A-Z][0-9]*\\.[0-9]", x= rownames(sce), value = T)
LOC <- grep('(^LOC|LINC)[1-9]*', x= rownames(sce), value = T)

# Hemoglobin genes
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes_total, rownames(sce@assays$RNA))
HB.genes <- rownames(sce@assays$RNA)[HB_m]
HB <- HB.genes[!is.na(HB.genes)]

sce <- PercentageFeatureSet(sce, pattern = "^ERCC-", col.name = "percent.ercc") %>%
  PercentageFeatureSet(pattern = "^RP[SL][[:digit:]]", col.name = "percent.RP") %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  PercentageFeatureSet(pattern = "^[A-Z][A-Z][0-9]*\\.[0-9]", col.name = "percent.ncRNA") %>%
  PercentageFeatureSet(pattern = "(^LOC|LINC)[1-9]*", col.name = "percent.LOC") %>%
  PercentageFeatureSet(features=HB, col.name = "percent.HB")

## Standard preprocessing pipeline ####
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2500)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(object = sce)

ElbowPlot(sce, ndims=50)

## Clustering before batch correction ####
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, dims = 1:15, n.neighbors = 30L, min.dist = 0.3)

## Save UMAP plots before Harmony ####
p1 <- UMAPPlot(sce, reduction = "umap", label = T, raster=FALSE) + NoLegend()
ggsave("./Result/01.L1/02.Ann/harmony_before_cluster.pdf", plot=p1, width=5, height=5)

p2 <- DimPlot(sce, reduction = "umap", group.by="Sample", shuffle = T, raster=FALSE)
ggsave("./Result/01.L1/02.Ann/harmony_before_Sample.pdf", plot=p2, width=5, height=5)

p3 <- DimPlot(sce, reduction = "umap", group.by="Class", shuffle = T, raster=FALSE)
ggsave("./Result/01.L1/02.Ann/harmony_before_Class.pdf", plot=p3, width=5, height=5)

p4 <- DimPlot(sce, reduction = "umap", group.by="Histologic_type", shuffle = T, raster=FALSE)
ggsave("./Result/01.L1/02.Ann/harmony_before_Histologic_type.pdf", plot=p4, width=5, height=5)

## Apply Harmony batch correction ####
sce <- RunHarmony(sce, "Sample", max.iter.harmony = 10)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:20)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:20)

## Rerun with optimized parameters ####
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(object = sce)
sce <- RunHarmony(sce, c("Sample"), max.iter.harmony = 10)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:20, k.param = 25)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:20)
sce <- FindClusters(sce, resolution = 1.5)

## Save UMAP plots after Harmony ####
p1 <- UMAPPlot(sce, reduction = "umap", label = T, raster=FALSE) + NoLegend()
ggsave("./Result/01.L1/02.Ann/harmony_after_cluster.pdf", plot=p1, width=5, height=5)

p2 <- DimPlot(sce, reduction = "umap", group.by="Sample", shuffle = T, raster=FALSE)
ggsave("./Result/01.L1/02.Ann/harmony_after_Sample.pdf", plot=p2, width=5, height=5)

p3 <- DimPlot(sce, reduction = "umap", group.by="Class", shuffle = T, raster=FALSE)
ggsave("./Result/01.L1/02.Ann/harmony_after_Class.pdf", plot=p3, width=5, height=5)

p4 <- DimPlot(sce, reduction = "umap", group.by="Histologic_type", shuffle = T, raster=FALSE)
ggsave("./Result/01.L1/02.Ann/harmony_after_Histologic_type.pdf", plot=p4, width=5, height=5)

## Test multiple resolutions ####
resolutions <- c(0.01, 0.05, 0.1, 0.5, 0.8, 1, 1.2, 1.5)
p <- list()

for (i in seq_along(resolutions)) {
  res <- resolutions[i]
  sce <- FindClusters(sce, resolution = res)
  p[[i]] <- UMAPPlot(sce, reduction = "umap", label = TRUE, raster = FALSE) + 
    ggtitle(paste("Resolution =", res)) + NoLegend()
}

fig3 <- plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]], ncol = 2)
save_plot("./Result/01.L1/02.Ann/UMAP_Clusters_by_Resolution.pdf", fig3,
          ncol = 2, nrow = 4, base_aspect_ratio = 2, 
          base_height = 5, base_width = 5)

## Clean up metadata columns ####
sce@meta.data <- sce@meta.data[, -c(15:24)]
sce@meta.data <- sce@meta.data[, -13]

saveRDS(sce, "./Data/01.All_QC_Rds/merge_harmony.rds")


# ==============================================================================
# Section 3: Harmony Visualization
# ==============================================================================

## Load data ####
library(Seurat)
library(ggplot2)
library(tidyverse)
library(scCustomize)

sce <- readRDS("project_dir/Data/01.All_QC_Rds/merge_harmony.rds")

## Generate UMAP plots by different groupings ####
# By People
p1 <- DimPlot(sce, reduction = "umap", group.by="People", shuffle = T, 
              raster=FALSE, pt.size = 0.1,
              cols = c("#D2EBC8","#7DBFA7","#AECDE1",
                       "#F5CFE4","#B383B9","#8FA4AE")) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"))
ggsave("./Result/01.L1/02.Ann/harmony_after_People.pdf", plot=p1, width=6.5, height=5)

# By Sample
p2 <- DimPlot(sce, reduction = "umap", group.by="Sample", shuffle = T, 
              raster=FALSE, pt.size = 0.1,
              cols = c('#E5D2DD', '#F1BB72', '#53A85F', '#F3B1A0', '#D6E7A3', '#57C3F3', 
                       '#E59CC4', '#6778AE', '#3A6963', '#BD956A', '#8C549C', '#C1E6F3')) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"))
ggsave("./Result/01.L1/02.Ann/harmony_after_Sample.pdf", plot=p2, width=7, height=5)

# By Class
p3 <- DimPlot(sce, reduction = "umap", group.by="Class", shuffle = T, 
              raster=FALSE, pt.size = 0.1,
              cols = c('#476D87','#E95C59')) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"))
ggsave("./Result/01.L1/02.Ann/harmony_after_Class.pdf", plot=p3, width=6, height=5)

# By Histologic type
p4 <- DimPlot(sce, reduction = "umap", group.by="Histologic_type", shuffle = T, 
              raster=FALSE, pt.size = 0.1,
              cols = c("#C381A8", "#447BAD")) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"))
ggsave("./Result/01.L1/02.Ann/harmony_after_Histologic_type.pdf", plot=p4, width=6, height=5)


# ==============================================================================
# Section 4: Level 1 (L1) Cell Type Annotation
# ==============================================================================

## Load data ####
library(Seurat)
sce <- readRDS("project_dir/Data/01.All_QC_Rds/merge_harmony.rds")

## Define cell type annotations ####
Idents(sce) <- sce$RNA_snn_res.1.5

new.cluster.ids <- c(
  "0" = "Mesenchymal cells",
  "1" = "Smooth muscle cells",
  "2" = "Mesenchymal cells",
  "3" = "Endometrial stroma cells",
  "4" = "Mesenchymal cells",
  "5" = "Smooth muscle cells",
  "6" = "Endometrial stroma cells",
  "7" = "Endometrial stroma cells", 
  "8" = "TNK cells",
  "9" = "Mesenchymal cells",
  "10" = "Mesenchymal cells", 
  "11" = "Smooth muscle cells", 
  "12" = "Endothelial cells", 
  "13" = "Myeloid cells", 
  "14" = "TNK cells",
  "15" = "Mesenchymal cells", 
  "16" = "Myeloid cells", 
  "17" = "TNK cells", 
  "18" = "Myeloid cells", 
  "19" = "Mesenchymal cells", 
  "20" = "Endometrial stroma cells", 
  "21" = "Mast cells",
  "22" = "Granulosa cells", 
  "23" = "Smooth muscle cells", 
  "24" = "B cells", 
  "25" = "Epithelial cells", 
  "26" = "B cells", 
  "27" = "Granulosa cells",
  "28" = "Myeloid cells", 
  "29" = "Smooth muscle cells", 
  "30" = "Mesenchymal cells", 
  "31" = "TNK cells", 
  "32" = "Endometrial stroma cells", 
  "33" = "Smooth muscle cells"
)

## Rename cluster identities ####
sce <- RenameIdents(sce, new.cluster.ids)
sce$L1_Celltype <- sce@active.ident

## Set factor levels ####
Celltype_order <- c("Endometrial stroma cells","Mesenchymal cells","Smooth muscle cells",
                    "TNK cells","Myeloid cells","Endothelial cells","Granulosa cells",
                    "Mast cells","B cells","Epithelial cells")

sce$L1_Celltype <- factor(sce$L1_Celltype, levels = Celltype_order)
Idents(sce) <- sce$L1_Celltype

sce$L1_Celltype_Class <- paste0(sce$L1_Celltype, "_", sce$Class)

saveRDS(sce, "./Data/02.Ann/annotation(L1).rds")


# ==============================================================================
# Section 5: L1 Cell Type Visualization
# ==============================================================================

## Load data and define colors ####
library(Seurat)
library(ggplot2)
library(scCustomize)

sce <- readRDS("project_dir/Data/02.Ann/annotation(L1).rds")

col <- c(
  `Endometrial stroma cells`='#476D77',
  `Mesenchymal cells`='#476D87',
  `Smooth muscle cells`='#E5D2DD',
  `TNK cells`='#F1BB72',
  `Myeloid cells`='#F3B1A0',
  `Endothelial cells`='#6778AE',
  `Granulosa cells`='#57C3F3',
  `Mast cells`='#53A85F',
  `B cells`='#E59CC4',
  `Epithelial cells`='#D6E7A3'
)

## Generate L1 cell type UMAP ####
p <- UMAPPlot(sce, reduction = "umap", label = T, raster=FALSE, 
              group.by="L1_Celltype", shuffle = T, pt.size = 0.1,
              cols = col, label.size = 5.5) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid")) + 
  NoLegend()
ggsave("./Result/01.L1/02.Ann/L1_Celltype.pdf", plot=p, width=10, height=10)


## Create dot plots for marker genes ####
marker_list <- list(
  Epithelial = c("KRT8", "KRT19"),
  B = c("CD79A", "MS4A1"),
  Mast = c("TPSB2"), 
  Granulosa = c("DSP","GATM","FOXL2"),
  Endothelial = c("CD34", "ESAM", "VWF"),
  Myeloid = c("LYZ", "CD14", "CD68"),
  TNK = c("CD2","CD3E","CCL5"),
  SMC = c("ACTA2", "TAGLN"),
  Endometrial = c("MME","ESR1","PGR"),
  Mesenchymal = c("DCN","PDGFRA","MMP2")
)

markers <- unlist(marker_list, use.names = FALSE)

## Dot plot with marker categories ####
p <- DotPlot(sce, features = marker_list, group.by = "L1_Celltype") +
  RotatedAxis() + 
  labs(x = 'Markers', y = 'L1_Celltypes') + 
  scale_color_gradientn(values = seq(0,1,0.2), 
                        colours = c("#5891BF","grey90","#AA3538"))
ggsave("./Result/01.L1/02.Ann/marker_Dotplot_L1_Celltype.pdf", p, 
       width = 14.5, height = 5)

## Dot plot without marker labels ####
p <- DotPlot(sce, features = markers, group.by = "L1_Celltype") +
  RotatedAxis() + 
  labs(x = 'Markers', y = 'L1_Celltypes') + 
  scale_color_gradientn(values = seq(0,1,0.2), 
                        colours = c("#5891BF","grey90","#AA3538"))
ggsave("./Result/01.L1/02.Ann/marker_Dotplot_L1_Celltype(nomarkerlabel).pdf", p, 
       width = 12, height = 5)

## Horizontal dot plot ####
plot <- DotPlot(sce, features = markers, group.by = "L1_Celltype") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle=60, hjust=1, vjust=1)) +
  labs(x=NULL, y=NULL) +
  guides(size=guide_legend(order=3)) +
  scale_color_gradientn(values = seq(0,1,0.2), 
                        colours = c("#5891BF","grey90","#AA3538"))
ggsave("./Result/01.L1/02.Ann/marker_Dotplot_L1_Celltype_1.pdf", plot, 
       width=5.5, height=8)

## Cell type distribution by people ####
pdf("./Result/01.L1/02.Ann/L1_Celltype_number_by_People.pdf", width = 20, height = 8)
gplots::balloonplot(table(sce$People, sce$L1_Celltype))
dev.off()


# ==============================================================================
# Section 6: Cell Proportion Analysis
# ==============================================================================

## Load libraries and data ####
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggalluvial)
library(dplyr)

sce <- readRDS("project_dir/Data/02.Ann/annotation(L1).rds")

## Prepare cell metadata ####
Idents(sce) <- sce$L1_Celltype
cell <- data.frame(
  Cell_ID = names(sce@active.ident),
  Cluster = sce@meta.data[["L1_Celltype"]],
  Class = sce@meta.data[["Class"]],
  Sample = sce@meta.data[["Sample"]],
  People = sce@meta.data[["People"]],
  Celltype = sce@active.ident,
  Histologic_type = sce@meta.data[["Histologic_type"]],
  Age = sce@meta.data[["Age"]]
)

## Cell numbers by cell type ####
cell_counts <- cell %>%
  group_by(Celltype) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

cell$Celltype <- factor(cell$Celltype, levels = as.character(cell_counts$Celltype))

p1 <- ggplot(cell, aes(y = Celltype, fill = Celltype)) +
  geom_bar(position = "identity") +
  labs(x = "Cell numbers", y = "Celltype") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  coord_flip() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = cell_counts, aes(y = Celltype, x = count, label = count),
            position = position_stack(vjust = 1), size = 3)
ggsave("./Result/01.L1/02.Ann/Cell_numbers_by_celltype.pdf", p1, width=8, height=4)

## Sample by cell type proportion ####
# Filter non-immune cell types
cell_filtered <- cell %>% 
  filter(Celltype %in% c("Mesenchymal cells","Endometrial stroma cells", 
                         "Smooth muscle cells","Endothelial cells", 
                         "Epithelial cells","Granulosa cells"))

cell_filtered$Sample <- factor(cell_filtered$Sample, levels = rev(levels(cell_filtered$Sample)))
cell_filtered$People <- factor(cell_filtered$People, levels = rev(levels(cell_filtered$People)))

# Calculate proportions by sample
cell_counts <- cell_filtered %>%
  group_by(Sample, Celltype) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Stacked bar plot
p4 <- ggplot(cell_counts, aes(x = proportion, y = Sample, fill = Celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
  labs(x = "Proportion", y = "Sample", fill = "Celltype", 
       title = "Sample_by_Celltype_proportion") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        panel.grid = element_blank())
ggsave("./Result/01.L1/02.Ann/Sample_by_Celltype_proportion.pdf", p4, width=8, height=5)

# Bubble plot by sample
cell_counts_bubble <- cell_filtered %>%
  group_by(Sample, Celltype) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p_bubble <- ggplot(cell_counts_bubble, aes(x = Celltype, y = Sample)) +
  geom_point(aes(size = proportion, color = count), stroke = 0.5) +
  scale_size(range = c(0.5, 10)) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  labs(x = "Cell Type", y = "Sample", size = "Proportion", color = "Cell Numbers",
       title = "Bubble plot of Cell Proportion (Sample)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8))
ggsave("./Result/01.L1/02.Ann/Sample_by_Celltype_bubbleplot.pdf", p_bubble, 
       width=7, height=5)

## People by cell type proportion ####
cell_counts <- cell_filtered %>%
  group_by(People, Celltype) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Stacked bar plot by people
p11 <- ggplot(cell_counts, aes(x = proportion, y = People, fill = Celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
  labs(x = "Proportion", y = "People", fill = "Celltype", 
       title = "People_by_Celltype_proportion") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        panel.grid = element_blank())
ggsave("./Result/01.L1/02.Ann/People_by_Celltype_proportion.pdf", p11, 
       width=8, height=5)

# Bubble plot by people
cell_counts_bubble <- cell_filtered %>%
  group_by(People, Celltype) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p_bubble <- ggplot(cell_counts_bubble, aes(x = Celltype, y = People)) +
  geom_point(aes(size = proportion, color = count), stroke = 0.5) +
  scale_size(range = c(0.5, 10)) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  labs(x = "Cell Type", y = "People", size = "Proportion", color = "Cell Numbers",
       title = "Bubble plot of Cell Proportion (People)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8))
ggsave("./Result/01.L1/02.Ann/People_by_Celltype_bubbleplot.pdf", p_bubble, 
       width=6, height=5)


# ==============================================================================
# Section 7: Differential Gene Expression Analysis (DEG by Class)
# ==============================================================================

## Load libraries and data ####
library(Seurat)
library(ggplot2)
library(ggpubr)

sce <- readRDS("project_dir/Data/02.Ann/annotation(L1).rds")

## Perform DEG analysis for each cell type ####
Idents(sce) <- sce$Class

degs <- lapply(unique(sce$L1_Celltype)[1:9], function(x){
  FindMarkers(
    sce[, sce$L1_Celltype == x],
    ident.1 = 'OEM',
    ident.2 = 'Normal',
    mean.fxn = function(x){return(log(x = rowMeans(x = expm1(x = x)) + 1, base = 2))}
  )
})
names(degs) <- unique(sce$L1_Celltype)[1:9]

## Merge DEG results ####
merged_degs <- do.call(rbind, lapply(names(degs), function(cell_type) {
  df <- degs[[cell_type]]
  df$gene <- rownames(df)
  df$CellType <- cell_type
  return(df)
}))

## Filter upregulated and downregulated genes ####
upregulated_degs <- merged_degs[merged_degs$avg_log2FC > 0.5 & merged_degs$p_val_adj < 0.05, ]
downregulated_degs <- merged_degs[merged_degs$avg_log2FC < -0.5 & merged_degs$p_val_adj < 0.05, ]

## Save DEG results ####
write.csv(upregulated_degs, file = "./Result/01.L1/03.ER/L1_Celltype_DEG_by_Class_UP.csv", 
          row.names = TRUE)
write.csv(downregulated_degs, file = "./Result/01.L1/03.ER/L1_Celltype_DEG_by_Class_DN.csv", 
          row.names = TRUE)
write.csv(merged_degs, file = "./Result/01.L1/03.ER/L1_Celltype_DEG_by_Class_all.csv", 
          row.names = TRUE)

## Create bar plot of DEG counts ####
data <- data.frame(
  Celltype = c("B cells", "Endothelial cells", "Epithelial cells", "Granulosa cells",
               "Mast cells", "Mesenchymal cells", "Myeloid cells", "Smooth muscle cells",
               "TNK cells"),
  Up = c(4, 152, 14, 883, 9, 266, 274, 202, 44),
  Down = c(255, 168, 159, 894, 128, 277, 349, 138, 66)
)

data_long <- reshape2::melt(data, id.vars = "Celltype", 
                            variable.name = "regulation", value.name = "count")

p1 <- ggplot(data_long, aes(x = Celltype, y = count, fill = regulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Up" = "#CC3333", "Down" = "#0099CC")) +
  labs(y = "Count", x = "Celltype", fill = "Regulation") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = 0.4) + 
  ggtitle("OEM VS Normal")
ggsave("./Result/01.L1/03.ER/L1_DEG_by_Class_Count_barplot.pdf", p1, 
       width = 12, height = 5)


# ==============================================================================
# Section 8: Sample Correlation Analysis
# ==============================================================================

## Load libraries and data ####
library(Seurat)
library(ComplexHeatmap)
library(ggplot2)
library(ggtree)
library(treeio)
library(dendextend)

sce <- readRDS("project_dir/Data/02.Ann/annotation(L1).rds")

## Aggregate expression by sample ####
av <- AggregateExpression(sce, group.by = c("Sample"), assays = "RNA")
av <- as.data.frame(av[[1]])

## Select top 2000 genes with highest variance ####
cg <- names(tail(sort(apply(log(av+1), 1, sd)), 2000))
log_av <- log(av[cg, ] + 1)

## Sample annotation ####
samples <- c("Normal.Y01.Cortex", "Normal.Y01.Medulla", "Normal.Y02.Cortex",
             "Normal.Y03.Cortex", "Normal.Y03.Medulla", 
             "OEM.Y01.Cortex","OEM.Y01.Medulla", "OEM.Y02.Cortex",
             "OEM.Y02.Medulla", "OEM.Y03.Cortex", "OEM.Y03.Medulla")

class_info <- data.frame(Class = ifelse(grepl("Normal", samples), "Normal", "OEM"))
row.names(class_info) <- samples

annotation_colors <- list(
  Class = c("Normal" = '#476D87', "OEM" = '#E95C59')
)

## Test normality and choose correlation method ####
normality_pvals <- apply(log_av, 2, function(x) shapiro.test(x)$p.value)
cor_method <- if (all(normality_pvals > 0.05)) "pearson" else "spearman"

## Calculate correlation matrix ####
df <- cor(as.matrix(log_av), method = cor_method)

## Define color palette ####
color_palette <- circlize::colorRamp2(
  breaks = c(0, 0.5, 1),
  colors = c("#67A9CF", "white","#B2182B")
)

## Generate correlation heatmap ####
p <- ComplexHeatmap::pheatmap(
  df,
  show_colnames = F, 
  show_rownames = T,
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = class_info,
  annotation_colors = annotation_colors,
  color = color_palette,
  row_split = class_info$Class,
  name = "Cor",
  column_split = class_info$Class
)

pdf(paste0('./Result/01.L1/02.Ann/cor_sample(', cor_method, '_euclidean_complete).pdf'),
    height = 5.5, width = 8)
p
dev.off()

## Generate dendrogram with ggtree ####
distance_matrix <- as.dist(1 - df)
hc <- hclust(distance_matrix, method = "ward.D2")
tree <- as.phylo(hc)

# Group samples
group <- cutree(hc, k = 11)
group_df <- data.frame(label = names(group), group = as.factor(group))
mycol <- c('#E5D2DD', '#F1BB72', '#53A85F', '#F3B1A0', '#D6E7A3', '#57C3F3', 
           '#E59CC4', '#6778AE', '#3A6963', '#BD956A', '#8C549C')

p <- ggtree(tree, layout = "rectangular", aes(color = group)) %<+% group_df +
  geom_tree(size = 1) +
  geom_tiplab(aes(color = group), size = 4, align = TRUE, vjust = 0.5, hjust = 0) +
  ggtitle("Clustering Dendrogram (Spearman Similarity)") +
  theme_tree() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = mycol) +
  NoLegend() +
  coord_cartesian(xlim = c(0, max(tree$edge.length) * 2.4))

pdf(paste0('./Result/01.L1/02.Ann/cor_sample(', cor_method, '_euclidean_ward.D2_ggtree).pdf'), 
    width = 6, height = 5)
p
dev.off()


# ==============================================================================
# Section 9: Cell Perturbation Analysis
# ==============================================================================

## 9.1 Augur Analysis ####
library(Seurat)
library(Augur)
library(ggplot2)
set.seed(123)

sce <- readRDS("./Data/02.Ann/annotation(L1).rds")

## Calculate AUC scores with Augur ####
augur <- calculate_auc(
  sce,
  cell_type_col = "L1_Celltype",
  label_col = "Class",
  n_threads = 8
)

saveRDS(augur, "./Result/01.L1/04.Perturbation/augur_result.rds")

## Create lollipop plot ####
aucs <- augur$AUC
range <- range(aucs$auc)
expand <- abs(diff(range)) * 0.1

p <- aucs %>% 
  ggplot(aes(x = reorder(cell_type, auc), y = auc)) + 
  geom_hline(aes(yintercept = 1), linetype = "dotted", size = 1) + 
  geom_text(aes(label = format(auc, digits = 4), 
                y = ifelse(auc < 0.5, 0.5, auc)), 
            size = 3.5, nudge_y = expand, hjust = 0.5) + 
  geom_segment(aes(xend = cell_type, yend = 0.5), size = 1) +
  geom_point(size = 3.5, aes(color = cell_type)) +
  scale_y_continuous("AUC Score", 
                     limits = c(min(range[1] - expand, 0.5), range[2] + expand * 1.5)) + 
  scale_color_manual(values = col) + 
  coord_flip() + 
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(title = "Cell Type AUC Score by Augur",
       x = "Cell Type", y = "AUC Score") +
  NoLegend()
ggsave("./Result/01.L1/04.Perturbation/Augur_lollipop_plot.pdf", p, 
       width = 7, height = 4.5)


## 9.2 scDist Analysis ####
library(dplyr)
library(Seurat)
library(scDist)
library(ggplot2)

## Prepare metadata for scDist (run on server) ####
sce <- readRDS("/data/liuruquan/Project/01.OEM/Analysis_V2/Data/01.Ann/annotation(L1).rds")

mtd <- sce@meta.data
level_mapping <- setNames(seq_along(levels(mtd$Sample)), levels(mtd$Sample))
mtd$batch <- level_mapping[mtd$Sample]
level_mapping <- setNames(seq_along(levels(mtd$Class)), levels(mtd$Class))
mtd$group <- level_mapping[mtd$Class]
sce@meta.data <- mtd 

sim <- list(
  Y = sce@assays$RNA@counts %>% as.data.frame(),
  meta.data = sce@meta.data %>% as.data.frame()
)

## Run scDist ####
out <- scDist(
  sim$Y,
  sim$meta.data,
  fixed.effects = "group",
  random.effects = "batch",
  clusters = "L1_Celltype"
)

saveRDS(out, "/data/liuruquan/Project/01.OEM/Analysis_V2/Result/01.L1_perturbation/scDist/scDist_out.rds")

## Visualize scDist results ####
sce <- readRDS("./Data/02.Ann/annotation(L1).rds")
out <- readRDS("./Result/01.L1/04.Perturbation/scDist_out.rds")

results <- out$results
results <- results[order(results$Dist., decreasing = FALSE), ]
results$Significance <- ifelse(results$p.val < 0.001, "***",
                               ifelse(results$p.val < 0.01, "**",
                                      ifelse(results$p.val < 0.05, "*", "ns")))
results$Cell_Type <- row.names(results)

## Create bar plot ####
p <- ggplot(results, aes(x = reorder(Cell_Type, Dist.), y = Dist., fill = Cell_Type)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.9) +
  geom_errorbar(aes(ymin = `95% CI (low)`, ymax = `95% CI (upper)`), 
                width = 0.2, color = "black") +
  scale_fill_manual(values = col) +
  geom_text(aes(label = Significance), hjust = -0.2, size = 4, color = "black") +
  coord_flip() +
  labs(title = "Cell Type Significance Distance by scDist",
       x = " ", y = "Distance") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
ggsave("./Result/01.L1/04.Perturbation/scDist_barplot2.pdf", p, 
       width = 5.5, height = 4)


## 9.3 Combined Augur and scDist Visualization ####
library(dplyr)
library(ggplot2)
library(tibble)
library(ggrepel)

augur <- readRDS("./Result/01.L1/04.Perturbation/augur_result.rds")
out   <- readRDS("./Result/01.L1/04.Perturbation/scDist_out.rds")
sce   <- readRDS("./Data/02.Ann/annotation(L1).rds")

## Extract and merge metrics ####
df_auc <- augur$AUC %>% dplyr::select(cell_type, auc)

df_dist <- out$results %>%
  tibble::rownames_to_column("cell_type") %>%
  dplyr::mutate(
    Significance = ifelse(p.val < 0.001, "***",
                          ifelse(p.val < 0.01,  "**",
                                 ifelse(p.val < 0.05,  "*", "ns")))
  ) %>%
  dplyr::select(cell_type, Dist. = Dist., p.val, Significance)

df_n <- sce@meta.data %>%
  count(L1_Celltype, name = "n_cells") %>%
  rename(cell_type = L1_Celltype)

df_plot <- df_auc %>%
  inner_join(df_dist, by = "cell_type") %>%
  left_join(df_n, by = "cell_type")

df_plot$cell_type <- factor(df_plot$cell_type, levels = rev(levels(sce$L1_Celltype)))

## Create scatter plot comparing Augur AUC and scDist ####
p <- ggplot(df_plot, aes(x = auc, y = Dist.)) +
  geom_point(aes(size = n_cells, color = cell_type), alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = Significance),
    size = 7,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = col) +
  scale_size(range = c(2.5, 12)) +
  labs(x = "Augur AUC score",
       y = "scDist Distance",
       size = "Cell count",
       color = "Cell type") +
  theme_classic(base_size = 18)

ggsave("./Result/01.L1/04.Perturbation/AUC_vs_scDist_scatter.pdf", p, 
       width = 7.5, height = 5)