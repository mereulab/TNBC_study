library(Seurat)
library(tidyverse)  
library(writexl)

dataset <- "Wuidart"
step <-  "02_norm_clust_annot"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/mouse/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

data_path <- file.path(projectRoot, "processed_data", dataset, "01_QC_metrics", "data.RDS")


data <- readRDS(file.path(data_path))

data <- NormalizeData(data)

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2509    3366    4496    4615    5564   11873 

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 4500)
data <- ScaleData(data, features = rownames(data))
data <- RunPCA(data, features = VariableFeatures(object = data))
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:12)
data <- FindClusters(data, resolution = 0.2)

plot <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = plot, filename = file.path(savePlots, "QC_after_clustering.pdf"), width = 8, height = 8)

data <- RunUMAP(data, reduction.use = "pca", dims = 1:12)

plot <- DimPlot(data)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_clusters.pdf"), width = 8, height = 8)

# Differential expression
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, min.diff.pct = 0.1)
write_xlsx(data.markers, path = file.path(saveResults, "markers_RNA_snn_res.0.2.xlsx"))
saveRDS(data.markers, file.path(saveResults, "markers_RNA_snn_res.0.2.RDS"))

genes <- c("Notch1", "Prom1", "Krt5", "Acta2"); names(genes) <- genes

plot <- FeaturePlot(data, features = genes, min.cutoff = "q10", order = T)
ggsave(plot, filename = file.path(savePlots, "fp_population_markers.pdf"), width = 8, height = 8)

top10 <- data.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
plot <- DoHeatmap(data, features = top10$gene) + NoLegend()
ggsave(plot, filename = file.path(savePlots, "heatmap_RNA_snn_res.0.2.png"))

# RNA_snn_res.0.2:
#···················
# Cluster 0 and 4 = EMP
# Cluster 1 = bc
# Cluster 2 = lc_pos
# Cluster 3 = lc_neg


data@meta.data <- data@meta.data %>% mutate(annotation = case_when(RNA_snn_res.0.2 %in% c("0", "4") ~ "EMP",
                                                                   RNA_snn_res.0.2 == "1" ~ "bc",
                                                                   RNA_snn_res.0.2 == "2" ~ "lc_pos",
                                                                   RNA_snn_res.0.2 == "3" ~ "lc_neg"))
plot <- DimPlot(data, group.by = "annotation")
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_with_EMP_cluster.pdf"), width = 8, height = 8)

saveRDS(data, file.path(saveProcessedData, "data_with_EMP_cluster.RDS"))


# Filter out EMP cluster
data = subset(data, subset = annotation == "EMP", invert = T)

plot <- DimPlot(data, group.by = "annotation")
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation.pdf"), width = 8, height = 8)

Idents(data) <- "annotation"
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, min.diff.pct = 0.1)
write_xlsx(data.markers, path = file.path(saveResults, "markers_annotation.xlsx"))
saveRDS(data.markers, file.path(saveResults, "markers_annotation.RDS"))

saveRDS(data, file.path(saveProcessedData, "data.RDS"))


