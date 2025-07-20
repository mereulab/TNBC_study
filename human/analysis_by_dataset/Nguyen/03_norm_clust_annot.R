library(Seurat)
library(tidyverse)  
library(writexl)
options(future.globals.maxSize = 8000 * 1024^2)

dataset <- "Nguyen"
step <-  "03_norm_clust_annot"

cols = c("ERα pos LC" = "#fccf03", 
         "ER-pos" = "#fccf03", 
         "ERα neg LC" = "#ff9b21", 
         "ER-neg" = "#ff9b21", 
         "BC" = "#37ad7a", 
         "Basal" = "#37ad7a", 
         "Immune" = "#6baddb", 
         "Vascular & Lymphatic" = "#ba2f64", 
         "VasLymph" = "#ba2f64", 
         "Prolif. cells" = "#a070b5", 
         "Fibroblasts" = "#824a30",
         "Fibroblast" = "#824a30"
)

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/human/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

data_path <- file.path(projectRoot, "processed_data", dataset, "02_doublets", "data.RDS")


data <- readRDS(file.path(data_path))
data <- JoinLayers(data)

plot <- DimPlot(data, group.by = "sample")
ggsave(plot, filename = file.path(savePlots, "umap_sample_before_integration.pdf"))

data[["RNA"]] <- split(data[["RNA"]], f = data$sample)
data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)

data <- IntegrateLayers(
  object = data, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "rpca",
  verbose = FALSE)

ElbowPlot(data)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:15)
data <- FindClusters(data, resolution = 0.05, cluster.name = "harmony_clusters")
data <- RunUMAP(data, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")

plot <- DimPlot(data, group.by = "sample", reduction = "umap.harmony")
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_sample.pdf"))

plot <- DimPlot(data, group.by = "harmony_clusters", reduction = "umap.harmony")
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_harmony_clusters.pdf"))

plot <- DimPlot(data, split.by = "sample", group.by = "harmony_clusters", reduction = "umap.harmony", pt.size = 2, ncol = 4)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_sample_split_by_annotation.pdf"), width = 24, height = 15)


genes <- c("NOTCH1", "PROM1", "KRT5", "ACTA2"); names(genes) <- genes

plot <- FeaturePlot(data, features = genes, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_population_markers.pdf"), width = 8, height = 8)

plot <- FeaturePlot(data, features = "MKI67", min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_MKI67.pdf"), width = 8, height = 8)

# paper markers
basal_markers_paper <- c("KRT5", "ACTA2", "MYLK", "SNAI2", "NOTCH4", "DKK3")
er_pos_markers_paper <- c("ESR1", "PGR", "FOXA1")
er_neg_markers_paper <- c("TNFRSF11A", "KIT", "SOX10")

plot <- FeaturePlot(data, features = basal_markers_paper, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_basal_markers_papers.pdf"))

plot <- FeaturePlot(data, features = er_pos_markers_paper, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_er_pos_markers_paper.pdf"))

plot <- FeaturePlot(data, features = er_neg_markers_paper, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_er_neg_markers_paper.pdf"))


data <- JoinLayers(data)

data@meta.data <- data@meta.data %>% mutate(annotation = case_when(harmony_clusters == "0" ~ "Basal",
                                                                   harmony_clusters == "1" ~ "ER-neg",
                                                                   harmony_clusters == "2" ~ "ER-pos",
                                                                   harmony_clusters == "3" ~ "Basal"
))

plot <- DimPlot(data, group.by = "annotation", label = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "UMAP_annotation.pdf"), width = 7, height = 6.5)
plot <- DimPlot(data, group.by = "annotation", label = F, reduction = "umap.harmony", cols = cols)
ggsave(plot, filename = file.path(savePlots, "UMAP_annotation_no_label.pdf"), width = 7, height = 6.5)

#saveRDS(data, file.path(saveProcessedData, "data.RDS"))
data <- readRDS(file.path(saveProcessedData, "data.RDS"))

table(data$common_doublets)

