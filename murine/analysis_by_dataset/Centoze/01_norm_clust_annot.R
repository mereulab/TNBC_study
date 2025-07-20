library(Seurat)
library(tidyverse)  
library(writexl)

dataset <- "Centoze"
step <-  "01_norm_clust_annot"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/mouse/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

raw_data_path <- file.path(projectRoot, "raw_data", dataset, "GSE148791_Centonze_et_al_Seurat.rds")


data <- readRDS(file.path(raw_data_path))
data <- Seurat::UpdateSeuratObject(data)

table(data$cell_type)
# BC_TOM       BC_WT      Hy_TOM      LC_TOM       LC_WT neg_control 
# 71          47          89          82          48           0 
data <- subset(data, subset = cell_type %in% c("BC_WT", "LC_WT"))


plot <- DimPlot(data, group.by = "sc3_8_clusters")
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_sc3_8_clusters.pdf"), width = 8, height = 8)

plot <- DimPlot(data, group.by = "cell_type")
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_cell_type.pdf"), width = 8, height = 8)

genes <- c("Notch1", "Prom1", "Krt5", "Acta2"); names(genes) <- genes

plot <- FeaturePlot(data, features = genes, min.cutoff = "q10", order = T)
ggsave(plot, filename = file.path(savePlots, "fp_population_markers.pdf"), width = 8, height = 8)

Idents(data) <- "sc3_8_clusters"
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, min.diff.pct = 0.1)
write_xlsx(data.markers, path = file.path(saveResults, "markers_sc3_8_clusters.xlsx"))
saveRDS(data.markers, file.path(saveResults, "markers_sc3_8_clusters.RDS"))


top10 <- data.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
plot <- DoHeatmap(data, features = top10$gene) + NoLegend()
ggsave(plot, filename = file.path(savePlots, "heatmap_sc3_8_clusters.png"))


data@meta.data <- data@meta.data %>% mutate(annotation = case_when(cell_type == "LC_WT" & sc3_8_clusters == "6" ~ "lc_neg",
                                                                   cell_type == "LC_WT" & sc3_8_clusters == "7" ~ "lc_pos",
                                                                   cell_type == "BC_WT" & sc3_8_clusters == "5" ~ "bc", .default = "other"
                                                                    ))

data <- subset(data, subset = annotation %in% c("lc_neg", "lc_pos", "bc"))

plot <- FeaturePlot(data, features = genes, min.cutoff = "q10", order = T)
ggsave(plot, filename = file.path(savePlots, "fp_population_markers_annotation.pdf"), width = 8, height = 8)

plot <- DimPlot(data, group.by = "annotation")
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation.pdf"), width = 8, height = 8)



Idents(data) <- "annotation"
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, min.diff.pct = 0.1)
write_xlsx(data.markers, path = file.path(saveResults, "markers_annotation.xlsx"))
saveRDS(data.markers, file.path(saveResults, "markers_annotation.RDS"))

saveRDS(data, file.path(saveProcessedData, "data.RDS"))




