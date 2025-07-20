library(Seurat)
library(tidyverse)  
library(writexl)


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

dataset <- "Gray"
step <-  "03_norm_clust_annot"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/human/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

data_path <- file.path(projectRoot, "processed_data", dataset, "02_doublets", "data.RDS")


data <- readRDS(file.path(data_path))

plot <- DimPlot(data, group.by = "Major.subtype", label = T)
ggsave(plot, filename = file.path(savePlots, "UMAP_Major.subtype.pdf"))

plot <- DimPlot(data, group.by = "Cell.subtype", label = T)
ggsave(plot, filename = file.path(savePlots, "UMAP_Cell.subtype.pdf"))


Idents(data) <- "Cell.subtype"
data <- subset(data, idents = c("F_doublets", "BA_doublets", "HS_doublets"), invert = T) 

plot <- DimPlot(data, group.by = "Major.subtype", label = T)
ggsave(plot, filename = file.path(savePlots, "UMAP_Major.subtype_without_doublets.pdf"))

plot <- DimPlot(data, group.by = "Cell.subtype", label = T)
ggsave(plot, filename = file.path(savePlots, "UMAP_Cell.subtype_without_doublets.pdf"))

data@meta.data <- data@meta.data %>% mutate(annotation = case_when(Major.subtype == "BA" ~ "Basal", 
                                                 Major.subtype == "AV" ~ "ER-neg", 
                                                 Major.subtype == "HS" ~ "ER-pos", 
                                                 .default = Major.subtype
                                                 ))

plot <- DimPlot(data, group.by = "annotation", label = T, cols = cols)
ggsave(plot, filename = file.path(savePlots, "UMAP_annotation_without_doublets.pdf"))
plot <- DimPlot(data, group.by = "annotation", label = F, cols = cols)
ggsave(plot, filename = file.path(savePlots, "UMAP_annotation_without_doublets_no_label.pdf"), width = 7, height = 6.5)

genes <- c("NOTCH1", "PROM1", "KRT5", "ACTA2"); names(genes) <- genes

plot <- FeaturePlot(data, features = genes, min.cutoff = "q10", order = T)
ggsave(plot, filename = file.path(savePlots, "fp_population_markers.pdf"), width = 8, height = 8)


#saveRDS(data, file.path(saveProcessedData, "data.RDS"))
data <- readRDS(file.path(saveProcessedData, "data.RDS"))


