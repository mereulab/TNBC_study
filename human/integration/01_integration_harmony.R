library(Seurat)
library(tidyverse)  
library(writexl)
library(harmony)
library(UCell)
library(RColorBrewer)

step <-  "01_integration_harmony"

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

saveProcessedData <- file.path(projectRoot, "processed_data", "integration", step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", "integration", step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", "integration", step); dir.create(saveResults, recursive = T)


data_path_Pal <- file.path(projectRoot, "processed_data", "Pal", "03_norm_clust_annot", "data.RDS")
data_path_Gray <- file.path(projectRoot, "processed_data", "Gray", "03_norm_clust_annot", "data.RDS")
data_path_Nguyen <- file.path(projectRoot, "processed_data", "Nguyen", "03_norm_clust_annot", "data.RDS")

data_Pal <- readRDS(data_path_Pal)
data_Gray <- readRDS(data_path_Gray)
data_Nguyen <- readRDS(data_path_Nguyen)

data_Pal$dataset <- "Pal"
data_Gray$dataset <- "Gray"
data_Nguyen$dataset <- "Nguyen"

data <- merge(x = data_Pal, y = c(data_Gray, data_Nguyen))
data <- JoinLayers(data)


table(data$sample)
table(data$common_doublets)

data[["RNA"]] <- split(data[["RNA"]], f = data$sample)
data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)

data <- IntegrateLayers(
  object = data, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = FALSE)

ElbowPlot(data)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:20)
data <- FindClusters(data, resolution = 0.2, cluster.name = "harmony_clusters")
data <- RunUMAP(data, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")


plot <- DimPlot(data, split.by = "dataset", group.by = "annotation", reduction = "umap.harmony", pt.size = 2)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_split_by_sample.pdf"), width = 16, height = 8)

plot <- DimPlot(data, split.by = "annotation", group.by = "dataset", reduction = "umap.harmony", pt.size = 2)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_dataset_split_by_annotation.pdf"), width = 16, height = 8)

plot <- DimPlot(data, split.by = "dataset", group.by = "harmony_clusters", reduction = "umap.harmony", pt.size = 2, label = T)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_harmony_clusters_split_by_annotation.pdf"), width = 16, height = 8)

plot <- DimPlot(data, group.by = "harmony_clusters", reduction = "umap.harmony", pt.size = 1, label = T)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_harmony_clusters.pdf"), width = 8, height = 8)


genes <- c("NOTCH1", "PROM1", "KRT5", "ACTA2"); names(genes) <- genes

plot <- FeaturePlot(data, features = genes, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_population_markers.pdf"), width = 8, height = 8)


basal_markers_paper <- c("KRT5", "ACTA2", "MYLK", "SNAI2", "NOTCH4", "DKK3")
er_pos_markers_paper <- c("ESR1", "PGR", "FOXA1")
er_neg_markers_paper <- c("TNFRSF11A", "KIT", "SOX10")

plot <- FeaturePlot(data, features = basal_markers_paper, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_basal_markers_papers.pdf"))

plot <- FeaturePlot(data, features = er_pos_markers_paper, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_er_pos_markers_paper.pdf"))

plot <- FeaturePlot(data, features = er_neg_markers_paper, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_er_neg_markers_paper.pdf"))


plot <- FeaturePlot(data, features = "MKI67", min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_MKI67.pdf"))


cell_type_markers_Gray <- readxl::read_excel(file.path(projectRoot, "misc", "NIHMS1807635-supplement-3_Gray.xlsx"), sheet = 1,)
cell_types <- as.vector(cell_type_markers_Gray[2, ]) %>% unlist() %>% unname()
cell_type_markers_Gray <- cell_type_markers_Gray[-c(1:2),]
cell_type_markers_Gray <- as.data.frame(cell_type_markers_Gray)
colnames(cell_type_markers_Gray) <- cell_types

gene_sig <- cell_types %>%
  setNames(cell_types) %>%
  purrr::map(~cell_type_markers_Gray %>%
               select(.x) %>%
               pull(.x))

gene_sig <- gene_sig %>%
  purrr::map(~ .x[!is.na(.x)])  


data <- data %>% UCell::AddModuleScore_UCell(features = gene_sig)

plots <- cell_types %>% 
  setNames(cell_types) %>%
  purrr::map(~FeaturePlot(data, features = paste0(.x, "_UCell"), min.cutoff = "q10", order = T, reduction = "umap.harmony") +
               theme_void() + 
               theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
               labs(title = paste0(.x, " signature")) +
               scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")) ))
               #scale_color_viridis_c(option = "rocket"))

MKI67_plot <- FeaturePlot(data, features = "MKI67", min.cutoff = "q10", order = T, reduction = "umap.harmony") +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
  labs(title = "MKI67") +
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")) ) 

all_plots <- c(plots, list("MKI67" = MKI67_plot))
names(all_plots)
#[1] "AV"                   "HS"                   "BA"                   "Fibroblast"           "Vascular & Lymphatic"
#[6] "Immune"               "MKI67"
names(all_plots) <- c("ERα neg LC","ERα pos LC","BC","Fibroblast","Vascular & Lymphatic", "Immune","MKI67")


all_plots %>% purrr::imap(~ggsave(.x, filename = file.path(savePlots, paste0(.y, ".pdf"))))
combined_plots <- cowplot::plot_grid(plotlist = all_plots %>% purrr::imap(~.x + NoLegend() + labs(title = paste0(.y, " signature"))), ncol = 4)
ggsave(combined_plots, filename = file.path(savePlots, paste0("fp_signatures_ct_pub.png")), height = 12, width = 26)
ggsave(plots[[1]], filename = file.path(savePlots, paste0("legend.pdf")), height = 8, width = 8)


data@meta.data <- data@meta.data %>% 
  mutate(annotation_refined = case_when(harmony_clusters == "0" ~ "ER-neg",
                                        harmony_clusters == "1" ~ "ER-pos",
                                        harmony_clusters == "2" ~ "Basal",
                                        harmony_clusters == "3" ~ "Fibroblasts",
                                        harmony_clusters == "4" ~ "VasLymph",
                                        harmony_clusters == "5" ~ "ER-neg",
                                        harmony_clusters == "6" ~ "Prolif. cells",
                                        harmony_clusters == "7" ~ "Immune",
                                        harmony_clusters == "8" ~ "Basal"
                                        ))



plot <- DimPlot(data, group.by = "annotation_refined", reduction = "umap.harmony", pt.size = 0.3, label = T)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_refined.pdf"), width = 8, height = 8)


genes_paper <- c("ACTA2", "TAGLN", "CXCL14", "APOE", "TPM2", 
                 "KRT14", "SPARC", "CALD1", "MYL9", "KRT17")

plot <- FeaturePlot(data, features = genes_paper, min.cutoff = "q10", order = T, reduction = "umap.harmony", ncol = 5)
ggsave(plot, filename = file.path(savePlots, "fp_genes_paper.pdf"), width = 16, height = 8)


plots <- genes_paper %>% 
  setNames(genes_paper) %>%
  purrr::map(~FeaturePlot(data, features = paste0(.x), min.cutoff = "q10", order = T, reduction = "umap.harmony") +
               theme_void() + 
               theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
               labs(title = paste0(.x)) +
               scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")) ))
combined_plots <- cowplot::plot_grid(plotlist = plots %>% purrr::imap(~.x + NoLegend()), ncol = 5)
ggsave(combined_plots, filename = file.path(savePlots, paste0("fp_genes_paper_pub.png")), height = 12, width = 30)
ggsave(plots[[1]], filename = file.path(savePlots, paste0("legend.pdf")), height = 8, width = 8)



# Differential expression
data <- JoinLayers(data)
Idents(data) <- "annotation_refined"
DEGs_annotation_refined <- FindAllMarkers(data, only.pos = TRUE)
saveRDS(DEGs_annotation_refined, file.path(saveResults, "markers_annotation_refined.RDS"))

markers_split <- list()
markers_split <- DEGs_annotation_refined %>%
  split(.$cluster) %>%
  map(~ .x %>% 
        arrange(p_val_adj, desc(avg_log2FC)))
# Save each cluster as a separate sheet in an Excel file
writexl::write_xlsx(markers_split, path = file.path(saveResults, paste0("markers_annotation_refined", ".xlsx")))
write_xlsx(DEGs_annotation_refined, path = file.path(saveResults, "markers_annotation_refined_one_sheet.xlsx"))


data$dataset <- factor(data$dataset, levels = c("Gray", "Pal", "Nguyen"))
plot <- DimPlot(data, group.by = "annotation_refined", reduction = "umap.harmony", pt.size = 0.1, label = F, cols = cols, split.by = "dataset")
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_refined_split_dataset.pdf"), width = 14.5, height = 6.5)

#saveRDS(data, file.path(saveProcessedData, "data.RDS"))
data <- readRDS(file.path(saveProcessedData, "data.RDS"))


