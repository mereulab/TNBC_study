library(Seurat)
library(tidyverse)  
library(writexl)
library(harmony)
library(RColorBrewer)

step <-  "02_integration_harmony"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/mouse/"

saveProcessedData <- file.path(projectRoot, "processed_data", "integration", step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", "integration", step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", "integration", step); dir.create(saveResults, recursive = T)


data_path_Pal <- file.path(projectRoot, "processed_data", "Pal", "03_norm_clust_annot", "data.RDS")
data_path_Wuidart <- file.path(projectRoot, "processed_data", "Wuidart", "02_norm_clust_annot", "data.RDS")
data_path_Centoze <- file.path(projectRoot, "processed_data", "Centoze", "01_norm_clust_annot", "data.RDS")

data_Pal <- readRDS(data_path_Pal)
data_Wuidart <- readRDS(data_path_Wuidart)
data_Centoze <- readRDS(data_path_Centoze)

data_Pal$dataset <- "Pal"
data_Wuidart$dataset <- "Wuidart"
data_Centoze$dataset <- "Centoze"

data <- merge(x = data_Pal, y = c(data_Wuidart, data_Centoze))
data <- JoinLayers(data)

data[["RNA"]] <- split(data[["RNA"]], f = data$dataset)
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
data <- FindNeighbors(data, reduction = "harmony", dims = 1:15)
data <- FindClusters(data, resolution = 0.1, cluster.name = "harmony_clusters")
data <- RunUMAP(data, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")




plot <- DimPlot(data, split.by = "dataset", group.by = "annotation", reduction = "umap.harmony", pt.size = 2)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_split_by_dataset.pdf"), width = 16, height = 8)

plot <- DimPlot(data, split.by = "annotation", group.by = "dataset", reduction = "umap.harmony", pt.size = 2)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_dataset_split_by_annotation.pdf"), width = 16, height = 8)

plot <- DimPlot(data, split.by = "dataset", group.by = "harmony_clusters", reduction = "umap.harmony", pt.size = 2)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_harmony_clusters_split_by_annotation.pdf"), width = 16, height = 8)


genes <- c("Notch1", "Prom1", "Krt5", "Acta2"); names(genes) <- genes

plot <- FeaturePlot(data, features = genes, min.cutoff = "q10", order = T, reduction = "umap.harmony")
ggsave(plot, filename = file.path(savePlots, "fp_population_markers.pdf"), width = 8, height = 8)

plots <- genes %>% purrr::imap(~FeaturePlot(data, features = .x, min.cutoff = "q10", order = T, pt.size = 1, reduction = "umap.harmony", cols = c("lightgray", "red")) + theme_void() + 
                                 theme(plot.title = element_text(hjust = 0.5, size = 40)) + 
                                 labs(title = .x) 
                                 #scale_color_viridis_c(option = "rocket")
                                 )
combined_plots <- cowplot::plot_grid(plotlist = plots %>% purrr::map(~.x + NoLegend()), ncol = 4)
ggsave(combined_plots, filename = file.path(savePlots, paste0("fp_population_markers_pub.pdf")), height = 8, width = 30)
ggsave(plots[[1]], filename = file.path(savePlots, paste0("legend.pdf")), height = 8, width = 8)



plot = DimPlot(data, reduction = "umap.harmony", group.by = "annotation")
data <- CellSelector(plot = plot, object = data, ident = 'bc')
data <- CellSelector(plot = plot, object = data, ident = 'lc_er_neg')
data <- CellSelector(plot = plot, object = data, ident = 'lc_er_pos')
data$annotation_refined = data@active.ident

plot <- DimPlot(data, split.by = "dataset", group.by = "annotation_refined", reduction = "umap.harmony", pt.size = 2)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_refined_split_by_dataset.pdf"), width = 16, height = 8)

plot <- DimPlot(data, split.by = "dataset", group.by = "annotation_refined", reduction = "umap.harmony", pt.size = 2)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_refined_split_by_dataset.pdf"), width = 16, height = 8)

data@meta.data <- data@meta.data  %>% mutate(annotation_refined_names_pub = case_when(annotation_refined == "lc_er_pos" ~ "ERα pos LC",
                                                                                      annotation_refined == "lc_er_neg" ~ "ERα neg LC",
                                                                                      annotation_refined == "bc" ~ "BC"))
plot <- DimPlot(data, group.by = "annotation_refined_names_pub", reduction = "umap.harmony", pt.size = 1, cols = c("ERα pos LC" = "#fccf03", "ERα neg LC" = "#ff9b21", "BC" = "#37ad7a")) & theme_void() &
    labs(title = "Annotation") &
  guides(color = guide_legend(override.aes = list(size = 8))) &
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent", color = NA),  # Transparent legend
        #panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)  ) 
  
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_annotation_refined.png"), bg = "transparent", width = 8, height = 6)

data@meta.data <- data@meta.data  %>% mutate(dataset_names_pub = case_when(dataset == "Pal" ~ "Pal et al.",
                                                                           dataset == "Wuidart" ~ "Wuidart et al.",
                                                                           dataset == "Centoze" ~ "Centoze et al."))
plot <- DimPlot(data, group.by = "dataset_names_pub", reduction = "umap.harmony", pt.size = 1, cols = c("Pal et al." = "#6bc9ff", "Wuidart et al." = "#c91e46", "Centoze et al." = "#497050")) + theme_void() + 
   theme(plot.title = element_text(hjust = 0.5, size = 25),
         legend.text = element_text(size = 16)) + 
   labs(title = "Dataset") +
   guides(color = guide_legend(override.aes = list(size = 8)))
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_dataset.png"), width = 8, height = 6)


saveRDS(data, file.path(saveProcessedData, "data.RDS"))

# Differential expression
data <- JoinLayers(data)
Idents(data) <- "annotation_refined"
data.markers <- FindAllMarkers(data, only.pos = TRUE)
saveRDS(data.markers, file.path(saveResults, "markers_annotation_refined.RDS"))

markers_split <- list()
markers_split <- data.markers %>%
  split(.$cluster) %>%
  map(~ .x %>% 
        arrange(p_val_adj, desc(avg_log2FC)))
# Save each cluster as a separate sheet in an Excel file
writexl::write_xlsx(markers_split, path = file.path(saveResults, paste0("markers_annotation_refined", ".xlsx")))
write_xlsx(data.markers, path = file.path(saveResults, "markers_annotation_refined_one_sheet.xlsx"))


# Luminal vs Basal
data@meta.data <- data@meta.data %>% mutate(basal_luminal_annotation = if_else(annotation_refined == "bc", "bc", "lc"))
plot <- DimPlot(data, split.by = "dataset", group.by = "basal_luminal_annotation", reduction = "umap.harmony", pt.size = 2)
ggsave(plot = plot, filename = file.path(savePlots, "UMAP_basal_luminal_annotation_split_by_annotation.pdf"), width = 16, height = 8)

# Differential expression
Idents(data) <- "basal_luminal_annotation"
DEGs_luminal_basal <- FindAllMarkers(data, only.pos = TRUE)
saveRDS(DEGs_luminal_basal, file.path(saveResults, "markers_basal_luminal_annotation.RDS"))

markers_split <- list()
markers_split <- DEGs_luminal_basal %>%
  split(.$cluster) %>%
  map(~ .x %>% 
        arrange(p_val_adj, desc(avg_log2FC)))
# Save each cluster as a separate sheet in an Excel file
writexl::write_xlsx(markers_split, path = file.path(saveResults, paste0("markers_basal_luminal_annotation", ".xlsx")))
write_xlsx(DEGs_luminal_basal, path = file.path(saveResults, "markers_basal_luminal_annotation_one_sheet.xlsx"))

#pct1 > 0.7
markers_split <- list()
markers_split <- DEGs_luminal_basal %>%
  split(.$cluster) %>%
  map(~ .x %>% 
        arrange(p_val_adj, desc(avg_log2FC)) %>% 
        filter(pct.1> 0.7))
# Save each cluster as a separate sheet in an Excel file
writexl::write_xlsx(markers_split, path = file.path(saveResults, paste0("markers_basal_luminal_annotation_pct1_0.7", ".xlsx")))
write_xlsx(DEGs_luminal_basal, path = file.path(saveResults, "markers_basal_luminal_annotation_one_sheet_pct1_0.7.xlsx"))

#pct1 > 0.9
markers_split <- list()
markers_split <- DEGs_luminal_basal %>%
  split(.$cluster) %>%
  map(~ .x %>% 
        arrange(p_val_adj, desc(avg_log2FC)) %>% 
        filter(pct.1> 0.9 & pct.2 <0.5))
# Save each cluster as a separate sheet in an Excel file
writexl::write_xlsx(markers_split, path = file.path(saveResults, paste0("markers_basal_luminal_annotation_pct1_0.9_pct2_0.5", ".xlsx")))
write_xlsx(DEGs_luminal_basal, path = file.path(saveResults, "markers_basal_luminal_annotation_one_sheet_pct1_0.9_pct2_0.5.xlsx"))



# Checking initial results
DEGs_luminal_basal %>% filter(p_val_adj < 0.01) %>% filter(avg_log2FC > 0.25) %>% filter(cluster == "bc") %>% nrow()
# 896, before = 219
#old markers
DEGs_basal_luminal <- readRDS("~/Documents/Labs/Cellular_systems_genomics/Projects/Veronica Rodilla/integration_pal_centoze_wuidart/14.07.23/results/DEGs_basal_luminal.RDS")
old <- DEGs_basal_luminal %>% filter(p_val_adj < 0.01) %>% filter(avg_log2FC > 0) %>% pull(gene) # old markers!!
new <- DEGs_luminal_basal %>% filter(p_val_adj < 0.01) %>% filter(avg_log2FC > 0) %>% filter(cluster == "bc") %>% pull(gene)

# all of them are there
intersect(new, old) %>% length()
#[1] 208
setdiff(new, old) %>% length()
#[1] 714
# NEW 12 MORE GENES SPECIFIC TO BASAL

top_20_genes_basal_manuscript <- c("Acta2", "Tagln", "Cxcl14", "Apoe", "Tpm2",
                                   "Krt14", "Sparc", "Cald1", "Myl9", "Krt17",
                                   "Cnn1", "Il17b", "1600014C10Rik", "1500015O10Rik", "Mylk",
                                   "Rbp1", "Myh11", "Igfbp2", "Col4a2", "Apoc1")

top_20_genes_luminal_manuscript <- c("Spp1", "Csn1s2a", "Csn3", "Trf", "Csn2", 
                                     "Wfdc18", "Lalba", "Ly6d", "Csn1s1", "Plet1", 
                                     "Krt18", "Prlr", "Gpx3", "Areg", "Krt19",
                                     "Wfdc2", "Stc2", "Krt8", "Ptn", "Tmed3"
)

intersect(c(top_20_genes_luminal_manuscript, top_20_genes_basal_manuscript), DEGs_luminal_basal$gene) %>% length()
setdiff(c(top_20_genes_luminal_manuscript, top_20_genes_basal_manuscript), DEGs_luminal_basal$gene) 

plot <- FeaturePlot(data, features = top_20_genes_basal_manuscript, reduction = "umap.harmony", ncol = 5)
ggsave(plot = plot, filename = file.path(savePlots, "fp_top_20_genes_basal_manuscript.pdf"), width = 18, height = 16)

plot <- FeaturePlot(data, features = top_20_genes_luminal_manuscript, reduction = "umap.harmony", ncol = 5)
ggsave(plot = plot, filename = file.path(savePlots, "fp_top_20_genes_luminal_manuscript.pdf"), width = 18, height = 16)

data@meta.data <- data@meta.data %>% mutate(basal_luminal_annotation = if_else(annotation_refined == "bc", "BC", "LC"))
plot <- DimPlot(data, reduction = "pca", group.by = "basal_luminal_annotation", cols = c("BC" = "#37ad7a", "LC" = "#FFAD28")) 
ggsave(plot = plot, filename = file.path(savePlots, "pca_plot.pdf"), width = 6, height = 6)


pca_heatmap_plot <- DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE, nfeatures = 100, fast = F)
ggsave(plot = pca_heatmap_plot, filename = file.path(savePlots, "pca_heatmap_plot.png"), width = 15, height = 12)

pca_heatmap_plot <- DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE, nfeatures = 40, fast = F)
ggsave(plot = pca_heatmap_plot, filename = file.path(savePlots, "pca_heatmap_plot_20genes.png"), width = 15, height = 12)




# Genes up in basal

top_10_genes_basal_manuscript <- c("Acta2", "Tagln", "Cxcl14", "Apoe", "Tpm2",
                                   "Krt14", "Sparc", "Cald1", "Myl9", "Krt17")

plots <- top_10_genes_basal_manuscript %>% 
  setNames(top_10_genes_basal_manuscript) %>% 
  purrr::map(~FeaturePlot(data, features = paste0(.x), min.cutoff = "q10", order = T, reduction = "umap.harmony", pt.size = 1) +
               theme_void() + 
               theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
               labs(title = paste0(.x)) +
               scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")) ))

combined_plots <- cowplot::plot_grid(plotlist = plots %>% purrr::imap(~.x + NoLegend()), ncol = 5)
ggsave(combined_plots, filename = file.path(savePlots, paste0("fp_top_10_genes_basal_pub.png")), height = 12, width = 30)


"#ff9b21", "BC" = "#37ad7a")) & theme_void() &



