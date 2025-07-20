library(Seurat)
library(tidyverse)  
library(writexl)

dataset <- "Pal"
step <-  "01_QC_metrics"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/mouse/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

raw_data_path <- file.path(projectRoot, "raw_data", dataset, "GSM4994966_Adult-FVB")

# Create Seurat object:
list.files(raw_data_path) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = raw_data_path)

data <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, project = "Pal")

data
# An object of class Seurat 
# 17543 features across 11936 samples within 1 assay 
# Active assay: RNA (17543 features, 0 variable features)
# 1 layer present: counts

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     39    1105    1574    1660    2134    5799 


# QC METRICS
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
# Visualize QC metrics before filtering
# number of cells before filtering: 384
plot <-  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = plot, filename = file.path(savePlots, "QC_before_filtering.pdf"), width = 8, height = 8)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(plot = plot1 + plot2, filename = file.path(savePlots, "QC_scatter_plot.pdf"), width = 12, height = 6)

pdf(file = file.path(savePlots, "n_feature_density.pdf"))
plot(density(data$nFeature_RNA))
abline(v = 500)
dev.off()

# Visualize the number UMIs/transcripts per cell
plot <- data@meta.data %>% 
  ggplot(aes(x=nFeature_RNA)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave(plot = plot, filename = file.path(savePlots, "n_feature_density_2.pdf"), width = 12, height = 6)


# Filtering
data <- subset(data, subset = nFeature_RNA > 500 & percent.mt < 20)

# Visualize QC metrics after filtering
# number of cells before filtering: 280
plot <-  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = plot, filename = file.path(savePlots, "QC_after_filtering.pdf"), width = 8, height = 8)

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 501    1214    1643    1763    2193    5799 


saveRDS(data, file.path(saveProcessedData, "data.RDS"))




