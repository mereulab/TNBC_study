library(Seurat)
library(tidyverse)  
library(writexl)

dataset <- "Nguyen"
step <-  "01_QC_metrics"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/human/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

raw_data_path <- file.path(projectRoot, "raw_data", dataset)

# Load data
counts_Ind4 <- read.delim(file.path(raw_data_path, "GSM3099846_Ind4_Expression_Matrix.txt.gz"), row.names=1)
counts_Ind5 <- read.delim(file.path(raw_data_path, "GSM3099847_Ind5_Expression_Matrix.txt.gz"), row.names=1)
counts_Ind6 <- read.delim(file.path(raw_data_path, "GSM3099848_Ind6_Expression_Matrix.txt.gz"), row.names=1)
counts_Ind7 <- read.delim(file.path(raw_data_path, "GSM3099849_Ind7_Expression_Matrix.txt.gz"), row.names=1)


Ind4 <- CreateSeuratObject(counts = counts_Ind4, min.cells=3)
Ind5 <- CreateSeuratObject(counts = counts_Ind5, min.cells=3)
Ind6 <- CreateSeuratObject(counts = counts_Ind6, min.cells=3)
Ind7 <- CreateSeuratObject(counts = counts_Ind7, min.cells=3)


# Add sample metadata
Ind4$sample <- "Ind4"
Ind5$sample <- "Ind5"
Ind6$sample <- "Ind6"
Ind7$sample <- "Ind7"

data = merge(Ind4, y = c(Ind5, Ind6, Ind7))


data
# An object of class Seurat 
# 21842 features across 24646 samples within 1 assay 
# Active assay: RNA (21842 features, 0 variable features)
# 4 layers present: counts.1, counts.2, counts.3, counts.4

data <- JoinLayers(data)

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  432    1930    2546    2710    3400    7198 


# QC METRICS
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
# Visualize QC metrics before filtering
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
data <- subset(data, subset = nFeature_RNA > 500 & percent.mt < 10)

# Visualize QC metrics after filtering
# number of cells before filtering: 24520
plot <-  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = plot, filename = file.path(savePlots, "QC_after_filtering.pdf"), width = 8, height = 8)

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 532    1937    2554    2717    3404    7198 


saveRDS(data, file.path(saveProcessedData, "data.RDS"))
