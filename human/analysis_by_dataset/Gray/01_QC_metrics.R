library(Seurat)
library(tidyverse)  
library(writexl)

dataset <- "Gray"
step <-  "01_QC_metrics"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/human/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

raw_data_path <- file.path(projectRoot, "raw_data", dataset)

# Load data
counts <- read.csv(file.path(raw_data_path, "GSE180878_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_human.csv.gz"), header = TRUE)
rownames(counts) <- counts$GeneSymbol
counts$GeneSymbol <- NULL
data <- CreateSeuratObject(counts = counts, min.cells=3)

metadata <- read.csv(file.path(raw_data_path, "GSE180878_Li_Brugge_10XscRNAseq_Metadata_human.csv.gz"), header = TRUE)
metadata <- metadata %>% mutate(cellID = str_replace(cellID, "-", "."))
rownames(metadata) <- metadata$cellID
metadata$cellID <- NULL
data <- data %>% AddMetaData(metadata)

#' Extract noncarriers patient samples (1-s2.0-S1534580722003318-mmc2.xlsx)
#' RM-A
#' RM-B
#' RM-C
data$sample <- data$orig.ident
Idents(data) <- "sample"
data <- subset(x = data, idents = c("RM.A", "RM.B", "RM.C"))
data$sample <- factor(data$sample)
table(data$sample)

data
# An object of class Seurat 
# 20437 features across 384 samples within 1 assay 
# Active assay: RNA (20437 features, 0 variable features)
# layer present: counts

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     500    1483    2016    2068    2582    7101 


# QC METRICS
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
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

print(summary(data$nFeature_RNA))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#500    1483    2016    2068    2582    7101 

#-----> Already filtered !!!!!!!!!!!!!

# Visualize the number UMIs/transcripts per cell
plot <- data@meta.data %>% 
  ggplot(aes(x=nFeature_RNA)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 2500)
ggsave(plot = plot, filename = file.path(savePlots, "n_feature_density_2.pdf"), width = 12, height = 6)



saveRDS(data, file.path(saveProcessedData, "data.RDS"))

