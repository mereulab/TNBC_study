library(Seurat)
library(tidyverse)  
library(writexl)

dataset <- "Pal"
step <-  "01_QC_metrics"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/human/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

raw_data_path <- file.path(projectRoot, "raw_data", dataset)

# Load data
samples <- c("N-0064-Epi", "N-0093-Epi", "N-0123-Epi", "N-0230.16-Epi",
             "N-0275-Epi", "N-0342-Epi", "N-0372-Epi", "N-0408-Epi",
             "N-1469-Epi", "N-N280-Epi", "N-N1105-Epi")
data.dirs <- samples %>% purrr::imap(~file.path(raw_data_path, .x)) %>% setNames(samples)
expression_matrixes <- data.dirs %>% purrr::imap(~Read10X(data.dir = .x))
se_objs <- expression_matrixes %>% purrr::imap(~CreateSeuratObject(counts = .x, min.cells = 3, project = .y))
se_objs <- se_objs %>% purrr::imap(~AddMetaData(
  object = .x,
  metadata = .y,
  col.name = 'sample'
))

# Merge
se_objs <- unname(se_objs)
data <- merge(x = se_objs[[1]], y = se_objs[2:length(se_objs)] %>% as.vector())


data
# An object of class Seurat 
# 2410 features across 58830 samples within 1 assay 
# Active assay: RNA (22410 features, 0 variable features)
# 11 layers present: counts.N-0064-Epi, counts.N-0093-Epi, counts.N-0123-Epi, counts.N-0230.16-Epi, counts.N-0275-Epi, counts.N-0342-Epi, counts.N-0372-Epi, counts.N-0408-Epi, counts.N-1469-Epi, counts.N-N280-Epi, counts.N-N1105-Epi
data <- JoinLayers(data)

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   21    1146    1778    1804    2346    7746 


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
abline(v = 400)
dev.off()

# Visualize the number UMIs/transcripts per cell
plot <- data@meta.data %>% 
  ggplot(aes(x=nFeature_RNA)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 400)
ggsave(plot = plot, filename = file.path(savePlots, "n_feature_density_2.pdf"), width = 12, height = 6)


# Filtering
data <- subset(data, subset = nFeature_RNA > 400 & percent.mt < 25)

# Visualize QC metrics after filtering
# number of cells before filtering: 280
plot <-  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = plot, filename = file.path(savePlots, "QC_after_filtering.pdf"), width = 8, height = 8)

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 401    1242    1828    1874    2380    7746 


saveRDS(data, file.path(saveProcessedData, "data.RDS"))
