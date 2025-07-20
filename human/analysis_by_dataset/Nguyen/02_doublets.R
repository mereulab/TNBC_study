library(Seurat)
library(tidyverse)  
library(writexl)
library(scDblFinder)
library(DoubletFinder)
library(BiocParallel)

dataset <- "Nguyen"
step <-  "02_doublets"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/human/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

data_path <- file.path(projectRoot, "processed_data", dataset, "01_QC_metrics", "data.RDS")


data_all_samples <- readRDS(data_path)
samples <- unique(data_all_samples$sample)

se_objs <- list()

for(sample in samples){
  #···············································································
  # DOUBLETFINDER ································································
  #···············································································
  print(sample)
  Idents(data_all_samples) <- "sample"
  data <- subset(data_all_samples, idents = sample)
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data <- ScaleData(data)
  data <- RunPCA(data)
  ElbowPlot(data)
  data <- RunUMAP(data, dims = 1:15)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep(data, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  nExp_poi <- round(0.075*nrow(data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  data <- doubletFinder(data, PCs = 1:15, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  
  #···············································································
  # SCDBLFINDER ··································································
  #···············································································
  sce <- as.SingleCellExperiment(data)
  sce <- scDblFinder(sce, BPPARAM=MulticoreParam(3))
  data2 <- as.Seurat(sce, counts = "counts", data = "logcounts")
  
  data <- data %>% AddMetaData(data2@meta.data %>% select(starts_with("scDblFinder")))
  
  data@meta.data <- data@meta.data %>% rename_with(~"DF.classifications", starts_with("DF.classifications"))
  
  data@meta.data <- data@meta.data %>% mutate(common_doublets = if_else(DF.classifications == "Doublet" & scDblFinder.class == "doublet", "Doublet", "Singlet"))
  
  se_objs[[sample]] <- data
}

data <- merge(x = se_objs[[1]], y = se_objs[2:length(se_objs)])
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data)
ElbowPlot(data)
data <- RunUMAP(data, dims = 1:15)

plot <- DimPlot(data, group.by = "common_doublets", cols = c("Doublet" = "red", "Singlet" = "lightblue"))
ggsave(plot, filename = file.path(savePlots, "UMAP_common_doublets.pdf"), height = 8, width = 8)

saveRDS(data, file.path(saveProcessedData, "data_with_doublets.RDS"))

data <- subset(data, subset = common_doublets == "Singlet")

saveRDS(data, file.path(saveProcessedData, "data.RDS"))

