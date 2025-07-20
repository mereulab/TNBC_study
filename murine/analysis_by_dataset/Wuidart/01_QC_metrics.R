library(Seurat)
library(tidyverse)  
library(biomaRt)
library(writexl)

dataset <- "Wuidart"
step <-  "01_QC_metrics"

projectRoot <- "~/Documents/PROJECTS/04.2_BreastCancer_ID/mouse/"

saveProcessedData <- file.path(projectRoot, "processed_data", dataset, step); dir.create(saveProcessedData, recursive = T)
savePlots <- file.path(projectRoot, "plots", dataset, step); dir.create(savePlots, recursive = T)
saveResults <- file.path(projectRoot, "results", dataset, step); dir.create(saveResults, recursive = T)

raw_data_path <- file.path(projectRoot, "raw_data", dataset, "GSE110351_MG_EMP_LC_BC_raw_counts.csv.gz")


# CONVERT ENSMBL IDS TO SYMBOL IDS
#-------------------------------------------------------------------------------
raw_counts <- read_csv(raw_data_path) 
colnames(raw_counts)[1] <- "ensembl_gene_id"

raw_counts <- raw_counts %>% mutate(ensembl_gene_id_v2 = gsub("\\.\\d+$", "", ensembl_gene_id)) # Regex to remove the version

#Select the Ensembl database and specify the dataset for Mus musculus (mouse)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Query for gene names, removing version numbers from the IDs
genes <- raw_counts$ensembl_gene_id_v2  
length(unique(genes))
result <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                filters = "ensembl_gene_id", 
                values = genes, 
                mart = mart)

print(result)
length(unique(result$mgi_symbol_unique))
result <- result %>% mutate(mgi_symbol_unique = make.unique(as.character(result$mgi_symbol), sep = "_"))

raw_counts_gene_symbol <- raw_counts %>% inner_join(result, by = c("ensembl_gene_id_v2"="ensembl_gene_id"))
raw_counts_gene_symbol <- raw_counts_gene_symbol %>% as_tibble() %>% 
  dplyr::select(-c(ensembl_gene_id_v2, mgi_symbol, ensembl_gene_id)) %>% as.data.frame()

rownames(raw_counts_gene_symbol) <- raw_counts_gene_symbol$mgi_symbol_unique
raw_counts_gene_symbol$mgi_symbol_unique <-  NULL

write.csv(raw_counts_gene_symbol, file.path(projectRoot, "raw_data", dataset, "raw_counts_gene_symbol.csv"), row.names = TRUE)

print(dim(raw_counts))
print(dim(raw_counts_gene_symbol))

#-------------------------------------------------------------------------------

# Find rownames labels as "" and filter it
empty_gene_rows <- which(rownames(raw_counts_gene_symbol) == "" | is.na(rownames(raw_counts_gene_symbol)))
print(empty_gene_rows)
print(raw_counts_gene_symbol[13861,] %>% rownames)
raw_counts_gene_symbol <- raw_counts_gene_symbol[-13861, ]

# Create Seurat object
data <- CreateSeuratObject(counts = raw_counts_gene_symbol, min.cells = 3, project = "Wuidart")

data
# An object of class Seurat 
# 19166 features across 384 samples within 1 assay 
# Active assay: RNA (19166 features, 0 variable features)
# layer present: counts

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     29    2516    3610    3820    5190   11873 


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
abline(v = 2500)
dev.off()

# Visualize the number UMIs/transcripts per cell
plot <- data@meta.data %>% 
  ggplot(aes(x=nFeature_RNA)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 2500)
ggsave(plot = plot, filename = file.path(savePlots, "n_feature_density_2.pdf"), width = 12, height = 6)


# Filtering
data <- subset(data, subset = nFeature_RNA > 2500 & percent.mt < 8 & nCount_RNA > 100000)

# Visualize QC metrics after filtering
# number of cells before filtering: 280
plot <-  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = plot, filename = file.path(savePlots, "QC_after_filtering.pdf"), width = 8, height = 8)

# Summary of the data
summary(Matrix::colSums(data@assays$RNA$counts[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2509    3366    4496    4615    5564   11873 


saveRDS(data, file.path(saveProcessedData, "data.RDS"))

