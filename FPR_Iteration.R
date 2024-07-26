#### Directories ####
# work directory
setwd("~/Dropbox (GaTech)/TLC/McGrath Multiome/Data Analysis/Harini")
source('CellBenderv3_FUNCTIONS.R')

#### Load Libraries ####
library(Seurat)
library(dplyr)
library(scran)
library(stringr)
library(BiocManager)
library(sctransform)
library(patchwork)
library(ggplot2)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(Signac)
library(metap)
library(DoubletFinder)
library(tidyverse)
library(qlcMatrix)
library(biomartr)
library(Biostrings)
library(reshape2)
library(hdf5r)
library(BSgenome)
library(ggpubr)
library(DropletUtils) # exporting 10x as H5
library(scCustomize)

library(dplyr)
library(Seurat)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(DropletQC)
library(harmony)
library(devtools)
library(GeneOverlap)
library(Rsamtools)
library(GenomicRanges)
library(rhdf5)
library(Matrix)

#### Read in the Count Matrix ####
# for our purposes, we'll just focus on RNA and ignore the ATAC
########################################FPR_0.01###################################################

CBcounts <- Read_CellBender_h5_Mat("cellbender_tel12_01_FPR_0.01_filtered.h5")

afterCB <- CreateSeuratObject(
  counts = CBcounts,
  assay = 'RNA',
  project = "tel12_01"
)

# add mito percent to meta data
afterCB[["percent.mt"]] <- PercentageFeatureSet(afterCB, pattern = "^APK84-")

# add intronic read ratio from original CellRanger (i.e. before ambient removal) to meta data
nucFrac <- readRDS("nucFractionList_tel12_ALL.RDS")
afterCB$preCB_intronRat = nucFrac[colnames(afterCB), 1]

# plot our metrics

VlnPlot(afterCB, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))

# test out some subsetting parameters, remember to separate by &

test <- subset(afterCB, subset = percent.mt < 25)

# plot the subsetted metrics

VlnPlot(test, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))


# using GG plot to visualize output

ggplot(afterCB@meta.data, aes(percent.mt,preCB_intronRat*100)) +
  theme_classic() +
  geom_point(alpha = 0.3)  +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  labs(title ="Intron Rat vs. Mito PCT", x = "Percent Mito", y = "Intron Ratio") +
  ylim(0,100)

########################################FPR_0.1###################################################
CBcounts1 <- Read_CellBender_h5_Mat("cellbender_tel12_01_FPR_0.1_filtered.h5")

afterCB1 <- CreateSeuratObject(
  counts = CBcounts1,
  assay = 'RNA',
  project = "tel12_01"
)

# add mito percent to meta data
afterCB1[["percent.mt"]] <- PercentageFeatureSet(afterCB1, pattern = "^APK84-")

# add intronic read ratio from original CellRanger (i.e. before ambient removal) to meta data
nucFrac <- readRDS("nucFractionList_tel12_ALL.RDS")
afterCB1$preCB_intronRat = nucFrac[colnames(afterCB1), 1]

# plot our metrics

VlnPlot(afterCB1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))

# test out some subsetting parameters, remember to separate by &

test1 <- subset(afterCB1, subset = percent.mt < 10 & nCount_RNA < 50000 & nFeature_RNA < 20000 & preCB_intronRat >= 0.1 & preCB_intronRat < 0.5)

# plot the subsetted metrics

VlnPlot(test1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))


# Count the number of cells in the 'test' subset
num_cells_in_test1 <- nrow(test1@meta.data)

# Print or use the count as needed
print(num_cells_in_test1)


# using GG plot to visualize output

ggplot(afterCB1@meta.data, aes(percent.mt,preCB_intronRat*100)) +
  theme_classic() +
  geom_point(alpha = 0.3)  +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  labs(title ="Intron Rat vs. Mito PCT", x = "Percent Mito", y = "Intron Ratio") +
  ylim(0,100)

########################################FPR_0.03###################################################
CBcounts2 <- Read_CellBender_h5_Mat("cellbender_tel12_01_FPR_0.03_filtered.h5")

afterCB2 <- CreateSeuratObject(
  counts = CBcounts2,
  assay = 'RNA',
  project = "tel12_01"
)

# add mito percent to meta data
afterCB2[["percent.mt"]] <- PercentageFeatureSet(afterCB2, pattern = "^APK84-")

# add intronic read ratio from original CellRanger (i.e. before ambient removal) to meta data
nucFrac <- readRDS("nucFractionList_tel12_ALL.RDS")
afterCB2$preCB_intronRat = nucFrac[colnames(afterCB2), 1]

# plot our metrics

VlnPlot(afterCB2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))

# test out some subsetting parameters, remember to separate by &

test2 <- subset(afterCB2, subset = percent.mt < 10 & nCount_RNA < 50000 & nFeature_RNA < 20000 & preCB_intronRat >= 0.1 & preCB_intronRat < 0.5)

# plot the subsetted metrics

VlnPlot(test2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))


# Count the number of cells in the 'test' subset
num_cells_in_test2 <- nrow(test2@meta.data)

# Print or use the count as needed
print(num_cells_in_test2)

########################################FPR_0.05###################################################
CBcounts3 <- Read_CellBender_h5_Mat("cellbender_tel12_01_FPR_0.05_filtered.h5")

afterCB3 <- CreateSeuratObject(
  counts = CBcounts3,
  assay = 'RNA',
  project = "tel12_01"
)

# add mito percent to meta data
afterCB3[["percent.mt"]] <- PercentageFeatureSet(afterCB3, pattern = "^APK84-")

# add intronic read ratio from original CellRanger (i.e. before ambient removal) to meta data
nucFrac <- readRDS("nucFractionList_tel12_ALL.RDS")
afterCB3$preCB_intronRat = nucFrac[colnames(afterCB3), 1]

# plot our metrics

VlnPlot(afterCB3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))

# test out some subsetting parameters, remember to separate by &

test3 <- subset(afterCB3, subset = percent.mt < 10 & nCount_RNA < 50000 & nFeature_RNA < 20000 & preCB_intronRat >= 0.1 & preCB_intronRat < 0.5)

# plot the subsetted metrics

VlnPlot(test3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))


# Count the number of cells in the 'test' subset
num_cells_in_test3 <- nrow(test3@meta.data)

# Print or use the count as needed
print(num_cells_in_test3)


# using GG plot to visualize output

ggplot(afterCB3@meta.data, aes(percent.mt,preCB_intronRat*100)) +
  theme_classic() +
  geom_point(alpha = 0.3)  +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  labs(title ="Intron Rat vs. Mito PCT", x = "Percent Mito", y = "Intron Ratio") +
  ylim(0,100)


########################################### Loop through different FPR values################################################################
fpr_values <- c(0.01, 0.03, 0.05, 0.1)
num_cells <- numeric(length(fpr_values))

for (i in seq_along(fpr_values)) {
  # Read CellBender matrix for the specific FPR value
  CBcounts <- Read_CellBender_h5_Mat(paste0("cellbender_tel12_01_FPR_", fpr_values[i], "_filtered.h5"))
  
  # Create Seurat object
  afterCB <- CreateSeuratObject(counts = CBcounts, assay = 'RNA', project = paste("tel12_01_FPR", fpr_values[i]))
  
  # Add metrics to meta data
  afterCB[["percent.mt"]] <- PercentageFeatureSet(afterCB, pattern = "^APK84-")
  nucFrac <- readRDS("nucFractionList_tel12_ALL.RDS")
  afterCB$preCB_intronRat <- nucFrac[colnames(afterCB), 1]
  
  # Subset based on criteria
  subset_data <- subset(afterCB, subset = percent.mt < 10 & nCount_RNA < 50000 & nFeature_RNA < 20000)# & preCB_intronRat >= 0.1 & preCB_intronRat < 0.5)
  
  # Count the number of cells
  num_cells[i] <- nrow(subset_data@meta.data)
  
  # Visualize metrics
  VlnPlot(subset_data, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "preCB_intronRat"))
}

# Compare the number of cells for different FPR values
num_cells_df <- data.frame(FPR = fpr_values, Num_Cells = num_cells)
print(num_cells_df)

# Visualize the results
ggplot(num_cells_df, aes(x = FPR, y = Num_Cells)) +
  geom_point() +
  labs(title = "Number of Cells vs. FPR", x = "False Positive Rate", y = "Number of Cells")