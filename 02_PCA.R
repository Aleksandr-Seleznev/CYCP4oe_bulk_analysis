library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(stringr)

# DESeq2 object ----
# Use the "Converted counts" button in the Pre-Process tab of iDEP website
# to download the low count filtered counts 
raw_counts <-  read.csv("data/bulk_RNA/03.Result_X204SC25084634-Z01-F003_Arabidopsis_thaliana/Result_X204SC25084634-Z01-F003_Arabidopsis_thaliana/3.Quant/1.Count/gene_count_TAIRs_for_pr_PCA_converted.csv")
row.names(raw_counts) <- raw_counts$User_ID # asign TAIR ids to be rownames of the count matrix
raw_counts <- raw_counts[, -(1:4)] # delete 4 columns of IDs

# Delete the CYCPs so they don't affect clustering!!!
raw_counts <- raw_counts[!(rownames(raw_counts) %in% c("AT2G44740", "AT5G61650", "AT5G07450")),]
#

col_data <- data.frame(
  "sample" = c("WT_1", "WT_3", "WT_5", "WT_6", "P41_22", "P41_23", "P41_24", "P41_26", "P42_28", "P42_30", "P42_31", "P42_32", "P43_34", "P43_35", "P43_37", "P43_38"), # sample
  "groups" = c("WT", "WT", "WT", "WT", "P41", "P41", "P41", "P41", "P42", "P42", "P42", "P42", "P43", "P43", "P43", "P43") # genotype
)

# Make a DESeq2 object ----
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = col_data,
  design = ~groups
)

# Set WT as a baseline level
dds$groups <- relevel(dds$groups, ref = "WT")

# Normalise the data ----
rld <- rlog(dds, blind = T)

# PCA ----
# custom PCA
pca_data <- plotPCA(rld, intgroup = "groups", ntop = 23000, pcsToUse = c(1,2), returnData = T)


percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +  # Draws the sample dots
  geom_text_repel(aes(label = sample), 
                  size = 3,                 # Controls text size
                  box.padding = 0.5,        # Distance from the point
                  point.padding = 0.3,      # Space around the dot
                  show.legend = FALSE) +
  labs(x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%"),
       color = "Genotype") +
  theme_classic()

# Sample exclusion:
dds <- dds[,!(colnames(dds) %in% c("WT_3", "P41_23"))]
# Repeat everything from data normalistion