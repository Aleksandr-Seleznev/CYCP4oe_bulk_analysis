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

# Delete the CYCPs!!!
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

# Sample exclusion:
dds <- dds[,!(colnames(dds) %in% c("WT_3", "P41_23"))]

# Set WT as a baseline level
dds$groups <- relevel(dds$groups, ref = "WT")

# Run DESeq2 ----
dds <- DESeq2::DESeq(dds)

# Normalise the data ----
rld <- rlog(dds, blind = T)

# Derive DEGs:
# Derive DEGs ----
resultsNames(dds)

# DEfine thresholds
FC <- 1.5 # set up Fold Change threshold
FDR <- 0.05 # Set up p-value threshold

# Comparison 1 of 3:  P41-WT
res <- results(dds, name = "groups_P41_vs_WT")

res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort


#### assign TAIRids with gene names
# convert it to the data frame
res_2 <- as.data.frame(res)

# make TAIR ids a separate column
res_2$gene_ids <- rownames(res_2)

#get then the common names from TAIR
#BiocManager::install("org.At.tair.db")
library(org.At.tair.db)

# Your gene IDs
gene_ids <- res_2$gene_ids

# Get symbols
mapped <- AnnotationDbi::select(
  org.At.tair.db,
  keys = gene_ids,
  keytype = "TAIR",
  columns = c("SYMBOL", "GENENAME")
)

#Data frame with the IDs, Gene names and their annotaions
#but with IDs duplicates

#Delete the duplicatess
mapped <- mapped %>%
  distinct(TAIR, .keep_all = TRUE) %>% 
  arrange(SYMBOL) 

#Merge with the data frame of DEGs
new <- merge(mapped, res_2, by.x = "TAIR", by.y = "gene_ids")

# save the DEGs with ATXGXXXXX ids and common gene names
write.csv(new, file = "data/DEGs_new/DEGs_P41_WT.csv")

# Comparison 2 of 3:  P42-WT
res <- results(dds, name = "groups_P42_vs_WT")

res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort

#### assign TAIRids with gene names
# convert it to the data frame
res_2 <- as.data.frame(res)

# make TAIR ids a separate column
res_2$gene_ids <- rownames(res_2)

#get then the common names from TAIR
#BiocManager::install("org.At.tair.db")
library(org.At.tair.db)

# Your gene IDs
gene_ids <- res_2$gene_ids

# Get symbols
mapped <- AnnotationDbi::select(
  org.At.tair.db,
  keys = gene_ids,
  keytype = "TAIR",
  columns = c("SYMBOL", "GENENAME")
)

#Data frame with the IDs, Gene names and their annotaions
#but with IDs duplicates
View(mapped)

#Delete the duplicatess
mapped <- mapped %>%
  distinct(TAIR, .keep_all = TRUE) %>% 
  arrange(SYMBOL) 

#Merge with the data frame of DEGs
new <- merge(mapped, res_2, by.x = "TAIR", by.y = "gene_ids")

# save the DEGs with ATXGXXXXX ids and common gene names
write.csv(new, file = "data/DEGs_new/DEGs_P42_WT.csv")

# Comparison 3 of 3:  P43-WT
res <- results(dds, name = "groups_P43_vs_WT")

res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort

#### assign TAIRids with gene names
# convert it to the data frame
res_2 <- as.data.frame(res)

# make TAIR ids a separate column
res_2$gene_ids <- rownames(res_2)

#get then the common names from TAIR
#BiocManager::install("org.At.tair.db")
library(org.At.tair.db)

# Your gene IDs
gene_ids <- res_2$gene_ids

# Get symbols
mapped <- AnnotationDbi::select(
  org.At.tair.db,
  keys = gene_ids,
  keytype = "TAIR",
  columns = c("SYMBOL", "GENENAME")
)


#Delete the duplicatess
mapped <- mapped %>%
  distinct(TAIR, .keep_all = TRUE) %>% 
  arrange(SYMBOL) 

#Merge with the data frame of DEGs
new <- merge(mapped, res_2, by.x = "TAIR", by.y = "gene_ids")

# save the DEGs with ATXGXXXXX ids and common gene names
write.csv(new, file = "data/DEGs_new/DEGs_P43_WT.csv")


# Venn Diagrams ----
# DEGs VennDiagramm ----
P41 <- read.csv("data/DEGs_new/DEGs_P41_WT.csv")
P42 <- read.csv("data/DEGs_new/DEGs_P42_WT.csv")
P43 <- read.csv("data/DEGs_new/DEGs_P43_WT.csv")

# Distinguish between up and down
l2FC_thr <- log2(2)
padj <- 0.05

CYCPs_up <- list(
  `CYCP4;1oe_UP` = P41[P41$log2FoldChange >= l2FC_thr & P41$padj < padj, "TAIR"],
  
  `CYCP4;2oe_UP` = P42[P42$log2FoldChange >= l2FC_thr & P42$padj < padj, "TAIR"],
  
  `CYCP4;3oe_UP` = P43[P43$log2FoldChange >= l2FC_thr & P43$padj < padj, "TAIR"]
)

CYCPs_down <- list(
  `CYCP4;1oe_DOWN` = P41[P41$log2FoldChange <= -l2FC_thr & P41$padj < padj, "TAIR"],
  
  `CYCP4;2oe_DOWN` = P42[P42$log2FoldChange <= -l2FC_thr & P42$padj < padj, "TAIR"],
  
  `CYCP4;3oe_DOWN` = P43[P43$log2FoldChange <= -l2FC_thr & P43$padj < padj, "TAIR"]
)


# Define a pretty, publication-ready color palette
pretty_colors <- c("#E94849", "#5F98C6", "#71BF6E") 

# Generate the plot
library(ggvenn)

# CHose the list for the plot
list_to_use <- CYCPs_down

ggvenn(
  list_to_use, 
  columns = names(list_to_use), # Choose which groups to show (up to 4)
  fill_color = pretty_colors,
  stroke_size = 0,         # Thins out the border lines for a cleaner look
  stroke_color = "grey40",   # Softens the border color from harsh black,
  stroke_alpha = 0.05,
  fill_alpha = 0.8,          # Keeps circles translucent so intersections look nice
  text_size = 3.5,           # Adjusts font size of the counts inside the circles
  set_name_size = 3.5,          # Adjusts font size of the outer group labels
  show_percentage = F
) 

# DEGs Volcano plots ----
library(EnhancedVolcano)
library(org.At.tair.db) # to add common gene names to the DEGs
library(tidyverse)
library(ggrepel)

# Custom function to add gene names to the results of the DESeq2
#### Add the gene names (TO BE OPTIMIZED)
get_IDs <- function(x){
  gene_ids <- data.frame(genes = rownames(x), order = c(1:nrow(x)))
  
  # Get symbols
  mapped <- AnnotationDbi::select(
    org.At.tair.db,
    keys = gene_ids$genes,
    keytype = "TAIR",
    columns = c("SYMBOL", "GENENAME")
  )
  
  
  
  #Delete the duplicatess
  mapped <- mapped %>%
    distinct(TAIR, .keep_all = TRUE) %>% 
    arrange(SYMBOL) 
  
  #Merge with the data frame of DEGs
  new <- merge(mapped, gene_ids, by.x = "TAIR", by.y = "genes")
  
  # Check if the order is maintained
  # if(any(new$order == c(1:nrow(x))) == F){
  #   print("order is not maintained")
  # } else {
  #   print("vse ok")
  # }
  
  # Input the gene names into DESeq2 xults
  x$gene_names <- new$SYMBOL
  return(x)
}
####

# Set up thresholds
l2FC_thr <- log2(2)
padj <- 0.05

# P41 VS WT
res <- results(dds, name = "groups_P41_vs_WT")

#Get the IDs
res <- get_IDs(x = res)

# Calculate counts 
n_down <- sum(res$log2FoldChange < -l2FC_thr & res$padj < padj, na.rm = TRUE)
n_up   <- sum(res$log2FoldChange > l2FC_thr & res$padj < padj, na.rm = TRUE)

# Create a dynamic subtitle
my_subtitle <- paste("Up:", n_up, "| Down:", n_down)


#generate the plot
EnhancedVolcano(res,
                #lab = res$gene_names,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = res$gene_names[which((res$log2FoldChange < -2.5 | res$log2FoldChange > 2.50) & res$padj < 0.001)],
                title = '',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 4.0,
                xlim = c(-7, 8),
                #ylim = c(-0.25, 25),
                subtitle = my_subtitle,
                captionLabSize = 10,
                legendLabSize = 10,
                legendIconSize = 3,
                legendPosition = "none",
                gridlines.major = F,
                gridlines.minor = F,
                axisLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                colConnectors = "black",
                boxedLabels = F,
                col=c('grey60', 'grey60', 'grey60', 'red3')
)




# P42 VS WT 
# Extract the results
res <- results(dds, name = "groups_P42_vs_WT")


#Get the IDs
res <- get_IDs(x = res)

# Calculate counts (same as above)
n_down <- sum(res$log2FoldChange < -l2FC_thr & res$padj < padj, na.rm = TRUE)
n_up   <- sum(res$log2FoldChange > l2FC_thr & res$padj < padj, na.rm = TRUE)

# Create a dynamic subtitle
my_subtitle <- paste("Up:", n_up, "| Down:", n_down)

#generate the plot
EnhancedVolcano(res,
                #lab = res$gene_names,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = res$gene_names[which(res$log2FoldChange < -3.1 | res$log2FoldChange > 5 | res$padj < 10**-18)],
                title = '',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 4.0,
                xlim = c(-8, 10),  # there is a CeQORH at around (-5,100)
                ylim = c(-0.25, 35),
                subtitle = my_subtitle,
                captionLabSize = 10,
                legendLabSize = 10,
                legendIconSize = 3,
                gridlines.major = F,
                gridlines.minor = F,
                axisLabSize = 10,
                legendPosition = "none",
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                colConnectors = "black",
                boxedLabels = F,
                col=c('grey60', 'grey60', 'grey60', 'red3')
)



# P43 VS WT 
res <- results(dds, name = "groups_P43_vs_WT")

#Get the IDs
res <- get_IDs(x = res)


# Calculate counts (same as above)
n_down <- sum(res$log2FoldChange < -l2FC_thr & res$padj < padj, na.rm = TRUE)
n_up   <- sum(res$log2FoldChange > l2FC_thr & res$padj < padj, na.rm = TRUE)

# Create a dynamic subtitle
my_subtitle <- paste("Up:", n_up, "| Down:", n_down)

res_df <- as.data.frame(res)
#generate the plot
EnhancedVolcano(res,
                #lab = res$gene_names,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = res$gene_names[which(res$log2FoldChange < -3 | res$log2FoldChange > 3.5 | res$padj < 10**-10)],
                title = '',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 4.0,
                xlim = c(-7, 9),
                ylim = c(-0.25, 30),
                subtitle = my_subtitle,
                captionLabSize = 10,
                legendLabSize = 10,
                legendIconSize = 3,
                gridlines.major = F,
                gridlines.minor = F,
                axisLabSize = 10,
                legendPosition = "none",
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                colConnectors = "black",
                boxedLabels = F,
)

# DEGs for k-mean clustering heatmap ----
P41 <- read.csv("data/DEGs_new/DEGs_P41_WT.csv")
P42 <- read.csv("data/DEGs_new/DEGs_P42_WT.csv")
P43 <- read.csv("data/DEGs_new/DEGs_P43_WT.csv")

FC_thr <- log2(2)
p_thr <- 0.05

CYCPs_all <- list(
  P41_all = P41[abs(P41$log2FoldChange) >= FC_thr & P41$padj < p_thr, ],
  
  P42_all = P42[abs(P42$log2FoldChange) >= FC_thr & P42$padj < p_thr,],
  
  P43_all = P43[abs(P43$log2FoldChange) >= FC_thr & P43$padj < p_thr,]
)

int <- merge(CYCPs_all[["P41_all"]], CYCPs_all[["P42_all"]], by = "TAIR", all = T)
DEGs_all_new <- merge(int, CYCPs_all[["P43_all"]], by = "TAIR", all = T)

deg_genes <- DEGs_all_new$TAIR

# DEGs k-mean clustering ----
deg_mat <- assay(rld)
deg_mat <- deg_mat[deg_genes[which(deg_genes %in% rownames(deg_mat))], ]

# 3. Z-score scale the rows (genes)
# scale() works on columns, so we transpose t(), scale, then transpose back t()
deg_mat_scaled <- t(scale(t(deg_mat)))

#install.packages("factoextra")
library(factoextra)

# Run the elbow method automatically
# 'my_data' should be a standardized/scaled matrix or data frame of your features
fviz_nbclust(
  deg_mat_scaled, 
  FUNcluster = kmeans, 
  method = "wss",   # "wss" stands for Within-Cluster Sum of Squares OR "silhouette"
  k.max = 10         # Maximum number of k to check
) + 
  theme_classic() + 
  labs(title = "")

# 1. Run K-means (let's assume 4 clusters for this example)
set.seed(13) # Set seed for reproducibility
k_chosen <- 4
km <- kmeans(deg_mat_scaled, centers = k_chosen)

# 2. Create a data frame for row annotations
annotation_row <- data.frame(Cluster = factor(km$cluster))
rownames(annotation_row) <- names(km$cluster)

# Create annotation for the columns
annotation_col <- data.frame(Genotype = factor(c(rep("WT",3), rep("CYCP4;1oe", 3), rep("CYCP4;2oe", 4), rep("CYCP4;3oe", 4)), levels = c("WT", "CYCP4;1oe", "CYCP4;2oe", "CYCP4;3oe")))
rownames(annotation_col) <- colnames(deg_mat_scaled)

# Create custom colors
my_color = list(
  Genotype = c(WT = "#984EA3", `CYCP4;1oe` = "#E94849", `CYCP4;2oe` = "#5F98C6", `CYCP4;3oe` = "#71BF6E"),
  Cluster = c(`1` = "#C54063", `2` = "#F4AE09", `3` = "#3BA799", `4` = "#CBA1E9")
)

# 3. Reorder the matrix rows so the clusters sit neatly together
deg_mat_scaled_sorted <- deg_mat_scaled[order(km$cluster), ]

max_val <- max(abs(deg_mat_scaled_sorted))
my_breaks <- seq(-max_val, max_val, length.out = 101)

# 4. Draw the heatmap
p <- pheatmap(deg_mat_scaled_sorted,
              annotation_row = annotation_row,
              annotation_col = annotation_col,
              annotation_colors = my_color,
              breaks = my_breaks,
              annotation_names_col = F,
              annotation_names_row = F,
              cluster_rows = FALSE,         # FALSE because we manually sorted by K-means
              cluster_cols = F,          # Let columns (samples) cluster hierarchically
              show_rownames = FALSE,        # Turn off if you have hundreds of DEGs
              show_colnames = F,
              cellwidth = 8,
              #main = "K-means Clustering of DEGs",
              color = colorRampPalette(c("blue", "white", "red"))(100))

# Extract the gene names from the clusters
cluster_gene_lists <- split(names(km$cluster), km$cluster)

# Use the following to insert the gene names to string.db for GO analysis
writeClipboard(cluster_gene_lists[["1"]])
