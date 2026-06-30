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

# Log2FC and padj matrix ----
# Extract the results per genotype
res_P41 <- results(dds, name = "groups_P41_vs_WT")
res_P41_df <- as.data.frame(res_P41)

res_P42 <- results(dds, name = "groups_P42_vs_WT")
res_P42_df <- as.data.frame(res_P42)

res_P43 <- results(dds, name = "groups_P43_vs_WT")
res_P43_df <- as.data.frame(res_P43)

# Change the names
colnames(res_P41_df) <- paste0("P41_", colnames(res_P41_df))
colnames(res_P42_df) <- paste0("P42_", colnames(res_P42_df))
colnames(res_P43_df) <- paste0("P43_", colnames(res_P43_df))

# add TAIR column
res_P41_df$TAIR <- rownames(res_P41_df)
res_P42_df$TAIR <- rownames(res_P42_df)
res_P43_df$TAIR <- rownames(res_P43_df)

# Merge
int <- merge(res_P41_df, res_P42_df, by = "TAIR")
final_l2FC <- merge(int, res_P43_df, by = "TAIR")

# Log2FC matrix
l2FC_mat <- final_l2FC[,c("TAIR", "P41_log2FoldChange", "P42_log2FoldChange", "P43_log2FoldChange")]

l2FC_mat <- l2FC_mat %>% 
  column_to_rownames(var = "TAIR") %>% 
  as.matrix()

# Rename the genotypes
colnames(l2FC_mat) <- c("CYCP4;1oe", "CYCP4;2oe", "CYCP4;3oe")

# Padj matrix
padj_mat <- final_l2FC[,c("TAIR", "P41_padj", "P42_padj", "P43_padj")]

padj_mat <- padj_mat %>% 
  column_to_rownames(var = "TAIR") %>% 
  as.matrix()

# Rename the genotypes
colnames(padj_mat) <- c("CYCP4;1oe", "CYCP4;2oe", "CYCP4;3oe")

# L2FC heatmap ----

##### chose on of processes of interest
# Flavonoid biosynthesis
flav <- read.csv("data/GO_terms/KEGG_flavonoid_synthesis.csv")

# Cell cycle
flav <- read.csv("data/Cell cycle related genes.csv")
colnames(flav)[3] <- "SYMBOL"

# cell wall
flav <- read.csv("data_output/mol_proc_gene_lists/cell_wall.csv")

# Root hairs
flav <- read.csv("data_output/mol_proc_gene_lists/Root_hair_development.csv")
flav <- flav[!duplicated(flav$TAIR),]

# Cell cycle for report
flav <- read.csv("data/Cell cycle related genes_report.csv")
colnames(flav)[3] <- "SYMBOL"

# cell wall for report
flav <- read.csv("data_output/mol_proc_gene_lists/cell_wall_report.csv")

# Root hairs for report
flav <- read.csv("data_output/mol_proc_gene_lists/Root_hair_development_report.csv")
flav <- flav[!duplicated(flav$TAIR),]
flav$SYMBOL <- str_sub(flav$SYMBOL, 1, 6)
####

mat <- l2FC_mat
flav$TAIR <- str_sub(flav$TAIR, 1, 9)
flav <- flav[which(flav$TAIR %in% rownames(mat)),]
panel_of_genes <- flav$TAIR


deg_mat <- mat[panel_of_genes[which(panel_of_genes %in% rownames(mat))], ]
deg_mat <- deg_mat[panel_of_genes,]

# Filter padj
final_mat <- padj_mat
final_mat <- final_mat[panel_of_genes,]

final_mat_display <- final_mat


final_mat_display[final_mat < 0.05] <- "*"
final_mat_display[final_mat < 0.001] <- "**"
final_mat_display[final_mat < 0.0001] <- "***"
final_mat_display[final_mat > 0.05 | is.na(final_mat)] <- ""

max_val <- max(abs(deg_mat))
#max_val <- 1
my_breaks <- seq(-max_val, max_val, length.out = 101)
# Make a heatmap
my_palette <- colorRampPalette(c("#67a9cf", "#f7f7f7", "#ef8a62"))(100)

# annotation <- as.data.frame(flav$Type)
# rownames(annotation) <- flav$TAIR

# Make the gene names and the genotypes in italic
italic_rows <- lapply(flav$SYMBOL, function(x) bquote(italic(.(x))))
italic_cols <- lapply(colnames(deg_mat), function(x) bquote(italic(.(x))))

p <- pheatmap(deg_mat,
              #annotation_row = annotation,
              cluster_rows = F,        
              cluster_cols = F,          
              show_rownames = T,
              labels_row = as.expression(italic_rows),
              labels_col = as.expression(italic_cols),
              color = my_palette,
              breaks = my_breaks,
              angle_col = 45,
              border_color = F,
              cellwidth = 30,
              display_numbers = final_mat_display
              #cutree_rows = 4
) 
