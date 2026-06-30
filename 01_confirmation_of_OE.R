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


# Barplots of CYCPs OE ----
# 1. Choose your gene of interest, CYCP4s
target_gene <- c("AT2G44740", "AT5G61650", "AT5G07450")

# 2. Extract the continuous rlog matrix 
rlog_matrix <- assay(rld)

# 3. Filter the matrix for your specific genes
subset_matrix <- rlog_matrix[rownames(rlog_matrix) %in% target_gene, ]

# 4. Convert the matrix into a data frame and keep gene names
matrix_df <- as.data.frame(subset_matrix) %>%
  mutate(Gene = rownames(.))

# 5. Reshape from Wide to Long format
long_df <- matrix_df %>%
  pivot_longer(
    cols = -Gene, 
    names_to = "Sample_ID", 
    values_to = "Rlog_Expression"
  )


# 3. Merge with your metadata mapping sheet (sample_info)
metadata_df <- as.data.frame(colData(dds)) %>% 
  mutate(Sample_ID = rownames(.))

plot_data <- left_join(long_df, metadata_df, by = "Sample_ID")

# Adjust the names for a proper figure
plot_data$Gene_name <- ifelse(plot_data$Gene == "AT2G44740", "CYCP4;1", ifelse(plot_data$Gene == "AT5G07450", "CYCP4;3", "CYCP4;2"))

# Change the genotype names for the proper line names
plot_data$Genotype <- ifelse(plot_data$groups == "WT", "WT", ifelse(plot_data$groups == "P41", "CYCP4;1oe", ifelse(plot_data$groups == "P42", "CYCP4;2oe", "CYCP4;3oe")))

# Make a summary for the barplot
plot_summary <- plot_data %>%
  group_by(Gene, Gene_name, Genotype, groups) %>% 
  summarise(
    mean_expr = mean(Rlog_Expression),
    sd_expr = sd(Rlog_Expression),
    n = n(),
    se_expr = sd_expr / sqrt(n),
    .groups = "drop"
  )

# Barplot with standard error and jitter plot
dodge_width <- 0.6

ggplot(plot_summary, aes(x = Gene_name, y = mean_expr, fill = Genotype)) +
  # 1. Draw the expression bars
  geom_bar(stat = "identity", 
           #color = "black", 
           width = 0.6, 
           alpha = 0.8,
           position = position_dodge(width = dodge_width)) +
  
  # 2. Overlay the individual sample rlog values as points
  geom_jitter(data = plot_data, aes(x = Gene_name, y = Rlog_Expression, group = Genotype),
              # width = 0.15,
              # size = 1.5,
              alpha = 0.5,
              color = "black",
              inherit.aes = F,
              position = position_dodge(width = dodge_width)) +
  
  # Add standard error bars
  geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                width = 0.2, 
                color = "black", 
                linewidth = 0.7,
                position = position_dodge(width = dodge_width)) +
  
  # Graph polishing
  theme_classic() + # Clean grid layout
  labs(
    x = "",
    y = "Regularized Log2 Expression") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    strip.text = element_text(face = "bold", size = 12), # Formats the gene labels
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10, color = "black"),
    legend.text = element_text(face = "italic"),
    legend.position = "right") +
  scale_fill_brewer(palette = "Set1",
                    breaks = c("WT", "CYCP4;1oe", "CYCP4;2oe", "CYCP4;3oe"))
