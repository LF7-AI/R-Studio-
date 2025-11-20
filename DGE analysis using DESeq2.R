# Differential Gene Expression Analysis using DESeq2

# Install all required packages which are not installed before!!!!!!!!

# Install Required Packages (if not already installed)
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}
#if (!requireNamespace("DESeq2", quietly = TRUE)) {
 # BiocManager::install("DESeq2")
#}
if (!requireNamespace("RUVSeq", quietly = TRUE)) {
  BiocManager::install("RUVSeq")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  BiocManager::install("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  BiocManager::install("RColorBrewer")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

# Load Packages
library(DESeq2)
library(RUVSeq)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# ## ----warning=FALSE, message=FALSE---------------------------------------------
# Import Gene Counts
COUNTS <- read.csv(file="C:/Users/steff/Documents/R-Projects/GSE227516_counts.csv", header=TRUE, row.names=1)
head(COUNTS)

## -----------------------------------------------------------------------------
# Import Sample Information
META <- read.csv(file="C:/Users/steff/Documents/R-Projects/sample_information.csv", header=TRUE)
head(META)

## -----------------------------------------------------------------------------
# Reorder columns
COUNTS <- COUNTS[, c("P1", paste0("P", 2:9), "P10")]
head(COUNTS)

## -----------------------------------------------------------------------------
# Rounding off and convert to matrix
COUNTS <- round(COUNTS)        #rounds the numbers
COUNTS <- as.matrix(COUNTS)   #transforms to matrixs if its not

head(COUNTS)
class(COUNTS) # to check which format

## -----------------------------------------------------------------------------
# List unique sample conditions
unique(META$condition)        #to check the condition (distinct sample groups) present in the Meta-data

#Ideally, used to check before using DESeq2

# Creating DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = COUNTS, 
                              colData=META, 
                              design=~condition)

## -----------------------------------------------------------------------------
dim(dds)

## -----------------------------------------------------------------------------
# Filtering low count genes
threshold <- 10
dds <- dds[ rowMeans(counts(dds)) >= threshold,]

dim(dds)   #To check after filtered results

## -----------------------------------------------------------------------------

# DESeq2 Analysis
prdds <- DESeq(dds)
prdds

## -----------------------------------------------------------------------------
#prdds is now a DESeqDataSet object with results attached, containing:
  
# 1. normalized counts

# 2. dispersion values

# 3. fitted model

# 4. Wald test statistics

# 5. FDR-adjusted p-values
## -----------------------------------------------------------------------------

# Normalization
norm_counts <- counts(prdds, normalized = TRUE)
norm_counts <- as.data.frame(norm_counts)
head(norm_counts)

#Why we convert to to DataFrame-
#as.data.frame() converts the norm_count to dataframe from matrix.
#we do this because if its in data.frame() , then:
# - use for Plotting, merging, exporting, tidyverse.
#Matrix - good for Calculations
## -----------------------------------------------------------------------------

# Transformation
mks <- estimateSizeFactors(dds)
rld <- rlogTransformation(prdds, blind = FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# mks- Size factors correct for differences in sequencing depth / library size.
# rld- Regularized log transformation (rlog). It does:
# It transforms your count data to stabilize:
# 1. variance
# 2. noise
# 3. big differences between low-count genes.
#blind = FALSE → uses experimental design in the model
# → This is recommended when analysing the same dataset you modeled (i.e., when performing visualization after running DESeq).
#vsd- Variance Stabilizing Transformation (VST) is another method to stabilize variance across genes. Useful if we are using PCA and clustering.

###---------------------------------------------------------------------------------------------------------------------------------------------------

## ----scatter_plots, dev='png', fig.show='hide'--------------------------------
# Scatter Plots Comparison
# Can use ggplot2 instead of the png mode which had issues like:  the margin of plot is large, doesn't save if after using dev.off()
# We had to do comparision in seperate ways... we can combine them but for demo, this is good.
library(tidyverse)

df <- data.frame(
  log2 = log2(counts(mks, normalized=TRUE)[,1] + 1),
  rlog = assay(rld)[,1],
  vst  = assay(vsd)[,1]
)

# Plot 1:
ggplot(df, aes(x = log2, y = rlog, color= rlog)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "log2(x + 1) vs rlog", color = "log2 intensity") +
  theme_minimal(base_size = 14)


# Plot 2:
ggplot(df, aes(x = rlog, y = vst, colour = rlog)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_gradient(low = "darkgreen", high = "yellow") +
  labs(title = "rlog vs VST",  color = "rlog intensity") +
  theme_minimal(base_size = 14)

########-----EXPLAINATION-----#########
# So we used ggplot to compare the differences between to factor values.
#All codes are pretty much self-explanatory.
#All additional codes are for the ggplot color scatter plotting.
#Unable to do png - scatter plot as it said "Margins are to large", 
########-------------------------------------------------------------------------------------------------------------------------------------

## ----s2s_heatmap_plot, dev='png', fig.show='hide'-----------------------------
# Sample-to-sample distances
sample_dist <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dist)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows=sample_dist,
         clustering_distance_cols=sample_dist,
         col=colors,)
########-------------------------------------------------------------------------------------------------------------------------------------

## ----pca_plot, dev='png', fig.show='hide', message=FALSE, warning=FALSE-------
# PCA Plot
pca_data <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, aes(color = condition)) +
  geom_text_repel(aes(label = rownames(pca_data)), nudge_x = 0, nudge_y = 0) +
  xlab(paste0("PC1: ", round(attr(pca_data, "percentVar")[1], 2) * 100, "% variance")) +
  ylab(paste0("PC2: ", round(attr(pca_data, "percentVar")[2], 2) * 100, "% variance"))

#####EXPLAINATION###########-------

# PC1 & PC2- the PC1 and PC2 values are created for you by DESeq2, 
#and they come directly from the rlog-transformed expression matrix. 
#When you run this- pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE),
#The DESeq2 takes the top 500 genes from rld transformef counts & run PCA mathematically.
#In xlab &ylab - we use past0 for PC1 & PC2 because- this is to convert the PC values numericall varience.
#PC1 explains the most variation
#PC2 explains the second-most.
###_____________________________________________________________________________________________________________

## ----dispersion_plot, dev='png', fig.show='hide'------------------------------
# Dispersion Plot
plotDispEsts(prdds, main = "Dispersion plot", 
             genecol="gray20", fitcol="red", 
             finalcol="dodgerblue3" 
) 
#####EXPLAINATION###########-------
#This plot is used for variance test. Check Online for its 3 color significance.
#-------------------------------------------------------------------------------

## -----------------------------------------------------------------------------
# DESeq2 Result
res05 <- results(prdds, alpha = 0.05)
res05 <- na.omit(res05)

## -----------------------------------------------------------------------------
# Order by adjusted p-value
res05ordered <- res05[order(res05$padj),]
head(as.data.frame(res05ordered))

## ----ma_plot, dev='png', fig.show='hide'--------------------------------------
# MA Plot
DESeq2::plotMA(
  res05, 
  main="Sedentary vs Exercise, alpha=0.05", 
  ylim=c(-5,10),
  cex=0.5, 
  colNonSig=adjustcolor("gray20", alpha.f=0.5), 
  colSig=adjustcolor("dodgerblue3", alpha.f=0.5) 
)
abline(h = 1, col = '#ff0000' , lwd = 1)
abline(h = -1, col= '#ff0000', lwd = 1)

##------------------------------------------------------
#EXPLAINATION:
#An MA plot visualizes log2 fold change (M) vs mean normalized expression (A) for all genes.
#X-axis (A): mean expression (normalized counts → log scale)
#Y-axis (M): log2FoldChange
#Red/Blue points: significant genes
#Gray points: non-significant genes
#This helps you see upregulated/downregulated genes clearly.
#No xlim - because it plots for mean normalized count for each gene.
#And DESeq2 does it automatically it sets it.( i.e. DESeq2 auto-scales log10(mean counts))

## ----volcano_plot, dev='png', fig.show='hide'---------------------------------
# Volcano Plot
res05$gene_status <- ifelse(
  res05$padj < 0.05, 
  ifelse(
    res05$log2FoldChange > 1, 
    "Up-Regulated",
    ifelse(
      res05$log2FoldChange < -1, 
      "Down-Regulated", 
      "Non-significant"
    )
  ), 
  "Non-significant"
)
ggplot(
  res05, 
  aes(x = log2FoldChange, y = -log10(padj), color = factor(gene_status))
) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = brewer.pal(3, "Set1")) +
  theme_minimal() +
  ggtitle("Volcano Plot of Differentially Expressed Genes") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")


## -----------------------------------------------------------------------------
sig_genes <- as.data.frame(res05[res05$padj < 0.05 & abs(res05$log2FoldChange) > 1, ])
head(sig_genes)


## -----------------------------------------------------------------------------
up_genes <- subset(sig_genes, log2FoldChange > 0)
head(up_genes)

## -----------------------------------------------------------------------------
down_genes <- subset(sig_genes,log2FoldChange < 0)
head(down_genes)

## ----upreg_heatmap, dev='png', fig.show='hide'--------------------------------
top_up <- head(up_genes[order(up_genes$log2FoldChange, decreasing = TRUE), ], 10)
top_up_exp <- assay(rld)[rownames(top_up), ]
pheatmap(top_up_exp,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         show_colnames = TRUE,
         col=brewer.pal(name="RdBu", n=11),
         main = "Top Up Regulated Genes Heatmap")


## ----downreg_heatmap, dev='png', fig.show='hide'------------------------------
top_down <- head(down_genes[order(down_genes$log2FoldChange, decreasing = TRUE), ], 10)
top_down_exp <- assay(rld)[rownames(top_down), ]
pheatmap(top_down_exp,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         show_colnames = TRUE,
         col=brewer.pal(name="RdBu", n=11),
         main = "Top Down Regulated Genes Heatmap")


## -----------------------------------------------------------------------------
n_total <- nrow(COUNTS)
n_de_genes <- nrow(sig_genes)
n_up_genes <- nrow(up_genes)
n_down_genes <- nrow(down_genes)

## -----------------------------------------------------------------------------
# Combine significant genes with their categories
sig_genes$gene_status <- ifelse(sig_genes$log2FoldChange > 0, "Up-Regulated", "Down-Regulated")
sig_genes$gene_id <- rownames(sig_genes)

# Write the data frame to a file
if (!file.exists("result")) {
  dir.create("result")
}

write.table(sig_genes, file = "result/significant_DE_genes.csv", sep = ",", row.names = FALSE)


## -----------------------------------------------------------------------------
sessionInfo()