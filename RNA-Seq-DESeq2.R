#Install packages (run once) (if required)
install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install(c("DESeq2", "pasilla", "pheatmap", "EnhancedVolcano"))

#===========================================================================

# Load Libraries
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)  # optional, but fine to load explicitly

#========================================================

#Step-2: Create demo count matrix (its important for finding RNA gene exp)
# Simulated gene expression counts for 6 samples (3 control, 3 treated)

#THIS MATRIX IS ROUGHLY OR DEMO, NOT REAL BIOLOGICAL DATA#

counts <- matrix(
  c(400,380,410,800,850,870,    #Gene1
    50,45,60,200,220,210,       #Gene2
    600,620,590,610,630,600,    #Gene3
    1000,950,970,200,180,190,   #Gene4
    300,310,290,100,120,110),   #Gene5
  nrow=5,
  byrow=TRUE)
rownames(counts) <- paste0("Gene", 1:5)
colnames(counts) <- paste0("Sample", 1:6)

#Check Matrix 
counts

#Step-3: Define sample conditions
condition <- factor(c("Control", "Control", "Control", "Treated", "Treated", "Treated"))
coldata <- data.frame(rownames= colnames(counts), condition)

#Step-4: Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

#Step-5: Run DESeq2 Analysis
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

#Step-6: Volcano plot

# Convert DESeq2 result to dataframe
res_df <- as.data.frame(res)

# Create threshold column for significance
res_df$threshold <- as.factor(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05)

# Volcano plot
ggplot(res_df, aes(x= log2FoldChange, y= -log10(padj), color= threshold))+
  geom_point(size= 3, alpha=0.7)+
  theme_minimal(base_size=14)+
  scale_color_manual(values= c("grey","red"))+
  labs(title= "Volcano Plot", x= "log2FoldChange", y="-log10(padj)")

# Heatmap of normalized counts----------------
norm_counts <- counts(dds, normalized= TRUE)
pheatmap(log2(norm_counts+ 1), annotation_col= coldata,
         main= "Normalized Expression Heatmap")

#Step:8- SUMMARY
cat("\nSignificant genes(padj < 0.05):", sum(res$padj< 0.05, na.rm= TRUE), "\n")

###################################################################
# Create a temporary annotation object (recommended)

annotation_col_tmp <- data.frame(condition = coldata$condition)
rownames(annotation_col_tmp) <- coldata$rownames

#Now run 
pheatmap(log2(norm_counts + 1),
         annotation_col = annotation_col_tmp,
         main = "Normalized Expression Heatmap")
