if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(
  c("DESeq2", "pheatmap", "EnhancedVolcano", "airway", "apeglm"),
  update = FALSE, ask = FALSE
)

library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(airway)
library(apeglm)
library(ggplot2)

data("airway")
countData <- assay(airway)
colData   <- as.data.frame(colData(airway))
# DESeq2 uses "dex" column (trt vs untrt)
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ dex
)
# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10,]

# Step 2: Run DESeq2
dds <- DESeq(dds)
# Check available coefficients
print(resultsNames(dds))
# Automatically choose correct coef
coef_name <- resultsNames(dds)[grep("dex", resultsNames(dds))]
cat("Using coefficient:", coef_name, "\n")
# Shrinkage with correct coef
resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")

# Step 3: Show Results
head(resLFC[order(resLFC$padj), ])

# Step 4: Volcano Plot
EnhancedVolcano(
  resLFC,
  lab = rownames(resLFC),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Volcano Plot: Treated vs Untreated"
)

# Step 5: Heatmap of Top 20 DE Genes
rld <- rlog(dds)
top20 <- rownames(resLFC[order(resLFC$padj), ])[1:20]
pheatmap(
  assay(rld)[top20, ],
  scale = "row",
  main = "Top 20 Differentially Expressed Genes"
)
# Step 6: PCA Plot
plotPCA(rld, intgroup = "dex")
# Step 7: Save Results
write.csv(as.data.frame(resLFC), "DESeq2_results_airway.csv")


###========================================================================#####

#Example 2 :DUAL RNA-SEQ: Human + Bacteria DESeq2 Tutorial


#(No need to reload BiocManager)

library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
set.seed(123)


# STEP 1 — SIMULATE DUAL RNA-SEQ DATA
human_genes    <- paste0("HUMAN_GENE_", 1:500)
bacteria_genes <- paste0("BAC_GENE_", 1:200)
all_genes <- c(human_genes, bacteria_genes)
samples <- paste0("S", 1:6)
condition <- factor(c("Control", "Control", "Control", "Infected", "Infected", "Infected"))
# Human gene counts (mild infection response)
human_counts <- matrix(
  rpois(500 * 6, lambda = rep(c(50, 55, 60, 200, 220, 240), each = 500)),
  nrow = 500
)

# Bacterial genes (appear only in infected)
microbe_counts <- matrix(
  rpois(200 * 6, lambda = rep(c(2, 3, 4, 150, 180, 200), each = 200)),
  nrow = 200
)

countData <- rbind(human_counts, microbe_counts)
rownames(countData) <- all_genes
colnames(countData) <- samples
colData <- data.frame(row.names = samples, condition = condition)


# STEP 2 — Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Infected", "Control"))


# Add organism annotation
res$organism <- ifelse(grepl("HUMAN", rownames(res)), "Human", "Bacteria")

# STEP 3 — Volcano Plot (CORRECTED)
color_map <- c(Human = "blue", Bacteria = "red")
EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.01,
  FCcutoff = 1.5,
  colCustom = color_map[res$organism],
  selectLab = rownames(res)[1:10],
  title = "Dual RNA-Seq Volcano Plot",
  subtitle = "Human (blue) vs Bacteria (red)"
)



# STEP 4 — Corrected VST + Heatmap
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsd_mat <- assay(vsd)
# Select top 30 variable genes
top_genes <- head(order(apply(vsd_mat, 1, var), decreasing = TRUE), 30)
pheatmap(
  vsd_mat[top_genes, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  main = "Top 30 Variable Genes (Dual RNA-Seq)"
)

# STEP 5 — PCA PLOT
plotPCA(vsd, intgroup = "condition") +
  ggplot2::ggtitle("PCA Plot – Human + Bacterial Transcriptomes")
# STEP 6 — Organism-wise Summary
summary_table <- data.frame(
  Organism = c("Human", "Bacteria"),
  Significant = c(
    sum(res$padj < 0.05 & res$organism == "Human", na.rm = TRUE),
    sum(res$padj < 0.05 & res$organism == "Bacteria", na.rm = TRUE)
  )
)
print(summary_table)


##########====================================#########











           

