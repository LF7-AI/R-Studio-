#Install packages (run once)
install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install(c("DESeq2", "pasilla", "pheatmap", "EnhancedVolcano"))

############

# Load Libraries
library(DESeq2)
library(pasilla)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)  # optional, but fine to load explicitly
###########################

# Load pasilla count data (Its a dataset for fruti fly)
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package= "pasilla")
raw_counts <- read.table(datafile, header= TRUE, row.names= 1, stringsAsFactors = FALSE)

head(raw_counts)

#View(raw_counts) - this will open completely a seperate tabular spot on your running field, same as head()

#Keep only sample columnbs (7 samples)
countData <- raw_counts[, c("untreated1", "untreated2", "untreated3", "untreated4", "treated1", "treated2", "treated3")]
# in Dataframes i.e. dataframe[rows x columns]. The part here before "," is rows and after that is "columns".
#In this above equations it means include all the rows (,) and only take the columns which I selected.

#Convert to numeric matrix safely (meaning don't miss row numeric values like gene IDs)

#Load and align metadata
annotation<- system.file("extdata/pasilla_sample_annotation.csv", package= "pasilla")
colData<- read.csv(annotation, row.names= 1)   #Sample metadata

head(colData)   #To check the progress (1st progress)
#So for DESeq2- countData must needs to exactly match and aligned with the row names in colData


#Create rownames from colData to match countData columns 
rownames(colData) <- c("untreated1", "untreated2", "untreated3", "untreated4", "treated1", "treated2", "treated3")

head(rownames(colData)) #To check the progress

#Reorder rows to match countData colummns
colData <- colData[colnames(countData),]

#print(colnames(colData)) - just cheeking

head(colData) #To check the progress

# Reorders rows in ColData to match column order in countData
#This is to align metadata with count matrix for DESeq2

# Verify ALIGNMENT
all(colnames(countData) == rownames(colData)) #Should return true if successful

#Now keep only condition column

colData <- colData["condition"]
colData$condition <- factor(colData$condition)
# "$" - means to modify, in above it means to access the "condition" from the dataframe of colData
head(colData) #checking the progress
#=======================================================

# Creating DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData= colData,
                              design= ~condition) #if gene expression depends on variable condition.

#Filter low-count genes
dds <- dds[rowSums(counts(dds))>10,]

#counts(dds)- Extract the raw count matrix stored inside in the DESeqDataset
#rowSums() gives you a vector of total counts per gene.
#rowSums(counts(dds)) > 10- This creates a logical filter — TRUE for genes whose total count is > 10, FALSE otherwise.
# If you want to see progress of your sets use- head() or View()
#============================================

# Run DESeq2

dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)

#resOrdered <- res[order(res$padj), ]- Sort the results table by the adjusted p-value (padj) from smallest to largest.
#This helps you bring the most statistically significant genes (lowest padj) to the top.
#res$padj → extracts the “padj” column (adjusted p-values)

#==============================================================================================================================


#Save results
write.csv(as.data.frame(resOrdered), "DESeq2_pasilla_results.csv")

#======================================================================

# 7. Visualization

# Variance stablizing transformation

vsd <- vst(dds, blind = FALSE)

#PCA plot
plotPCA(vsd, intgroup= "condition")

#Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows= sampleDists,
         clustering_distance_cols= sampleDists,
         main= "Sample-to-Sample Distances")

#Volcano plot
EnhancedVolcano(res,
                lab= rownames(res),
                x= 'log2FoldChange',
                y= 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = 'Pasilla Knockdown vs Control')

#EXPLAINATION#=======================================================================================================================

#blind= FALSE- When blind = FALSE, it uses your experimental design (e.g., “condition”) to keep biological differences meaningful.

                