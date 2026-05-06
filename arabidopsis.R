###########################################################
# Step 1. Install & Load Packages
###########################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# install.packages("RcppArmadillo", repos = "https://cloud.r-project.org") # Dependency for DESeq2
# 
# BiocManager::install(c("pheatmap", "EnhancedVolcano", 
#                        "org.At.tair.db", "AnnotationDbi",
#                        "apeglm"))
# BiocManager::install(version = "3.21") # DESeq2 won't cooperate with R version 4.5.1 naturally
# BiocManager::install("DESeq2", force = TRUE)

library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(org.At.tair.db)
library(AnnotationDbi)

###########################################################
# Step 2. Import Data
###########################################################
source("data_setup.R",
       local = FALSE,   # if TRUE, variables stay local to the source; if FALSE, go to global env
       echo = TRUE,     # print each line as it is executed (useful for debugging)
       max.deparse.length = Inf)

###########################################################
# Step 3. Create DESeq2 Object
###########################################################
dds <- DESeqDataSetFromMatrix(countData = counts_int,
                              colData = metadata,
                              design = ~ condition)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

###########################################################
# Step 4. Run DESeq2
###########################################################
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition","drought","control"))

# Shrink fold changes for better visualization
resLFC <- lfcShrink(dds, coef=2, type="apeglm")

#####
# Examining Warning Message
resultsNames(dds)
suppressWarnings(
  resLFC <- lfcShrink(dds, coef="condition_drought_vs_control", type="apeglm")
)
summary(resLFC)

sum(is.na(resLFC$log2FoldChange)) # 0
sum(is.infinite(resLFC$log2FoldChange)) # 0 
#####

# Order by adjusted p-value
#resOrdered <- resLFC[order(resLFC$padj), ] 

resOrdered <- as.data.frame(resLFC[order(resLFC$padj), ])
resOrdered$gene_id <- rownames(resLFC)[order(resLFC$padj)]
rownames(resOrdered) <- resOrdered$gene_id

head(resOrdered, 10)  # top genes

###########################################################
# Step 5. Annotation (TAIR gene IDs -> gene names)
###########################################################
if (nrow(resOrdered) > 0 && !is.null(rownames(resOrdered))) {
  possible_keytypes <- keytypes(org.At.tair.db)
  
  keytype <- if ("TAIR" %in% possible_keytypes) "TAIR" else "ENTREZID"
  rownames_fixed <- sub("\\.\\d+$", "", rownames(resOrdered))  # remove .1, .2, etc.
  
  resOrdered$symbol <- mapIds(org.At.tair.db,
                              keys = rownames_fixed,
                              column = "SYMBOL",
                              keytype = keytype,
                              multiVals = "first")
} else {
  message("No genes to annotate — resOrdered is empty or has no rownames")
}

keytypes(org.At.tair.db)


###########################################################
# Step 6. Visualization
###########################################################

# (A) MA plot
plotMA(resLFC, ylim=c(-5,5))

# (B) Volcano plot
EnhancedVolcano(resOrdered,
                lab = resOrdered$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Arabidopsis Roots: Drought vs Control',
                pCutoff = 0.05,
                FCcutoff = 2)

# (C) Heatmap of top 30 DEGs
vsd <- vst(dds, blind=FALSE)
topGenes <- head(rownames(resOrdered), 30)
pheatmap(assay(vsd)[topGenes, ],
         annotation_col = metadata,
         show_rownames = TRUE,
         fontsize_row = 8,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

# (D) PCA plot
plotPCA(vsd, intgroup="condition")


# Suggestions from Dr. Brady
# tpm instead of estimated counts
  # DESeq2 is built to model raw count data using a negative binomial distribution
  # The normalization (for sequencing depth and RNA composition) is done inside DESeq2 via size factors
  # If you give it data that’s already normalized (like TPM), you break this assumption — DESeq2 can’t correctly model variance or dispersion anymore.
# tmm
  # It uses log-fold changes (M-values) between samples, trims the extremes, and computes a weighted mean.
  # Median of ratios method is used for DESeq
  # Both have the same goal: Make counts comparable across samples without distorting relative expression.
# t-test, ggplot2 instead for plots potentially
# something down the road, e-plant website, computable arabidopsis plant

