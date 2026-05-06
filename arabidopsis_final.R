###########################################################
# Differential Gene Expression Analysis: Arabidopsis Drought Stress
# Complete Analysis Pipeline with DESeq2
###########################################################

###########################################################
# Step 1. Install & Load Packages
###########################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# install.packages("RcppArmadillo", repos = "https://cloud.r-project.org")
# 
# BiocManager::install(c("pheatmap", "EnhancedVolcano", 
#                        "org.At.tair.db", "AnnotationDbi",
#                        "apeglm", "clusterProfiler", "pathview"))
# BiocManager::install(version = "3.21")
# BiocManager::install("DESeq2", force = TRUE)
# BiocManager::install("clusterProfiler") # NEW!

library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(org.At.tair.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ggplot2)

###########################################################
# Step 2. Import Data
###########################################################
source("data_setup.R",
       local = FALSE,
       echo = TRUE,
       max.deparse.length = Inf)

###########################################################
# Step 3. Quality Control - Pre-filtering
###########################################################

# Visualize gene count distribution before filtering
pdf("QC_gene_count_distribution.pdf", width=8, height=6)
hist(log10(rowSums(counts_int) + 1), breaks=50,
     main="Gene Count Distribution (Before Filtering)", 
     xlab="log10(counts + 1)",
     col="skyblue")
abline(v=log10(10), col="red", lty=2, lwd=2)
legend("topright", legend="Filtering threshold (10 counts)", 
       col="red", lty=2, lwd=2)
dev.off()

# Summary statistics
cat("\n=== Pre-filtering Summary ===\n")
cat("Total genes:", nrow(counts_int), "\n")
cat("Genes with >10 total counts:", sum(rowSums(counts_int) > 10), "\n")
cat("Genes to be filtered out:", sum(rowSums(counts_int) <= 10), "\n\n")

###########################################################
# Step 4. Create DESeq2 Object
###########################################################
dds <- DESeqDataSetFromMatrix(countData = counts_int,
                              colData = metadata,
                              design = ~ condition)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

cat("Genes retained after filtering:", nrow(dds), "\n\n")

###########################################################
# Step 5. Sample Quality Control
###########################################################

# Transform data for visualizaztion (before running DESeq)
vsd <- vst(dds, blind=TRUE)

# (A) Sample correlation heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pdf("QC_sample_correlation.pdf", width=8, height=7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = metadata,
         main="Sample-to-Sample Distances")
dev.off()

# (B) PCA plot - initial quality check
pdf("QC_PCA_plot.pdf", width=8, height=6)
pca_plot <- plotPCA(vsd, intgroup="condition")
print(pca_plot + 
        ggtitle("PCA: Sample Clustering by Condition") +
        theme_bw())
dev.off()

###########################################################
# Step 6. Run DESeq2 Differential Expression Analysis
###########################################################
dds <- DESeq(dds)

# Check dispersion estimates
pdf("QC_dispersion_plot.pdf", width=8, height=6)
plotDispEsts(dds, main="Dispersion Estimates")
dev.off()

# Extract results
res <- results(dds, contrast = c("condition","drought","control"))

# Apply log-fold change shrinkage for better estimates
resLFC <- lfcShrink(dds, coef="condition_drought_vs_control", type="apeglm")

# Verify no NA or infinite values
cat("=== Checking result quality ===\n")
cat("NA log2FoldChange values:", sum(is.na(resLFC$log2FoldChange)), "\n")
cat("Infinite log2FoldChange values:", sum(is.infinite(resLFC$log2FoldChange)), "\n\n")

###########################################################
# Step 7. Statistical Quality Control
###########################################################

# P-value distribution
pdf("QC_pvalue_distribution.pdf", width=8, height=6)
hist(res$pvalue, breaks=50, col="lightblue",
     main="P-value Distribution", xlab="P-value")
dev.off()

# Independent filtering check
pdf("QC_independent_filtering.pdf", width=8, height=6)
plot(metadata(res)$filterNumRej, 
     type="b", ylab="Number of rejections",
     xlab="Quantiles of filter",
     main="Independent Filtering")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.off()

###########################################################
# Step 8. Prepare Results Table
###########################################################

# Order by adjusted p-value
resOrdered <- as.data.frame(resLFC[order(resLFC$padj), ])
resOrdered$gene_id <- rownames(resLFC)[order(resLFC$padj)]
rownames(resOrdered) <- resOrdered$gene_id

# Annotation: Add gene symbols
if (nrow(resOrdered) > 0 && !is.null(rownames(resOrdered))) {
  rownames_fixed <- sub("\\.\\d+$", "", rownames(resOrdered))
  
  resOrdered$symbol <- mapIds(org.At.tair.db,
                              keys = rownames_fixed,
                              column = "SYMBOL",
                              keytype = "TAIR",
                              multiVals = "first")
  
  resOrdered$genename <- mapIds(org.At.tair.db,
                                keys = rownames_fixed,
                                column = "GENENAME",
                                keytype = "TAIR",
                                multiVals = "first")
} else {
  message("No genes to annotate — resOrdered is empty")
}

# Summary statistics
cat("\n=== Differential Expression Summary ===\n")
summary_stats <- data.frame(
  Total_Genes_Tested = nrow(resOrdered),
  Significant_padj005 = sum(resOrdered$padj < 0.05, na.rm=TRUE),
  Upregulated_FC2 = sum(resOrdered$padj < 0.05 & resOrdered$log2FoldChange > 1, na.rm=TRUE),
  Downregulated_FC2 = sum(resOrdered$padj < 0.05 & resOrdered$log2FoldChange < -1, na.rm=TRUE),
  Upregulated_FC4 = sum(resOrdered$padj < 0.05 & resOrdered$log2FoldChange > 2, na.rm=TRUE),
  Downregulated_FC4 = sum(resOrdered$padj < 0.05 & resOrdered$log2FoldChange < -2, na.rm=TRUE)
)
print(summary_stats)
cat("\n")

# Display top results
cat("Top 10 differentially expressed genes:\n")
print(head(resOrdered[, c("symbol", "log2FoldChange", "padj")], 10))
cat("\n")

###########################################################
# Step 9. Core Visualizations
###########################################################

# Transform data for visualization (using final model)
vsd <- vst(dds, blind=FALSE)

# (A) MA plot
pdf("Results_MA_plot.pdf", width=8, height=6)
plotMA(resLFC, ylim=c(-5,5), 
       main="MA Plot: Drought vs Control")
dev.off()

# (B) Volcano plot
pdf("Results_volcano_plot.pdf", width=10, height=8)
EnhancedVolcano(resOrdered,
                lab = resOrdered$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Arabidopsis Roots: Drought vs Control',
                subtitle = 'Differential Expression Analysis',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 4.0)
dev.off()

# (C) Heatmap of top 30 DEGs
pdf("Results_heatmap_top30.pdf", width=8, height=10)
topGenes <- head(rownames(resOrdered), 30)
pheatmap(assay(vsd)[topGenes, ],
         annotation_col = metadata,
         show_rownames = TRUE,
         fontsize_row = 8,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main="Top 30 Differentially Expressed Genes")
dev.off()

# (D) PCA plot (final)
pdf("Results_PCA_final.pdf", width=8, height=6)
pca_final <- plotPCA(vsd, intgroup="condition")
print(pca_final + 
        ggtitle("PCA: Drought vs Control Separation") +
        theme_bw())
dev.off()

###########################################################
# Step 10. Expression Pattern Analysis
###########################################################

# K-means clustering of significant genes
sig_genes <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)

if (nrow(sig_genes) > 10) {
  mat <- assay(vsd)[rownames(sig_genes), ]
  mat_scaled <- t(scale(t(mat)))
  
  # Determine optimal number of clusters (using 4 as default)
  set.seed(42)
  k <- kmeans(mat_scaled, centers=4, nstart=25) # Talk about more on Thurs
  
  # Create cluster annotation
  cluster_anno <- data.frame(Cluster = factor(k$cluster))
  rownames(cluster_anno) <- rownames(mat_scaled)
  
  pdf("Results_expression_clusters.pdf", width=8, height=12)
  pheatmap(mat_scaled[order(k$cluster), ],
           cluster_rows=FALSE,
           show_rownames=FALSE,
           annotation_row=cluster_anno,
           annotation_col=metadata,
           main="Expression Patterns of Significant DEGs")
  dev.off()
  
  # Save cluster assignments
  sig_genes$cluster <- k$cluster[match(rownames(sig_genes), names(k$cluster))]
  
  cat("Expression clusters identified:\n")
  print(table(k$cluster))
  cat("\n")
}

###########################################################
# Step 11. Gene Ontology Enrichment Analysis
###########################################################

# Prepare gene lists
sig_up <- subset(resOrdered, padj < 0.05 & log2FoldChange > 1)
sig_down <- subset(resOrdered, padj < 0.05 & log2FoldChange < -1)

gene_list_up <- sub("\\.\\d+$", "", rownames(sig_up))
gene_list_down <- sub("\\.\\d+$", "", rownames(sig_down))

# GO enrichment for up-regulated genes
if (length(gene_list_up) > 5) {
  cat("Running GO enrichment for UP-regulated genes...\n")
  ego_up <- enrichGO(gene = gene_list_up,
                     OrgDb = org.At.tair.db,
                     keyType = "TAIR",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     readable = FALSE)
  
  if (!is.null(ego_up) && nrow(ego_up@result) > 0) {
    pdf("Enrichment_GO_upregulated.pdf", width=10, height=8)
    print(dotplot(ego_up, showCategory=20, title="GO Enrichment: Up-regulated Genes"))
    dev.off()
    
    write.csv(as.data.frame(ego_up), "GO_enrichment_upregulated.csv", row.names=FALSE)
    cat("Found", nrow(ego_up@result), "enriched GO terms in up-regulated genes\n")
  } else {
    cat("No significant GO enrichment for up-regulated genes\n")
  }
}

# GO enrichment for down-regulated genes
if (length(gene_list_down) > 5) {
  cat("Running GO enrichment for DOWN-regulated genes...\n")
  ego_down <- enrichGO(gene = gene_list_down,
                       OrgDb = org.At.tair.db,
                       keyType = "TAIR",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = FALSE)
  
  if (!is.null(ego_down) && nrow(ego_down@result) > 0) {
    pdf("Enrichment_GO_downregulated.pdf", width=10, height=8)
    print(dotplot(ego_down, showCategory=20, title="GO Enrichment: Down-regulated Genes"))
    dev.off()
    
    write.csv(as.data.frame(ego_down), "GO_enrichment_downregulated.csv", row.names=FALSE)
    cat("Found", nrow(ego_down@result), "enriched GO terms in down-regulated genes\n")
  } else {
    cat("No significant GO enrichment for down-regulated genes\n")
  }
}

###########################################################
# Step 12. KEGG Pathway Analysis
###########################################################

# Note: KEGG requires Entrez gene IDs - convert if needed
if (length(gene_list_up) > 5) {
  cat("\nRunning KEGG pathway analysis...\n")
  
  # Try KEGG enrichment (may need additional setup for Arabidopsis)
  tryCatch({
    kk_up <- enrichKEGG(gene = gene_list_up,
                        organism = 'ath',
                        keyType = 'kegg',
                        pvalueCutoff = 0.05)
    
    if (!is.null(kk_up) && nrow(kk_up@result) > 0) {
      pdf("Enrichment_KEGG_upregulated.pdf", width=10, height=8)
      print(dotplot(kk_up, showCategory=15, title="KEGG Pathways: Up-regulated Genes"))
      dev.off()
      
      write.csv(as.data.frame(kk_up), "KEGG_enrichment_upregulated.csv", row.names=FALSE)
      cat("Found", nrow(kk_up@result), "enriched KEGG pathways\n")
    }
  }, error = function(e) {
    cat("KEGG analysis skipped (may require additional gene ID conversion)\n")
  })
}

###########################################################
# Step 13. Export Results
###########################################################

# Export complete results table
write.csv(resOrdered, "DESeq2_results_complete.csv", row.names=TRUE)

# Export significant genes only
write.csv(subset(resOrdered, padj < 0.05), 
          "DESeq2_results_significant_padj005.csv", row.names=TRUE)

# Export up-regulated genes
write.csv(sig_up, "DESeq2_upregulated_genes.csv", row.names=TRUE)

# Export down-regulated genes
write.csv(sig_down, "DESeq2_downregulated_genes.csv", row.names=TRUE)

# Export gene lists for external tools (e.g., ePlant)
writeLines(rownames(subset(resOrdered, padj < 0.05)), 
           "significant_gene_list.txt")
writeLines(gene_list_up, "upregulated_gene_list.txt")
writeLines(gene_list_down, "downregulated_gene_list.txt")

# Export summary statistics
write.csv(summary_stats, "analysis_summary_statistics.csv", row.names=FALSE)

###########################################################
# Step 14. Session Info
###########################################################

cat("\n=== Analysis Complete ===\n")
cat("Results and plots saved to working directory\n\n")

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")

cat("Session info saved. Analysis complete!\n")
cat("\nKey output files:\n")
cat("  - DESeq2_results_complete.csv\n")
cat("  - DESeq2_results_significant_padj005.csv\n")
cat("  - Multiple QC and results PDF plots\n")
cat("  - GO and KEGG enrichment results (if available)\n")
cat("  - Gene lists for further analysis\n")

