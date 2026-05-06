###########################################################
# Differential Gene Expression Analysis: Comparing Drought/Salt Stresses
# SAMPLE CODE
###########################################################

###########################################################
# Venn Diagram Analysis: Drought vs Salt DEGs
###########################################################

library(VennDiagram)
library(ggvenn)

# Get significant DEGs for each condition
drought_up <- rownames(subset(drought_results, padj < 0.05 & log2FoldChange > 1))
drought_down <- rownames(subset(drought_results, padj < 0.05 & log2FoldChange < -1))

salt_up <- rownames(subset(salt_results, padj < 0.05 & log2FoldChange > 1))
salt_down <- rownames(subset(salt_results, padj < 0.05 & log2FoldChange < -1))

# Create Venn diagram for UP-regulated genes
venn_up_list <- list(
  Drought = drought_up,
  Salt = salt_up
)

pdf("Venn_Upregulated_Drought_vs_Salt.pdf", width=8, height=8)
ggvenn(venn_up_list, 
       fill_color = c("#E41A1C", "#377EB8"),
       stroke_size = 1,
       set_name_size = 6,
       text_size = 5) +
  ggtitle("Up-regulated DEGs: Drought vs Salt Stress")
dev.off()

# Create Venn diagram for DOWN-regulated genes
venn_down_list <- list(
  Drought = drought_down,
  Salt = salt_down
)

pdf("Venn_Downregulated_Drought_vs_Salt.pdf", width=8, height=8)
ggvenn(venn_down_list, 
       fill_color = c("#E41A1C", "#377EB8"),
       stroke_size = 1,
       set_name_size = 6,
       text_size = 5) +
  ggtitle("Down-regulated DEGs: Drought vs Salt Stress")
dev.off()

# Calculate overlap statistics
overlap_up <- intersect(drought_up, salt_up)
overlap_down <- intersect(drought_down, salt_down)

drought_only_up <- setdiff(drought_up, salt_up)
salt_only_up <- setdiff(salt_up, drought_up)

drought_only_down <- setdiff(drought_down, salt_down)
salt_only_down <- setdiff(salt_down, drought_down)

# Print summary statistics
cat("=== UP-REGULATED GENES ===\n")
cat("Drought-responsive:", length(drought_up), "\n")
cat("Salt-responsive:", length(salt_up), "\n")
cat("Shared (both stresses):", length(overlap_up), 
    sprintf("(%.1f%%)", 100*length(overlap_up)/length(drought_up)), "\n")
cat("Drought-specific:", length(drought_only_up), "\n")
cat("Salt-specific:", length(salt_only_up), "\n\n")

cat("=== DOWN-REGULATED GENES ===\n")
cat("Drought-responsive:", length(drought_down), "\n")
cat("Salt-responsive:", length(salt_down), "\n")
cat("Shared (both stresses):", length(overlap_down), 
    sprintf("(%.1f%%)", 100*length(overlap_down)/length(drought_down)), "\n")
cat("Drought-specific:", length(drought_only_down), "\n")
cat("Salt-specific:", length(salt_only_down), "\n\n")

# Export gene lists for each category
write.csv(data.frame(gene_id = overlap_up), "Shared_Upregulated_Genes.csv", row.names=FALSE)
write.csv(data.frame(gene_id = overlap_down), "Shared_Downregulated_Genes.csv", row.names=FALSE)
write.csv(data.frame(gene_id = drought_only_up), "Drought_Specific_Upregulated.csv", row.names=FALSE)
write.csv(data.frame(gene_id = salt_only_up), "Salt_Specific_Upregulated.csv", row.names=FALSE)

###########################################################
# Scatter Plot: Drought vs Salt Log2FC Comparison
###########################################################

library(ggplot2)

# Merge drought and salt results
merged_results <- merge(
  drought_results[, c("log2FoldChange", "padj")],
  salt_results[, c("log2FoldChange", "padj")],
  by = "row.names",
  suffixes = c("_drought", "_salt")
)
rownames(merged_results) <- merged_results$Row.names
merged_results$Row.names <- NULL

# Classify genes
merged_results$category <- "Not Significant"
merged_results$category[merged_results$padj_drought < 0.05 & 
                          merged_results$padj_salt < 0.05] <- "Both Significant"
merged_results$category[merged_results$padj_drought < 0.05 & 
                          merged_results$padj_salt >= 0.05] <- "Drought Only"
merged_results$category[merged_results$padj_drought >= 0.05 & 
                          merged_results$padj_salt < 0.05] <- "Salt Only"

# Create scatter plot
pdf("Scatterplot_Drought_vs_Salt_FC.pdf", width=10, height=8)
ggplot(merged_results, aes(x = log2FoldChange_drought, 
                           y = log2FoldChange_salt,
                           color = category)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = c(
    "Both Significant" = "#984EA3",
    "Drought Only" = "#E41A1C", 
    "Salt Only" = "#377EB8",
    "Not Significant" = "gray80"
  )) +
  xlim(-8, 8) + ylim(-8, 8) +
  labs(
    title = "Fold Change Comparison: Drought vs Salt Stress",
    x = "Log2 Fold Change (Drought)",
    y = "Log2 Fold Change (Salt)",
    color = "Significance Category"
  ) +
  theme_bw() +
  theme(legend.position = "right")
dev.off()

# Calculate correlation
cor_all <- cor(merged_results$log2FoldChange_drought, 
               merged_results$log2FoldChange_salt,
               use = "complete.obs")

sig_both <- merged_results[merged_results$category == "Both Significant", ]
cor_sig <- cor(sig_both$log2FoldChange_drought, 
               sig_both$log2FoldChange_salt,
               use = "complete.obs")

cat("Correlation (all genes):", round(cor_all, 3), "\n")
cat("Correlation (significant in both):", round(cor_sig, 3), "\n")

###########################################################
# Side-by-Side Volcano Plots
###########################################################

library(EnhancedVolcano)
library(gridExtra)

# Drought volcano plot
p1 <- EnhancedVolcano(drought_results,
                      lab = drought_results$symbol,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = 'Drought Stress',
                      subtitle = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize = 1.5,
                      labSize = 3.0,
                      col = c('gray30', 'forestgreen', 'royalblue', 'red2'),
                      legendPosition = 'none')

# Salt volcano plot
p2 <- EnhancedVolcano(salt_results,
                      lab = salt_results$symbol,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = 'Salt Stress',
                      subtitle = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      pointSize = 1.5,
                      labSize = 3.0,
                      col = c('gray30', 'forestgreen', 'royalblue', 'red2'))

# Combine side by side
pdf("Volcano_Comparison_Drought_vs_Salt.pdf", width=16, height=8)
grid.arrange(p1, p2, ncol=2)
dev.off()

###########################################################
# GO Enrichment Comparison: Shared vs Specific
###########################################################

library(clusterProfiler)
library(ggplot2)

# Run GO enrichment for each category
# Shared up-regulated
shared_up_ids <- sub("\\.\\d+$", "", overlap_up)
ego_shared_up <- enrichGO(gene = shared_up_ids,
                          OrgDb = org.At.tair.db,
                          keyType = "TAIR",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# Drought-specific up-regulated
drought_only_ids <- sub("\\.\\d+$", "", drought_only_up)
ego_drought_up <- enrichGO(gene = drought_only_ids,
                           OrgDb = org.At.tair.db,
                           keyType = "TAIR",
                           ont = "BP",
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05)

# Salt-specific up-regulated
salt_only_ids <- sub("\\.\\d+$", "", salt_only_up)
ego_salt_up <- enrichGO(gene = salt_only_ids,
                        OrgDb = org.At.tair.db,
                        keyType = "TAIR",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

# Create comparison list
ego_list <- list(
  "Shared\nUp-regulated" = ego_shared_up,
  "Drought-Specific\nUp-regulated" = ego_drought_up,
  "Salt-Specific\nUp-regulated" = ego_salt_up
)

# Compare clusters
library(enrichplot)
pdf("GO_Comparison_Shared_vs_Specific.pdf", width=14, height=10)

# Method 1: Side-by-side dotplots
p1 <- dotplot(ego_shared_up, showCategory=10, title="Shared: Both Stresses")
p2 <- dotplot(ego_drought_up, showCategory=10, title="Drought-Specific")
p3 <- dotplot(ego_salt_up, showCategory=10, title="Salt-Specific")

grid.arrange(p1, p2, p3, ncol=3)
dev.off()

# Method 2: Compare using compareCluster (if same gene universe)
# This creates a unified comparison plot
gene_categories <- list(
  Shared = shared_up_ids,
  Drought_Only = drought_only_ids,
  Salt_Only = salt_only_ids
)

comparison <- compareCluster(geneClusters = gene_categories,
                             fun = "enrichGO",
                             OrgDb = org.At.tair.db,
                             keyType = "TAIR",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05)

pdf("GO_CompareCluster_Analysis.pdf", width=12, height=10)
dotplot(comparison, showCategory=15) + 
  ggtitle("GO Enrichment: Shared vs Stress-Specific Genes")
dev.off()

###########################################################
# Summary Bar Chart: DEG Counts Across Categories
###########################################################

library(ggplot2)

# Create summary data frame
summary_data <- data.frame(
  Category = rep(c("Shared", "Drought-Specific", "Salt-Specific"), 2),
  Direction = rep(c("Up-regulated", "Down-regulated"), each = 3),
  Count = c(
    length(overlap_up), length(drought_only_up), length(salt_only_up),
    length(overlap_down), length(drought_only_down), length(salt_only_down)
  )
)

# Create stacked bar chart
pdf("Summary_DEG_Counts.pdf", width=10, height=6)
ggplot(summary_data, aes(x = Category, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Up-regulated" = "#E41A1C", 
                               "Down-regulated" = "#377EB8")) +
  labs(
    title = "Differentially Expressed Genes: Drought vs Salt Comparison",
    x = "Gene Category",
    y = "Number of DEGs",
    fill = "Regulation Direction"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "top")
dev.off()

# Create percentage breakdown
total_drought <- length(drought_up) + length(drought_down)
total_salt <- length(salt_up) + length(salt_down)

pct_data <- data.frame(
  Stress = rep(c("Drought", "Salt"), each = 2),
  Category = rep(c("Shared", "Specific"), 2),
  Percentage = c(
    100 * length(overlap_up) / length(drought_up),
    100 * length(drought_only_up) / length(drought_up),
    100 * length(overlap_up) / length(salt_up),
    100 * length(salt_only_up) / length(salt_up)
  )
)

pdf("Percentage_Overlap.pdf", width=8, height=6)
ggplot(pct_data, aes(x = Stress, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_stack(vjust = 0.5),
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Shared" = "#984EA3", "Specific" = "#FF7F00")) +
  labs(
    title = "Percentage of Shared vs Specific Responses",
    subtitle = "Up-regulated genes only",
    x = "Stress Type",
    y = "Percentage of DEGs (%)"
  ) +
  theme_bw() +
  theme(legend.position = "top")
dev.off()

###########################################################
# Summary Statistics Table
###########################################################

summary_table <- data.frame(
  Metric = c(
    "Total DEGs (up)",
    "Total DEGs (down)",
    "Total DEGs (both directions)",
    "Shared DEGs (up)",
    "Shared DEGs (down)",
    "Stress-specific DEGs (up)",
    "Stress-specific DEGs (down)",
    "Percentage overlap (up)",
    "Percentage overlap (down)",
    "Top enriched GO term (shared)",
    "Top enriched GO term (drought-specific)",
    "Top enriched GO term (salt-specific)"
  ),
  Drought = c(
    length(drought_up),
    length(drought_down),
    length(drought_up) + length(drought_down),
    length(overlap_up),
    length(overlap_down),
    length(drought_only_up),
    length(drought_only_down),
    sprintf("%.1f%%", 100*length(overlap_up)/length(drought_up)),
    sprintf("%.1f%%", 100*length(overlap_down)/length(drought_down)),
    "Response to ABA",
    "Water transport",
    "-"
  ),
  Salt = c(
    length(salt_up),
    length(salt_down),
    length(salt_up) + length(salt_down),
    length(overlap_up),
    length(overlap_down),
    length(salt_only_up),
    length(salt_only_down),
    sprintf("%.1f%%", 100*length(overlap_up)/length(salt_up)),
    sprintf("%.1f%%", 100*length(overlap_down)/length(salt_down)),
    "Response to ABA",
    "-",
    "Sodium ion transport"
  )
)

write.csv(summary_table, "Comparative_Summary_Table.csv", row.names=FALSE)
print(summary_table)

###########################################################
# Fisher's Exact Test: Is overlap significant?
###########################################################

# Create contingency table
# Example for up-regulated genes

#                  In Salt    Not in Salt
# In Drought         a            b
# Not in Drought     c            d

total_genes <- nrow(drought_results)  # Genes tested in both

a <- length(overlap_up)  # In both
b <- length(drought_only_up)  # Drought only
c <- length(salt_only_up)  # Salt only
d <- total_genes - a - b - c  # Neither

contingency <- matrix(c(a, b, c, d), nrow = 2)

fisher_result <- fisher.test(contingency)
# p-value = P(seeing >= 2,400 shared genes | null hypothesis is true)

cat("Fisher's Exact Test for Overlap:\n")
cat("Odds Ratio:", fisher_result$estimate, "\n")
cat("P-value:", fisher_result$p.value, "\n")

if (fisher_result$p.value < 0.001) {
  cat("Conclusion: Overlap is HIGHLY significant\n")
  cat("Drought and salt responses are significantly correlated\n")
}

#########  Example Contingency Table: ######## 
#                   Salt-Responsive    Salt Non-Responsive    Total
# Drought-Responsive      2,400              2,100            4,500
# Not Drought-Resp        1,400             14,100           15,500
# Total                   3,800             16,200           20,000
# 
#### Random Chance: #### 
# Expected_a = (Total_Drought × Total_Salt) / Total_Genes
# Expected_a = (4,500 × 3,800) / 20,000
#
# (4500/20000) * (3800/20000)
# Expected_a = 855 genes
#
#### Idea we're testing: #### 
# p-value = P(seeing ≥ 2,400 shared genes | null hypothesis is true)
# 
####  Sample output intepretation: #### 
# "Of 4,500 drought-responsive genes, 2,400 (53.3%) also responded to salt stress. 
# Fisher's exact test revealed this overlap was highly significant
# (p < 2.2e-16, odds ratio = 7.85, 95% CI: 7.32-8.41), representing a 2.8-fold enrichment
# over the 855 genes expected by random chance alone. This indicates that drought and
# salinity stress responses share substantial molecular mechanisms in Arabidopsis roots."

