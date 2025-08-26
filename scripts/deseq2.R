#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)

# ----------------------------
# Read data
# ----------------------------
counts <- read.delim("../results/counts_matrix.tsv", check.names=FALSE, stringsAsFactors=FALSE)
rownames(counts) <- counts[,1]; counts <- counts[,-1]; counts <- round(as.matrix(counts))

coldata <- read.csv("../results/coldata.csv", check.names=FALSE, stringsAsFactors=FALSE)
if(any(duplicated(coldata$sample_id))) coldata$sample_id <- make.unique(coldata$sample_id)
rownames(coldata) <- coldata$sample_id

# Ensure matching
stopifnot(all(colnames(counts) == rownames(coldata)))

# ----------------------------
# DESeq2
# ----------------------------
coldata$group <- factor(paste(coldata$condition, coldata$timepoint, sep="_"))
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~group)
counts(dds) <- round(counts(dds))
dds <- DESeq(dds)

# ----------------------------
# Extract results
# ----------------------------
results_dir <- "../results/"
dir.create(results_dir, showWarnings=FALSE)

contrasts <- combn(levels(coldata$group), 2, simplify=FALSE)
for (c in contrasts) {
  res <- results(dds, contrast=c("group", c[1], c[2]))
  write.csv(as.data.frame(res), file=paste0(results_dir, "DESeq2_", c[1], "_vs_", c[2], ".csv"))
}

# ----------------------------
# PCA
# ----------------------------
vsd <- vst(dds, blind=FALSE)
pca <- prcomp(t(assay(vsd)))
pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=coldata$group)
ggplot(pca_df, aes(PC1, PC2, color=group)) + geom_point(size=4) + theme_minimal()
ggsave(paste0(results_dir, "PCA_plot.png"))

# ----------------------------
# Volcano plot (example)
# ----------------------------
res <- results(dds, contrast=c("group","PD_d120","HC_d120"))
res$gene <- rownames(res)
res$significant <- ifelse(res$padj<0.05, "yes","no")
ggplot(res, aes(log2FoldChange, -log10(padj), color=significant)) +
  geom_point(alpha=0.5) + scale_color_manual(values=c("grey","red")) +
  theme_minimal() + ggtitle("Volcano Plot")
ggsave(paste0(results_dir,"Volcano_PD_d120_vs_HC_d120.png"))

# ----------------------------
# GO enrichment (example)
# ----------------------------
sig_genes <- rownames(subset(res, padj<0.05))
ego <- enrichGO(gene=sig_genes, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="BP", pAdjustMethod="BH")
write.csv(as.data.frame(ego), file=paste0(results_dir,"GO_BP_PD_d120_vs_HC_d120.csv"))
