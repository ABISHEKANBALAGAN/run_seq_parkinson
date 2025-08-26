# ============================
# DESeq2 Workflow with Visualization
# ============================

# Load libraries
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

# ----------------------------
# 1. Set paths
# ----------------------------
counts_file  <- "counts_matrix.tsv"
coldata_file <- "coldata.csv"
results_dir  <- "E:/run_seq_parkinson/results/"

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ----------------------------
# 2. Load data
# ----------------------------
counts <- read.delim(counts_file, row.names=1, check.names=FALSE)
coldata <- read.csv(coldata_file, check.names=FALSE, stringsAsFactors=FALSE)

# Use sample_id as rownames
rownames(coldata) <- coldata$sample_id

# Ensure sample names match
if(!all(colnames(counts) == rownames(coldata))) stop("Column names of counts do not match row names of coldata")

# ----------------------------
# 3. Create DESeq2 dataset
# ----------------------------
coldata$group <- factor(paste(coldata$condition, coldata$timepoint, sep="_"))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = coldata,
  design    = ~ group
)

# ----------------------------
# 4. Run DESeq2
# ----------------------------
dds <- DESeq(dds)

# ----------------------------
# 5. Generate PCA plot
# ----------------------------
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA of samples")
ggsave(filename=file.path(results_dir, "PCA_plot.png"), plot=p, width=6, height=5)

# ----------------------------
# 6. Loop through all contrasts
# ----------------------------
groups <- levels(coldata$group)
combs <- combn(groups, 2, simplify=FALSE)

for (contrast_pair in combs) {
  res <- results(dds, contrast=c("group", contrast_pair[1], contrast_pair[2]))
  
  # Save all results
  all_file <- paste0(results_dir, "DESeq2_", contrast_pair[1], "_vs_", contrast_pair[2], ".csv")
  write.csv(as.data.frame(res), all_file)
  
  # Save significant genes
  sig_res <- subset(res, padj < 0.05 & !is.na(padj))
  sig_file <- paste0(results_dir, "DESeq2_significant_", contrast_pair[1], "_vs_", contrast_pair[2], ".csv")
  write.csv(as.data.frame(sig_res), sig_file)
  
  # MA-plot
  png(filename=paste0(results_dir, "MAplot_", contrast_pair[1], "_vs_", contrast_pair[2], ".png"), width=800, height=600)
  plotMA(res, main=paste0(contrast_pair[1], " vs ", contrast_pair[2]), ylim=c(-5,5))
  dev.off()

cat("DESeq2 analysis complete. All results saved in", results_dir, "\n")
