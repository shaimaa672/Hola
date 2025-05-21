# miRNA pipeline


#Before running, ensure you have all required packages installed. Run this in the R console
#install.packages(c("ggplot2", "pheatmap", "ggrepel"))
# Install Bioconductor manager (if not already installed)
#if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# Install Bioconductor packages
#BiocManager::install(c("edgeR", "EnhancedVolcano"))

# Set working directory and load libraries
setwd("D:/Shaimaa/amg_dea")
library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(EnhancedVolcano)

# Helper function to create design matrix
create_design_matrix <- function(metadata, ref_level = "control") {
  metadata$treatment <- factor(metadata$treatment)
  metadata$treatment <- relevel(metadata$treatment, ref = ref_level)
  design <- model.matrix(~ treatment, data = metadata)
  rownames(design) <- rownames(metadata)
  return(design)
}

# Helper function to process DGE data
process_dge <- function(count_data, metadata, min_cpm = 20, min_samples = 1) {
  # Create DGEList object
  dge <- DGEList(counts = count_data, group = metadata$treatment)
  
  # Filter low counts
  keep <- rowSums(cpm(dge) > min_cpm) >= min_samples
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalize
  dge <- calcNormFactors(dge)
  return(dge)
}

# Helper function to perform DE analysis
perform_de_analysis <- function(dge, design, contrast_name, output_prefix) {
  # Fit GLM model
  fit <- glmFit(dge, design)
  
  # Create contrast
  contrast <- makeContrasts(contrasts = contrast_name, levels = design)
  
  # Perform LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Get results
  top_genes <- topTags(lrt, n = Inf)
  
  # Filter significant genes (FDR < 0.05 and |logFC| > 1.5)
  significant_genes <- top_genes$table[top_genes$table$FDR < 0.05 & abs(top_genes$table$logFC) > 1.5, ]
  
  # Save results
  write.csv(significant_genes, file = paste0(output_prefix, "_significant_genes.csv"))
  
  return(list(
    all_results = as.data.frame(lrt$table),
    significant_genes = significant_genes,
    fit = fit,
    lrt = lrt
  ))
}

# Helper function to create heatmap
create_heatmap <- function(dge, significant_genes, title, filename = NULL) {
  # Get significant gene names
  sig_genes <- rownames(significant_genes)
  
  # Subset and transform counts
  sig_counts <- dge$counts[sig_genes, , drop = FALSE]
  log_counts <- log2(sig_counts + 1)
  scaled_counts <- t(scale(t(log_counts)))
  
  # Create heatmap
  p <- pheatmap(
    scaled_counts,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = title,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    filename = filename
  )
  
  return(p)
}

# Helper function to create lollipop plot
create_lollipop_plot <- function(significant_genes, title) {
  plot_data <- data.frame(
    Gene = rownames(significant_genes),
    LogFC = significant_genes$logFC,
    FDR = significant_genes$FDR,
    NegLogFDR = -log10(significant_genes$FDR),
    Regulation = ifelse(significant_genes$logFC > 0, "Upregulated", "Downregulated"),
    stringsAsFactors = FALSE
  )
  
  # Order genes by LogFC
  plot_data$Gene <- factor(plot_data$Gene, levels = plot_data$Gene[order(plot_data$LogFC)])
  
  ggplot(plot_data, aes(x = Gene, y = LogFC, color = Regulation, size = NegLogFDR)) +
    geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = LogFC), size = 1) +
    geom_point() +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
    scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
    theme_minimal() +
    labs(title = title, x = "Gene", y = "Log2 Fold Change") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Helper function to create volcano plot
create_volcano_plot <- function(de_results, title) {
  # Calculate FDR if not present
  if (!"FDR" %in% colnames(de_results)) {
    de_results$FDR <- p.adjust(de_results$PValue, method = "BH")
  }
  
  # Create significance column
  de_results$significance <- ifelse(
    de_results$FDR < 0.05 & abs(de_results$logFC) > 1.5,
    "Significant",
    "Not Significant"
  )
  
  ggplot(de_results, aes(x = logFC, y = -log10(FDR), color = significance)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-log10 FDR", color = "Significance") +
    theme(legend.position = "right") +
    geom_text_repel(
      aes(label = ifelse(significance == "Significant", rownames(de_results), "")),
      size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 50
    )
}

# Main analysis function
run_analysis <- function(cell_line, data, metadata) {
  # Subset data for cell line
  cell_metadata <- metadata[metadata$cell_line == cell_line, ]
  cell_data <- data[, metadata$cell_line == cell_line]
  
  # Process DGE data
  dge <- process_dge(cell_data, cell_metadata)
  
  # Create design matrix
  design <- create_design_matrix(cell_metadata)
  
  # Ensure sample order matches
  design <- design[match(colnames(dge$counts), rownames(design)), , drop = FALSE]
  
  # Calculate normalized counts (CPM)
  normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
  write.csv(normalized_counts, file = paste0("normalized_counts_", cell_line, ".csv"))
  
  # Perform DE analysis
  contrast_name <- paste0("treatment", levels(cell_metadata$treatment)[2])
  de_results <- perform_de_analysis(dge, design, contrast_name, cell_line)
  
  # Create visualizations
  create_heatmap(
    dge, de_results$significant_genes,
    paste("Heatmap of Significant Genes in", cell_line),
    paste0("heatmap_", cell_line, ".pdf")
  )
  
  print(create_lollipop_plot(
    de_results$significant_genes,
    paste("Lollipop Plot of Significant Genes in", cell_line)
  ))
  
  print(create_volcano_plot(
    de_results$all_results,
    paste("Volcano Plot of Genes in", cell_line)
  ))
  
  return(de_results)
}

# Load data
data <- as.matrix(read.csv("amg_final_counts.csv", row.names = 1))
metadata <- read.csv("metadata2.csv", row.names = 1)

# Run analysis for both cell lines
mcf7_results <- run_analysis("MCF7", data, metadata)
hepg2_results <- run_analysis("HepG2", data, metadata)

