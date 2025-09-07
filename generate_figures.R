# Title: R script to generate figures for 'Identifying Alzheimer’s disease-related pathways based on whole-genome sequencing data'
# Author: Yongheng Wang
# Date: 2025-09-06
# Description: This script loads processed data to generate Figure 2C, Figure 3A, 
#              and the SEM path diagram for the manuscript.

# --- 1. Environment Setup ---

# List of required packages
packages <- c("ggplot2", "ggpubr", "semPlot", "lisrelToR")
bioc_packages <- c("DESeq2")

# Install CRAN packages if they are not already installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Install Bioconductor packages if they are not already installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
installed_bioc <- bioc_packages %in% rownames(installed.packages())
if (any(installed_bioc == FALSE)) {
  BiocManager::install(bioc_packages[!installed_bioc])
}

# Load all required packages
invisible(lapply(c(packages, bioc_packages), library, character.only = TRUE))

cat("All packages installed and loaded successfully.\n")


# --- 2. Define Paths and Create Directories ---

# Assumes the script is run from the root of the project directory.
# Project structure should be:
# ./
# |- R/generate_figures.R
# |- data/
# |- figures/

data_dir <- "data"
output_dir <- "figures"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# --- 3. Generate Figure 2C ---
# Description: Bar plot of Spearman correlation between predicted and observed RNA-seq expression.

cat("Generating Figure 2C...\n")

# Load data
expr_matrix_counts <- read.table(file.path(data_dir, "MayoRNAseq_RNAseq_CBE_geneCounts.tsv"), header = TRUE, row.names = 1, quote = "\t", check.names = FALSE)
expr_matrix_predict <- read.table(file.path(data_dir, "UTMOST_Brain_Cerebellum.MayoRNA.predict.txt"), header = TRUE, sep = "\t", row.names = 1)

# Normalize RNA-seq counts using VST from DESeq2
dds <- DESeqDataSetFromMatrix(countData = expr_matrix_counts, colData = data.frame(row.names = colnames(expr_matrix_counts)), design = ~1)
vst_data <- vst(dds, blind = TRUE)
expr_matrix_normalized_vst <- assay(vst_data)

# Process predicted expression matrix
expr_matrix_predict_t <- t(expr_matrix_predict[, -1]) # Genes in rows, samples in columns
samples_name <- sapply(colnames(expr_matrix_predict_t), function(x) strsplit(x, "_")[[1]][1])
colnames(expr_matrix_predict_t) <- paste(samples_name, "CER", sep = "_")

# Find common samples and genes
common_samples <- intersect(colnames(expr_matrix_predict_t), colnames(expr_matrix_normalized_vst))
common_genes <- intersect(rownames(expr_matrix_predict_t), rownames(expr_matrix_normalized_vst))

expr_matrix_RNA <- expr_matrix_normalized_vst[common_genes, common_samples]
expr_matrix_UTMOST <- expr_matrix_predict_t[common_genes, common_samples]

# Calculate per-sample Spearman correlation
if (!all(colnames(expr_matrix_RNA) == colnames(expr_matrix_UTMOST))) {
  stop("Sample names or order do not match between matrices!")
}

sample_correlation_list <- sapply(colnames(expr_matrix_RNA), function(sample_name) {
  cor(expr_matrix_RNA[, sample_name], expr_matrix_UTMOST[, sample_name], method = "spearman")
}, USE.NAMES = TRUE)

correlation_df <- data.frame(
  Sample = names(sample_correlation_list),
  Spearman_Correlation = sample_correlation_list
)

# Create and save plot
p_bar <- ggplot(correlation_df, aes(x = Sample, y = Spearman_Correlation)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", fill = "skyblue") +
  labs(title = "Spearman Correlation: Predicted vs. RNA-seq Expression (Cerebellum)",
       x = "Sample Name", y = "Spearman Correlation Coefficient (r)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")

ggsave(filename = file.path(output_dir, "Fig2C_cerebellum_correlation_barplot.pdf"),
       plot = p_bar, width = 12, height = 7)

cat("Figure 2C saved successfully.\n")


# --- 4. Generate Figure 3A ---
# Description: Boxplots showing correlation between Polygenic Risk Score (PRS) and clinical phenotypes.

cat("Generating Figure 3A...\n")

# Load data
PRS.score.pheno <- read.table(file.path(data_dir, "MayoRNA_PRS_cov.txt"), sep = "\t", header = TRUE, row.names = NULL)

# Data wrangling
PRS.score.pheno$group <- as.factor(ifelse(PRS.score.pheno$PHENO == "1", "Control", "Case"))
PRS.score.pheno$Thal <- as.factor(PRS.score.pheno$Thal)
PRS.score.pheno$Braak <- as.factor(PRS.score.pheno$Braak)

# Create the three individual plots
p_pheno <- ggplot(PRS.score.pheno, aes(x = group, y = PRScore, color = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 0.7) +
  labs(x = "Diagnosis", y = "Polygenic Risk Score") +
  theme_bw() +
  theme(legend.position = "none") +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top")

p_thal <- ggplot(PRS.score.pheno, aes(x = Thal, y = PRScore, color = Thal)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 0.7) +
  labs(x = "Thal Stage", y = "Polygenic Risk Score") +
  theme_bw() +
  theme(legend.position = "none") +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top")

p_braak <- ggplot(PRS.score.pheno, aes(x = Braak, y = PRScore, color = Braak)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(x = "Braak Stage", y = "Polygenic Risk Score") +
  theme_bw() +
  theme(legend.position = "none") +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top")

# Combine plots into a single figure
fig3a_combined <- ggarrange(p_pheno, p_thal, p_braak, ncol = 3, nrow = 1, labels = "AUTO")

# Save the combined plot
ggsave(filename = file.path(output_dir, "Fig3A_PRS_phenotype_correlations.pdf"),
       plot = fig3a_combined, width = 12, height = 5)

cat("Figure 3A saved successfully.\n")


# --- 5. Generate SEM Figure (Cerebellum) ---
# Description: Structural Equation Modeling (SEM) path diagram for the cerebellum region.

cat("Generating SEM figure...\n")

region <- "cerebellum"
lisrel_output_file <- file.path(data_dir, paste0("SEM_lisrel(", region, ").OUT"))

# Check if the input file exists
if (!file.exists(lisrel_output_file)) {
  warning(paste("LISREL output file not found:", lisrel_output_file, "- Skipping SEM figure generation."))
} else {
  # Save plot to PDF
  pdf(file = file.path(output_dir, paste0("SEM_diagram_", region, ".pdf")), width = 8, height = 8)
  semPaths(lisrel_output_file,
           what = "std", # Use standardized estimates
           whatLabels = "est",
           layout = "tree2",  #tree,circle,spring,tree2,circle2,
           rotation = 2,  #1, 2, 3 and 4 : top, left side, bottom and right side
           edge.label.cex = 0.8,
           nCharNodes = 0, # Allow for longer node names
           sizeMan = 8,    # Manifest variable size
           sizeLat = 10,   # Latent variable size
           style = "lisrel",
           color = "lightblue",
           border.color = "black",
           intAtSide = TRUE,
           ask = FALSE)
  dev.off()
  
  cat("SEM figure saved successfully.\n")
}

cat("\n--- Script finished. All figures have been generated in the '", output_dir, "' directory. ---\n")