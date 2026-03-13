
library(pheatmap)

df <- read.table("k562_features_to_TPMs", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set row names and keep only numeric columns
rownames(df) <- df$Feature
mat <- as.matrix(df[, c("TPM_K562", "TPM_HepG2")])

# Optional: log-transform because TPM values span a wide range
#mat_log <- log2(mat + 1)

# transpose
mat_t <- t(mat)

# keep first 30 cols
mat_t <- mat_t[, 1:50]

# Plot heatmap, keep values as they are (not log-transformed) to show actual expression levels
pheatmap(log2(mat_t + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Log2 TPM Heatmap (K562 features)")


df <- read.table("hepg2_features_to_TPMs", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set row names and keep only numeric columns
rownames(df) <- df$Feature
mat <- as.matrix(df[, c("TPM_K562", "TPM_HepG2")])

# Optional: log-transform because TPM values span a wide range
#mat_log <- log2(mat + 1)

# transpose
mat_t <- t(mat)

# keep first 30 cols
mat_t <- mat_t[, 1:50]

# Plot heatmap, keep values as they are (not log-transformed) to show actual expression levels
pheatmap(log2(mat_t + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Log2 TPM Heatmap (HepG2 features)")



df <- read.table("sknsh_features_to_TPMs", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set row names and keep only numeric columns
rownames(df) <- df$Feature
mat <- as.matrix(df[, c("TPM_K562", "TPM_HepG2")])

# Optional: log-transform because TPM values span a wide range
#mat_log <- log2(mat + 1)

# transpose
mat_t <- t(mat)

# keep first 30 cols
mat_t <- mat_t[, 1:50]

# Plot heatmap, keep values as they are (not log-transformed) to show actual expression levels
pheatmap(log2(mat_t + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Log2 TPM Heatmap (SK-N-SH features)")


df <- read.table("features_in_k562_and_hepg2_and_sknsh_TPMs", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set row names and keep only numeric columns
rownames(df) <- df$Feature
mat <- as.matrix(df[, c("TPM_K562", "TPM_HepG2")])

# Optional: log-transform because TPM values span a wide range
#mat_log <- log2(mat + 1)

# transpose
mat_t <- t(mat)

# keep first 30 cols
mat_t <- mat_t[, 1:50]

# Plot heatmap, keep values as they are (not log-transformed) to show actual expression levels
pheatmap(log2(mat_t + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Log2 TPM Heatmap (Features in all three)")


df <- read.table("features_unique_to_k562_TPMs", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set row names and keep only numeric columns
rownames(df) <- df$Feature
mat <- as.matrix(df[, c("TPM_K562", "TPM_HepG2")])

# Optional: log-transform because TPM values span a wide range
#mat_log <- log2(mat + 1)

# transpose
mat_t <- t(mat)

# Plot heatmap, keep values as they are (not log-transformed) to show actual expression levels
pheatmap(log2(mat_t + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Log2 TPM Heatmap (Features only in K562)")



df <- read.table("features_unique_to_hepg2_TPMs", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set row names and keep only numeric columns
rownames(df) <- df$Feature
mat <- as.matrix(df[, c("TPM_K562", "TPM_HepG2")])

# Optional: log-transform because TPM values span a wide range
#mat_log <- log2(mat + 1)

# transpose
mat_t <- t(mat)

# Plot heatmap, keep values as they are (not log-transformed) to show actual expression levels
pheatmap(log2(mat_t + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Log2 TPM Heatmap (Features only in HepG2)")


df <- read.table("features_unique_to_sknsh_TPMs", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set row names and keep only numeric columns
rownames(df) <- df$Feature
mat <- as.matrix(df[, c("TPM_K562", "TPM_HepG2")])

# Optional: log-transform because TPM values span a wide range
#mat_log <- log2(mat + 1)

# transpose
mat_t <- t(mat)

# Plot heatmap, keep values as they are (not log-transformed) to show actual expression levels
pheatmap(log2(mat_t + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Log2 TPM Heatmap (Features only in SK-N-SH)")