args <- commandArgs(trailingOnly = TRUE)

expr_file <- args[1]
meta_file <- args[2]
out_file  <- args[3]

suppressPackageStartupMessages(library(limma))

expr_df <- read.csv(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
meta_df <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)

gene_col <- if ("gene_symbol" %in% names(expr_df)) {
  "gene_symbol"
} else if ("Gene_symbol" %in% names(expr_df)) {
  "Gene_symbol"
} else {
  stop("Expression file must contain a gene_symbol or Gene_symbol column.")
}

gene_symbols <- expr_df[[gene_col]]

expr_mat <- as.matrix(expr_df[, !(names(expr_df) %in% gene_col), drop = FALSE])
rownames(expr_mat) <- gene_symbols

meta_df <- meta_df[match(colnames(expr_mat), meta_df$sample_id), , drop = FALSE]

if (!all(colnames(expr_mat) == meta_df$sample_id)) {
  stop("Sample order mismatch between expression matrix and metadata.")
}

meta_df$response <- factor(meta_df$response, levels = c(0, 1))
meta_df$gi_involvement <- factor(meta_df$gi_involvement)

# Choose one design first
# simpler:
design <- model.matrix(~response, data = meta_df)
#design <- model.matrix(~ response + ses_cd, data = meta_df)
#design <- model.matrix(~ response + gi_involvement, data = meta_df)
colnames(design)
head(design)

fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

tt <- topTable(fit, coef = "response1", number = Inf, sort.by = "P")
tt_sig <- subset(tt, adj.P.Val < 0.05 & abs(logFC) > 0.5)

tt$Gene_symbol <- if ("ID" %in% names(tt)) tt$ID else rownames(tt)

write.csv(tt, out_file, row.names = FALSE)
