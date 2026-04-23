# ==========================================
# BayesPrism: Downstream analysis
# ==========================================

suppressWarnings(library(BayesPrism))
suppressWarnings(library(DESeq2))

# ---- 1. Load BayesPrism result ----
bp.res <- readRDS("bp.res.rds")

# ---------------------------------
# 2. Extract theta
# ---------------------------------
get_theta_matrix <- function(bp.res) {
  possible <- list()

  try(possible[[length(possible) + 1]] <- bp.res@posterior.initial.cellType@theta, silent = TRUE)
  try(possible[[length(possible) + 1]] <- bp.res@posterior.theta_f, silent = TRUE)
  try(possible[[length(possible) + 1]] <- bp.res@theta, silent = TRUE)

  for (x in possible) {
    if (!is.null(x)) return(x)
  }

  stop("Could not find theta in bp.res. Run slotNames(bp.res) and str(bp.res, max.level = 2).")
}

theta <- get_theta_matrix(bp.res)

write.table(theta, file = "theta.tsv", sep = "\t", quote = FALSE, col.names = NA)
cat("Saved theta.tsv\n")

# ---------------------------------
# 3. Cluster bulk samples by theta
# ---------------------------------
theta.dist <- dist(theta)
theta.hc <- hclust(theta.dist, method = "ward.D2")

pdf("theta_clustering.pdf", width = 8, height = 6)
plot(theta.hc, main = "Bulk sample clustering by theta", xlab = "", sub = "")
dev.off()

cat("Saved theta_clustering.pdf\n")

# ---------------------------------
# 4. Extract Z
# ---------------------------------
get_Z_matrix <- function(bp.res) {
  possible <- list()

  try(possible[[length(possible) + 1]] <- bp.res@posterior.initial.cellState@Z, silent = TRUE)
  try(possible[[length(possible) + 1]] <- bp.res@Z, silent = TRUE)
  try(possible[[length(possible) + 1]] <- bp.res@posterior.Z, silent = TRUE)

  for (x in possible) {
    if (!is.null(x)) return(x)
  }

  return(NULL)
}

Z.raw <- get_Z_matrix(bp.res)

if (!is.null(Z.raw)) {
  Z.raw <- as.matrix(Z.raw)
  saveRDS(Z.raw, file = "Z_raw.rds")
  cat("Saved Z_raw.rds\n")
} else {
  cat("Z matrix not found in this BayesPrism result object\n")
}

# ---------------------------------
# 5. Normalize Z with vst
# ---------------------------------
Z.vst <- NULL

if (!is.null(Z.raw)) {
  Z.mat <- as.matrix(Z.raw)

  # Heuristic orientation check
  if (!is.null(rownames(Z.mat)) && nrow(Z.mat) < ncol(Z.mat)) {
    Z.sample_gene <- Z.mat
  } else {
    Z.sample_gene <- t(Z.mat)
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(t(Z.sample_gene)),
    colData = data.frame(row.names = rownames(Z.sample_gene), condition = "all"),
    design = ~ 1
  )

  Z.vst.obj <- vst(dds, blind = TRUE)
  Z.vst <- t(assay(Z.vst.obj))

  saveRDS(Z.vst, file = "Z_vst.rds")
  write.table(Z.vst, file = "Z_vst.tsv", sep = "\t", quote = FALSE, col.names = NA)

  cat("Saved Z_vst.rds and Z_vst.tsv\n")
}

# ---------------------------------
# 6. Cluster samples by Z
# ---------------------------------
if (!is.null(Z.vst)) {
  z.dist <- dist(Z.vst)
  z.hc <- hclust(z.dist, method = "ward.D2")

  pdf("Z_clustering.pdf", width = 8, height = 6)
  plot(z.hc, main = "Bulk sample clustering by Z (vst)", xlab = "", sub = "")
  dev.off()

  cat("Saved Z_clustering.pdf\n")
}

# ---------------------------------
# 7. Signature gene z-scores
# ---------------------------------
# Edit this list for your genes of interest
signature.genes <- c(
  # "GENE1", "GENE2", "GENE3"
)

if (!is.null(Z.vst) && length(signature.genes) > 0) {
  sig.genes.present <- intersect(signature.genes, colnames(Z.vst))

  if (length(sig.genes.present) > 0) {
    sig.mat <- Z.vst[, sig.genes.present, drop = FALSE]
    sig.z <- scale(sig.mat)
    signature.score <- rowMeans(sig.z, na.rm = TRUE)

    write.table(
      data.frame(sample = rownames(Z.vst), signature_score = signature.score),
      file = "signature_scores.tsv",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    cat("Saved signature_scores.tsv\n")
  } else {
    cat("No signature genes were found in Z_vst\n")
  }
}

# ---------------------------------
# 8. Correlate Z with theta
# ---------------------------------
# Edit this if your malignant label is different
tumor.label <- "tumor"

if (!is.null(Z.vst)) {
  nonmal.theta <- theta

  if (tumor.label %in% colnames(nonmal.theta)) {
    nonmal.theta <- nonmal.theta[, setdiff(colnames(nonmal.theta), tumor.label), drop = FALSE]
  }

  for (ct in colnames(nonmal.theta)) {
    ct.vec <- nonmal.theta[, ct]

    cor.df <- data.frame(
      gene = colnames(Z.vst),
      cor = apply(Z.vst, 2, function(g) cor(g, ct.vec, method = "spearman", use = "pairwise.complete.obs"))
    )

    write.table(
      cor.df,
      file = paste0("cor_Z_vs_", ct, ".tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  cat("Saved correlation tables between Z and theta\n")
}

# ---------------------------------
# 9. Export table for clinical analysis
# ---------------------------------
theta.df <- data.frame(sample = rownames(theta), theta, check.names = FALSE)

write.table(
  theta.df,
  file = "theta_for_clinical_analysis.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Saved theta_for_clinical_analysis.tsv\n")
cat("Downstream analysis completed.\n")
