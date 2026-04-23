# ==========================================
# BayesPrism: Filter outlier genes script
# ==========================================

suppressWarnings(library(BayesPrism))

# ---- 1. Load your data ----
# Edit these file names/paths as needed
bulk_file <- "bulk_counts.tsv"
sc_file   <- "sc_counts.tsv"
meta_file <- "sc_metadata.tsv"

bk.dat <- read.table(bulk_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
sc.dat <- read.table(sc_file,   header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta   <- read.table(meta_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

bk.dat <- as.matrix(bk.dat)
sc.dat <- as.matrix(sc.dat)

# ---- 2. Align metadata to sc.dat ----
meta <- meta[rownames(sc.dat), , drop=FALSE]

# Edit these to match your metadata column names
cell.type.labels  <- meta$cell_type
cell.state.labels <- meta$cell_state

# Optional fallback if you do not have cell states yet
# cell.state.labels <- cell.type.labels

# ---- 3. Basic checks ----
if (nrow(sc.dat) != length(cell.type.labels)) {
  stop("Mismatch: nrow(sc.dat) != length(cell.type.labels)")
}
if (nrow(sc.dat) != length(cell.state.labels)) {
  stop("Mismatch: nrow(sc.dat) != length(cell.state.labels)")
}

cat("sc.dat dimensions:", dim(sc.dat), "\n")
cat("bk.dat dimensions:", dim(bk.dat), "\n")
cat("Shared genes:", length(intersect(colnames(sc.dat), colnames(bk.dat))), "\n\n")

# ---- 4. Visualize outlier genes from scRNA-seq ----
# species = "hs" for human, "mm" for mouse
sc.stat <- plot.scRNA.outlier(
  input = sc.dat,
  cell.type.labels = cell.type.labels,
  species = "hs",
  return.raw = TRUE
  # pdf.prefix = "mydata.sc.stat"
)

cat("\nHead of sc.stat:\n")
print(head(sc.stat))

# ---- 5. Visualize outlier genes from bulk RNA-seq ----
bk.stat <- plot.bulk.outlier(
  bulk.input = bk.dat,
  sc.input = sc.dat,
  cell.type.labels = cell.type.labels,
  species = "hs",
  return.raw = TRUE
  # pdf.prefix = "mydata.bk.stat"
)

cat("\nHead of bk.stat:\n")
print(head(bk.stat))

# ---- 6. Filter outlier genes from scRNA-seq ----
# Common groups from the tutorial:
# Rb, Mrp, other_Rb, chrM, MALAT1, chrX, chrY
# exp.cells = 5 removes genes expressed in fewer than 5 cells
sc.dat.filtered <- cleanup.genes(
  input = sc.dat,
  input.type = "count.matrix",
  species = "hs",
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells = 5
)

cat("\nFiltered sc.dat dimensions:\n")
print(dim(sc.dat.filtered))

# ---- 7. Check bulk vs scRNA concordance ----
# BayesPrism notes this function is for human data
plot.bulk.vs.sc(
  sc.input = sc.dat.filtered,
  bulk.input = bk.dat
  # pdf.prefix = "mydata.bk.vs.sc"
)

# ---- 8. Keep protein-coding genes ----
sc.dat.filtered.pc <- select.gene.type(
  sc.dat.filtered,
  gene.type = "protein_coding"
)

cat("\nProtein-coding filtered dimensions:\n")
print(dim(sc.dat.filtered.pc))

# ---- 9. Optional: select signature genes ----
# This can help if cell types are very similar or batch effects are strong.
# The tutorial filters genes first to reduce memory use.
sc.dat.for.de <- sc.dat[, colSums(sc.dat > 0) > 3, drop = FALSE]

diff.exp.stat <- get.exp.stat(
  sc.dat = sc.dat.for.de,
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  pseudo.count = 0.1,      # 0.1 for 10x-like data; larger for Smart-seq
  cell.count.cutoff = 50,
  n.cores = 1
)

sc.dat.filtered.pc.sig <- select.marker(
  sc.dat = sc.dat.filtered.pc,
  stat = diff.exp.stat,
  pval.max = 0.01,
  lfc.min = 0.1
)

cat("\nSignature-gene filtered dimensions:\n")
print(dim(sc.dat.filtered.pc.sig))

# ---- 10. Save outputs if desired ----
# saveRDS(sc.stat, "sc_outlier_stats.rds")
# saveRDS(bk.stat, "bulk_outlier_stats.rds")
# saveRDS(sc.dat.filtered, "sc_filtered.rds")
# saveRDS(sc.dat.filtered.pc, "sc_filtered_protein_coding.rds")
# saveRDS(sc.dat.filtered.pc.sig, "sc_filtered_protein_coding_signature.rds")

cat("\nOutlier filtering workflow completed.\n")
