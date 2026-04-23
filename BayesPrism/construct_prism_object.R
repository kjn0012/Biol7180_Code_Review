# ==========================================
# BayesPrism: Construct prism object script
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

# Edit these column names to match your metadata
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

cat("Cell type counts:\n")
print(sort(table(cell.type.labels)))

cat("\nCell state counts:\n")
print(sort(table(cell.state.labels)))

# ---- 4. Optional preprocessing before prism construction ----
# If you already created sc.dat.filtered.pc in your previous filtering step,
# use that instead of raw sc.dat.
#
# Example:
# reference.dat <- sc.dat.filtered.pc
#
# Otherwise, this script uses raw sc.dat as a fallback:
reference.dat <- sc.dat

# ---- 5. Choose the malignant key ----
# Set key to the malignant cell type label in cell.type.labels, e.g. "tumor"
# Set to NULL if:
# - there is no malignant population, or
# - malignant cells are matched between reference and bulk
key.label <- "tumor"
# key.label <- NULL

# Optional safety check
if (!is.null(key.label) && !(key.label %in% unique(cell.type.labels))) {
  stop("The chosen key.label was not found in cell.type.labels")
}

# ---- 6. Construct the prism object ----
myPrism <- new.prism(
  reference = reference.dat,
  mixture = bk.dat,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = key.label,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

# ---- 7. Save prism object if desired ----
# saveRDS(myPrism, file = "myPrism.rds")

cat("\nPrism object created successfully.\n")
