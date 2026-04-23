# ==============================

# BayesPrism: Load Data Script

# ==============================

# ---- 0. Load libraries ----

suppressWarnings(library(BayesPrism))

# ---- 1. File paths (EDIT THESE) ----

bulk_file <- "bulk_counts.tsv"
sc_file   <- "sc_counts.tsv"
meta_file <- "sc_metadata.tsv"

# ---- 2. Read data ----

bk.dat <- read.table(bulk_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
sc.dat <- read.table(sc_file,   header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta   <- read.table(meta_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

# ---- 3. Convert to matrices ----

bk.dat <- as.matrix(bk.dat)
sc.dat <- as.matrix(sc.dat)

# ---- 4. Align metadata with sc.dat ----

meta <- meta[rownames(sc.dat), , drop=FALSE]

# ---- 5. Extract labels (EDIT COLUMN NAMES IF NEEDED) ----

cell.type.labels  <- meta$cell_type
cell.state.labels <- meta$cell_state

# ---- 6. Basic sanity checks ----

if (nrow(sc.dat) != length(cell.type.labels)) {
stop("Mismatch: sc.dat rows != cell.type.labels length")
}
if (nrow(sc.dat) != length(cell.state.labels)) {
stop("Mismatch: sc.dat rows != cell.state.labels length")
}

# ---- 7. Check gene overlap ----

common_genes <- intersect(colnames(bk.dat), colnames(sc.dat))
cat("Number of shared genes:", length(common_genes), "\n")

if (length(common_genes) < 1000) {
warning("Low gene overlap — check gene naming consistency")
}

# ---- 8. Optional: subset to shared genes ----

bk.dat <- bk.dat[, common_genes, drop=FALSE]
sc.dat <- sc.dat[, common_genes, drop=FALSE]

# ---- 9. Optional: fallback if no cell states ----

if (is.null(cell.state.labels) || all(is.na(cell.state.labels))) {
cell.state.labels <- cell.type.labels
message("Using cell.type.labels as cell.state.labels")
}

# ---- 10. Create BayesPrism object ----

bp <- new.prism(
reference = sc.dat,
mixture = bk.dat,
input.type = "count.matrix",
cell.type.labels = cell.type.labels,
cell.state.labels = cell.state.labels,
key = NULL
)

# ---- 11. Done ----

cat("BayesPrism object created successfully\n")

# Optional: save object

# saveRDS(bp, file = "bp_object.rds")
