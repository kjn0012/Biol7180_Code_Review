# ==========================================
# BayesPrism: Run script
# ==========================================

suppressWarnings(library(BayesPrism))

# ---- 1. Load your data ----
bulk_file <- "bulk_counts.tsv"
sc_file   <- "sc_counts.tsv"
meta_file <- "sc_metadata.tsv"

bk.dat <- read.table(bulk_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
sc.dat <- read.table(sc_file,   header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta   <- read.table(meta_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

bk.dat <- as.matrix(bk.dat)
sc.dat <- as.matrix(sc.dat)

# ---- 2. Align metadata ----
meta <- meta[rownames(sc.dat), , drop=FALSE]

# Edit these column names for your metadata
cell.type.labels  <- meta$cell_type
cell.state.labels <- meta$cell_state

# Optional fallback
# cell.state.labels <- cell.type.labels

# ---- 3. Choose reference matrix ----
# Replace this with your filtered matrix if available:
# reference.dat <- sc.dat.filtered.pc
# or:
# reference.dat <- sc.dat.filtered.pc.sig
reference.dat <- sc.dat

# ---- 4. Set malignant cell type ----
key.label <- "tumor"
# key.label <- NULL

# ---- 5. Construct prism object ----
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

# ---- 6. Run BayesPrism ----
bp.res <- run.prism(prism = myPrism)

# ---- 7. Save result ----
saveRDS(bp.res, file = "bp.res.rds")

cat("BayesPrism run completed and saved to bp.res.rds\n")
