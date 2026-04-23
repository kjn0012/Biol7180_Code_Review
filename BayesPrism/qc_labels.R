# ==========================================
# BayesPrism QC of cell type and state labels
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

# ---- 2. Align metadata to scRNA-seq matrix ----
meta <- meta[rownames(sc.dat), , drop=FALSE]

# Edit these column names to match your metadata
cell.type.labels  <- meta$cell_type
cell.state.labels <- meta$cell_state

# ---- 3. Basic sanity checks ----
if (nrow(sc.dat) != length(cell.type.labels)) {
  stop("Mismatch: nrow(sc.dat) != length(cell.type.labels)")
}
if (nrow(sc.dat) != length(cell.state.labels)) {
  stop("Mismatch: nrow(sc.dat) != length(cell.state.labels)")
}
if (any(is.na(cell.type.labels))) {
  stop("cell.type.labels contains NA values")
}
if (any(is.na(cell.state.labels))) {
  stop("cell.state.labels contains NA values")
}

cat("sc.dat dimensions:", dim(sc.dat), "\n")
cat("bulk dimensions:", dim(bk.dat), "\n")
cat("Shared genes:", length(intersect(colnames(sc.dat), colnames(bk.dat))), "\n\n")

# ---- 4. Summarize label counts ----
cat("Cell type counts:\n")
print(sort(table(cell.type.labels)))

cat("\nCell state counts:\n")
print(sort(table(cell.state.labels)))

cat("\nCell state x cell type table:\n")
print(table(cbind.data.frame(cell.state.labels, cell.type.labels)))

# Optional warnings for small groups
small.types  <- names(which(table(cell.type.labels) < 20))
small.states <- names(which(table(cell.state.labels) < 20))

if (length(small.types) > 0) {
  warning("Some cell types have <20 cells: ", paste(small.types, collapse=", "))
}
if (length(small.states) > 0) {
  warning("Some cell states have <20 cells: ", paste(small.states, collapse=", "))
}

# ---- 5. Plot QC correlation heatmaps ----
# The tutorial uses plot.cor.phi on sc.dat with type/state labels. :contentReference[oaicite:1]{index=1}
# Their example uses larger text for cell types and smaller text/margins for cell states. :contentReference[oaicite:2]{index=2}

# Cell type correlation
plot.cor.phi(
  input = sc.dat,
  input.labels = cell.type.labels,
  title = "cell type correlation",
  cexRow = 0.5,
  cexCol = 0.5
  # pdf.prefix = "mydata.cor.ct"
)

# Cell state correlation
plot.cor.phi(
  input = sc.dat,
  input.labels = cell.state.labels,
  title = "cell state correlation",
  cexRow = 0.2,
  cexCol = 0.2,
  margins = c(2, 2)
  # pdf.prefix = "mydata.cor.cs"
)

# ---- 6. Optional: save plots directly to PDF ----
# Uncomment this section if you want files written automatically

# plot.cor.phi(
#   input = sc.dat,
#   input.labels = cell.type.labels,
#   title = "cell type correlation",
#   pdf.prefix = "mydata.cor.ct",
#   cexRow = 0.5,
#   cexCol = 0.5
# )
#
# plot.cor.phi(
#   input = sc.dat,
#   input.labels = cell.state.labels,
#   title = "cell state correlation",
#   pdf.prefix = "mydata.cor.cs",
#   cexRow = 0.2,
#   cexCol = 0.2,
#   margins = c(2, 2)
# )

cat("\nQC script completed.\n")
