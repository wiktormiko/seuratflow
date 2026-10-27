# Small simulated dataset written to a temp directory
example_dir <- function(n_samples = 3, n_cells = 150, n_genes = 300) {
  d <- tempfile("example_10x_")
  write_example_10x(d, n_samples = n_samples, n_cells = n_cells, n_genes = n_genes)
  d
}

quiet <- function(x) suppressWarnings(suppressMessages(x))

# Write every sample of the example directory as a 10X-style h5 file (same layout that
# Read10X_h5() and Read_CellBender_h5_Mat() read), at the paths each reader expects.
write_example_h5 <- function(src, outs_dir, cellbender_dir) {
  staged <- seuratflow:::stage_10x_dirs(src)
  on.exit(unlink(staged, recursive = TRUE))
  write_h5 <- function(m, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    f <- hdf5r::H5File$new(path, mode = "w")
    on.exit(f$close_all())
    g <- f$create_group("matrix")
    g[["data"]] <- as.integer(m@x)
    g[["indices"]] <- as.integer(m@i)
    g[["indptr"]] <- as.integer(m@p)
    g[["shape"]] <- as.integer(dim(m))
    g[["barcodes"]] <- colnames(m)
    fg <- g$create_group("features")
    fg[["id"]] <- rownames(m)
    fg[["name"]] <- rownames(m)
    fg[["feature_type"]] <- rep("Gene Expression", nrow(m))
    fg[["genome"]] <- rep("GRCh38", nrow(m))
  }
  for (s in basename(list.dirs(src, recursive = FALSE))) {
    m <- methods::as(Seurat::Read10X(file.path(staged, s)), "CsparseMatrix")
    write_h5(m, file.path(outs_dir, s, "outs", "filtered_feature_bc_matrix.h5"))
    write_h5(m, file.path(cellbender_dir, s, paste0(s, "_out_filtered.h5")))
  }
  file.copy(file.path(src, "metadata.csv"), file.path(outs_dir, "metadata.csv"))
  file.copy(file.path(src, "metadata.csv"), file.path(cellbender_dir, "metadata.csv"))
}
