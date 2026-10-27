#' Write a small simulated 10X dataset
#'
#' Creates one directory per sample in the old 10X mtx layout (`matrix.mtx.gz`,
#' `barcodes.tsv.gz`, `genes.tsv.gz`) plus a `metadata.csv` sample sheet with `sample` and
#' `condition` columns. Cells come from `n_groups` simulated populations and the last sample
#' has a higher library size, so Harmony has something to correct. Gene names include `MT-`,
#' `RPL`, lncRNA and hemoglobin genes so the QC metrics have something to find.
#'
#' It is meant for trying the pipeline and for the package tests, not for biology.
#'
#' @param path Output directory (created if needed).
#' @param n_samples Number of samples (`s1`, `s2`, ...).
#' @param n_cells Cells per sample.
#' @param n_genes Number of genes (at least 60).
#' @param n_groups Number of simulated cell populations.
#' @param seed Random seed.
#' @return `path`, invisibly.
#' @examples
#' d <- file.path(tempdir(), "example_10x")
#' write_example_10x(d, n_samples = 2, n_cells = 100, n_genes = 200)
#' list.files(d, recursive = TRUE)
#' \dontrun{
#' so <- read10X_multi_directory_add_QC(
#'   d, file.path(d, "metadata.csv"),
#'   default_10X_path = FALSE, detect_doublets = FALSE
#' )
#' }
#' @export
write_example_10x <- function(path, n_samples = 3, n_cells = 700, n_genes = 1500,
                              n_groups = 3, seed = 1) {
  stopifnot(n_genes >= 60, n_samples >= 1, n_cells >= 1)
  set.seed(seed)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  special <- c(paste0("MT-ND", 1:13), paste0("RPL", 1:30),
               "DDX11L16", "WASH7P", "MIR1302-2HG", "FAM138A", "LINC00115", "LINC01409",
               "HBB", "HBA1")
  genes <- c(special, paste0("GENE", seq_len(n_genes - length(special))))
  profiles <- matrix(stats::rexp(n_genes * n_groups, 20), n_genes, n_groups)
  samples <- paste0("s", seq_len(n_samples))

  gzip_copy <- function(from, to) {
    con_in <- file(from, "rb")
    con_out <- gzfile(to, "wb")
    on.exit({ close(con_in); close(con_out) })
    writeBin(readBin(con_in, "raw", file.size(from)), con_out)
  }

  for (i in seq_along(samples)) {
    s <- samples[i]
    grp <- sample.int(n_groups, n_cells, replace = TRUE)
    lib <- stats::rlnorm(n_cells, 0, 0.4) * (if (i == n_samples && n_samples > 1) 1.5 else 1)
    mu <- profiles[, grp, drop = FALSE] * rep(lib, each = n_genes) * 15
    m <- matrix(stats::rpois(n_genes * n_cells, mu), n_genes, n_cells)
    m[1:13, ] <- m[1:13, ] + matrix(stats::rpois(13 * n_cells, 0.3), 13, n_cells)
    m <- Matrix::Matrix(m, sparse = TRUE)

    d <- file.path(path, s)
    dir.create(d, showWarnings = FALSE)
    tmp <- tempfile(fileext = ".mtx")
    Matrix::writeMM(m, tmp)
    gzip_copy(tmp, file.path(d, "matrix.mtx.gz"))
    unlink(tmp)
    writeLines(paste0(s, "_BC", seq_len(n_cells)), gzfile(file.path(d, "barcodes.tsv.gz")))
    utils::write.table(data.frame(genes, genes), gzfile(file.path(d, "genes.tsv.gz")),
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  utils::write.csv(
    data.frame(sample = samples, condition = rep(c("a", "b"), length.out = n_samples)),
    file.path(path, "metadata.csv"), row.names = FALSE)
  invisible(path)
}
