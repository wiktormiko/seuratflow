#' Mirror a 10X directory with symlinks so old `genes.tsv.gz` files are readable
#'
#' Seurat expects `features.tsv.gz`. Instead of renaming files in the raw data, this builds a
#' temporary copy of `base_path` made of symlinks (`genes.tsv.gz` is linked as
#' `features.tsv.gz`). `base_path` is returned untouched when no `genes.tsv.gz` exists.
#' @noRd
stage_10x_dirs <- function(base_path) {
  files <- list.files(base_path, recursive = TRUE)
  if (!any(basename(files) == "genes.tsv.gz")) return(base_path)

  tmp <- tempfile("seurat_staged_")
  dir.create(tmp)
  for (f in files) {
    target <- file.path(tmp, sub("genes\\.tsv\\.gz$", "features.tsv.gz", f))
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    file.symlink(normalizePath(file.path(base_path, f)), target)
  }
  tmp
}

#' Parameter record stored in `@misc$pipeline`
#' @noRd
pipeline_info <- function(params) {
  list(params   = params)
}

#' Named vector of sample sub-directories
#' @noRd
.sample_dirs <- function(base_path) {
  d <- list.dirs(base_path, recursive = FALSE)
  stats::setNames(d, basename(d))
}

#' file.path() that tolerates an empty sub-directory
#' @noRd
.join_path <- function(dir, sub, file) {
  if (is.null(sub) || !nzchar(sub)) file.path(dir, file) else file.path(dir, sub, file)
}

#' Stop early, naming the samples, when input files are missing
#' @noRd
.check_inputs <- function(inputs, what) {
  nf <- inputs[!file.exists(inputs)]
  if (length(nf)) stop("Missing ", what, " for: ", paste(names(nf), collapse = ", "),
                       "\n  e.g. ", nf[[1]], call. = FALSE)
}

#' Run one output step; a failure becomes a warning
#' @noRd
.qc_step <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    warning("[", label, "] skipped: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

#' Write a csv and return a one-row manifest
#' @noRd
.csv_out <- function(df, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, file, row.names = FALSE)
  data.frame(file = file, width = NA_real_, height = NA_real_)
}
