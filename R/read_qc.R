#' Shared engine of the three readers
#'
#' `inputs` is a named character vector (names = sample names, values = what `read_fun()`
#' receives); `read_fun(input)` returns a genes x cells count matrix.
#' @noRd
.read_add_qc <- function(inputs, read_fun, sample_sheet_path, sample_sheet_column, species,
                         detect_doublets, nmads, mito_floor, seed, workers) {
  sample_sheet <- if (tools::file_ext(sample_sheet_path) == "xlsx") {
    read_xlsx(sample_sheet_path)
  } else {
    read_csv(sample_sheet_path, show_col_types = FALSE)
  }
  stopifnot(sample_sheet_column %in% names(sample_sheet),
            !anyDuplicated(sample_sheet[[sample_sheet_column]]))

  missing <- setdiff(names(inputs), sample_sheet[[sample_sheet_column]])
  if (length(missing)) stop("Samples not in sample sheet: ", paste(missing, collapse = ", "))

  if (!is.null(workers)) {
    old_plan <- if (workers > 1) future::plan(future::multisession, workers = workers) else
      future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  # Read + QC for one sample (runs in a worker)
  qc_one <- function(input, s) {
    m <- read_fun(input)

    so <- CreateSeuratObject(m, project = s, min.features = 1)
    so$sample <- s

    # sample-level metadata repeated for every cell
    md <- sample_sheet[sample_sheet[[sample_sheet_column]] == s, ] |>
      dplyr::select(-dplyr::all_of(sample_sheet_column)) |> as.data.frame()
    md <- md[rep(1, ncol(so)), , drop = FALSE]
    rownames(md) <- colnames(so)
    so <- AddMetaData(so, md)

    # QC metrics
    so$log10_genes <- log10(so$nFeature_RNA)
    so$log10_UMI   <- log10(so$nCount_RNA)
    # cell cycle normalizes each sample separately and errors when the marker genes are absent;
    # if that happens, warn and rerun without it instead of failing the whole run
    qc_args <- list(species = species, add_MSigDB = FALSE,
                    add_IEG = FALSE, add_IEG_module_score = FALSE,
                    add_cell_cycle = TRUE, add_lncRNA = FALSE, add_top_pct = TRUE)
    so <- tryCatch(
      do.call(Add_Cell_QC_Metrics, c(list(so), qc_args)),
      error = function(e) {
        warning("[", s, "] cell cycle scoring failed (", conditionMessage(e),
                "); continuing without it", call. = FALSE)
        do.call(Add_Cell_QC_Metrics, c(list(so), modifyList(qc_args, list(add_cell_cycle = FALSE))))
      })

    # MAD outliers (keep thresholds so they can be reported / plotted)
    umi <- isOutlier(so$log10_UMI,    type = "both",   nmads = nmads)
    gen <- isOutlier(so$log10_genes,  type = "both",   nmads = nmads)
    mt  <- isOutlier(so$percent_mito, type = "higher", nmads = nmads)
    so$umi_outlier    <- as.logical(umi)
    so$ngenes_outlier <- as.logical(gen)
    so$mt_outlier     <- as.logical(mt) & so$percent_mito > mito_floor
    so@misc$qc_thresholds <- data.frame(
      sample       = s,
      UMI_low      = round(10^attr(umi, "thresholds")[["lower"]]),
      UMI_high     = round(10^attr(umi, "thresholds")[["higher"]]),
      genes_low    = round(10^attr(gen, "thresholds")[["lower"]]),
      genes_high   = round(10^attr(gen, "thresholds")[["higher"]]),
      mito_cutoff  = round(max(attr(mt, "thresholds")[["higher"]], mito_floor), 2)
    )

    # Doublets (per sample, on raw counts); RNG is handled by future.seed
    if (detect_doublets) {
      sce <- SingleCellExperiment(list(counts = GetAssayData(so, layer = "counts")))
      sce <- scDblFinder(sce, clusters = TRUE)
      so$doublet_score <- sce$scDblFinder.score
      so$doublet       <- sce$scDblFinder.class == "doublet"
    }

    so$is_outlier <- so$umi_outlier | so$ngenes_outlier | so$mt_outlier |
      (if (detect_doublets) so$doublet else FALSE)
    so
  }

  seurat_objs <- future.apply::future_Map(
    qc_one, inputs, names(inputs),
    future.seed = seed,
    future.packages = c("Seurat", "scCustomize", "scuttle", "scDblFinder", "SingleCellExperiment")
  )

  thresholds <- dplyr::bind_rows(lapply(seurat_objs, function(x) x@misc$qc_thresholds))

  merged <- if (length(seurat_objs) == 1) seurat_objs[[1]] else
    merge(seurat_objs[[1]], y = seurat_objs[-1], add.cell.ids = names(seurat_objs))
  merged@misc$qc_thresholds <- thresholds
  # only the filtering parameters; the outcome is qc_thresholds + qc_flagged
  merged@misc$pipeline$read_qc <- pipeline_info(list(
    nmads = nmads, mito_floor = mito_floor, detect_doublets = detect_doublets))

  flagged <- as.data.frame(
    merged[[]] |>
      dplyr::group_by(sample) |>
      dplyr::summarise(n_cells = dplyr::n(), n_flagged = sum(is_outlier),
                       pct_flagged = round(100 * mean(is_outlier), 1),
                       umi = sum(umi_outlier), genes = sum(ngenes_outlier),
                       mito = sum(mt_outlier),
                       doublet = if (detect_doublets) sum(doublet) else NA_integer_)
  )
  if (!detect_doublets) flagged$doublet <- NULL
  merged@misc$qc_flagged <- flagged   # counts before filtering; survives subset()

  message("\nMAD cutoffs (nmads = ", nmads, "; mito also needs > ", mito_floor, "%):")
  print(as.data.frame(thresholds), row.names = FALSE)

  message("\nCells flagged per sample:")
  print(flagged, row.names = FALSE)

  merged
}

#' Read 10X matrix directories and add per-cell QC
#'
#' Reads one 10X matrix directory per sample (in parallel), builds a Seurat object per sample,
#' attaches the sample sheet metadata, computes QC metrics with
#' [scCustomize::Add_Cell_QC_Metrics()], flags outliers by MAD (UMI, genes, mitochondrial
#' percentage) and optionally doublets with scDblFinder, then merges all samples.
#'
#' Old-style `genes.tsv.gz` files are read through a temporary symlinked copy, the raw data is
#' never modified. Cell cycle scoring is attempted per sample; if it fails (for example because
#' the marker genes are absent) a warning is raised and the QC metrics are computed without it.
#'
#' @section What is stored in `@misc`:
#' * `qc_thresholds`: the MAD cutoffs per sample (UMI, genes, mitochondrial percentage).
#' * `qc_flagged`: cells and flags per sample, counted before filtering.
#' * `pipeline$read_qc$params`: `nmads`, `mito_floor` and `detect_doublets`.
#'
#' @param base_path Directory with one sub-directory per sample (the sub-directory names are the
#'   sample names).
#' @param sample_sheet_path Path to a `.xlsx` or `.csv` sample sheet with one row per sample.
#' @param sample_sheet_column Column of the sample sheet that matches the sample names.
#' @param species Passed to [scCustomize::Add_Cell_QC_Metrics()] (for example `"human"`,
#'   `"mouse"`).
#' @param default_10X_path If `TRUE`, matrices are read from
#'   `<sample>/outs/filtered_feature_bc_matrix`; if `FALSE`, from `<sample>` directly.
#' @param nmads Number of median absolute deviations used for the outlier cutoffs.
#' @param mito_floor Cells are only flagged for mitochondrial content above this percentage,
#'   even if the MAD cutoff is lower.
#' @param seed Seed for the random steps (doublet detection); handled per worker.
#' @param workers `NULL` uses the current [future::plan()]; `1` runs sequentially; larger values
#'   start that many parallel R sessions for the duration of the call.
#' @param detect_doublets Run scDblFinder per sample. If `FALSE` no doublet columns are added
#'   and `is_outlier` ignores doublets.
#' @param ... Passed to [Seurat::Read10X()].
#'
#' @return A merged Seurat object. Metadata gains `sample`, the sample sheet columns, the
#'   scCustomize QC metrics, `log10_genes`, `log10_UMI`, `umi_outlier`, `ngenes_outlier`,
#'   `mt_outlier`, (`doublet_score` and `doublet`) and `is_outlier`.
#'
#' @seealso [read10X_h5_multi_directory_add_QC()], [read_cellbender_multi_directory_add_QC()],
#'   [plot_qc_ridges()], [save_read_qc()]
#' @examples
#' \dontrun{
#' so <- read10X_multi_directory_add_QC(
#'   "data/10x", "data/metadata.xlsx",
#'   sample_sheet_column = "Sample ID", species = "mouse",
#'   default_10X_path = FALSE, workers = 4
#' )
#' so_f <- subset(so, subset = !is_outlier)
#' }
#' @export
read10X_multi_directory_add_QC <- function(base_path, sample_sheet_path,
                                           sample_sheet_column = "sample",
                                           species = "human",
                                           default_10X_path = TRUE,
                                           nmads = 3, mito_floor = 5,
                                           seed = 1, workers = NULL,
                                           detect_doublets = TRUE, ...) {
  read_args <- list(...)

  # symlink-staged copy of the data if genes.tsv.gz needs to be features.tsv.gz
  staged <- stage_10x_dirs(base_path)
  if (!identical(staged, base_path)) on.exit(unlink(staged, recursive = TRUE), add = TRUE)

  dirs <- .sample_dirs(staged)
  inputs <- if (default_10X_path) file.path(dirs, "outs", "filtered_feature_bc_matrix") else dirs
  names(inputs) <- names(dirs)
  .check_inputs(inputs, "matrix directories")

  read_fun <- function(path) {
    m <- do.call(Read10X, c(list(data.dir = path), read_args))
    if (is.list(m)) m[["Gene Expression"]] else m
  }

  .read_add_qc(inputs, read_fun, sample_sheet_path, sample_sheet_column, species,
               detect_doublets, nmads, mito_floor, seed, workers)
}

#' Read 10X h5 files and add per-cell QC
#'
#' Same workflow and output as [read10X_multi_directory_add_QC()], reading
#' `<base_path>/<sample>/<secondary_path>/<h5_filename>` with [Seurat::Read10X_h5()]. File
#' locations follow [scCustomize::Read10X_h5_Multi_Directory()]; cellranger multi layouts are
#' not supported.
#'
#' @inheritParams read10X_multi_directory_add_QC
#' @param secondary_path Directory between the sample directory and the h5 file. Must be `NULL`
#'   when `default_10X_path = TRUE`.
#' @param default_10X_path If `TRUE`, the h5 file is read from `<sample>/outs/`.
#' @param h5_filename Name of the h5 file in every sample directory.
#' @param ... Passed to [Seurat::Read10X_h5()].
#' @inherit read10X_multi_directory_add_QC return
#' @export
read10X_h5_multi_directory_add_QC <- function(base_path, sample_sheet_path,
                                              sample_sheet_column = "sample",
                                              species = "human",
                                              secondary_path = NULL,
                                              default_10X_path = TRUE,
                                              h5_filename = "filtered_feature_bc_matrix.h5",
                                              nmads = 3, mito_floor = 5,
                                              seed = 1, workers = NULL,
                                              detect_doublets = TRUE, ...) {
  read_args <- list(...)
  if (default_10X_path) {
    if (!is.null(secondary_path)) stop("default_10X_path = TRUE requires secondary_path = NULL")
    secondary_path <- "outs"
  }

  dirs <- .sample_dirs(base_path)
  inputs <- vapply(dirs, .join_path, character(1), sub = secondary_path, file = h5_filename)
  .check_inputs(inputs, "h5 files")

  read_fun <- function(path) {
    m <- do.call(Read10X_h5, c(list(filename = path), read_args))
    if (is.list(m)) m[["Gene Expression"]] else m
  }

  .read_add_qc(inputs, read_fun, sample_sheet_path, sample_sheet_column, species,
               detect_doublets, nmads, mito_floor, seed, workers)
}

#' Read CellBender output and add per-cell QC
#'
#' Same workflow and output as [read10X_multi_directory_add_QC()], reading CellBender h5
#' files with [scCustomize::Read_CellBender_h5_Mat()]. File naming follows
#' [scCustomize::Read_CellBender_h5_Multi_Directory()]: by default
#' `<base_path>/<sample>/<sample>_out_filtered.h5`.
#'
#' @inheritParams read10X_multi_directory_add_QC
#' @param secondary_path Directory between the sample directory and the h5 file.
#' @param filtered_h5 `TRUE` reads `<sample>_out_filtered.h5`, `FALSE` reads `<sample>_out.h5`.
#' @param custom_name Replaces the file suffix (must end in `.h5`).
#' @param no_file_prefix If `TRUE` the file name is just the suffix, without the sample name.
#' @param ... Passed to [scCustomize::Read_CellBender_h5_Mat()] (for example `h5_group_name`,
#'   `feature_slot_name`).
#' @inherit read10X_multi_directory_add_QC return
#' @export
read_cellbender_multi_directory_add_QC <- function(base_path, sample_sheet_path,
                                                   sample_sheet_column = "sample",
                                                   species = "human",
                                                   secondary_path = NULL,
                                                   filtered_h5 = TRUE,
                                                   custom_name = NULL,
                                                   no_file_prefix = FALSE,
                                                   nmads = 3, mito_floor = 5,
                                                   seed = 1, workers = NULL,
                                                   detect_doublets = TRUE, ...) {
  read_args <- list(...)

  suffix <- custom_name %||% if (filtered_h5) "_out_filtered.h5" else "_out.h5"
  dirs <- .sample_dirs(base_path)
  fname <- if (no_file_prefix) rep(suffix, length(dirs)) else paste0(names(dirs), suffix)
  inputs <- stats::setNames(
    vapply(seq_along(dirs), function(i) .join_path(dirs[[i]], secondary_path, fname[i]), character(1)),
    names(dirs))
  .check_inputs(inputs, "CellBender h5 files")

  read_fun <- function(path) do.call(Read_CellBender_h5_Mat, c(list(file_name = path), read_args))

  .read_add_qc(inputs, read_fun, sample_sheet_path, sample_sheet_column, species,
               detect_doublets, nmads, mito_floor, seed, workers)
}
