#' Attach a size hint (inches) to a plot
#' @noRd
.set_size <- function(p, width, height) {
  attr(p, "size_in") <- c(width = width, height = height)
  p
}

#' Rows and columns a patchwork uses for n panels
#' @noRd
.panel_grid <- function(n, ncol = NULL, nrow = NULL) {
  d <- ggplot2::wrap_dims(n, nrow = nrow, ncol = ncol)
  c(nrow = d[[1]], ncol = d[[2]])
}

#' Fallback size for plain ggplots without a hint
#' @noRd
.facet_size <- function(p) {
  lay <- tryCatch(ggplot2::ggplot_build(p)$layout$layout, error = function(e) NULL)
  if (is.null(lay) || !all(c("ROW", "COL") %in% names(lay))) return(NULL)
  c(width = 3.5 * max(lay$COL) + 1.5, height = 3 * max(lay$ROW) + 1)
}

#' Save a plot in several formats with an automatic size
#'
#' The plot functions of this package attach a size hint (inches) computed from what they draw:
#' number of samples, panels, clusters, genes and legend entries. `save_plot()` uses it, so you
#' do not have to guess a width and height. Adding layers to a returned plot (`p + theme(...)`)
#' drops the hint; a plain ggplot is then sized from its facet layout (3.5 x 3 inches per panel
#' plus margins) and a patchwork gets 8 x 6 inches.
#'
#' @param p A ggplot or patchwork.
#' @param file Path without extension (an image extension is stripped if present).
#' @param width,height Override the automatic size (inches).
#' @param formats File formats to write.
#' @param dpi Resolution of raster formats.
#' @param scale Multiplier applied to the size (automatic or given).
#' @return A data frame with the files written and the size used (invisibly useful for logging).
#' @export
save_plot <- function(p, file, width = NULL, height = NULL, formats = c("png", "pdf"),
                      dpi = 300, scale = 1) {
  hint <- attr(p, "size_in")
  # a patchwork without a hint has no single facet layout to read, so it gets the default
  if (is.null(hint) && inherits(p, "ggplot") && !inherits(p, "patchwork")) hint <- .facet_size(p)
  width  <- (width  %||% hint[["width"]]  %||% 8) * scale
  height <- (height %||% hint[["height"]] %||% 6) * scale

  # only strip a real image extension: names like "dotplot_leiden_0.5" contain a dot
  base <- sub("\\.(png|pdf|svg|jpe?g|tiff?)$", "", file, ignore.case = TRUE)
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  files <- paste0(base, ".", formats)
  for (f in files) {
    ggplot2::ggsave(f, p, width = width, height = height, dpi = dpi, limitsize = FALSE)
  }
  data.frame(file = files, width = round(width, 1), height = round(height, 1))
}

#' Write pipeline parameters and session info next to the plots
#' @noRd
.write_params <- function(seurat_obj, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- character()
  p <- seurat_obj@misc$pipeline
  if (!is.null(p)) {
    f <- file.path(out_dir, "pipeline_parameters.yml")
    txt <- tryCatch(yaml::as.yaml(p),
                    error = function(e) paste(utils::capture.output(utils::str(p)), collapse = "\n"))
    writeLines(txt, f)
    files <- c(files, f)
  }
  f <- file.path(out_dir, "sessionInfo.txt")
  writeLines(utils::capture.output(utils::sessionInfo()), f)
  data.frame(file = c(files, f), width = NA_real_, height = NA_real_)
}

#' Save the read-stage QC outputs to the QualityControl folder
#'
#' Writes to `<out_dir>/read_qc/`: `qc_thresholds.csv`, `qc_flagged.csv`, the ridge and scatter
#' plots, and (in `<out_dir>`) `pipeline_parameters.yml` and `sessionInfo.txt`. Run it on the
#' object returned by a `*_add_QC()` reader; call it again on the filtered object with
#' `tag = "postfilter"` to save the after-filtering plots next to the originals.
#'
#' @param seurat_obj A Seurat object from one of the readers.
#' @param out_dir Output folder.
#' @param tag `"prefilter"` draws the MAD cutoffs on the ridges; any other tag does not.
#' @param sample_col Metadata column with the sample.
#' @param formats,scale Passed to [save_plot()].
#' @return A data frame of the files written (invisibly).
#' @seealso [save_preprocess_qc()]
#' @export
save_read_qc <- function(seurat_obj, out_dir = "QualityControl", tag = "prefilter",
                         sample_col = "sample", formats = c("png", "pdf"), scale = 1) {
  d <- file.path(out_dir, "read_qc")
  m <- list()
  if (!is.null(seurat_obj@misc$qc_thresholds)) m$thr <- .csv_out(seurat_obj@misc$qc_thresholds, file.path(d, "qc_thresholds.csv"))
  if (!is.null(seurat_obj@misc$qc_flagged))    m$flg <- .csv_out(seurat_obj@misc$qc_flagged,    file.path(d, "qc_flagged.csv"))

  m$ridges <- .qc_step("ridges", save_plot(
    plot_qc_ridges(seurat_obj, sample_col = sample_col, show_thresholds = identical(tag, "prefilter")),
    file.path(d, paste0("qc_ridges_", tag)), formats = formats, scale = scale))
  m$scatter <- .qc_step("scatter", save_plot(
    plot_qc_scatter(seurat_obj, sample_col = sample_col),
    file.path(d, paste0("qc_scatter_", tag)), formats = formats, scale = scale))
  m$params <- .write_params(seurat_obj, out_dir)

  manifest <- do.call(rbind, m)
  rownames(manifest) <- NULL
  message("Saved ", nrow(manifest), " files to ", normalizePath(d, mustWork = FALSE))
  invisible(manifest)
}

#' Save the preprocessing UMAPs to the QualityControl folder
#'
#' Writes `umap_integration` (before and after Harmony), `umap_qc_metrics` and `umap_clusters`
#' to `<out_dir>/preprocessing/`. [save_preprocess_qc()] calls this; use it alone to refresh the
#' UMAPs without the resolution and marker steps.
#'
#' @param seurat_obj A preprocessed Seurat object.
#' @param out_dir Output folder.
#' @param group_cols Metadata columns for the integration UMAPs.
#' @param prefix Prefix of the cluster columns.
#' @param fixed_coord `FALSE` gives square panels, which also makes the automatic size more
#'   accurate.
#' @param formats,scale Passed to [save_plot()].
#' @return A data frame of the files written (invisibly).
#' @export
save_preprocess_umaps <- function(seurat_obj, out_dir = "QualityControl",
                                  group_cols = c("sample", "condition"), prefix = "leiden_",
                                  fixed_coord = FALSE, formats = c("png", "pdf"), scale = 1) {
  d <- file.path(out_dir, "preprocessing")
  sv <- function(p, name) save_plot(p, file.path(d, name), formats = formats, scale = scale)
  m <- list()

  m$integration <- .qc_step("integration UMAPs", sv(
    plot_integration_umaps(seurat_obj, group_cols = group_cols, fixed_coord = fixed_coord),
    "umap_integration"))
  m$qc_umaps <- .qc_step("QC UMAPs", sv(
    plot_qc_umaps(seurat_obj, fixed_coord = fixed_coord), "umap_qc_metrics"))
  m$clusters <- .qc_step("cluster UMAPs", sv(
    plot_cluster_umaps(seurat_obj, prefix = prefix, fixed_coord = fixed_coord), "umap_clusters"))

  manifest <- do.call(rbind, m)
  if (!is.null(manifest)) rownames(manifest) <- NULL
  message("Saved ", NROW(manifest), " files to ", normalizePath(d, mustWork = FALSE))
  invisible(manifest)
}

#' Save the preprocessing QC outputs to the QualityControl folder
#'
#' Writes to `<out_dir>/preprocessing/`: the elbow plot, the UMAPs
#' ([save_preprocess_umaps()]), the resolution choice (clustree, silhouette plots and tables)
#' and, with `markers = TRUE`, marker tables and a DotPlot per resolution in `markers/`. It also
#' writes `pipeline_parameters.yml` and `sessionInfo.txt` in `<out_dir>`. Every output is its
#' own step: a failure becomes a warning and the remaining outputs are still written.
#'
#' @inheritParams save_preprocess_umaps
#' @param markers Also run [find_markers_all_resolutions()] and save its results.
#' @param top_n Markers per cluster for the DotPlots.
#' @return A data frame of the files written (invisibly).
#' @seealso [save_read_qc()]
#' @export
save_preprocess_qc <- function(seurat_obj, out_dir = "QualityControl",
                               group_cols = c("sample", "condition"),
                               markers = FALSE, top_n = 5, fixed_coord = FALSE,
                               formats = c("png", "pdf"), scale = 1) {
  d <- file.path(out_dir, "preprocessing")
  sv <- function(p, name, sub = NULL) {
    save_plot(p, do.call(file.path, as.list(c(d, sub, name))), formats = formats, scale = scale)
  }
  m <- list()
  info <- seurat_obj@misc$pipeline$preprocess$params

  if (!is.null(info$npcs)) m$elbow <- .qc_step("elbow", sv(
    .set_size(ElbowPlot(seurat_obj, ndims = 30) +
                geom_vline(xintercept = info$npcs, colour = "darkorange", linetype = "dashed"),
              width = 6, height = 4), "elbow"))

  m$umaps <- suppressMessages(save_preprocess_umaps(
    seurat_obj, out_dir = out_dir, group_cols = group_cols, fixed_coord = fixed_coord,
    formats = formats, scale = scale))

  res <- .qc_step("resolution choice", choose_resolution(seurat_obj))
  if (!is.null(res)) {
    m$res_tab <- .csv_out(res$table,       file.path(d, "resolution_summary.csv"))
    m$cl_tab  <- .csv_out(res$per_cluster, file.path(d, "silhouette_per_cluster.csv"))
    m$clustree   <- .qc_step("clustree",   sv(res$clustree, "clustree"))
    m$silhouette <- .qc_step("silhouette", sv(res$silhouette, "silhouette_by_resolution"))
    m$sil_cl     <- .qc_step("silhouette per cluster", sv(res$silhouette_by_cluster, "silhouette_by_cluster"))
  }

  if (markers) {
    mk <- .qc_step("markers", find_markers_all_resolutions(seurat_obj, top_n = top_n))
    for (r in names(mk$dotplots)) {
      m[[paste0("mk_", r)]] <- .csv_out(mk$markers[[r]],     file.path(d, "markers", paste0("all_markers_", r, ".csv")))
      m[[paste0("top_", r)]] <- .csv_out(mk$top_markers[[r]], file.path(d, "markers", paste0("top_markers_", r, ".csv")))
      if (!is.null(mk$dotplots[[r]]))
        m[[paste0("dp_", r)]] <- .qc_step(paste("dotplot", r), sv(mk$dotplots[[r]], paste0("dotplot_", r), "markers"))
    }
  }
  m$params <- .write_params(seurat_obj, out_dir)

  manifest <- do.call(rbind, m)
  rownames(manifest) <- NULL
  message("Saved ", nrow(manifest), " files to ", normalizePath(d, mustWork = FALSE))
  invisible(manifest)
}
