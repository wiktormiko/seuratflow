#' UMAPs before and after Harmony
#'
#' One row per grouping column and one panel per reduction. Columns of `group_cols` that are not
#' in the metadata are skipped. Uses [umap_style()]; the returned plot carries a size hint (see
#' [save_plot()]).
#'
#' @param seurat_obj A preprocessed Seurat object.
#' @param group_cols Metadata columns to colour by.
#' @param reductions Named character vector of reductions; the names become the panel titles.
#' @param seed Seed for shuffling the cells.
#' @param fixed_coord `TRUE` uses `coord_fixed()`, so the panel shape follows the UMAP ranges;
#'   `FALSE` gives square panels (see [umap_style()]).
#' @return A patchwork object.
#' @export
plot_integration_umaps <- function(seurat_obj, group_cols = c("sample", "condition"),
                                   reductions = c("Before Harmony" = "umap_unintegrated",
                                                  "After Harmony"  = "umap_harmony"),
                                   seed = 1, fixed_coord=T) {
  reductions <- reductions[reductions %in% Reductions(seurat_obj)]
  stopifnot(length(reductions) > 0)
  group_cols <- intersect(group_cols, colnames(seurat_obj[[]]))
  stopifnot(length(group_cols) > 0)

  rows <- lapply(group_cols, function(g) {
    ps <- Map(function(r, title) {
      DimPlot(seurat_obj, reduction = r, group.by = g, shuffle = TRUE, seed = seed, raster = T) +
        ggtitle(title) +
        umap_style(fixed=fixed_coord)
    }, reductions, names(reductions))
    patchwork::wrap_plots(ps, nrow = 1) + patchwork::plot_layout(guides = "collect") 
  })
  # size hint: one row per grouping column; a row is as tall as its legend needs
  n_lv  <- vapply(group_cols, function(g) dplyr::n_distinct(seurat_obj[[g]][[1]]), integer(1))
  lab_w <- vapply(group_cols, function(g) max(nchar(as.character(unique(seurat_obj[[g]][[1]])))), numeric(1))
  .set_size(patchwork::wrap_plots(rows, ncol = 1) & umap_style(fixed=fixed_coord),
            width  = length(reductions) * 4.2 + 1.2 + 0.09 * max(lab_w),
            height = sum(pmax(4, 0.8 + 0.24 * n_lv)))
}

#' UMAPs coloured by QC metrics
#'
#' One panel per metric found in the metadata. Numeric metrics are drawn with
#' [Seurat::FeaturePlot()] and clipped to the 1st to 99th percentile (cutoffs are dropped when
#' both are equal, for example no hemoglobin reads); categorical metrics such as `Phase` use
#' [Seurat::DimPlot()].
#'
#' @inheritParams plot_integration_umaps
#' @param reduction Reduction to plot (default: the default dimensional reduction).
#' @param metrics Metadata columns. `NULL` uses `log10_UMI`, `log10_genes`, `percent_mito`,
#'   `percent_ribo`, `percent_mito_ribo`, `percent_hemo`, `log10GenesPerUMI`, `percent_top50`,
#'   `S.Score`, `G2M.Score`, `Phase` and `doublet_score` (the ones present).
#' @param ncol Number of panel columns.
#' @return A patchwork object.
#' @export
plot_qc_umaps <- function(seurat_obj, reduction = DefaultDimReduc(seurat_obj),
                          metrics = NULL, ncol = 3, seed = 1, fixed_coord=T) {
  candidates <- c("log10_UMI", "log10_genes", "percent_mito", "percent_ribo", "percent_mito_ribo",
                  "percent_hemo", "log10GenesPerUMI", "percent_top50",
                  "S.Score", "G2M.Score", "Phase", "doublet_score")
  md <- seurat_obj[[]]
  metrics <- intersect(metrics %||% candidates, colnames(md))
  stopifnot(length(metrics) > 0)

  # DimPlot turns group.by into a factor, so a numeric metric would get one colour level per
  # distinct value (tens of thousands). Numeric -> FeaturePlot, categorical (Phase) -> DimPlot.
  ps <- lapply(metrics, function(m) {
    if (is.numeric(md[[m]])) {
      # q01 == q99 (e.g. almost no hemoglobin reads) would give an empty colour scale -> no cutoff
      q <- quantile(md[[m]], c(0.01, 0.99), na.rm = TRUE)
      cut <- if (diff(q) == 0) c(NA, NA) else c("q01", "q99")
      FeaturePlot(seurat_obj, features = m, reduction = reduction, raster = TRUE,
                  min.cutoff = cut[1], max.cutoff = cut[2]) +
        umap_style(fixed = fixed_coord, type = "continuous")
    } else {
      DimPlot(seurat_obj, reduction = reduction, group.by = m, shuffle = TRUE, seed = seed,
              raster = TRUE) +
        ggtitle(m) +
        umap_style(fixed = fixed_coord, type = "discrete")
    }
  })
  # size hint: ncol panels per row, each with its own colour bar / legend
  g <- .panel_grid(length(ps), ncol = ncol)
  .set_size(patchwork::wrap_plots(ps, ncol = ncol),
            width = g[["ncol"]] * 4.4, height = g[["nrow"]] * 3.8)
}

#' UMAPs of the clusters at every resolution
#'
#' One panel per `<prefix><resolution>` column, with cluster labels on the plot instead of a
#' legend.
#'
#' @inheritParams plot_qc_umaps
#' @param prefix Prefix of the cluster columns.
#' @return A patchwork object.
#' @export
plot_cluster_umaps <- function(seurat_obj, prefix = "leiden_",
                               reduction = DefaultDimReduc(seurat_obj),
                               ncol = 2, seed = 1, fixed_coord = TRUE) {
  cols <- grep(paste0("^", prefix), colnames(seurat_obj[[]]), value = TRUE)
  stopifnot(length(cols) > 0)
  cols <- cols[order(as.numeric(sub(prefix, "", cols)))]

  ps <- lapply(cols, function(cl) {
    DimPlot(seurat_obj, reduction = reduction, group.by = cl, label = TRUE, repel = TRUE,
            shuffle = TRUE, seed = seed, raster = TRUE) +
      ggtitle(paste0(cl, " (", dplyr::n_distinct(seurat_obj[[cl]][[1]]), " clusters)")) +
      umap_style(fixed = fixed_coord, type = "discrete") +
      NoLegend()
  })
  g <- .panel_grid(length(ps), ncol = ncol)
  .set_size(patchwork::wrap_plots(ps, ncol = ncol),
            width = g[["ncol"]] * 4.4, height = g[["nrow"]] * 4.2)
}
