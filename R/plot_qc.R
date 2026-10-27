#' Ridge plots of the QC metrics per sample
#'
#' One ridge panel per metric with one ridge per sample. The MAD cutoffs stored by the readers
#' (`@misc$qc_thresholds`) are drawn as red dashed lines. The returned plot carries a size hint
#' (see [save_plot()]) that grows with the number of samples.
#'
#' @param seurat_obj A Seurat object from one of the `*_add_QC()` readers.
#' @param sample_col Metadata column with the sample.
#' @param metrics Metadata columns to plot.
#' @param log_metrics Metrics shown on a log10 axis.
#' @param zoom_quantile For metrics that are not logged, the x axis is cut at this quantile so a
#'   long tail does not flatten the distributions. `NULL` disables it.
#' @param show_thresholds Draw the stored MAD cutoffs.
#' @return A patchwork object.
#' @seealso [plot_qc_scatter()], [save_read_qc()]
#' @export
plot_qc_ridges <- function(seurat_obj, sample_col = "sample",
                           metrics = c("nCount_RNA", "nFeature_RNA", "percent_mito"),
                           log_metrics = c("nCount_RNA", "nFeature_RNA"),
                           zoom_quantile = 0.95, show_thresholds = TRUE) {
  df <- seurat_obj[[]]
  df[[sample_col]] <- factor(df[[sample_col]])
  thr <- seurat_obj@misc$qc_thresholds
  thr_cols <- list(nCount_RNA = c("UMI_low", "UMI_high"),
                   nFeature_RNA = c("genes_low", "genes_high"),
                   percent_mito = "mito_cutoff")
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel2"))(nlevels(df[[sample_col]]))

  plots <- lapply(metrics, function(metric) {
    df$value <- df[[metric]]
    p <- ggplot(df, aes(x = value, y = .data[[sample_col]], fill = .data[[sample_col]])) +
      # ndensity scales each sample to its own peak, trim stops the density at the data range
      ggridges::geom_density_ridges(stat = "density", aes(height = after_stat(ndensity)), trim = TRUE) +
      scale_fill_manual(values = pal) +
      labs(x = metric, y = NULL) +
      ggridges::theme_ridges() +
      theme(legend.position = "none")

    if (metric %in% log_metrics) {
      p <- p + scale_x_log10()
    } else if (!is.null(zoom_quantile)) {
      # zoom on the bulk of the cells so a long tail does not squish the distributions
      p <- p + coord_cartesian(xlim = c(min(df$value), quantile(df$value, zoom_quantile)))
    }

    if (show_thresholds && !is.null(thr) && metric %in% names(thr_cols)) {
      seg <- thr |>
        tidyr::pivot_longer(dplyr::all_of(thr_cols[[metric]]), values_to = "x") |>
        dplyr::mutate(y = as.numeric(factor(sample, levels = levels(df[[sample_col]]))))
      p <- p + geom_segment(data = seg, aes(x = x, xend = x, y = y, yend = y + 0.9),
                            inherit.aes = FALSE, colour = "red", linetype = "dashed")
    }
    p
  })

  # size hint: one panel per metric (patchwork's default grid), height grows with the samples
  g <- .panel_grid(length(plots))
  .set_size(patchwork::wrap_plots(plots),
            width  = g[["ncol"]] * 4.5 + 0.5,
            height = g[["nrow"]] * max(2.5, 1 + 0.4 * nlevels(df[[sample_col]])) + 0.3)
}

#' Scatter plots of genes versus UMI and genes versus mitochondrial percentage
#'
#' Both plots are faceted by sample and coloured by `flag_col`. The returned plot carries a size
#' hint (see [save_plot()]).
#'
#' @inheritParams plot_qc_ridges
#' @param flag_col Logical metadata column used for the colour (default `is_outlier`).
#' @param ncol Number of facet columns.
#' @return A patchwork object.
#' @export
plot_qc_scatter <- function(seurat_obj, sample_col = "sample", flag_col = "is_outlier", ncol = 2) {
  df <- seurat_obj[[]]
  if (!flag_col %in% names(df)) df[[flag_col]] <- FALSE

  common <- list(
    geom_point(size = 0.3, alpha = 0.5),
    facet_wrap(sample_col, ncol = ncol),
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red"), name = flag_col),
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    theme_bw()
  )
  p1 <- ggplot(df, aes(nFeature_RNA, nCount_RNA, colour = .data[[flag_col]])) +
    common + scale_x_log10() + scale_y_log10()
  p2 <- ggplot(df, aes(nFeature_RNA, percent_mito, colour = .data[[flag_col]])) +
    common + scale_x_log10()

  # size hint: two plots side by side, each faceted by sample in `ncol` columns
  n_s <- dplyr::n_distinct(df[[sample_col]])
  .set_size(patchwork::wrap_plots(p1, p2) + patchwork::plot_layout(guides = "collect"),
            width  = 2 * (min(ncol, n_s) * 2.6 + 0.8) + 1,
            height = ceiling(n_s / ncol) * 2.6 + 1)
}
