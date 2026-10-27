#' Markers and DotPlots at every resolution
#'
#' Runs [Seurat::FindAllMarkers()] on each `<prefix><resolution>` column, keeps the `top_n`
#' significant markers per cluster and draws one DotPlot per resolution. Works on copies of the
#' object: layers are joined for log-normalized data, and
#' [Seurat::PrepSCTFindMarkers()] is run for SCT data.
#'
#' Install presto for a much faster Wilcoxon test; `max_cells_per_ident` downsamples each
#' cluster; a [future::plan()] parallelises across clusters.
#'
#' @param seurat_obj A clustered Seurat object.
#' @param prefix Prefix of the cluster columns.
#' @param resolutions Optional subset of resolutions (the numbers after the prefix).
#' @param top_n Markers kept per cluster.
#' @param rank_by Column of the marker table used to rank (default `avg_log2FC`).
#' @param p_adj_cutoff Adjusted p-value cutoff for the top markers.
#' @param only_pos,min_pct,logfc_threshold,max_cells_per_ident Passed to
#'   [Seurat::FindAllMarkers()].
#' @param assay Assay to use (default: the default assay).
#' @param seed Random seed.
#' @param verbose Print progress of Seurat.
#' @param ... Passed to [Seurat::FindAllMarkers()] (for example `test.use`).
#' @return A list with `markers` (all markers), `top_markers` and `dotplots`, each a list named
#'   by the cluster column.
#' @export
find_markers_all_resolutions <- function(seurat_obj, prefix = "leiden_", resolutions = NULL,
                                         top_n = 5, rank_by = "avg_log2FC", p_adj_cutoff = 0.05,
                                         only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0.25,
                                         max_cells_per_ident = Inf, assay = NULL,
                                         seed = 1, verbose = FALSE, ...) {
  cols <- grep(paste0("^", prefix), colnames(seurat_obj[[]]), value = TRUE)
  if (!is.null(resolutions)) cols <- intersect(cols, paste0(prefix, resolutions))
  stopifnot(length(cols) > 0)
  cols <- cols[order(as.numeric(sub(prefix, "", cols)))]

  assay <- assay %||% DefaultAssay(seurat_obj)
  DefaultAssay(seurat_obj) <- assay
  if (assay == "SCT") {
    seurat_obj <- PrepSCTFindMarkers(seurat_obj, verbose = verbose)
  } else if (length(Layers(seurat_obj[[assay]], search = "data")) > 1) {
    seurat_obj <- JoinLayers(seurat_obj, assay = assay)
  }

  out <- lapply(cols, function(cl) {
    # clusters in numeric order so the DotPlot reads as a diagonal
    x <- seurat_obj[[cl]][[1]]
    lv <- levels(factor(x))
    seurat_obj[[cl]] <- factor(x, levels = lv[order(suppressWarnings(as.numeric(lv)))])
    Idents(seurat_obj) <- cl

    message("[", cl, "] FindAllMarkers on ", nlevels(Idents(seurat_obj)), " clusters")
    mk <- FindAllMarkers(seurat_obj, only.pos = only_pos, min.pct = min_pct,
                         logfc.threshold = logfc_threshold,
                         max.cells.per.ident = max_cells_per_ident,
                         random.seed = seed, verbose = verbose, ...)

    top <- mk |>
      dplyr::filter(p_val_adj < p_adj_cutoff) |>
      dplyr::group_by(cluster) |>
      dplyr::slice_max(order_by = .data[[rank_by]], n = top_n, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::arrange(cluster)

    genes <- unique(top$gene)
    p <- if (length(genes) == 0) NULL else
      DotPlot(seurat_obj, features = genes, group.by = cl, assay = assay) +
        scale_colour_viridis_c(option = "rocket", direction = -1) +
        RotatedAxis() +
        labs(title = paste0(cl, " (top ", top_n, " per cluster)"), x = NULL, y = "cluster")

    # size hint: width follows the number of genes, height the number of clusters
    if (!is.null(p)) p <- .set_size(p, width = 0.24 * length(genes) + 2.8,
                                    height = 0.26 * nlevels(Idents(seurat_obj)) + 2.2)

    list(markers = mk, top = top, dotplot = p)
  })
  names(out) <- cols

  list(markers     = lapply(out, `[[`, "markers"),
       top_markers = lapply(out, `[[`, "top"),
       dotplots    = lapply(out, `[[`, "dotplot"))
}
