#' Compare clustering resolutions with clustree and silhouette scores
#'
#' Silhouette widths are computed on the embedding the clusters were built on (Harmony if it
#' was run, otherwise PCA) using the first `npcs` dimensions, on a random subsample because the
#' distance matrix grows with the square of the number of cells. Silhouette tends to favour few
#' clusters, so read it together with the clustree and the biology.
#'
#' @param seurat_obj A Seurat object processed by [seurat_preprocess()] (needs
#'   `@misc$pipeline$preprocess$params$npcs`).
#' @param prefix Prefix of the cluster columns.
#' @param n_cells Maximum number of cells used for the silhouette.
#' @param seed Seed for the subsample.
#' @return A list with `table` (one row per resolution), `per_cluster` (one row per cluster and
#'   resolution), and the plots `clustree`, `silhouette` and `silhouette_by_cluster`. The plots
#'   carry size hints (see [save_plot()]).
#' @export
choose_resolution <- function(seurat_obj, prefix = "leiden_", n_cells = 5000, seed = 1) {
  info <- seurat_obj@misc$pipeline$preprocess$params
  stopifnot(!is.null(info$npcs))
  # same space the neighbors/clusters were built in: harmony if it was run, else pca
  red <- if ("harmony" %in% Reductions(seurat_obj)) "harmony" else "pca"
  emb <- Embeddings(seurat_obj, red)[, 1:info$npcs]

  cols <- grep(paste0("^", prefix), colnames(seurat_obj[[]]), value = TRUE)
  cols <- cols[order(as.numeric(sub(prefix, "", cols)))]

  set.seed(seed)
  cells <- if (nrow(emb) > n_cells) sample(rownames(emb), n_cells) else rownames(emb)
  d <- dist(emb[cells, ])

  # silhouette width per (subsampled) cell and resolution
  cell_sil <- do.call(rbind, lapply(cols, function(cl) {
    lab <- as.character(seurat_obj[[cl]][cells, 1])
    if (length(unique(lab)) < 2) return(NULL)
    data.frame(column = cl,
               resolution = as.numeric(sub(prefix, "", cl)),
               cluster = lab,
               sil_width = cluster::silhouette(as.integer(factor(lab)), d)[, "sil_width"])
  }))
  stopifnot(!is.null(cell_sil))

  # one row per resolution
  res <- data.frame(column = cols,
                    resolution = as.numeric(sub(prefix, "", cols)),
                    n_clusters = vapply(cols, function(cl) dplyr::n_distinct(seurat_obj[[cl]][[1]]),
                                        integer(1), USE.NAMES = FALSE))
  res$mean_silhouette <- vapply(res$column, function(cl)
    if (cl %in% cell_sil$column) mean(cell_sil$sil_width[cell_sil$column == cl]) else NA_real_,
    numeric(1), USE.NAMES = FALSE)
  print(res, row.names = FALSE)

  # one row per cluster within each resolution (n_cells is within the subsample)
  per_cluster <- cell_sil |>
    dplyr::group_by(resolution, cluster) |>
    dplyr::summarise(n_cells = dplyr::n(), mean_silhouette = mean(sil_width), .groups = "drop") |>
    dplyr::left_join(dplyr::select(res, resolution, resolution_mean = mean_silhouette), by = "resolution")
  cl_levels <- unique(per_cluster$cluster[order(suppressWarnings(as.numeric(per_cluster$cluster)))])
  per_cluster$cluster <- factor(per_cluster$cluster, levels = cl_levels)

  p_sil <- ggplot(res, aes(resolution, mean_silhouette)) +
    geom_line() + geom_point(aes(size = n_clusters)) +
    theme_bw() + labs(y = "mean silhouette", size = "clusters")

  # bars below the resolution's own mean (orange) are the weakest-separated clusters
  p_cl <- ggplot(per_cluster, aes(cluster, mean_silhouette,
                                  fill = mean_silhouette < resolution_mean)) +
    geom_col() +
    geom_hline(aes(yintercept = resolution_mean), linetype = "dashed") +
    geom_hline(yintercept = 0) +
    facet_wrap(~resolution, scales = "free_x", labeller = label_both) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "darkorange"), guide = "none") +
    theme_bw() + labs(x = "cluster", y = "mean silhouette")

  # size hints: panel width follows the number of clusters, clustree height the number of resolutions
  n_res <- nrow(res)
  max_k <- max(res$n_clusters)
  g <- .panel_grid(n_res)

  list(table = res,
       per_cluster = per_cluster,
       # ggplot2 >= 4.0 looks the "edge_colourbar" guide up by name among attached packages, so
       # without library(ggraph) the plot fails to draw; passing the guide object avoids that
       clustree = .set_size(clustree::clustree(seurat_obj, prefix = prefix) +
                              guides(edge_colour = ggraph::guide_edge_colourbar()),
                            width = max(6, 0.45 * max_k + 3), height = max(5, 1.5 * n_res + 1)),
       silhouette = .set_size(p_sil, width = 5, height = 4),
       silhouette_by_cluster = .set_size(p_cl,
                                         width  = g[["ncol"]] * max(3, 0.28 * max_k + 1) + 0.5,
                                         height = g[["nrow"]] * 3 + 0.5))
}
