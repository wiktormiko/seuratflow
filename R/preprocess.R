#' Normalize, reduce, integrate and cluster
#'
#' Runs normalization (log-normalization or SCTransform), PCA, a UMAP of the unintegrated data,
#' optional Harmony integration with a second UMAP, and Leiden clustering at several
#' resolutions. The number of PCs is chosen with [findPC::findPC()] unless `npcs` is given. An
#' elbow plot with the chosen cutoff is printed.
#'
#' Reductions created: `pca`, `umap_unintegrated`, and with Harmony `harmony` and
#' `umap_harmony` (the default reduction). Clusters are stored in `leiden_<resolution>`.
#'
#' With `normalization = "sct"` on a merged object that still has one layer per sample,
#' [Seurat::SCTransform()] fits one model per sample. The glmGamPoi package makes it much
#' faster. Marker detection on SCT data needs [Seurat::PrepSCTFindMarkers()], which
#' [find_markers_all_resolutions()] runs for you.
#'
#' `@misc$pipeline$preprocess$params` stores `npcs` and `vars_to_regress`.
#'
#' @param seurat_obj A Seurat object, typically filtered with
#'   `subset(so, subset = !is_outlier)`.
#' @param run_harmony Run Harmony. Skipped with a warning if `batch` has a single level.
#' @param batch Metadata column(s) to integrate over.
#' @param resolutions Clustering resolutions.
#' @param npcs Number of PCs. `NULL` chooses it automatically.
#' @param pc_buffer PCs added to the automatic choice.
#' @param max_pcs Number of PCs computed (and upper bound for `npcs`).
#' @param n_hvg Number of variable features.
#' @param vars_to_regress Metadata columns to regress out (for example `"percent_mito"`).
#' @param normalization `"lognorm"` or `"sct"`.
#' @param seed Seed for PCA, UMAP and clustering.
#' @param verbose Print progress of the Seurat steps.
#' @return The Seurat object with the new reductions and cluster columns.
#' @seealso [choose_resolution()], [plot_integration_umaps()], [save_preprocess_qc()]
#' @export
seurat_preprocess <- function(seurat_obj, run_harmony = TRUE, batch = "sample",
                              resolutions = c(0.2, 0.5, 0.8, 1.2),
                              npcs = NULL, pc_buffer = 2, max_pcs = 50,
                              n_hvg = 2000, vars_to_regress = NULL,
                              normalization = c("lognorm", "sct"),
                              seed = 1, verbose = FALSE) {
  normalization <- match.arg(normalization)
  if (run_harmony) {
    stopifnot(all(batch %in% colnames(seurat_obj[[]])))
    if (nrow(unique(seurat_obj[[batch]])) < 2) {
      warning("Only one batch level; skipping Harmony.")
      run_harmony <- FALSE
    }
  }

  if (normalization == "sct") {
    seurat_obj <- seurat_obj |>
      SCTransform(variable.features.n = n_hvg, vars.to.regress = vars_to_regress,
                  verbose = verbose) |>
      RunPCA(npcs = max_pcs, seed.use = seed, verbose = verbose)
  } else {
    seurat_obj <- seurat_obj |>
      NormalizeData(verbose = verbose) |>
      FindVariableFeatures(nfeatures = n_hvg, verbose = verbose) |>
      ScaleData(vars.to.regress = vars_to_regress, verbose = verbose) |>
      RunPCA(npcs = max_pcs, seed.use = seed, verbose = verbose)
  }

  if (is.null(npcs)) {
    sdev <- Stdev(seurat_obj, reduction = "pca")
    npcs <- as.integer(findPC(sdev = sdev, number = length(sdev),
                              method = "perpendicular line", figure = FALSE)) + pc_buffer
  }
  npcs <- min(npcs, max_pcs)
  message("Using ", npcs, " PCs")

  seurat_obj <- RunUMAP(seurat_obj, dims = 1:npcs, reduction = "pca",
                        reduction.name = "umap_unintegrated", reduction.key = "UMAPunint_",
                        seed.use = seed, verbose = verbose)
  red <- "pca"
  if (run_harmony) {
    seurat_obj <- seurat_obj |>
      RunHarmony(group.by.vars = batch, reduction.use = "pca",
                 reduction.save = "harmony", verbose = verbose) |>
      RunUMAP(dims = 1:npcs, reduction = "harmony", reduction.name = "umap_harmony", reduction.key = "UMAP_", 
              verbose = verbose)
    red <- "harmony"
    DefaultDimReduc(seurat_obj) <- 'umap_harmony'
  }

  cl_names <- paste0("leiden_", resolutions)
  seurat_obj <- seurat_obj |>
    FindNeighbors(dims = 1:npcs, reduction = red,
                  verbose = verbose) |>
    FindClusters(resolution = resolutions, algorithm = 4, cluster.name = cl_names,
                 random.seed = seed, verbose = verbose)

  # npcs is reused by choose_resolution() and save_preprocess_qc()
  seurat_obj@misc$pipeline$preprocess <- pipeline_info(list(
    npcs = npcs, vars_to_regress = vars_to_regress))

  print(ElbowPlot(seurat_obj, ndims=30) + geom_vline(xintercept = npcs, color='darkorange', linetype = 'dashed'))
  seurat_obj
}
