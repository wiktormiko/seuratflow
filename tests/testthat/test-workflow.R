skip_on_cran()

# The whole workflow on the simulated dataset: built once, checked step by step.
src <- example_dir()
sheet <- file.path(src, "metadata.csv")
so <- quiet(read10X_multi_directory_add_QC(src, sheet, default_10X_path = FALSE,
                                           detect_doublets = FALSE, workers = 1))

test_that("write_example_10x writes the old 10X layout and a sample sheet", {
  expect_true(all(file.exists(file.path(src, "s1",
                                        c("matrix.mtx.gz", "barcodes.tsv.gz", "genes.tsv.gz")))))
  expect_equal(readr::read_csv(sheet, show_col_types = FALSE)$sample, c("s1", "s2", "s3"))
})

test_that("the mtx reader returns a merged object with QC columns and a lean @misc", {
  expect_equal(ncol(so), 3 * 150)
  expect_true(all(c("sample", "condition", "percent_mito", "umi_outlier", "ngenes_outlier",
                    "mt_outlier", "is_outlier") %in% names(so[[]])))
  expect_false(any(c("doublet", "doublet_score") %in% names(so[[]])))
  expect_setequal(names(so@misc), c("qc_thresholds", "qc_flagged", "pipeline"))
  expect_equal(nrow(so@misc$qc_thresholds), 3)
  expect_false("doublet" %in% names(so@misc$qc_flagged))
  expect_setequal(names(so@misc$pipeline$read_qc$params),
                  c("nmads", "mito_floor", "detect_doublets"))
  # raw data untouched: genes.tsv.gz was read through a temporary symlinked copy
  expect_false(file.exists(file.path(src, "s1", "features.tsv.gz")))
})

test_that("readers fail early with useful messages", {
  bad <- example_dir(n_samples = 2, n_cells = 50, n_genes = 100)
  file.rename(file.path(bad, "s2"), file.path(bad, "s9"))
  expect_error(quiet(read10X_multi_directory_add_QC(bad, file.path(bad, "metadata.csv"),
                                                    default_10X_path = FALSE)),
               "not in sample sheet")
  expect_error(read10X_h5_multi_directory_add_QC(src, sheet, secondary_path = "x"),
               "secondary_path = NULL")
  expect_error(read10X_h5_multi_directory_add_QC(src, sheet), "Missing h5 files")
})

test_that("detect_doublets = TRUE adds doublet columns", {
  skip_if_not_installed("scDblFinder")
  d <- example_dir(n_samples = 2, n_cells = 300, n_genes = 400)
  x <- quiet(read10X_multi_directory_add_QC(d, file.path(d, "metadata.csv"),
                                            default_10X_path = FALSE, detect_doublets = TRUE))
  expect_true(all(c("doublet", "doublet_score") %in% names(x[[]])))
  expect_true("doublet" %in% names(x@misc$qc_flagged))
})

test_that("h5 and CellBender readers give the same cells and flags as the mtx reader", {
  skip_if_not_installed("hdf5r")
  h5 <- tempfile("h5_"); cb <- tempfile("cb_")
  write_example_h5(src, h5, cb)
  x <- quiet(read10X_h5_multi_directory_add_QC(h5, file.path(h5, "metadata.csv"),
                                               detect_doublets = FALSE))
  y <- quiet(read_cellbender_multi_directory_add_QC(cb, file.path(cb, "metadata.csv"),
                                                    detect_doublets = FALSE))
  expect_equal(dim(x), dim(so))
  expect_equal(dim(y), dim(so))
  expect_equal(sum(x$nCount_RNA), sum(so$nCount_RNA))
  expect_identical(x$is_outlier, so$is_outlier)
  expect_identical(y$is_outlier, so$is_outlier)
})

test_that("parallel reading matches sequential reading", {
  skip_if_not_installed("pkgload")
  skip_if(pkgload::is_dev_package("seuratflow"), "workers need the installed package")
  par <- quiet(read10X_multi_directory_add_QC(src, sheet, default_10X_path = FALSE,
                                              detect_doublets = FALSE, workers = 2))
  expect_identical(par$is_outlier, so$is_outlier)
})

test_that("QC plots carry a size hint that grows with the number of samples", {
  r <- plot_qc_ridges(so)
  s <- plot_qc_scatter(so)
  expect_s3_class(r, "patchwork")
  hint <- attr(r, "size_in")
  expect_named(hint, c("width", "height"))

  many <- so
  many$sample <- rep(paste0("S", 1:12), length.out = ncol(many))
  expect_gt(attr(plot_qc_ridges(many), "size_in")[["height"]], hint[["height"]])
  expect_gt(attr(plot_qc_scatter(many), "size_in")[["height"]],
            attr(s, "size_in")[["height"]])
})

so_f <- subset(so, subset = !is_outlier)
pp <- quiet(seurat_preprocess(so_f, resolutions = c(0.2, 0.5)))

test_that("seurat_preprocess (log-normalization) adds reductions, clusters and a lean @misc", {
  expect_true(all(c("pca", "umap_unintegrated", "harmony", "umap_harmony") %in%
                    SeuratObject::Reductions(pp)))
  expect_equal(SeuratObject::DefaultDimReduc(pp), "umap_harmony")
  expect_true(all(c("leiden_0.2", "leiden_0.5") %in% names(pp[[]])))
  expect_setequal(names(pp@misc$pipeline$preprocess$params), c("npcs", "vars_to_regress"))
  expect_equal(SeuratObject::DefaultAssay(pp), "RNA")
})

test_that("seurat_preprocess works with SCT and without Harmony", {
  skip_if_not_installed("sctransform")
  sct <- quiet(seurat_preprocess(so_f, resolutions = 0.3, normalization = "sct",
                                 run_harmony = FALSE, npcs = 5))
  expect_equal(SeuratObject::DefaultAssay(sct), "SCT")
  expect_false("harmony" %in% SeuratObject::Reductions(sct))
  expect_equal(sct@misc$pipeline$preprocess$params$npcs, 5)
  expect_true("leiden_0.3" %in% names(sct[[]]))
})

# Drawing a plot (not just building the object) is what exposed a missing-guide bug once,
# so every returned plot is rendered to a file.
renders <- function(p) {
  f <- tempfile(fileext = ".png")
  ggplot2::ggsave(f, p, width = 6, height = 5, dpi = 40, limitsize = FALSE)
  file.exists(f)
}

test_that("UMAP plots return patchworks with size hints and draw", {
  for (p in list(plot_integration_umaps(pp), plot_qc_umaps(pp), plot_cluster_umaps(pp))) {
    expect_s3_class(p, "patchwork")
    expect_true(all(attr(p, "size_in") > 0))
    expect_no_error(renders(p))
  }
})

test_that("QC plots draw", {
  expect_no_error(renders(plot_qc_ridges(so)))
  expect_no_error(renders(plot_qc_scatter(so)))
})

test_that("choose_resolution returns tables and plots", {
  res <- quiet(choose_resolution(pp))
  expect_equal(res$table$resolution, c(0.2, 0.5))
  expect_true(all(c("table", "per_cluster", "clustree", "silhouette",
                    "silhouette_by_cluster") %in% names(res)))
  expect_named(attr(res$clustree, "size_in"), c("width", "height"))
  # clustree must draw even when ggraph is not attached
  expect_no_error(renders(res$clustree))
  expect_no_error(renders(res$silhouette))
  expect_no_error(renders(res$silhouette_by_cluster))
})

test_that("find_markers_all_resolutions returns markers, top markers and DotPlots", {
  mk <- quiet(find_markers_all_resolutions(pp, top_n = 3))
  expect_named(mk, c("markers", "top_markers", "dotplots"))
  expect_named(mk$dotplots, c("leiden_0.2", "leiden_0.5"))
  expect_lte(max(table(mk$top_markers[["leiden_0.5"]]$cluster)), 3)
  expect_named(attr(mk$dotplots[[1]], "size_in"), c("width", "height"))
})

test_that("the QualityControl folder gets tables, plots and parameters", {
  out <- tempfile("QualityControl_")
  m1 <- quiet(save_read_qc(so, out_dir = out, formats = "png", scale = 0.5))
  m1b <- quiet(save_read_qc(so_f, out_dir = out, tag = "postfilter", formats = "png", scale = 0.5))
  m2 <- quiet(save_preprocess_qc(pp, out_dir = out, markers = TRUE, top_n = 3,
                                 formats = "png", scale = 0.5))
  expected <- c("read_qc/qc_thresholds.csv", "read_qc/qc_flagged.csv",
                "read_qc/qc_ridges_prefilter.png", "read_qc/qc_scatter_prefilter.png",
                "read_qc/qc_ridges_postfilter.png",
                "preprocessing/elbow.png", "preprocessing/umap_integration.png",
                "preprocessing/umap_qc_metrics.png", "preprocessing/umap_clusters.png",
                "preprocessing/clustree.png", "preprocessing/silhouette_by_resolution.png",
                "preprocessing/silhouette_by_cluster.png",
                "preprocessing/resolution_summary.csv", "preprocessing/silhouette_per_cluster.csv",
                "preprocessing/markers/dotplot_leiden_0.2.png",
                "preprocessing/markers/dotplot_leiden_0.5.png",
                "preprocessing/markers/top_markers_leiden_0.5.csv",
                "pipeline_parameters.yml", "sessionInfo.txt")
  expect_true(all(file.exists(file.path(out, expected))),
              info = paste(expected[!file.exists(file.path(out, expected))], collapse = ", "))
  expect_true(all(c("file", "width", "height") %in% names(m2)))
})
