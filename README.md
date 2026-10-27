# seuratflow

Helpers that automate a **Seurat v5** single-cell RNA-seq preprocessing workflow and write every
quality control table and plot to a `QualityControl/` folder, with figure sizes chosen
automatically.

- **Read + QC in one step.** 10X matrix directories, 10X h5 files and CellBender output are read
  in parallel; each sample gets QC metrics, MAD-based outlier flags and (optionally) scDblFinder
  doublet calls, then everything is merged.
- **Preprocess.** Log-normalization or SCTransform, automatic number of PCs, Harmony, and Leiden
  clustering at several resolutions.
- **Choose a resolution.** clustree plus overall and per-cluster silhouette scores.
- **Markers.** Top markers and a DotPlot for every resolution.
- **Reports.** QC plots, UMAPs, tables, parameters and `sessionInfo()` saved to one folder.
- **ggplot2 style for UMAPs.** `theme_umap()` and `umap_style()`.

## Installation

Requires R >= 4.1 and Seurat v5. seuratflow uses packages from CRAN, Bioconductor and GitHub.
`remotes::install_github()` installs the CRAN ones for you, but installing everything first
(below) avoids the usual failures from the Bioconductor and GitHub packages and from packages
that some minimal R installations lack (for example `cluster` and `ggraph`).

```r
# 1. CRAN
install.packages(c(
  "remotes", "BiocManager",
  "Seurat", "SeuratObject", "scCustomize", "harmony",            # single-cell
  "ggplot2", "patchwork", "ggraph", "clustree", "ggridges", "RColorBrewer",   # plots
  "dplyr", "tidyr", "readr", "readxl", "yaml", "rlang", "Matrix",
  "future", "future.apply", "cluster"
))

# 2. Bioconductor
BiocManager::install(c("scDblFinder", "scuttle", "SingleCellExperiment"))

# 3. GitHub
remotes::install_github("haotian-zhuang/findPC")   # picks the number of PCs

# 4. seuratflow
remotes::install_github("wiktormiko/seuratflow")
```

Check that nothing is missing (the result should be `character(0)`):

```r
required <- c("Seurat", "SeuratObject", "scCustomize", "harmony", "ggplot2", "patchwork",
              "ggraph", "clustree", "ggridges", "RColorBrewer", "dplyr", "tidyr", "readr",
              "readxl", "yaml", "rlang", "Matrix", "future", "future.apply", "cluster",
              "scDblFinder", "scuttle", "SingleCellExperiment", "findPC")
setdiff(required, rownames(installed.packages()))
```

Optional packages:

| package | where | used for |
|---|---|---|
| `hdf5r` | CRAN | `read10X_h5_multi_directory_add_QC()` and `read_cellbender_multi_directory_add_QC()` |
| `presto` | GitHub: `immunogenomics/presto` | much faster marker detection |
| `glmGamPoi` | Bioconductor | much faster `seurat_preprocess(normalization = "sct")` |

If you reinstall seuratflow while R is running, restart R before using it again; otherwise R
can report a corrupt lazy-load database.

## Quick start

`write_example_10x()` simulates a small dataset so you can try everything first.

```r
library(seuratflow)

d <- file.path(tempdir(), "example_10x")
write_example_10x(d)                      # 3 samples, sample sheet in d/metadata.csv

# 1. Read + QC (parallel over samples)
so <- read10X_multi_directory_add_QC(
  d, file.path(d, "metadata.csv"),
  default_10X_path = FALSE,               # matrices sit directly in each sample folder
  workers = 3,
  detect_doublets = TRUE                  # FALSE skips scDblFinder
)

# 2. Look at the cutoffs, then filter
plot_qc_ridges(so)                        # red dashed lines = the MAD cutoffs
plot_qc_scatter(so)
so_f <- subset(so, subset = !is_outlier)

# 3. Preprocess: normalize, PCA, UMAP, Harmony, Leiden at several resolutions
so_f <- seurat_preprocess(so_f, resolutions = c(0.2, 0.5, 0.8),
                          normalization = "lognorm")        # or "sct"

# 4. Diagnostics
plot_integration_umaps(so_f)              # before / after Harmony
plot_qc_umaps(so_f)                       # QC metrics on the UMAP
plot_cluster_umaps(so_f)                  # clusters at every resolution
res <- choose_resolution(so_f)            # clustree + silhouette
res$clustree; res$silhouette_by_cluster

# 5. Markers + DotPlot for every resolution
mk <- find_markers_all_resolutions(so_f, top_n = 5)
mk$dotplots[["leiden_0.5"]]

# 6. Save everything
save_read_qc(so)                          # tables + QC plots before filtering
save_read_qc(so_f, tag = "postfilter")    # the same plots after filtering
save_preprocess_qc(so_f, markers = TRUE)  # UMAPs, resolution choice, markers
```

### Other inputs

```r
# 10X h5: <base>/<sample>/outs/filtered_feature_bc_matrix.h5
so <- read10X_h5_multi_directory_add_QC(base_path, sample_sheet_path)

# CellBender: <base>/<sample>/<sample>_out_filtered.h5
so <- read_cellbender_multi_directory_add_QC(base_path, sample_sheet_path)
```

All three readers accept the same QC arguments (`species`, `nmads`, `mito_floor`,
`detect_doublets`, `workers`, `seed`). Old `genes.tsv.gz` files are read through a temporary
symlinked copy, so the raw data is never renamed.

## The `QualityControl/` folder

```
QualityControl/
├── read_qc/
│   ├── qc_thresholds.csv            MAD cutoffs per sample
│   ├── qc_flagged.csv               cells flagged per sample and why
│   ├── qc_ridges_prefilter.png/pdf
│   ├── qc_scatter_prefilter.png/pdf
│   └── qc_ridges_postfilter.png/pdf ...
├── preprocessing/
│   ├── elbow.png/pdf
│   ├── umap_integration.png/pdf     before / after Harmony
│   ├── umap_qc_metrics.png/pdf
│   ├── umap_clusters.png/pdf
│   ├── clustree.png/pdf
│   ├── silhouette_by_resolution.png/pdf
│   ├── silhouette_by_cluster.png/pdf
│   ├── resolution_summary.csv
│   ├── silhouette_per_cluster.csv
│   └── markers/                     all_markers_*.csv, top_markers_*.csv, dotplot_*.png/pdf
├── pipeline_parameters.yml
└── sessionInfo.txt
```

Each output is its own step: if one fails you get a warning and the rest is still written.

### Automatic figure size

Every plot function attaches a size hint computed from what it draws (samples, panels,
clusters, genes, legend entries), and `save_plot()` reads it:

```r
save_plot(plot_qc_ridges(so), "figures/ridges")             # png + pdf, size chosen for you
save_plot(plot_qc_ridges(so), "figures/ridges", scale = 1.2) # or scale / override width, height
```

Adding layers to a returned plot drops the hint; plain ggplots are then sized from their
facets and patchworks get 8 x 6 inches.

## What is stored in `@misc`

Only what is needed to interpret the object later:

| slot | content |
|---|---|
| `qc_thresholds` | MAD cutoffs per sample |
| `qc_flagged` | cells and flags per sample, counted before filtering |
| `pipeline$read_qc$params` | `nmads`, `mito_floor`, `detect_doublets` |
| `pipeline$preprocess$params` | `npcs`, `vars_to_regress` |

## UMAP styling

```r
library(Seurat)
DimPlot(so_f, label = TRUE) & umap_style()                  # scCustomize palette for clusters
FeaturePlot(so_f, c("CD3E", "MS4A1")) & umap_style()        # reversed viridis "rocket"
DimPlot(so_f) & umap_style(fixed = FALSE)                   # square panels instead of coord_fixed()

ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point() + theme_umap()
```

`umap_style(fixed = TRUE)` (the default) uses `coord_fixed()`, so each panel's shape follows its
embedding; use `fixed = FALSE` when comparing panels side by side.

## Functions

| Step | Functions |
|---|---|
| Read + QC | `read10X_multi_directory_add_QC()`, `read10X_h5_multi_directory_add_QC()`, `read_cellbender_multi_directory_add_QC()` |
| QC plots | `plot_qc_ridges()`, `plot_qc_scatter()` |
| Preprocess | `seurat_preprocess()` |
| UMAP plots | `plot_integration_umaps()`, `plot_qc_umaps()`, `plot_cluster_umaps()` |
| Resolution + markers | `choose_resolution()`, `find_markers_all_resolutions()` |
| Saving | `save_plot()`, `save_read_qc()`, `save_preprocess_umaps()`, `save_preprocess_qc()` |
| Styling | `theme_umap()`, `umap_style()` |
| Example data | `write_example_10x()` |

## License

MIT © Wiktor Mikolajczak
