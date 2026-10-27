#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom stats quantile setNames dist
#' @importFrom utils modifyList
#' @importFrom Seurat Read10X Read10X_h5 NormalizeData FindVariableFeatures ScaleData RunPCA
#' @importFrom Seurat RunUMAP SCTransform PrepSCTFindMarkers FindNeighbors FindClusters
#' @importFrom Seurat FindAllMarkers DotPlot DimPlot FeaturePlot ElbowPlot RotatedAxis NoLegend
#' @importFrom SeuratObject CreateSeuratObject AddMetaData GetAssayData DefaultAssay
#' @importFrom SeuratObject "DefaultAssay<-" Idents "Idents<-" Layers JoinLayers Embeddings
#' @importFrom SeuratObject Reductions Stdev DefaultDimReduc "DefaultDimReduc<-"
#' @importFrom scuttle isOutlier
#' @importFrom scDblFinder scDblFinder
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scCustomize Add_Cell_QC_Metrics Read_CellBender_h5_Mat
#' @importFrom harmony RunHarmony
#' @importFrom findPC findPC
#' @importFrom readxl read_xlsx
#' @importFrom readr read_csv
## usethis namespace: end
NULL

# Column names used inside dplyr / ggplot2 expressions
utils::globalVariables(c(
  "value", "sample", "x", "y", "cluster", "p_val_adj", "sil_width", "resolution",
  "mean_silhouette", "resolution_mean", "is_outlier", "umi_outlier", "ngenes_outlier",
  "mt_outlier", "doublet", "n_clusters", "ndensity", "nFeature_RNA", "nCount_RNA",
  "percent_mito"
))

`%||%` <- function(x, y) if (is.null(x)) y else x
