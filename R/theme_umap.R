#' Minimal theme for UMAP and other embedding plots
#'
#' [ggplot2::theme_light()] without axis text, ticks and grid.
#'
#' @param base_size,base_family Passed to [ggplot2::theme_light()].
#' @return A ggplot2 theme.
#' @seealso [umap_style()]
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_umap()
#' @export
theme_umap <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}

#' One-call style for UMAP plots
#'
#' Add it to a Seurat plot (or to a patchwork with `&`) to set the colour scale, coordinate
#' system and theme in one go. Continuous colours (FeaturePlot) use reversed viridis "rocket";
#' discrete colours (DimPlot) use the scCustomize rule: `"polychrome"` for up to 36 levels,
#' shuffled `"varibow"` above that.
#'
#' With `fixed = TRUE` the panels use `coord_fixed()`, so their shape follows the x and y range
#' of each embedding (panels of different embeddings differ in size). With `fixed = FALSE` they
#' are square (`aspect.ratio = 1`), which is what you want when comparing panels side by side.
#'
#' @param continuous_scale Scale used when colour is mapped to a numeric variable.
#' @param discrete_scale Scale used when colour is mapped to a factor. `NULL` uses the
#'   scCustomize rule above.
#' @param type `"auto"` detects the scale type from the plot; or force `"continuous"` or
#'   `"discrete"`.
#' @param fixed `TRUE` for `coord_fixed()`, `FALSE` for square panels.
#' @param base_size,base_family Passed to [theme_umap()].
#' @return An object of class `umap_style` that is added to a plot with `+` (or `&`).
#' @seealso [theme_umap()]
#' @examples
#' \dontrun{
#' Seurat::DimPlot(so, label = TRUE) & umap_style()
#' Seurat::FeaturePlot(so, c("CD3E", "MS4A1")) & umap_style(fixed = FALSE)
#' }
#' @export
umap_style <- function(
    continuous_scale = scale_colour_viridis_c(option = "rocket", direction = -1),
    discrete_scale = NULL,
    type = c("auto", "continuous", "discrete"),
    fixed = TRUE,
    base_size = 11,
    base_family = ""
) {
  th <- theme_umap(base_size = base_size, base_family = base_family)
  if (!fixed) th <- th + theme(aspect.ratio = 1)   # adding to a complete theme keeps it complete
  
  structure(
    list(
      continuous_scale = continuous_scale,
      discrete_scale = discrete_scale,
      type = match.arg(type),
      coord = if (fixed) coord_fixed(),            # NULL when fixed = FALSE
      theme = th
    ),
    class = "umap_style"
  )
}

#' Scale type and levels of a plot's colour scale
#' @noRd
.umap_colour_info <- function(plot) {
  sc <- tryCatch(
    ggplot_build(plot)$plot$scales$get_scales("colour"),
    error = function(e) NULL
  )
  if (is.null(sc)) {
    return(list(type = "none", levels = NULL))
  }
  if (inherits(sc, "ScaleDiscrete")) {
    list(type = "discrete", levels = sc$get_limits())
  } else {
    list(type = "continuous", levels = NULL)
  }
}

#' scCustomize default palette rule applied to the plotted levels
#' @noRd
.sccustom_discrete_scale <- function(levels) {
  n <- length(levels)
  pal <- if (n <= 36) {
    scCustomize::DiscretePalette_scCustomize(
      num_colors = n,
      palette = "polychrome"
    )
  } else {
    scCustomize::DiscretePalette_scCustomize(
      num_colors = n,
      palette = "varibow",
      shuffle_pal = TRUE
    )
  }
  # name colours by level so mapping is explicit (scCustomize may return
  # names that don't match your cluster labels)
  scale_colour_manual(values = setNames(unname(pal), levels))
}

#' @exportS3Method ggplot2::ggplot_add
#' @noRd
ggplot_add.umap_style <- function(object, plot, object_name, ...) {
  info <- .umap_colour_info(plot)
  type <- if (object$type == "auto") info$type else object$type

  scale <- switch(
    type,
    continuous = object$continuous_scale,
    discrete = if (!is.null(object$discrete_scale)) {
      object$discrete_scale
    } else if (!is.null(info$levels)) {
      .sccustom_discrete_scale(info$levels)
    },
    none = NULL
  )
  suppressMessages(plot + list(scale, object$coord, object$theme))
}
