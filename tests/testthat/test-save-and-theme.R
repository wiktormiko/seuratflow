test_that(".panel_grid follows ggplot2's wrap rule", {
  expect_equal(seuratflow:::.panel_grid(5), c(nrow = 2, ncol = 3))
  expect_equal(seuratflow:::.panel_grid(3, ncol = 3), c(nrow = 1, ncol = 3))
})

test_that("save_plot uses the size hint, scale and keeps dots in file names", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  p <- seuratflow:::.set_size(p, width = 6, height = 3)
  out <- file.path(tempfile("plots_"), "dotplot_leiden_0.5")

  res <- save_plot(p, out, formats = c("png", "pdf"), dpi = 72)
  expect_equal(res$width, c(6, 6))
  expect_equal(res$height, c(3, 3))
  expect_true(all(file.exists(paste0(out, c(".png", ".pdf")))))   # "0.5" was not cut off

  res2 <- save_plot(p, out, formats = "png", dpi = 72, scale = 2)
  expect_equal(c(res2$width, res2$height), c(12, 6))

  res3 <- save_plot(p, out, width = 4, height = 4, formats = "png", dpi = 72)
  expect_equal(c(res3$width, res3$height), c(4, 4))
})

test_that("save_plot without a hint sizes plain ggplots from their facets, patchworks 8 x 6", {
  plain <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  faceted <- plain + ggplot2::facet_wrap(~cyl, ncol = 3)

  r <- save_plot(faceted, file.path(tempfile(), "faceted"), formats = "png", dpi = 72)
  expect_equal(c(r$width, r$height), c(3.5 * 3 + 1.5, 3 * 1 + 1))

  r1 <- save_plot(plain, file.path(tempfile(), "plain"), formats = "png", dpi = 72)
  expect_equal(c(r1$width, r1$height), c(3.5 + 1.5, 3 + 1))

  r2 <- save_plot(patchwork::wrap_plots(plain, plain), file.path(tempfile(), "pw"),
                  formats = "png", dpi = 72)
  expect_equal(c(r2$width, r2$height), c(8, 6))
})

test_that("umap_style sets the colour scale for discrete and continuous colours", {
  base <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()

  p_disc <- base + ggplot2::aes(colour = factor(cyl)) + umap_style()
  expect_no_error(ggplot2::ggplot_build(p_disc))

  p_cont <- base + ggplot2::aes(colour = hp) + umap_style()
  expect_no_error(ggplot2::ggplot_build(p_cont))

  # forcing the type works
  expect_no_error(ggplot2::ggplot_build(
    base + ggplot2::aes(colour = hp) + umap_style(type = "continuous")))
})

test_that("umap_style(fixed = FALSE) keeps a square aspect ratio", {
  base <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg, colour = hp)) + ggplot2::geom_point()
  expect_equal(unname((base + umap_style(fixed = FALSE))$theme$aspect.ratio), 1)
  expect_null((base + umap_style(fixed = TRUE))$theme$aspect.ratio)
})

test_that("theme_umap removes axis text, ticks and grid", {
  th <- theme_umap()
  expect_s3_class(th$axis.text, "element_blank")
  expect_s3_class(th$axis.ticks, "element_blank")
  expect_s3_class(th$panel.grid, "element_blank")
})
