# =============================================================================
# Unit Tests: Base Graphics Device Safety
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "plotting.R"))

test_that("with_png_device closes failed devices and removes partial files", {
  root <- tempfile("plotting_device_")
  dir.create(root)
  failed <- file.path(root, "failed.png")
  original_device <- grDevices::dev.cur()
  expect_error(with_png_device(failed, function() stop("draw failure")), "draw failure")
  expect_equal(grDevices::dev.cur(), original_device)
  expect_false(file.exists(failed))

  complete <- file.path(root, "complete.png")
  expect_invisible(with_png_device(complete, function() graphics::plot(1, 1)))
  expect_true(file.exists(complete))
  expect_gt(file.info(complete)$size, 0)
})
