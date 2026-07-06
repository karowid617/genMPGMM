# ============================================================
# Tests for plot.MPObject and report_MPObject
# ============================================================

library(genMPGMM)

# Skip all tests silently if visualisation deps are absent
skip_if_no_plot_deps <- function() {
  for (pkg in c("ggplot2", "gridExtra", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste("Package", pkg, "not available"))
    }
  }
}

# A small, fast fixture used by every test
small_dat <- suppressWarnings(genMPGMM(
  n_feature_patterns        = c(2, 2),
  n_components              = c(2, 2),
  feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
  mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
  dist_mahalanobis          = 3,
  target_ari                = c(1, 0.3),
  M = 20, N = 30,
  covariance_spec           = list(type = "diagonal"),
  noise_feature_fraction    = 0,
  seed = 1
))

# ============================================================
# plot.MPObject -- return value
# ============================================================

test_that("plot.MPObject returns invisible list of ggplots", {
  skip_if_no_plot_deps()
  plots <- suppressWarnings(
    plot(small_dat, ask = FALSE)
  )
  expect_type(plots, "list")
  expect_length(plots, 7)
  expect_named(plots, c("heatmap", "feature_partitions",
                         "observation_partitions", "pca",
                         "proportions", "mahalanobis",
                         "component_means"))
})

test_that("plot.MPObject: each top-level entry is a ggplot or list of ggplots", {
  skip_if_no_plot_deps()
  plots <- suppressWarnings(plot(small_dat, ask = FALSE))
  for (nm in names(plots)) {
    entry <- plots[[nm]]
    ok <- inherits(entry, "gg") ||
      (is.list(entry) && all(vapply(entry, inherits, logical(1), "gg")))
    expect_true(ok, label = paste("Entry", nm, "is ggplot or list of ggplots"))
  }
})

# ============================================================
# plot.MPObject -- which argument
# ============================================================

test_that("plot.MPObject: which = integer subset returns correct plots", {
  skip_if_no_plot_deps()
  plots <- suppressWarnings(plot(small_dat, which = c(1, 6), ask = FALSE))
  expect_named(plots, c("heatmap", "mahalanobis"))
})

test_that("plot.MPObject: which = character names works", {
  skip_if_no_plot_deps()
  plots <- suppressWarnings(
    plot(small_dat, which = c("pca", "proportions"), ask = FALSE)
  )
  expect_named(plots, c("pca", "proportions"))
})

test_that("plot.MPObject: invalid integer index raises error", {
  skip_if_no_plot_deps()
  expect_error(
    plot(small_dat, which = 99, ask = FALSE),
    regexp = "out of range"
  )
})

test_that("plot.MPObject: unknown name raises error", {
  skip_if_no_plot_deps()
  expect_error(
    plot(small_dat, which = "banana", ask = FALSE),
    regexp = "Unknown plot name"
  )
})

# ============================================================
# plot.MPObject -- K_p = 1 edge case
# ============================================================

test_that("plot.MPObject: K_p=1 completes without error", {
  skip_if_no_plot_deps()
  dat_k1 <- suppressWarnings(genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(1, 1),
    feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
    mixing_proportions        = list(c(1), c(1)),
    dist_mahalanobis          = 3,
    target_ari                = c(1, 0.3),
    M = 20, N = 30,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = 0, seed = 2
  ))
  expect_no_error(suppressWarnings(plot(dat_k1, ask = FALSE)))
})

# ============================================================
# plot.MPObject -- noise features present
# ============================================================

test_that("plot.MPObject: noise features do not break heatmap", {
  skip_if_no_plot_deps()
  dat_noise <- suppressWarnings(genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(2, 2),
    feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
    mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
    dist_mahalanobis          = 3,
    target_ari                = c(1, 0.3),
    M = 30, N = 40,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = 0.2, seed = 3
  ))
  expect_no_error(suppressWarnings(
    plot(dat_noise, which = "heatmap", ask = FALSE)
  ))
})

# ============================================================
# report_MPObject -- saves file
# ============================================================

test_that("report_MPObject creates a PDF file", {
  skip_if_no_plot_deps()
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp))
  suppressWarnings(
    report_MPObject(small_dat, file = tmp, open = FALSE)
  )
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 1000L)
})

test_that("report_MPObject: which subset writes smaller PDF than full report", {
  skip_if_no_plot_deps()
  tmp_full    <- tempfile(fileext = ".pdf")
  tmp_partial <- tempfile(fileext = ".pdf")
  on.exit({ unlink(tmp_full); unlink(tmp_partial) })

  suppressWarnings(report_MPObject(small_dat, file = tmp_full,    which = "all",   open = FALSE))
  suppressWarnings(report_MPObject(small_dat, file = tmp_partial, which = c(1, 6), open = FALSE))

  expect_gt(file.size(tmp_full), file.size(tmp_partial))
})

test_that("report_MPObject returns path invisibly", {
  skip_if_no_plot_deps()
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp))
  result <- suppressWarnings(report_MPObject(small_dat, file = tmp, open = FALSE))
  expect_equal(normalizePath(result), normalizePath(tmp))
})

# ============================================================
# report_MPObject -- multi-profile fixture (P=5, mixed K)
# ============================================================

test_that("report_MPObject handles complex multi-profile fixture without error", {
  skip_if_no_plot_deps()
  complex_dat <- suppressWarnings(genMPGMM(
    n_feature_patterns        = c(2, 5, 3, 4, 2),
    n_components              = c(2, 5, 2, 3, 3),
    feature_group_proportions = list(
      c(0.5, 0.5),
      c(0.3, 0.2, 0.2, 0.15, 0.15),
      c(0.4, 0.3, 0.3),
      c(0.25, 0.25, 0.25, 0.25),
      c(0.6, 0.4)
    ),
    mixing_proportions = list(
      c(0.5, 0.5),
      c(0.2, 0.2, 0.2, 0.2, 0.2),
      c(0.5, 0.5),
      c(1/3, 1/3, 1/3),
      c(1/3, 1/3, 1/3)
    ),
    dist_mahalanobis       = c(3, 4, 3, 5, 4),
    target_ari             = c(1, 0.5, 0.8, 0.2, 0.6),
    M = 40, N = 60,
    covariance_spec        = list(type = "diagonal"),
    noise_feature_fraction = 0,
    seed = 42
  ))
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp))
  expect_no_error(suppressWarnings(
    report_MPObject(complex_dat, file = tmp, open = FALSE)
  ))
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 1000L)
})

# ============================================================
# Dependency guard: informative error when packages missing
# ============================================================

test_that("plot.MPObject fails informatively if deps missing (mocked)", {
  skip_if_no_plot_deps()
  fn_body <- deparse(body(genMPGMM:::.mp_check_deps))
  expect_true(any(grepl("install.packages", fn_body)))
})
