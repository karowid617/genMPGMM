# ============================================================
# Tests for genMPGMM
# ============================================================

library(genMPGMM)

# --------------- helpers ---------------

quick_call <- function(seed = 42, M = 40, N = 60,
                       noise_feature_fraction = 0,
                       target_ari     = c(1, 0.3),
                       covariance_spec = list(type = "diagonal"),
                       ...) {
  genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(2, 2),
    feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
    mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
    dist_mahalanobis          = 3,
    target_ari                = target_ari,
    M                         = M,
    N                         = N,
    covariance_spec           = covariance_spec,
    noise_feature_fraction    = noise_feature_fraction,
    seed                      = seed,
    ...
  )
}

# ============================================================
# I1–I2: Matrix dimensions
# ============================================================

test_that("I1/I2: output matrix dimensions are correct", {
  res <- quick_call(M = 50, N = 80)
  expect_equal(dim(res$X),        c(50L, 80L))
  expect_equal(dim(res$X_signal), c(50L, 80L))
})

# ============================================================
# I3–I4: Feature partitions are well-formed
# ============================================================

test_that("I3/I4: feature partitions have correct length and labels", {
  res <- quick_call(M = 40, N = 60)
  for (p in seq_len(res$settings$P)) {
    L_p <- res$settings$n_feature_patterns[p]
    expect_length(res$s[[p]], 40)
    expect_true(all(res$s[[p]] %in% seq_len(L_p)),
                info = sprintf("profile %d: labels out of range", p))
    expect_true(all(seq_len(L_p) %in% res$s[[p]]),
                info = sprintf("profile %d: some labels missing entirely", p))
  }
})

# ============================================================
# I5–I6: Observation partitions are well-formed
# ============================================================

test_that("I5/I6: observation partitions have correct length and labels", {
  res <- quick_call(M = 40, N = 60)
  for (p in seq_len(res$settings$P)) {
    K_p <- res$settings$n_components[p]
    expect_length(res$z[[p]], 60)
    expect_true(all(res$z[[p]] %in% seq_len(K_p)))
  }
})

# ============================================================
# I7: Mixing proportions respected (within 1/N rounding error)
# ============================================================

test_that("I7: mixing proportions are approximately correct", {
  mp <- list(c(0.2, 0.8), c(0.6, 0.4))
  res <- genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(2, 2),
    feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
    mixing_proportions        = mp,
    dist_mahalanobis          = 3,
    target_ari                = c(1, 0.3),
    M = 40, N = 200,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = 0,
    seed = 1
  )
  for (p in 1:2) {
    achieved <- res$achieved_mixing[[p]]
    expect_equal(achieved, mp[[p]], tolerance = 1 / 200 + 1e-9,
                 label = sprintf("mixing proportions profile %d", p))
  }
})

# ============================================================
# I8: Feature-group proportions respected
# ============================================================

test_that("I8: feature-group proportions are approximately correct", {
  fp <- list(c(0.3, 0.7), c(0.4, 0.6))
  res <- genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(2, 2),
    feature_group_proportions = fp,
    mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
    dist_mahalanobis          = 3,
    target_ari                = c(1, 0.3),
    M = 200, N = 60,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = 0,
    seed = 1
  )
  for (p in 1:2) {
    achieved <- res$achieved_feature_group_props[[p]]
    expect_equal(achieved, fp[[p]], tolerance = 1 / 200 + 1e-9,
                 label = sprintf("feature proportions profile %d", p))
  }
})

# ============================================================
# I9: ARI vs_reference — target achieved within tolerance
# ============================================================

test_that("I9: ARI vs_reference is within tolerance for feasible targets", {
  targets <- c(0.0, 0.3, 0.6)
  for (tgt in targets) {
    res <- suppressWarnings(genMPGMM(
      n_feature_patterns        = c(2, 2),
      n_components              = c(2, 2),
      feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
      mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
      dist_mahalanobis          = 3,
      target_ari                = c(1, tgt),
      ari_tol                   = 0.02,
      ari_max_iter              = 10000,
      M = 80, N = 60,
      covariance_spec           = list(type = "diagonal"),
      noise_feature_fraction    = 0,
      seed = 7
    ))
    achieved <- res$achieved_ari_features[1, 2]
    # Absolute tolerance: ARI lives in [-1, 1], so relative tolerance is not meaningful
    expect_true(
      abs(achieved - tgt) < 0.06,
      label = sprintf("ARI target %.2f: achieved %.4f (diff = %.4f)", tgt, achieved, abs(achieved - tgt))
    )
  }
})

# ============================================================
# I10: Mahalanobis distance — K_p = 2 case is exact
# ============================================================

test_that("I10: Mahalanobis distance K_p=2 matches target exactly", {
  for (tgt in c(2, 5, 10)) {
    res <- genMPGMM(
      n_feature_patterns        = c(2, 2),
      n_components              = c(2, 2),
      feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
      mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
      dist_mahalanobis          = tgt,
      target_ari                = c(1, 0.3),
      M = 40, N = 60,
      covariance_spec           = list(type = "diagonal"),
      noise_feature_fraction    = 0,
      seed = 99
    )
    D <- res$achieved_mahalanobis[[1]]
    expect_equal(D[1, 2], tgt, tolerance = 1e-6,
                 label = sprintf("Mahalanobis distance target %g", tgt))
  }
})

# ============================================================
# I10b: Mahalanobis distance — K_p > 2: median achieves target
# ============================================================

test_that("I10b: median Mahalanobis distance K_p>2 matches target", {
  tgt <- 4
  res <- genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(3, 3),
    feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
    mixing_proportions        = list(rep(1/3, 3), rep(1/3, 3)),
    dist_mahalanobis          = tgt,
    target_ari                = c(1, 0.3),
    M = 60, N = 90,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = 0,
    seed = 11
  )
  D     <- res$achieved_mahalanobis[[1]]
  dvals <- D[upper.tri(D)]
  expect_equal(median(dvals), tgt, tolerance = 1e-6)
})

# ============================================================
# I11: Covariance matrices are positive definite
# ============================================================

test_that("I11: all profile covariance matrices are positive definite", {
  res <- genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(2, 2),
    feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
    mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
    dist_mahalanobis          = 3,
    target_ari                = c(1, 0.3),
    M = 10, N = 40,
    covariance_spec           = list(type = "full_shared"),
    noise_feature_fraction    = 0,
    seed = 5
  )
  for (p in seq_len(res$settings$P)) {
    evals <- eigen(res$cov_mtx[[p]], symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(evals > 0),
                info = sprintf("profile %d covariance is not positive definite", p))
  }
})

# ============================================================
# I12: Diagonal mode produces diagonal covariance
# ============================================================

test_that("I12: diagonal covariance spec produces diagonal Sigma", {
  res <- quick_call(M = 10, N = 40,
                    covariance_spec = list(type = "diagonal", diag_values = rep(2, 10)))
  for (p in seq_len(res$settings$P)) {
    S <- res$cov_mtx[[p]]
    off <- S[lower.tri(S)]
    expect_true(all(off == 0), info = sprintf("profile %d: off-diagonal not zero", p))
    expect_equal(diag(S), rep(2, 10))
  }
})

# ============================================================
# I13: K_p = 1 edge case
# ============================================================

test_that("I13: K_p=1 gives zero Mahalanobis distance and one component", {
  res <- genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(1, 1),
    feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
    mixing_proportions        = list(c(1), c(1)),
    dist_mahalanobis          = 5,
    target_ari                = c(1, 0.3),
    M = 30, N = 40,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = 0,
    seed = 3
  )
  expect_equal(res$achieved_mahalanobis[[1]], matrix(0, 1, 1))
  expect_equal(res$achieved_mahalanobis[[2]], matrix(0, 1, 1))
  expect_equal(res$achieved_mixing[[1]], c(1))
})

# ============================================================
# I14: Reproducibility via set.seed
# ============================================================

test_that("I14: same seed gives identical results", {
  r1 <- quick_call(seed = 17)
  r2 <- quick_call(seed = 17)
  expect_identical(r1$X, r2$X)
  expect_identical(r1$s, r2$s)
  expect_identical(r1$z, r2$z)
})

test_that("I14b: different seeds give different X", {
  r1 <- quick_call(seed = 1)
  r2 <- quick_call(seed = 2)
  expect_false(identical(r1$X, r2$X))
})

# ============================================================
# I15: Noise features are unrelated to observation partition
# ============================================================

test_that("I15: noise features have no strong correlation with z", {
  res <- genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(2, 2),
    feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
    mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
    dist_mahalanobis          = 10,
    target_ari                = c(1, 0.3),
    M = 100, N = 200,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = 0.2,
    seed = 8
  )
  nf <- res$noise_feature_indices
  if (length(nf) > 0) {
    z_bin <- as.numeric(res$z[[1]])
    cors  <- apply(res$X[nf, , drop = FALSE], 1, function(row) abs(cor(row, z_bin)))
    expect_true(all(cors < 0.3),
                info = paste("Some noise features show unexpectedly high correlation with z:",
                             paste(round(cors[cors >= 0.3], 3), collapse = ", ")))
  }
})

# ============================================================
# I16: Invalid inputs raise informative errors
# ============================================================

test_that("I16a: P < 2 raises an error", {
  expect_error(
    genMPGMM(n_feature_patterns=c(2), n_components=c(2),
             feature_group_proportions=list(c(.5,.5)),
             mixing_proportions=list(c(.5,.5)),
             dist_mahalanobis=3, target_ari=c(1),
             M=20, N=20, covariance_spec=list(type="diagonal"),
             noise_feature_fraction=0),
    regexp = "P.*>=.*2|>= 2"
  )
})

test_that("I16b: mismatched feature_group_proportions raises error", {
  expect_error(
    genMPGMM(n_feature_patterns=c(2,2), n_components=c(2,2),
             feature_group_proportions=list(c(.5,.5,.5), c(.5,.5)),  # wrong length
             mixing_proportions=list(c(.5,.5),c(.5,.5)),
             dist_mahalanobis=3, target_ari=c(1,.3),
             M=20, N=20, covariance_spec=list(type="diagonal"),
             noise_feature_fraction=0),
    regexp = "feature_group_proportions"
  )
})

test_that("I16c: pairwise_matrix with non-symmetric target raises error", {
  asym <- matrix(c(1, 0.3, 0.7, 1), nrow = 2)
  expect_error(
    genMPGMM(n_feature_patterns=c(2,2), n_components=c(2,2),
             feature_group_proportions=list(c(.5,.5),c(.5,.5)),
             mixing_proportions=list(c(.5,.5),c(.5,.5)),
             dist_mahalanobis=3, target_ari=asym,
             ari_mode="pairwise_matrix",
             M=20, N=20, covariance_spec=list(type="diagonal"),
             noise_feature_fraction=0),
    regexp = "symmetric"
  )
})

test_that("I16d: pairwise_matrix with diagonal != 1 raises error", {
  bad_diag <- matrix(c(0, 0.5, 0.5, 0), nrow = 2)
  expect_error(
    genMPGMM(n_feature_patterns=c(2,2), n_components=c(2,2),
             feature_group_proportions=list(c(.5,.5),c(.5,.5)),
             mixing_proportions=list(c(.5,.5),c(.5,.5)),
             dist_mahalanobis=3, target_ari=bad_diag,
             ari_mode="pairwise_matrix",
             M=20, N=20, covariance_spec=list(type="diagonal"),
             noise_feature_fraction=0),
    regexp = "diagonal"
  )
})

# ============================================================
# I17: X_signal == X when add_noise=FALSE and noise_feature_fraction=0
# ============================================================

test_that("I17: X equals X_signal when no noise is applied", {
  res <- quick_call(add_noise = FALSE, noise_feature_fraction = 0)
  expect_identical(res$X, res$X_signal)
})

# ============================================================
# Fix A: vs_reference target_ari[1] != 1 triggers a warning
# ============================================================

test_that("Fix-A: vs_reference target_ari[1] != 1 emits a warning", {
  # target_ari of length P with [1] != 1 should warn (it is silently ignored)
  expect_warning(
    genMPGMM(
      n_feature_patterns        = c(2, 2),
      n_components              = c(2, 2),
      feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
      mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
      dist_mahalanobis          = 3,
      target_ari                = c(0.5, 0.3),   # [1] is not 1
      M = 40, N = 60,
      covariance_spec           = list(type = "diagonal"),
      noise_feature_fraction    = 0, seed = 42
    ),
    regexp = "ignored|always 1"
  )
})

test_that("Fix-A: vs_reference target_ari[1] = 1 produces no spurious warning", {
  # Only warnings from ARI non-convergence are acceptable; no "ignored" warning
  warns <- character(0)
  withCallingHandlers(
    genMPGMM(
      n_feature_patterns        = c(2, 2),
      n_components              = c(2, 2),
      feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
      mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
      dist_mahalanobis          = 3,
      target_ari                = c(1, 0.3),  # correct value
      M = 40, N = 60,
      covariance_spec           = list(type = "diagonal"),
      noise_feature_fraction    = 0, seed = 42
    ),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  ignored_warns <- grepl("ignored|always 1", warns)
  expect_false(any(ignored_warns))
})

# ============================================================
# Fix C: ARI convergence warning when target is unreachable
# ============================================================

test_that("Fix-C: convergence warning when ARI target is geometrically impossible", {
  # 50/50 vs 90/10 proportions → ARI = 1 is impossible
  expect_warning(
    genMPGMM(
      n_feature_patterns        = c(2, 2),
      n_components              = c(2, 2),
      feature_group_proportions = list(c(0.5, 0.5), c(0.9, 0.1)),
      mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
      dist_mahalanobis          = 3,
      target_ari                = c(1, 1.0),
      ari_max_iter              = 500,
      ari_tol                   = 0.01,
      M = 60, N = 60,
      covariance_spec           = list(type = "diagonal"),
      noise_feature_fraction    = 0,
      seed = 1
    ),
    regexp = "converge|ARI"
  )
})

# ============================================================
# Additive noise sanity check
# ============================================================

test_that("add_noise increases per-cell variance", {
  r_no  <- quick_call(add_noise = FALSE, noise_feature_fraction = 0)
  r_yes <- quick_call(add_noise = TRUE,  noise_sd = 1, noise_feature_fraction = 0)
  expect_gt(var(as.vector(r_yes$X)), var(as.vector(r_no$X)))
})

# ============================================================
# noise_feature_indices: correct count and valid range
# ============================================================

test_that("noise_feature_indices length and range are correct", {
  M   <- 100
  frc <- 0.15
  res <- genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(2, 2),
    feature_group_proportions = list(c(.5,.5),c(.5,.5)),
    mixing_proportions        = list(c(.5,.5),c(.5,.5)),
    dist_mahalanobis          = 3,
    target_ari                = c(1, 0.3),
    M = M, N = 60,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = frc,
    seed = 6
  )
  expect_equal(length(res$noise_feature_indices), floor(M * frc))
  expect_true(all(res$noise_feature_indices %in% seq_len(M)))
})

# ============================================================
# pairwise_matrix mode: correct target achieved
# ============================================================

test_that("pairwise_matrix mode achieves correct ARI for P=2", {
  tgt <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  res <- genMPGMM(
    n_feature_patterns        = c(2, 2),
    n_components              = c(2, 2),
    feature_group_proportions = list(c(.5,.5),c(.5,.5)),
    mixing_proportions        = list(c(.5,.5),c(.5,.5)),
    dist_mahalanobis          = 3,
    target_ari                = tgt,
    ari_mode                  = "pairwise_matrix",
    ari_tol                   = 0.02,
    ari_max_iter              = 5000,
    M = 80, N = 60,
    covariance_spec           = list(type = "diagonal"),
    noise_feature_fraction    = 0,
    seed = 55
  )
  expect_equal(res$achieved_ari_features[1, 2], 0.4, tolerance = 0.05, scale = 1)
})

# ============================================================
# profile_weights: different weights change output
# ============================================================

test_that("profile_weights change the combined X", {
  r1 <- quick_call(profile_weights = c(1, 1))
  r2 <- quick_call(profile_weights = c(1, 0))
  expect_false(identical(r1$X, r2$X))
})
