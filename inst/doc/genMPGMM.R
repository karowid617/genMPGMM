## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----logo, echo=FALSE, out.width="160px", out.extra='style="float:right; padding:10px"'----
knitr::include_graphics("../man/figures/logo.png")

## ----load---------------------------------------------------------------------
library(genMPGMM)

## ----basic--------------------------------------------------------------------
set.seed(42)

dat <- genMPGMM(
  n_feature_patterns        = c(2, 2),
  n_components              = c(2, 2),
  feature_group_proportions = list(c(0.5, 0.5), c(0.6, 0.4)),
  mixing_proportions        = list(c(0.5, 0.5), c(0.4, 0.6)),
  dist_mahalanobis          = c(3, 4),
  target_ari                = c(1, 0.3),
  M                         = 60,
  N                         = 80,
  covariance_spec           = list(type = "diagonal"),
  noise_feature_fraction    = 0,
  seed                      = 42
)

class(dat)

## ----summary------------------------------------------------------------------
summary(dat)

## ----Xdim---------------------------------------------------------------------
dim(dat$X)          # M x N

## ----s------------------------------------------------------------------------
# s[[p]] is an integer vector of length M
# giving the feature-group label for each feature under profile p
length(dat$s)          # P
table(dat$s[[1]])      # group sizes in profile 1
table(dat$s[[2]])      # group sizes in profile 2

## ----z------------------------------------------------------------------------
# z[[p]] is an integer vector of length N
# giving the component label for each observation under profile p
length(dat$z)
table(dat$z[[1]])
table(dat$z[[2]])

## ----ari----------------------------------------------------------------------
# P x P matrix; entry [p, q] = ARI between s[[p]] and s[[q]]
round(dat$achieved_ari_features, 3)

## ----mahal--------------------------------------------------------------------
# achieved_mahalanobis[[p]] is a K_p x K_p matrix of pairwise distances
dat$achieved_mahalanobis[[1]]   # profile 1, K=2 components
dat$achieved_mahalanobis[[2]]   # profile 2, K=2 components

## ----ari_example, eval = FALSE------------------------------------------------
# # Profile 2 feature partition almost identical to profile 1
# target_ari = c(1, 0.9)
# 
# # Profile 2 feature partition completely independent of profile 1
# target_ari = c(1, 0.0)

## ----mahal_example------------------------------------------------------------
set.seed(1)
well_separated <- genMPGMM(
  n_feature_patterns        = c(2, 2),
  n_components              = c(3, 2),
  feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
  mixing_proportions        = list(c(1/3, 1/3, 1/3), c(0.5, 0.5)),
  dist_mahalanobis          = c(5, 3),
  target_ari                = c(1, 0.5),
  M = 40, N = 60,
  covariance_spec           = list(type = "diagonal"),
  noise_feature_fraction    = 0,
  seed = 1
)

# Profile 1 has K=3 components: three off-diagonal pairwise distances; median ≈ 5
D <- well_separated$achieved_mahalanobis[[1]]
round(D[upper.tri(D)], 2)

## ----cov_diag, eval = FALSE---------------------------------------------------
# covariance_spec = list(type = "diagonal")                   # unit variances
# covariance_spec = list(type = "diagonal", diag_values = 2)  # variance = 2 everywhere
# covariance_spec = list(type = "diagonal",
#                        diag_values = runif(M, 0.5, 1.5))    # feature-specific variances

## ----cov_full, eval = FALSE---------------------------------------------------
# covariance_spec = list(type = "full_shared")                 # random positive-definite matrix
# covariance_spec = list(type = "full_shared", cov_mtx = Sigma) # supply your own M x M matrix

## ----noise--------------------------------------------------------------------
noisy <- genMPGMM(
  n_feature_patterns        = c(2, 2),
  n_components              = c(2, 2),
  feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
  mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
  dist_mahalanobis          = 3,
  target_ari                = c(1, 0.5),
  M = 60, N = 80,
  covariance_spec           = list(type = "diagonal"),
  noise_feature_fraction    = 0.2,   # 20 % of features are pure noise
  seed = 7
)

cat("Noise features:", length(noisy$noise_feature_indices), "out of", nrow(noisy$X), "\n")
cat("Noise feature indices:", head(noisy$noise_feature_indices), "...\n")

## ----addnoise, eval = FALSE---------------------------------------------------
# genMPGMM(..., add_noise = TRUE, noise_sd = 0.3)

## ----multilabel---------------------------------------------------------------
set.seed(10)
mp <- genMPGMM(
  n_feature_patterns        = c(3, 2),
  n_components              = c(3, 2),
  feature_group_proportions = list(c(0.4, 0.3, 0.3), c(0.6, 0.4)),
  mixing_proportions        = list(c(0.4, 0.3, 0.3), c(0.5, 0.5)),
  dist_mahalanobis          = c(4, 3),
  target_ari                = c(1, 0.2),
  M = 60, N = 90,
  covariance_spec           = list(type = "diagonal"),
  noise_feature_fraction    = 0,
  seed = 10
)

# Profile 1 sees 3 groups; profile 2 sees 2
cat("Profile 1 labels (first 10 obs):", head(mp$z[[1]], 10), "\n")
cat("Profile 2 labels (first 10 obs):", head(mp$z[[2]], 10), "\n")

## ----report, eval = FALSE-----------------------------------------------------
# # Landscape A4 report -- all plots
# report_MPObject(dat, file = "my_report.pdf")
# 
# # Selected plots only
# report_MPObject(dat, file = "my_report_partial.pdf", which = c(1, 4, 6))

## ----repro--------------------------------------------------------------------
args <- list(
  n_feature_patterns        = c(2, 2),
  n_components              = c(2, 2),
  feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
  mixing_proportions        = list(c(0.5, 0.5), c(0.5, 0.5)),
  dist_mahalanobis          = c(3, 3),
  target_ari                = c(1, 0.5),
  M = 30, N = 40,
  covariance_spec           = list(type = "diagonal"),
  noise_feature_fraction    = 0,
  seed = 99
)

a <- do.call(genMPGMM, args)
b <- do.call(genMPGMM, args)
identical(a$X, b$X)

