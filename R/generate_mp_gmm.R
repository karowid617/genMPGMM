#' Generate multi-profile Gaussian mixture data
#'
#' Generate a synthetic data matrix with multiple latent profiles.
#'
#' @param n_feature_patterns Number of feature groups in each profile.
#' @param n_components Number of observation components in each profile.
#' @param feature_group_proportions List of feature-group proportions.
#' @param mixing_proportions List of component proportions.
#' @param dist_mahalanobis Target pairwise Mahalanobis separation.
#' @param target_ari ARI control for feature partitions.
#' @param M Number of features.
#' @param N Number of observations.
#' @param covariance_spec A list specifying the covariance structure within each
#'   profile. Supported types are \code{"diagonal"} and \code{"full_shared"}.
#'   For \code{"diagonal"}, use \code{diag_values} to provide either a single
#'   variance value, a vector of feature-specific variances of length M, or
#'   a list of such vectors for profile-specific covariance matrices. For
#'   \code{"full_shared"}, use \code{cov_mtx} to provide a full M by M
#'   covariance matrix (or a list of such matrices, one per profile). If
#'   \code{cov_mtx} is not supplied, a random positive-definite covariance matrix
#'   is generated using \code{scale}.
#' @param add_noise Should additive Gaussian noise be added? By default \code{add_noise = FALSE}.
#' @param noise_sd Standard deviation of additive noise. By default \code{noise_sd = 0.1}.
#' @param seed Random seed (optional).
#' @param noise_feature_fraction Fraction of features replaced by noise features. By default \code{noise_feature_fraction = 0.1}.
#' @param ari_mode ARI control mode. This should be \code{"vs_reference"} (similarity specified relative to the first profile) or \code{"pairwise_matrix"} (similarity specified by user for all profile pairs).
#' @param ari_tol Tolerance for ARI optimization. By default \code{ari_tol = 0.01}.
#' @param ari_max_iter Maximum number of ARI optimization iterations. By default \code{ari_max_iter = 10000}.
#' @param ari_swap_frac Fraction of labels swapped in ARI perturbation. By default \code{ari_swap_frac = 0.1}
#' @param profile_weights Profile weights (optional). Vector of length P,
#'        controlling how strongly each profile contributes to the final data matrix.
#'        If \code{NULL}, all profiles are weighted equally.
#'
#' @details
#' The argument \code{covariance_spec} determines how within-component
#' variability is represented in each profile.
#'
#' With \code{type = "diagonal"}, features are assumed to be uncorrelated, so
#' the covariance matrix contains non-zero values only on the diagonal. This is
#' a convenient choice when different features should have different variances,
#' but no feature-to-feature dependence is required.
#'
#' With \code{type = "full_shared"}, a full covariance matrix is used within
#' each profile, allowing arbitrary covariance between features. The covariance
#' matrix is shared across all mixture components within a profile, while
#' component separation is controlled through the component means.
#'
#' @examples
#' data <- genMPGMM(
#' n_feature_patterns = c(2, 2),
#' n_components = c(2, 2),
#' feature_group_proportions = list(c(0.5, 0.5),c(0.7, 0.3)),
#' mixing_proportions = list(c(0.2, 0.8),c(0.6, 0.4)),
#' dist_mahalanobis = c(3, 4),
#' target_ari = c(1, 0.2),
#' M = 80,
#' N = 60,
#' covariance_spec = list(type = "diagonal",diag_values = list(rep(1, 80),rep(1.5, 80))),
#' seed = 123)
#'
#' @return An object of class \code{"MPObject"}.
#'
#' @export

genMPGMM <- function(
  # ------------------------------
  # REQUIRED
  n_feature_patterns,                          # number of feature groups in each profile
  n_components,                          # number of observation components in each profile
  feature_group_proportions,      # list of length P; each sums to 1; lengths = L_p
  mixing_proportions,             # list of length P; each sums to 1; lengths = K_p
  dist_mahalanobis,               # scalar or length P
  target_ari,            # ARI control for feature partitions s[[p]]
  M,                              # number of features (rows)
  N,                              # number of observations (cols)
  covariance_spec,                # list(type=...)
  add_noise = FALSE,              # add global Gaussian noise?

  # ------------------------------
  # STRONGLY RECOMMENDED
  noise_sd = 0.1,                 # additive iid noise sd
  seed = NULL,
  noise_feature_fraction = 0.1,     # overwrite a fraction of rows with pure noise
  ari_mode = "vs_reference",      # 'vs_reference' or 'pairwise_matrix'

  # ------------------------------
  # TECHNICAL
  ari_tol = 0.01,
  ari_max_iter = 10000,
  ari_swap_frac = 0.1,
  profile_weights = NULL         # optional weights when summing profile signals
) {
  # ------------------------------
  # PROFILE MEAN STRUCTURE
  template_sd = 1
  feature_sd_within_group = 0.05
  baseline_sd = 0.5
  P = length(n_feature_patterns)

  validate_generator_inputs(
    P = P,
    n_feature_patterns = n_feature_patterns,
    n_components = n_components,
    M = M,
    N = N,
    feature_group_proportions = feature_group_proportions,
    mixing_proportions = mixing_proportions,
    dist_mahalanobis = dist_mahalanobis,
    covariance_spec = covariance_spec,
    ari_mode = ari_mode,
    target_ari = target_ari,
    add_noise = add_noise,
    noise_sd = noise_sd,
    noise_feature_fraction = noise_feature_fraction
  )

  if (!is.null(seed)) set.seed(seed)

  feature_group_proportions <- lapply(feature_group_proportions, normalize_probs)
  mixing_proportions <- lapply(mixing_proportions, normalize_probs)
  dist_vec <- coerce_profile_vector(dist_mahalanobis, P, "dist_mahalanobis")

  if (is.null(profile_weights)) {
    profile_weights <- rep(1, P)
  } else {
    profile_weights <- coerce_profile_vector(profile_weights, P, "profile_weights")
  }

  # 1) Feature partitions s[[p]] with ARI control
  s_list <- generate_feature_partitions(
    M = M,
    P = P,
    n_feature_patterns = n_feature_patterns,
    feature_group_proportions = feature_group_proportions,
    target_ari = target_ari,
    ari_mode = ari_mode,
    ari_tol = ari_tol,
    ari_max_iter = ari_max_iter,
    ari_swap_frac = ari_swap_frac
  )

  # 2) Observation partitions z[[p]]
  z_list <- generate_observation_partitions(
    N = N,
    P = P,
    n_components = n_components,
    mixing_proportions = mixing_proportions
  )

  # 3) Profile-specific covariance + centroid structure
  cov_mtx_list <- vector("list", P)
  templates_list <- vector("list", P)
  baselines_list <- vector("list", P)
  delta_list <- vector("list", P)
  mu_list <- vector("list", P)
  D_list <- vector("list", P)

  for (p in seq_len(P)) {
    cov_mtx_p <- build_profile_covariance(
      M = M,
      p = p,
      P = P,
      covariance_spec = covariance_spec
    )

    means_obj <- generate_profile_mean_structure(
      s_p = s_list[[p]],
      L_p = n_feature_patterns[p],
      K_p = n_components[p],
      cov_mtx_p = cov_mtx_p,
      target_mahalanobis = dist_vec[p],
      template_sd = template_sd,
      feature_sd_within_group = feature_sd_within_group,
      baseline_sd = baseline_sd
    )

    cov_mtx_list[[p]] <- cov_mtx_p
    templates_list[[p]] <- means_obj$templates
    baselines_list[[p]] <- means_obj$feature_baselines
    delta_list[[p]] <- means_obj$delta
    mu_list[[p]] <- means_obj$mu
    D_list[[p]] <- means_obj$D
  }

  # 4) Generate each profile signal
  X_profiles <- vector("list", P)
  X_signal <- matrix(0, nrow = M, ncol = N)

  for (p in seq_len(P)) {
    X_p <- generate_profile_signal(
      M = M,
      N = N,
      z_p = z_list[[p]],
      mu_p = mu_list[[p]],
      cov_mtx_p = cov_mtx_list[[p]]
    )

    X_profiles[[p]] <- X_p
    X_signal <- X_signal + profile_weights[p] * X_p
  }

  # 5) Optional additive noise
  X <- X_signal
  if (add_noise && noise_sd > 0) {
    X <- X + matrix(rnorm(M * N, mean = 0, sd = noise_sd), nrow = M, ncol = N)
  }

  # 6) Optional pure-noise rows
  noise_feature_indices <- integer(0)
  if (noise_feature_fraction > 0) {
    n_noise <- floor(M * noise_feature_fraction)
    if (n_noise > 0) {
      noise_feature_indices <- sort(sample(seq_len(M), size = n_noise, replace = FALSE))
      X[noise_feature_indices, ] <- matrix(
        rnorm(n_noise * N, mean = 0, sd = 1),
        nrow = n_noise,
        ncol = N
      )
    }
  }

  # 7) Reports
  achieved_ari_features <- ari_matrix_from_list(s_list)

  achieved_mixing <- lapply(z_list, function(zp) {
    tab <- table(factor(zp, levels = seq_len(max(zp))))
    as.numeric(tab) / sum(tab)
  })

  achieved_feature_group_props <- lapply(s_list, function(sp) {
    tab <- table(factor(sp, levels = seq_len(max(sp))))
    as.numeric(tab) / sum(tab)
  })

  achieved_mahalanobis <- lapply(D_list, function(D) round(D, 4))

  out <- list(
    X = X,
    X_signal = X_signal,
    X_profiles = X_profiles,

    s = s_list,                        # feature partitions
    z = z_list,                        # observation partitions

    templates = templates_list,        # L_p x K_p group templates
    feature_baselines = baselines_list,# length M per profile
    delta = delta_list,                # K_p x M component-separating part
    mu = mu_list,                      # K_p x M final component centroids
    cov_mtx = cov_mtx_list,                # M x M covariance per profile

    achieved_ari_features = achieved_ari_features,
    achieved_mixing = achieved_mixing,
    achieved_feature_group_props = achieved_feature_group_props,
    achieved_mahalanobis = achieved_mahalanobis,

    noise_feature_indices = noise_feature_indices,

    settings = list(
      P = P,
      n_feature_patterns = n_feature_patterns,
      n_components = n_components,
      feature_group_proportions = feature_group_proportions,
      mixing_proportions = mixing_proportions,
      dist_mahalanobis = dist_mahalanobis,
      target_ari = target_ari,
      M = M,
      N = N,
      covariance_spec = covariance_spec,
      add_noise = add_noise,
      noise_sd = noise_sd,
      seed = seed,
      noise_feature_fraction = noise_feature_fraction,
      ari_mode = ari_mode,
      ari_tol = ari_tol,
      ari_max_iter = ari_max_iter,
      ari_swap_frac = ari_swap_frac,
      profile_weights = profile_weights,
      template_sd = template_sd,
      feature_sd_within_group = feature_sd_within_group,
      baseline_sd = baseline_sd
    )
  )

  class(out) <- "MPObject"
  out
}
