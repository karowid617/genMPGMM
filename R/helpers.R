normalize_probs <- function(x, name = "probabilities") {
  if (is.null(x)) stop(sprintf("`%s` cannot be NULL.", name))
  if (any(x < 0)) stop(sprintf("`%s` must be non-negative.", name))
  s <- sum(x)
  if (s <= 0) stop(sprintf("`%s` must sum to a positive value.", name))
  x / s
}

counts_from_proportions <- function(n, probs) {
  probs <- normalize_probs(probs)
  raw <- n * probs
  counts <- floor(raw)
  remainder <- n - sum(counts)

  if (remainder > 0) {
    ord <- order(raw - counts, decreasing = TRUE)
    counts[ord[seq_len(remainder)]] <- counts[ord[seq_len(remainder)]] + 1
  }

  counts
}

make_exact_labels <- function(n, probs, labels = NULL) {
  probs <- normalize_probs(probs)
  K <- length(probs)

  if (is.null(labels)) labels <- seq_len(K)
  if (length(labels) != K) stop("`labels` must have same length as `probs`.")

  counts <- counts_from_proportions(n, probs)
  out <- rep(labels, times = counts)
  sample(out, length(out), replace = FALSE)
}

ari_pair <- function(a, b) {
  mclust::adjustedRandIndex(a, b)
}

ari_matrix_from_list <- function(labels_list) {
  P <- length(labels_list)
  A <- matrix(1, nrow = P, ncol = P)
  for (p in seq_len(P)) {
    for (q in seq_len(P)) {
      if (p < q) {
        A[p, q] <- ari_pair(labels_list[[p]], labels_list[[q]])
        A[q, p] <- A[p, q]
      }
    }
  }
  A
}

coerce_profile_vector <- function(x, P, name) {
  if (length(x) == 1) {
    rep(x, P)
  } else if (length(x) == P) {
    x
  } else {
    stop(sprintf("`%s` must have length 1 or P.", name))
  }
}

get_profile_element <- function(x, p, P, name = "value") {
  if (is.null(x)) return(NULL)

  if (is.list(x)) {
    if (length(x) == 1) return(x[[1]])
    if (length(x) == P) return(x[[p]])
    stop(sprintf("List `%s` must have length 1 or P.", name))
  }

  if (length(x) == 1) return(x)
  if (length(x) == P) return(x[p])

  x
}

random_spd <- function(d, scale = 1, jitter = 1e-6) {
  A <- matrix(rnorm(d * d), nrow = d)
  S <- crossprod(A) / d
  S <- scale * S + diag(jitter, d)
  S
}

# ============================================================
# 2. Input validation
# ============================================================

validate_generator_inputs <- function(
    P,
    n_feature_patterns,
    n_components,
    M,
    N,
    feature_group_proportions,
    mixing_proportions,
    dist_mahalanobis,
    covariance_spec,
    ari_mode,
    target_ari,
    add_noise,
    noise_sd,
    noise_feature_fraction
) {
  if (P < 2) stop("`P` must be >= 2.")
  if (length(n_feature_patterns) != P) stop("`n_feature_patterns` must have length P.")
  if (length(n_components) != P) stop("`n_components` must have length P.")
  if (any(n_feature_patterns < 1)) stop("All values in `n_feature_patterns` must be >= 1.")
  if (any(n_components < 1)) stop("All values in `n_components` must be >= 1.")
  if (M < 2) stop("`M` must be >= 2.")
  if (N < 2) stop("`N` must be >= 2.")

  if (!is.list(feature_group_proportions) || length(feature_group_proportions) != P) {
    stop("`feature_group_proportions` must be a list of length P.")
  }
  if (!is.list(mixing_proportions) || length(mixing_proportions) != P) {
    stop("`mixing_proportions` must be a list of length P.")
  }

  for (p in seq_len(P)) {
    if (length(feature_group_proportions[[p]]) != n_feature_patterns[p]) {
      stop(sprintf("`feature_group_proportions[[%d]]` must have length n_feature_patterns[%d] = %d.",
                   p, p, n_feature_patterns[p]))
    }
    if (length(mixing_proportions[[p]]) != n_components[p]) {
      stop(sprintf("`mixing_proportions[[%d]]` must have length n_components[%d] = %d.",
                   p, p, n_components[p]))
    }
  }

  if (!(length(dist_mahalanobis) %in% c(1, P))) {
    stop("`dist_mahalanobis` must have length 1 or P.")
  }

  if (!is.list(covariance_spec) || is.null(covariance_spec$type)) {
    stop("`covariance_spec` must be a list with at least `$type`.")
  }

  if (!ari_mode %in% c("vs_reference", "pairwise_matrix")) {
    stop("`ari_mode` must be 'vs_reference' or 'pairwise_matrix'.")
  }

  if (ari_mode == "vs_reference") {
    if (!(length(target_ari) %in% c(1, P - 1, P))) {
      stop("For `ari_mode = 'vs_reference'`, `target_ari` must have length 1, P-1, or P.")
    }
  }

  if (ari_mode == "pairwise_matrix") {
    if (!is.matrix(target_ari) || any(dim(target_ari) != c(P, P))) {
      stop("For `ari_mode = 'pairwise_matrix'`, `target_ari` must be a P x P matrix.")
    }
    if (!isSymmetric(target_ari)) {
      stop("For `ari_mode = 'pairwise_matrix'`, `target_ari` must be a symmetric matrix (ARI is symmetric).")
    }
    if (any(abs(diag(target_ari) - 1) > 1e-9)) {
      stop("For `ari_mode = 'pairwise_matrix'`, all diagonal entries of `target_ari` must equal 1 (ARI of a partition with itself is always 1).")
    }
    if (any(target_ari < -1 - 1e-9) || any(target_ari > 1 + 1e-9)) {
      stop("All entries of `target_ari` must be in [-1, 1].")
    }
  }

  if (!is.logical(add_noise) || length(add_noise) != 1) {
    stop("`add_noise` must be TRUE/FALSE.")
  }

  if (!is.numeric(noise_sd) || length(noise_sd) != 1 || noise_sd < 0) {
    stop("`noise_sd` must be a single non-negative number.")
  }

  if (!is.numeric(noise_feature_fraction) ||
      length(noise_feature_fraction) != 1 ||
      noise_feature_fraction < 0 || noise_feature_fraction >= 1) {
    stop("`noise_feature_fraction` must be in [0, 1).")
  }
}

# ============================================================
# 3. Feature partitions s[[p]] with ARI control
# ============================================================

perturb_labels_to_target_ari <- function(
    ref_labels,
    initial_labels,
    target_ari,
    max_iter = 10000,
    tol = 0.01,
    swap_frac = 0.1
) {
  ref_labels <- as.vector(ref_labels)
  best_labels <- as.vector(initial_labels)

  best_diff <- abs(ari_pair(ref_labels, best_labels) - target_ari)
  if (best_diff <= tol) return(best_labels)

  n <- length(best_labels)
  n_swap <- max(2, floor(swap_frac * n))
  if (n_swap %% 2 == 1) n_swap <- n_swap + 1
  n_swap <- min(n_swap, n)

  for (iter in seq_len(max_iter)) {
    cand <- best_labels
    idx <- sample(seq_len(n), size = n_swap, replace = FALSE)
    cand[idx] <- sample(cand[idx], length(idx), replace = FALSE)

    cand_diff <- abs(ari_pair(ref_labels, cand) - target_ari)

    if (cand_diff < best_diff) {
      best_labels <- cand
      best_diff <- cand_diff
    }
    if (best_diff <= tol) break
  }

  if (best_diff > tol) {
    warning(sprintf(
      paste0("ARI optimization did not converge: |achieved ARI - target ARI| = %.4f > tol = %.4f ",
             "after %d iterations. Consider increasing `ari_max_iter` or relaxing `ari_tol`. ",
             "Note: target ARI = %.4f may not be achievable if feature-group proportions differ ",
             "between profiles."),
      best_diff, tol, max_iter, target_ari
    ))
  }

  best_labels
}

make_target_ari_vs_reference <- function(target_ari, P) {
  if (length(target_ari) == 1) {
    c(1, rep(target_ari, P - 1))
  } else if (length(target_ari) == P - 1) {
    c(1, target_ari)
  } else if (length(target_ari) == P) {
    if (abs(target_ari[1] - 1) > 1e-9) {
      warning(sprintf(
        paste0("In `ari_mode = 'vs_reference'`, `target_ari[1]` is ignored: ",
               "ARI of profile 1 with itself is always 1. Got %.4f."),
        target_ari[1]
      ))
    }
    target_ari
  } else {
    stop("Invalid `target_ari` length for `vs_reference` mode.")
  }
}

generate_feature_partitions_vs_reference <- function(
    M,
    P,
    n_feature_patterns,
    feature_group_proportions,
    target_ari,
    ari_tol = 0.01,
    ari_max_iter = 10000,
    ari_swap_frac = 0.1
) {
  target_full <- make_target_ari_vs_reference(target_ari, P)
  s_list <- vector("list", P)

  s_list[[1]] <- make_exact_labels(
    n = M,
    probs = feature_group_proportions[[1]],
    labels = seq_len(n_feature_patterns[1])
  )

  for (p in 2:P) {
    init_labels <- make_exact_labels(
      n = M,
      probs = feature_group_proportions[[p]],
      labels = seq_len(n_feature_patterns[p])
    )

    s_list[[p]] <- perturb_labels_to_target_ari(
      ref_labels = s_list[[1]],
      initial_labels = init_labels,
      target_ari = target_full[p],
      max_iter = ari_max_iter,
      tol = ari_tol,
      swap_frac = ari_swap_frac
    )
  }

  s_list
}

pairwise_ari_objective_for_profile <- function(labels_list, p, target_ari_matrix) {
  P <- length(labels_list)
  obj <- 0
  for (q in seq_len(P)) {
    if (q != p) {
      obj <- obj + (ari_pair(labels_list[[p]], labels_list[[q]]) - target_ari_matrix[p, q])^2
    }
  }
  obj
}

generate_feature_partitions_pairwise <- function(
    M,
    P,
    n_feature_patterns,
    feature_group_proportions,
    target_ari_matrix,
    ari_tol = 0.01,
    ari_max_iter = 20000
) {
  s_list <- vector("list", P)

  for (p in seq_len(P)) {
    s_list[[p]] <- make_exact_labels(
      n = M,
      probs = feature_group_proportions[[p]],
      labels = seq_len(n_feature_patterns[p])
    )
  }

  for (iter in seq_len(ari_max_iter)) {
    improved <- FALSE
    p <- sample(seq_len(P), 1)
    current <- s_list[[p]]

    idx_by_label <- split(seq_len(M), current)
    present_labels <- names(idx_by_label)[lengths(idx_by_label) > 0]

    if (length(present_labels) >= 2) {
      labs <- sample(present_labels, 2, replace = FALSE)
      i <- sample(idx_by_label[[labs[1]]], 1)
      j <- sample(idx_by_label[[labs[2]]], 1)

      cand <- current
      tmp <- cand[i]
      cand[i] <- cand[j]
      cand[j] <- tmp

      old_obj <- pairwise_ari_objective_for_profile(s_list, p, target_ari_matrix)
      s_cand <- s_list
      s_cand[[p]] <- cand
      new_obj <- pairwise_ari_objective_for_profile(s_cand, p, target_ari_matrix)

      if (new_obj < old_obj) {
        s_list[[p]] <- cand
        improved <- TRUE
      }
    }

    if (!improved && iter %% 1000 == 0) {
      current_ari <- ari_matrix_from_list(s_list)
      off_diag    <- lower.tri(target_ari_matrix)
      max_abs_diff <- max(abs(current_ari[off_diag] - target_ari_matrix[off_diag]))
      if (max_abs_diff <= ari_tol) break
    }
  }

  current_ari  <- ari_matrix_from_list(s_list)
  off_diag     <- lower.tri(target_ari_matrix)
  final_max_diff <- max(abs(current_ari[off_diag] - target_ari_matrix[off_diag]))
  if (final_max_diff > ari_tol) {
    warning(sprintf(
      paste0("Pairwise ARI optimization did not converge: max |achieved - target| = %.4f > tol = %.4f ",
             "after %d iterations."),
      final_max_diff, ari_tol, ari_max_iter
    ))
  }

  s_list
}

generate_feature_partitions <- function(
    M,
    P,
    n_feature_patterns,
    feature_group_proportions,
    target_ari,
    ari_mode = "vs_reference",
    ari_tol = 0.01,
    ari_max_iter = 10000,
    ari_swap_frac = 0.1
) {
  if (ari_mode == "vs_reference") {
    generate_feature_partitions_vs_reference(
      M = M,
      P = P,
      n_feature_patterns = n_feature_patterns,
      feature_group_proportions = feature_group_proportions,
      target_ari = target_ari,
      ari_tol = ari_tol,
      ari_max_iter = ari_max_iter,
      ari_swap_frac = ari_swap_frac
    )
  } else {
    generate_feature_partitions_pairwise(
      M = M,
      P = P,
      n_feature_patterns = n_feature_patterns,
      feature_group_proportions = feature_group_proportions,
      target_ari_matrix = target_ari,
      ari_tol = ari_tol,
      ari_max_iter = ari_max_iter
    )
  }
}

# ============================================================
# 4. Observation partitions z[[p]]
# ============================================================

generate_observation_partitions <- function(
    N,
    P,
    n_components,
    mixing_proportions
) {
  z_list <- vector("list", P)

  for (p in seq_len(P)) {
    z_list[[p]] <- make_exact_labels(
      n = N,
      probs = mixing_proportions[[p]],
      labels = seq_len(n_components[p])
    )
  }

  z_list
}

# ============================================================
# 5. Covariance builder
# ============================================================

build_profile_covariance <- function(
    M,
    p,
    P,
    covariance_spec
) {
  type <- covariance_spec$type

  if (type == "diagonal") {
    diag_p <- get_profile_element(covariance_spec$diag_values, p, P, "covariance_spec$diag_values")
    if (is.null(diag_p)) {
      diag_p <- rep(1, M)
    } else if (length(diag_p) == 1) {
      diag_p <- rep(diag_p, M)
    } else if (length(diag_p) != M) {
      stop(sprintf("Diagonal covariance for profile %d must have length 1 or M=%d.", p, M))
    }
    return(diag(as.numeric(diag_p), M))
  }

  if (type == "full_shared") {
    cov_mtx_p <- get_profile_element(covariance_spec$cov_mtx, p, P, "covariance_spec$cov_mtx")
    if (is.null(cov_mtx_p)) {
      scale_p <- get_profile_element(covariance_spec$scale, p, P, "covariance_spec$scale")
      if (is.null(scale_p)) scale_p <- 1
      cov_mtx_p <- random_spd(M, scale = scale_p)
    } else {
      if (!is.matrix(cov_mtx_p) || any(dim(cov_mtx_p) != c(M, M))) {
        stop(sprintf("Supplied cov_mtx for profile %d must be %d x %d.", p, M, M))
      }
    }
    return(cov_mtx_p)
  }

  stop("Unsupported covariance type. Use 'diagonal', or 'full_shared'.")
}

# ============================================================
# 6. Profile mean structure + Mahalanobis control
# ============================================================

mahalanobis_distance_matrix <- function(mu_mat, cov_mtx) {
  # mu_mat: K x M
  K <- nrow(mu_mat)
  if (K == 1) return(matrix(0, 1, 1))

  cov_mtx_inv <- solve(cov_mtx)
  D <- matrix(0, K, K)

  for (k1 in seq_len(K)) {
    for (k2 in seq_len(K)) {
      d <- mu_mat[k1, ] - mu_mat[k2, ]
      D[k1, k2] <- sqrt(drop(t(d) %*% cov_mtx_inv %*% d))
    }
  }
  D
}

make_regular_simplex <- function(K) {
  # Returns a K x (K-1) matrix whose rows are the K vertices of a regular
  # simplex: all pairs of rows have equal Euclidean distance (= sqrt(2)).

  if (K == 1) {
    return(matrix(0, nrow = 1, ncol = 1))
  }

  # Centering matrix I - (1/K) J has K-1 non-zero eigenvalues all equal to 1.
  # Its left singular vectors (scaled by singular values) embed the K simplex
  # vertices into (K-1)-dimensional Euclidean space.
  A <- diag(K) - matrix(1 / K, nrow = K, ncol = K)

  sv     <- svd(A)
  coords <- sv$u[, 1:(K - 1), drop = FALSE] %*% diag(sv$d[1:(K - 1)], nrow = K - 1)

  # Numerical re-centering (theoretically already centered)
  coords <- scale(coords, center = TRUE, scale = FALSE)

  coords
}

generate_profile_mean_structure <- function(
    s_p,
    L_p,
    K_p,
    cov_mtx_p,
    target_mahalanobis,
    template_sd = 1,
    feature_sd_within_group = 0.05,
    baseline_sd = 0.5
) {
  M <- length(s_p)

  # ---------------------------------------
  # 1. K_p vertices of a regular simplex in (K_p-1) dimensions.
  #    All pairwise Euclidean distances between vertices are equal.
  # ---------------------------------------
  simplex <- make_regular_simplex(K_p)   # K_p x (K_p - 1)

  if (K_p == 1) {
    simplex <- matrix(0, nrow = 1, ncol = 1)
  }

  latent_dim <- ncol(simplex)

  # ---------------------------------------
  # 2. Per feature-group: a random loading vector of dimension latent_dim.
  #    Features within the same group share the same base loading, so their
  #    component-mean profiles are similar (but not identical due to per-feature
  #    noise added below).
  # ---------------------------------------
  group_loadings <- matrix(
    rnorm(L_p * latent_dim, mean = 0, sd = template_sd),
    nrow = L_p,
    ncol = latent_dim
  )

  # ---------------------------------------
  # 3. Build delta (K_p x M): component-mean offsets for each feature.
  #    delta[k, m] = simplex[k, ] %*% loading_m, where loading_m is the
  #    group loading for feature m plus a small per-feature perturbation.
  # ---------------------------------------
  delta <- matrix(0, nrow = K_p, ncol = M)

  for (m in seq_len(M)) {
    l <- s_p[m]

    loading_m <- group_loadings[l, ]

    # Small per-feature deviation from the group loading
    if (latent_dim == 1) {
      loading_m <- loading_m + rnorm(1, mean = 0, sd = feature_sd_within_group)
    } else {
      loading_m <- loading_m + rnorm(latent_dim, mean = 0, sd = feature_sd_within_group)
    }

    delta[, m] <- simplex %*% loading_m
  }

  # ---------------------------------------
  # 4. Scale delta so that the MEDIAN pairwise Mahalanobis distance between
  #    component means equals target_mahalanobis.
  #    For K_p = 2 this is exact; for K_p > 2 individual pairs may deviate.
  # ---------------------------------------
  if (K_p > 1) {
    D0 <- mahalanobis_distance_matrix(delta, cov_mtx_p)
    dvals <- D0[upper.tri(D0)]
    dvals <- dvals[is.finite(dvals) & dvals > 0]

    current_median <- if (length(dvals) == 0) 1 else median(dvals)
    scale_factor <- target_mahalanobis / current_median
    delta <- delta * scale_factor
  }

  # ---------------------------------------
  # 5. Per-feature baseline shift (same for all K_p components, so it cancels
  #    in pairwise Mahalanobis distances between means).
  # ---------------------------------------
  feature_baselines <- rnorm(M, mean = 0, sd = baseline_sd)

  mu <- delta
  for (k in seq_len(K_p)) {
    mu[k, ] <- mu[k, ] + feature_baselines
  }

  D <- mahalanobis_distance_matrix(mu, cov_mtx_p)

  list(
    templates = group_loadings,         # L_p x latent_dim
    feature_baselines = feature_baselines,
    delta = delta,
    mu = mu,
    D = D
  )
}

# ============================================================
# 7. Generate profile-specific signal
# ============================================================

generate_profile_signal <- function(
    M,
    N,
    z_p,
    mu_p,
    cov_mtx_p
) {
  K_p <- nrow(mu_p)
  X_p <- matrix(0, nrow = M, ncol = N)

  for (k in seq_len(K_p)) {
    cols_k <- which(z_p == k)
    if (length(cols_k) == 0) next

    if (M == 1) {
      X_p[1, cols_k] <- rnorm(
        n = length(cols_k),
        mean = mu_p[k, 1],
        sd = sqrt(cov_mtx_p[1, 1])
      )
    } else {
      block <- MASS::mvrnorm(
        n = length(cols_k),
        mu = mu_p[k, ],
        Sigma = cov_mtx_p
      )
      X_p[, cols_k] <- t(block)
    }
  }

  X_p
}
