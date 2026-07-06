# ============================================================
# Visualisation for MPObject
# ============================================================
#
# Public API
# ----------
#   plot.MPObject(x, which, ask, umap_*, alluvial_colour_profile, ...)
#   report_MPObject(x, file, which, open, umap_*, alluvial_colour_profile, ...)
#
# Internal helpers are prefixed with .mp_
# ============================================================

# ---- dependency guard -------------------------------------

.mp_check_deps <- function() {
  missing_pkgs <- Filter(
    function(p) !requireNamespace(p, quietly = TRUE),
    c("ggplot2", "gridExtra", "scales")
  )
  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages must be installed to use plotting functions:\n",
      "  ", paste(missing_pkgs, collapse = ", "), "\n",
      "Install them with:\n",
      '  install.packages(c(', paste0('"', missing_pkgs, '"', collapse = ", "), '))',
      call. = FALSE
    )
  }
}

# ---- shared theme & colours --------------------------------

.mp_theme <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle     = ggplot2::element_text(colour = "grey40", size = base_size - 1,
                                                lineheight = 1.3),
      legend.key.size   = ggplot2::unit(0.4, "cm"),
      strip.text        = ggplot2::element_text(face = "bold"),
      strip.background  = ggplot2::element_rect(fill = "white")
    )
}

.mp_group_colours <- function(n) {
  all_cols <- c("#4F7FC7", "#D97A32", "#B84A82", "#7A6BC7", "#D9A441")
  all_cols[seq_len(min(n, length(all_cols)))]
}

.mp_comp_colours <- function(n) {
  all_cols <- c("#4F7FC7", "#D97A32", "#B84A82", "#7A6BC7", "#D9A441")
  all_cols[seq_len(min(n, length(all_cols)))]
}

# ---- profile-aware subspace helpers ------------------------
#
# Used by silhouette and within/between-distance plots.
# Observations are projected onto the mean-contrast subspace of each
# profile -- the directions in feature space along which that profile's
# component means vary most.

# Signal feature row indices for profile p (non-noise, valid group assignment).
.mp_profile_signal_rows <- function(dat, p) {
  M          <- dat$settings$M
  noise_rows <- dat$noise_feature_indices
  if (is.null(noise_rows)) noise_rows <- integer(0)
  rows <- setdiff(seq_len(M), noise_rows)

  sp <- dat$s[[p]]
  if (!is.null(sp) && length(sp) == M) {
    rows <- intersect(rows, which(is.finite(as.numeric(sp)) & as.numeric(sp) > 0))
  }
  rows
}

# Returns a list: X_scaled (N x |rows_keep|), rows_keep, center, scale.
.mp_scaled_profile_matrix <- function(dat, p) {
  rows <- .mp_profile_signal_rows(dat, p)
  if (length(rows) < 2L)
    stop("Profile ", p, ": fewer than 2 usable signal features.", call. = FALSE)

  X_obs <- t(dat$X[rows, , drop = FALSE])
  vars  <- apply(X_obs, 2, stats::var, na.rm = TRUE)
  keep  <- is.finite(vars) & vars > 0

  if (sum(keep) < 2L)
    stop("Profile ", p, ": fewer than 2 non-constant signal features.", call. = FALSE)

  X_obs    <- X_obs[, keep, drop = FALSE]
  X_scaled <- scale(X_obs)
  list(
    X_scaled  = X_scaled,
    rows_keep = rows[keep],
    center    = attr(X_scaled, "scaled:center"),
    scale     = attr(X_scaled, "scaled:scale")
  )
}

# Returns N x r score matrix in the mean-contrast subspace of profile p.
.mp_profile_contrast_scores <- function(dat, p) {
  obj <- .mp_scaled_profile_matrix(dat, p)
  if (is.null(dat$mu) || length(dat$mu) < p)
    stop("dat$mu[[", p, "]] is not available.", call. = FALSE)

  mu_p <- dat$mu[[p]][, obj$rows_keep, drop = FALSE]   # K_p x |rows_keep|
  if (nrow(mu_p) < 2L)
    stop("Profile ", p, ": fewer than 2 component means.", call. = FALSE)

  mu_scaled   <- sweep(sweep(mu_p, 2, obj$center, "-"), 2, obj$scale, "/")
  mu_centered <- scale(mu_scaled, center = TRUE, scale = FALSE)
  sv          <- svd(mu_centered)
  tol         <- max(dim(mu_centered)) * max(sv$d) * .Machine$double.eps
  r           <- max(1L, sum(sv$d > tol))
  B           <- sv$v[, seq_len(r), drop = FALSE]
  scores      <- obj$X_scaled %*% B
  colnames(scores) <- paste0("contrast_", seq_len(r))
  scores
}

# Pairwise observation distances for profile p.
.mp_profile_distances <- function(dat, p,
                                   space = c("mean_contrast", "full_scaled")) {
  space <- match.arg(space)
  if (space == "mean_contrast")
    return(stats::dist(.mp_profile_contrast_scores(dat, p)))
  stats::dist(.mp_scaled_profile_matrix(dat, p)$X_scaled)
}

# ---- individual plot builders ------------------------------

# 1. Feature partitions: one coloured bar per profile
.mp_plot_feature_parts <- function(dat) {
  P      <- dat$settings$P
  M      <- dat$settings$M
  L_vals <- dat$settings$n_feature_patterns

  fp_df <- do.call(rbind, lapply(seq_len(P), function(p) {
    data.frame(
      Feature = seq_len(M),
      Profile = sprintf("Profile %d  (L=%d)", p, L_vals[p]),
      Group   = factor(dat$s[[p]])
    )
  }))

  ari_mat  <- dat$achieved_ari_features
  all_levels <- as.character(seq_len(max(L_vals)))
  pal        <- stats::setNames(.mp_group_colours(length(all_levels)), all_levels)

  ggplot2::ggplot(fp_df,
                  ggplot2::aes(.data$Feature, 1, fill = .data$Group)) +
    ggplot2::geom_tile(height = 0.8) +
    ggplot2::scale_fill_manual(values = pal, name = "Feature group",
                               drop = FALSE) +
    ggplot2::facet_wrap(~Profile, ncol = 1) +
    ggplot2::labs(
      title    = "Feature partitions per profile",
      subtitle = paste0(
        "Each bar = M features coloured by group assignment.\n",
        "ARI close to 1 -> groups align across profiles.  ",
        "ARI close to 0 -> groups independent.\n",
        "Target ARI: ", paste(
          vapply(seq_len(P), function(p) {
            if (p == 1) "1.000 (ref)"
            else sprintf("%.3f", .mp_scalar_target_ari(dat, p))
          }, character(1)),
          collapse = "  |  "
        )
      ),
      x = "Feature index", y = NULL
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    .mp_theme() +
    ggplot2::theme(axis.text.y  = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid   = ggplot2::element_blank())
}

# Helper: extract scalar ARI target for profile p (vs reference)
.mp_scalar_target_ari <- function(dat, p) {
  ta <- dat$settings$target_ari
  if (is.matrix(ta)) return(ta[1, p])
  full <- if (length(ta) == 1) c(1, rep(ta, dat$settings$P - 1))
          else if (length(ta) == dat$settings$P - 1) c(1, ta)
          else ta
  full[p]
}

# 2. Observation (component) partitions
.mp_plot_obs_parts <- function(dat) {
  P      <- dat$settings$P
  N      <- dat$settings$N
  K_vals <- dat$settings$n_components

  obs_df <- do.call(rbind, lapply(seq_len(P), function(p) {
    data.frame(
      Obs     = seq_len(N),
      Profile = sprintf("Profile %d  (K=%d)", p, K_vals[p]),
      Comp    = factor(dat$z[[p]])
    )
  }))

  K_max  <- max(K_vals)
  all_lv <- as.character(seq_len(K_max))
  pal    <- stats::setNames(.mp_comp_colours(K_max), all_lv)

  prop_strs <- vapply(seq_len(P), function(p) {
    vals <- round(dat$settings$mixing_proportions[[p]], 2)
    sprintf("P%d: (%s)", p, paste(vals, collapse = ", "))
  }, character(1))

  ggplot2::ggplot(obs_df,
                  ggplot2::aes(.data$Obs, 1, fill = .data$Comp)) +
    ggplot2::geom_tile(height = 0.8) +
    ggplot2::scale_fill_manual(values = pal, name = "Component", drop = FALSE) +
    ggplot2::facet_wrap(~Profile, ncol = 1) +
    ggplot2::labs(
      title    = "Observation (component) partitions",
      subtitle = paste0(
        "Each bar = N observations coloured by component assignment.\n",
        "Bar widths reflect mixing proportions: ",
        paste(prop_strs, collapse = "  |  ")
      ),
      x = "Observation index", y = NULL
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    .mp_theme() +
    ggplot2::theme(axis.text.y  = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid   = ggplot2::element_blank())
}

# 3. PCA: list of P ggplot objects (one per profile)
.mp_plot_pca <- function(dat) {
  P           <- dat$settings$P
  noise_rows  <- dat$noise_feature_indices
  signal_rows <- setdiff(seq_len(dat$settings$M), noise_rows)

  if (length(signal_rows) < 2) {
    return(list(ggplot2::ggplot() +
                  ggplot2::annotate("text", x = 0.5, y = 0.5,
                                    label = "Not enough signal features for PCA.") +
                  ggplot2::theme_void()))
  }

  pca_fit <- stats::prcomp(t(dat$X[signal_rows, ]), scale. = TRUE)
  var_pct <- round(100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2), 1)
  x_lab   <- sprintf("PC1  (%.1f %%)", var_pct[1])
  y_lab   <- sprintf("PC2  (%.1f %%)", var_pct[2])

  K_max <- max(dat$settings$n_components)
  pal   <- stats::setNames(.mp_comp_colours(K_max), as.character(seq_len(K_max)))

  lapply(seq_len(P), function(p) {
    pca_df <- data.frame(
      PC1  = pca_fit$x[, 1],
      PC2  = pca_fit$x[, 2],
      Comp = factor(dat$z[[p]])
    )
    ggplot2::ggplot(pca_df,
                    ggplot2::aes(.data$PC1, .data$PC2, colour = .data$Comp)) +
      ggplot2::geom_point(size = 1.8, alpha = 0.7) +
      ggplot2::scale_colour_manual(values = pal, name = "Component", drop = FALSE) +
      ggplot2::labs(
        title    = sprintf("Profile %d  (K=%d)", p, dat$settings$n_components[p]),
        subtitle = sprintf(
          "Mahalanobis target: %g  |  achieved median: %.3f",
          dat$settings$dist_mahalanobis[p],
          {
            D  <- dat$achieved_mahalanobis[[p]]
            dv <- D[upper.tri(D)]
            if (length(dv) == 0) 0 else median(dv)
          }
        ),
        x = x_lab, y = y_lab
      ) +
      .mp_theme()
  })
}

# 4. UMAP: list of P ggplot objects (one per profile)
.mp_plot_umap <- function(dat,
                          n_neighbors    = 15,
                          min_dist       = 0.1,
                          metric         = "euclidean",
                          n_epochs       = NULL,
                          nn_method = "annoy",
                          seed           = 42,
                          scale_features = TRUE) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    msg <- paste0("Package 'uwot' is required for UMAP.\n",
                  "Install it with: install.packages('uwot')")
    return(lapply(seq_len(dat$settings$P), function(p)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg) +
        ggplot2::theme_void()))
  }

  P          <- dat$settings$P
  noise_rows <- dat$noise_feature_indices
  if (is.null(noise_rows)) noise_rows <- integer(0)
  signal_rows <- setdiff(seq_len(dat$settings$M), noise_rows)

  fallback <- function(msg)
    lapply(seq_len(P), function(p)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg) +
        ggplot2::theme_void())

  if (length(signal_rows) < 2)
    return(fallback("Not enough signal features for UMAP."))
  if (dat$settings$N < 3)
    return(fallback("Not enough observations for UMAP."))

  X_umap <- t(dat$X[signal_rows, , drop = FALSE])
  vars   <- apply(X_umap, 2, stats::var, na.rm = TRUE)
  keep   <- is.finite(vars) & vars > 0

  if (sum(keep) < 2)
    return(fallback("Not enough non-constant signal features for UMAP."))

  X_umap <- X_umap[, keep, drop = FALSE]
  if (scale_features) X_umap <- scale(X_umap)

  n_neighbors <- min(n_neighbors, nrow(X_umap) - 1L)
  set.seed(seed)

  umap_fit <- uwot::umap(
    X_umap,
    n_neighbors  = n_neighbors,
    min_dist     = min_dist,
    metric       = metric,
    n_components = 2L,
    n_epochs     = n_epochs,
    verbose      = FALSE,
    ret_model    = FALSE,
    init = "pca"
  )

  K_max <- max(dat$settings$n_components)
  pal   <- stats::setNames(.mp_comp_colours(K_max), as.character(seq_len(K_max)))

  lapply(seq_len(P), function(p) {
    umap_df <- data.frame(
      UMAP1 = umap_fit[, 1],
      UMAP2 = umap_fit[, 2],
      Comp  = factor(dat$z[[p]])
    )
    ggplot2::ggplot(umap_df,
                    ggplot2::aes(.data$UMAP1, .data$UMAP2, colour = .data$Comp)) +
      ggplot2::geom_point(size = 1.8, alpha = 0.7) +
      ggplot2::scale_colour_manual(values = pal, name = "Component", drop = FALSE) +
      ggplot2::labs(
        title    = sprintf("Profile %d  (K=%d)", p, dat$settings$n_components[p]),
        subtitle = sprintf(
          "n_neighbors = %d  |  min_dist = %.2f  |  metric = %s",
          n_neighbors, min_dist, metric
        ),
        x = "UMAP 1", y = "UMAP 2"
      ) +
      .mp_theme()
  })
}

# 5. Mahalanobis tile matrix
.mp_plot_mahal <- function(dat) {
  P <- dat$settings$P

  mah_df <- do.call(rbind, lapply(seq_len(P), function(p) {
    D   <- dat$achieved_mahalanobis[[p]]
    tgt <- dat$settings$dist_mahalanobis[p]
    K   <- nrow(D)
    g   <- expand.grid(k1 = seq_len(K), k2 = seq_len(K))
    data.frame(
      Profile = sprintf("Profile %d  (target = %g)", p, tgt),
      k1      = factor(g$k1),
      k2      = factor(g$k2),
      Dist    = as.vector(D)
    )
  }))

  max_d <- max(mah_df$Dist)

  ggplot2::ggplot(mah_df,
                  ggplot2::aes(.data$k1, .data$k2, fill = .data$Dist)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$Dist)),
      size = 3.8, colour = "white", fontface = "bold"
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#4F7FC7", "#A9C7E8", "#F7F7F2", "#E8B58B", "#D97A32"),
      limits  = c(0, max_d),
      name    = "Distance"
    ) +
    ggplot2::facet_wrap(~Profile) +
    ggplot2::coord_equal() +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      title    = "Pairwise Mahalanobis distances between component means",
      subtitle = paste0(
        "Diagonal = 0 (same component).  Off-diagonal = separation between components.\n",
        "For K=2: target is met exactly.  For K>2: target is the MEDIAN of all off-diagonal values."
      ),
      x = "Component", y = "Component"
    ) +
    .mp_theme()
}

# 6. Silhouette widths: list of P ggplot objects
#    Distances are computed in the profile mean-contrast subspace.
.mp_plot_silhouette <- function(dat) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    msg <- paste0("Package 'cluster' is required for silhouette plots.\n",
                  "Install it with: install.packages('cluster')")
    return(lapply(seq_len(dat$settings$P), function(p)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                          hjust = 0.5) +
        ggplot2::theme_void()))
  }

  P      <- dat$settings$P
  K_vals <- dat$settings$n_components
  K_max  <- max(K_vals)
  pal    <- stats::setNames(.mp_comp_colours(K_max), as.character(seq_len(K_max)))

  lapply(seq_len(P), function(p) {
    K   <- K_vals[p]

    if (K < 2L) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5,
                                 label = sprintf("Profile %d: K=1 — silhouette not defined.", p)) +
               ggplot2::theme_void())
    }

    z_p <- as.integer(dat$z[[p]])
    d_p <- tryCatch(
      .mp_profile_distances(dat, p, space = "mean_contrast"),
      error = function(e) {
        message("Profile ", p, ": silhouette skipped — ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(d_p)) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5,
                                 label = sprintf("Profile %d: distances unavailable.", p)) +
               ggplot2::theme_void())
    }

    sil_mat <- as.matrix(cluster::silhouette(z_p, d_p))
    sil_df  <- data.frame(
      cluster   = as.integer(sil_mat[, 1]),
      sil_width = sil_mat[, 3]
    )
    sil_df$Comp  <- factor(sil_df$cluster)
    sil_df$Facet <- factor(
      paste0("Cluster ", sil_df$cluster),
      levels = paste0("Cluster ", seq_len(K))
    )
    sil_df <- sil_df[order(sil_df$cluster, sil_df$sil_width), , drop = FALSE]
    sil_df$ObsNum <- ave(sil_df$sil_width, sil_df$cluster, FUN = seq_along)
    sil_df$xmin   <- pmin(0, sil_df$sil_width)
    sil_df$xmax   <- pmax(0, sil_df$sil_width)

    avg_sil  <- mean(sil_df$sil_width, na.rm = TRUE)
    neg_frac <- mean(sil_df$sil_width < 0, na.rm = TRUE)
    x_min    <- max(-1, min(-0.25,
                            floor(min(sil_df$sil_width, na.rm = TRUE) * 4) / 4))

    ggplot2::ggplot(sil_df) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.5) +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$ObsNum - 0.45, ymax = .data$ObsNum + 0.45,
                     fill = .data$Comp, colour = .data$Comp),
        alpha = 0.88, linewidth = 0.12
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$Facet),
                          scales = "free_y", space = "free_y") +
      ggplot2::scale_fill_manual(
        values = pal[seq_len(K)], name = "Component"
      ) +
      ggplot2::scale_colour_manual(
        values = pal[seq_len(K)], name = "Component"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(x_min, 1),
        breaks = seq(x_min, 1, by = 0.25),
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        breaks = NULL,
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::labs(
        title    = sprintf("Silhouette — Profile %d  (K=%d)", p, K),
        subtitle = sprintf(
          "Mean width = %.3f  |  %.1f%% negative  |  profile mean-contrast subspace",
          avg_sil, 100 * neg_frac
        ),
        x = "Silhouette width", y = NULL
      ) +
      .mp_theme() +
      ggplot2::theme(
        panel.spacing.y    = ggplot2::unit(0, "lines"),
        panel.grid.major.y = ggplot2::element_blank(),
        strip.text.y       = ggplot2::element_text(angle = 90),
        legend.position    = "bottom"
      )
  })
}

# 7. Within/between-component distances: list of P ggplot objects
#    Distances are computed in the profile mean-contrast subspace.
.mp_plot_distance_dist <- function(dat) {
  P      <- dat$settings$P
  K_vals <- dat$settings$n_components

  pal_wb <- c(
    "Within component"   = "#D97A32",
    "Between components" = "#4F7FC7"
  )

  lapply(seq_len(P), function(p) {
    K   <- K_vals[p]

    if (K < 2L) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5,
                                 label = sprintf("Profile %d: K=1 — distance plot not meaningful.", p)) +
               ggplot2::theme_void())
    }

    z_p <- as.integer(dat$z[[p]])
    d_p <- tryCatch(
      .mp_profile_distances(dat, p, space = "mean_contrast"),
      error = function(e) {
        message("Profile ", p, ": distance plot skipped — ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(d_p)) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5,
                                 label = sprintf("Profile %d: distances unavailable.", p)) +
               ggplot2::theme_void())
    }

    dist_mat  <- as.matrix(d_p)
    upper_idx <- which(upper.tri(dist_mat), arr.ind = TRUE)

    pair_type <- ifelse(
      z_p[upper_idx[, 1]] == z_p[upper_idx[, 2]],
      "Within component", "Between components"
    )

    dist_df <- data.frame(
      Distance = dist_mat[upper_idx],
      Type     = factor(pair_type,
                        levels = c("Within component", "Between components"))
    )

    type_means <- tapply(dist_df$Distance, dist_df$Type, mean, na.rm = TRUE)
    sep_ratio  <- unname(type_means["Between components"] /
                           type_means["Within component"])

    ggplot2::ggplot(dist_df,
                    ggplot2::aes(.data$Distance,
                                 fill = .data$Type, colour = .data$Type)) +
      ggplot2::geom_density(alpha = 0.35, linewidth = 0.8, adjust = 1.1) +
      ggplot2::geom_vline(xintercept = type_means["Within component"],
                          colour = pal_wb["Within component"],
                          linetype = "dashed", linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = type_means["Between components"],
                          colour = pal_wb["Between components"],
                          linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = pal_wb, name = NULL) +
      ggplot2::scale_colour_manual(values = pal_wb, name = NULL) +
      ggplot2::labs(
        title    = sprintf(
          "Within / between-component distances — Profile %d  (K=%d)", p, K),
        subtitle = sprintf(
          "Within mean = %.2f  |  Between mean = %.2f  |  Separation ratio = %.2f\nProfile mean-contrast subspace",
          type_means["Within component"],
          type_means["Between components"],
          sep_ratio
        ),
        x = "Pairwise distance", y = "Density"
      ) +
      .mp_theme() +
      ggplot2::theme(legend.position = "top")
  })
}

# 8. Alluvial / Sankey: component assignment flows across profiles.
#    Returns a single ggplot coloured by colour_by_profile's assignments.
.mp_plot_alluvial <- function(dat, colour_by_profile = 1L) {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    return(ggplot2::ggplot() +
             ggplot2::annotate(
               "text", x = 0.5, y = 0.5,
               label = paste0(
                 "Package 'ggalluvial' is required for the alluvial plot.\n",
                 "Install it with: install.packages('ggalluvial')"
               ),
               hjust = 0.5
             ) +
             ggplot2::theme_void())
  }

  P      <- dat$settings$P
  K_vals <- dat$settings$n_components
  cbp    <- as.integer(colour_by_profile)

  if (cbp < 1L || cbp > P)
    stop("`colour_by_profile` must be an integer between 1 and P (", P, ").",
         call. = FALSE)

  obs_wide <- do.call(data.frame, lapply(seq_len(P), function(p) {
    factor(as.integer(dat$z[[p]]), levels = seq_len(K_vals[p]))
  }))
  colnames(obs_wide) <- paste0("Profile_", seq_len(P))

  agg_wide <- as.data.frame(table(obs_wide), stringsAsFactors = FALSE)
  agg_wide <- agg_wide[agg_wide$Freq > 0, , drop = FALSE]
  agg_wide$alluvium_id <- seq_len(nrow(agg_wide))

  lodes <- ggalluvial::to_lodes_form(
    agg_wide,
    axes  = seq_len(P),
    id    = "alluvium_id",
    key   = "Profile",
    value = "Component"
  )

  colour_col    <- paste0("Profile_", cbp)
  colour_lookup <- stats::setNames(
    as.character(agg_wide[[colour_col]]),
    as.character(agg_wide$alluvium_id)
  )
  K_c <- K_vals[cbp]
  pal <- stats::setNames(.mp_comp_colours(K_c), as.character(seq_len(K_c)))

  lodes$ColourComp <- factor(
    colour_lookup[as.character(lodes$alluvium_id)],
    levels = as.character(seq_len(K_c))
  )

  x_labels <- stats::setNames(
    paste0("Profile ", seq_len(P)),
    paste0("Profile_", seq_len(P))
  )

  ggplot2::ggplot(
    lodes,
    ggplot2::aes(x        = .data$Profile,
                 y        = .data$Freq,
                 stratum  = .data$Component,
                 alluvium = .data$alluvium_id)
  ) +
    ggalluvial::geom_flow(
      ggplot2::aes(fill = .data$ColourComp, colour = .data$ColourComp),
      alpha     = 0.70,
      width     = 0.35,
      knot.pos  = 0.4
    ) +
    ggalluvial::geom_stratum(
      width     = 0.35,
      fill      = "white",
      colour    = "grey45",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      stat      = ggalluvial::StatStratum,
      ggplot2::aes(label = ggplot2::after_stat(stratum)),
      size      = 3.5,
      fontface  = "bold",
      colour    = "grey20"
    ) +
    ggplot2::scale_x_discrete(labels = x_labels) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_manual(
      values = pal,
      name   = sprintf("Profile %d\ncomponent", cbp)
    ) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::labs(
      title    = sprintf(
        "Component assignments across profiles  (colour = Profile %d)", cbp),
      subtitle = paste0(
        "Each stream = observations sharing the same component-label pattern across profiles.\n",
        sprintf(
          "Width proportional to count.  Change with alluvial_colour_profile = 1..%d.", P)
      ),
      x = NULL, y = "Count"
    ) +
    .mp_theme() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "right"
    )
}

# ---- build catalogue (lazy: only builds what is selected) --

.mp_build_all <- function(dat,
                          selected                = .mp_resolve_which("all"),
                          umap_args               = list(),
                          alluvial_colour_profile = 1L) {
  plots <- list()

  if ("feature_partitions"     %in% selected)
    plots$feature_partitions     <- .mp_plot_feature_parts(dat)
  if ("observation_partitions" %in% selected)
    plots$observation_partitions <- .mp_plot_obs_parts(dat)
  if ("pca"                    %in% selected)
    plots$pca                    <- .mp_plot_pca(dat)
  if ("umap"                   %in% selected)
    plots$umap                   <- do.call(.mp_plot_umap,
                                            c(list(dat), umap_args))
  if ("mahalanobis"            %in% selected)
    plots$mahalanobis            <- .mp_plot_mahal(dat)
  if ("silhouette"             %in% selected)
    plots$silhouette             <- .mp_plot_silhouette(dat)
  if ("distance_distributions" %in% selected)
    plots$distance_distributions <- .mp_plot_distance_dist(dat)
  if ("alluvial"               %in% selected)
    plots$alluvial               <- .mp_plot_alluvial(dat, alluvial_colour_profile)

  plots
}

# Resolve `which` argument to a character vector of plot names
.mp_resolve_which <- function(which_arg) {
  plot_names <- c(
    "feature_partitions", "observation_partitions",
    "pca", "umap", "mahalanobis",
    "silhouette", "distance_distributions", "alluvial"
  )
  if (identical(which_arg, "all")) return(plot_names)
  if (is.numeric(which_arg)) {
    idx <- as.integer(which_arg)
    bad <- idx[idx < 1L | idx > length(plot_names)]
    if (length(bad))
      stop("Plot indices out of range [1, ", length(plot_names), "]: ",
           paste(bad, collapse = ", "), call. = FALSE)
    return(plot_names[idx])
  }
  if (is.character(which_arg)) {
    bad <- setdiff(which_arg, plot_names)
    if (length(bad))
      stop("Unknown plot name(s): ", paste(bad, collapse = ", "),
           "\nValid names: ", paste(plot_names, collapse = ", "),
           call. = FALSE)
    return(which_arg)
  }
  stop("`which` must be \"all\", an integer vector 1-", length(plot_names),
       ", or a character vector of plot names.", call. = FALSE)
}

# ---- interactive rendering (used by plot.MPObject) ---------

.mp_list_title <- function(name) {
  switch(name,
    pca                    = "PCA of X  (noise features excluded)",
    umap                   = "UMAP of X  (noise features excluded)",
    silhouette             = "Silhouette widths  (profile mean-contrast subspace)",
    distance_distributions = "Within / between-component distances",
    name
  )
}

.mp_print_entry <- function(entry, name) {
  if (inherits(entry, "gg")) {
    print(entry)
    return(invisible(NULL))
  }

  if (is.list(entry) && all(vapply(entry, inherits, logical(1), "gg"))) {
    n   <- length(entry)
    ttl <- .mp_list_title(name)

    gridExtra::grid.arrange(
      grid::textGrob(ttl, gp = grid::gpar(fontface = "bold", fontsize = 13)),
      gridExtra::arrangeGrob(grobs = entry, ncol = min(n, 2L)),
      ncol    = 1,
      heights = grid::unit(c(0.8, 1), c("cm", "null"))
    )
    return(invisible(NULL))
  }

  warning("Could not render plot '", name, "'.")
}

# ---- grob conversion (single ggplot or list -> one grob) ---

.mp_entry_to_grob <- function(entry, name) {
  if (inherits(entry, "gg"))
    return(ggplot2::ggplotGrob(entry))

  if (is.list(entry) && all(vapply(entry, inherits, logical(1), "gg"))) {
    n   <- length(entry)
    ttl <- .mp_list_title(name)
    return(gridExtra::arrangeGrob(
      grid::textGrob(ttl, x = 0.02, hjust = 0,
                     gp = grid::gpar(fontface = "bold", fontsize = 12)),
      gridExtra::arrangeGrob(grobs = entry, ncol = min(n, 2L)),
      ncol    = 1,
      heights = grid::unit(c(0.8, 1), c("cm", "null"))
    ))
  }
  grid::nullGrob()
}

# ---- margined page draw ------------------------------------

.mp_draw_margined <- function(grob, page_w, page_h, margin = 1) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x      = grid::unit(margin, "inches"),
    y      = grid::unit(margin, "inches"),
    width  = grid::unit(page_w - 2 * margin, "inches"),
    height = grid::unit(page_h - 2 * margin, "inches"),
    just   = c("left", "bottom")
  ))
  grid::grid.draw(grob)
  grid::popViewport()
}

# ---- pagination helper (adaptive layout) -------------------

.mp_paginate_plots <- function(plot_list, per_page, layout_ncol,
                               title, subtitle = NULL) {
  n     <- length(plot_list)
  n_pgs <- ceiling(n / per_page)

  lapply(seq_len(n_pgs), function(i) {
    idx   <- seq((i - 1L) * per_page + 1L, min(i * per_page, n))
    grobs <- plot_list[idx]

    while (length(grobs) < per_page)
      grobs <- c(grobs, list(grid::nullGrob()))

    pg_title <- if (n_pgs > 1L)
      sprintf("%s  (%d / %d)", title, i, n_pgs)
    else
      title

    plots_grob <- gridExtra::arrangeGrob(grobs = grobs, ncol = layout_ncol)

    if (!is.null(subtitle)) {
      header <- gridExtra::arrangeGrob(
        grid::textGrob(pg_title, x = 0.02, hjust = 0,
                       gp = grid::gpar(fontface = "bold", fontsize = 12)),
        grid::textGrob(subtitle, x = 0.02, hjust = 0,
                       gp = grid::gpar(fontsize = 9, lineheight = 1.1,
                                       col = "grey40")),
        ncol    = 1,
        heights = grid::unit(c(0.7, 1.5), "cm")
      )
      gridExtra::arrangeGrob(header, plots_grob,
                             ncol    = 1,
                             heights = grid::unit(c(2.5, 1), c("cm", "null")))
    } else {
      gridExtra::arrangeGrob(
        grid::textGrob(pg_title, x = 0.02, hjust = 0,
                       gp = grid::gpar(fontface = "bold", fontsize = 12)),
        plots_grob,
        ncol    = 1,
        heights = grid::unit(c(0.8, 1), c("cm", "null"))
      )
    }
  })
}

# ---- summary table grob ------------------------------------

.mp_summary_table_grob <- function(dat) {
  s <- dat$settings
  P <- s$P

  rows <- lapply(seq_len(P), function(p) {
    D   <- dat$achieved_mahalanobis[[p]]
    dv  <- D[upper.tri(D)]
    mah <- if (length(dv) == 0) "n/a (K=1)"
            else sprintf("%.3f / %.1f", median(dv), s$dist_mahalanobis[p])
    ari <- if (p == 1) "reference"
            else sprintf("%.3f / %.3f",
                         dat$achieved_ari_features[1, p],
                         .mp_scalar_target_ari(dat, p))
    mix_str  <- paste(round(s$mixing_proportions[[p]], 2),        collapse = ", ")
    feat_str <- paste(round(s$feature_group_proportions[[p]], 2), collapse = ", ")
    data.frame(
      "Profile"                        = sprintf("P%d", p),
      "K"                              = s$n_components[p],
      "Mixing proportions"             = mix_str,
      "L"                              = s$n_feature_patterns[p],
      "Feature group proportions"      = feat_str,
      "Mahalanobis\n(achieved/target)" = mah,
      "ARI\n(achieved/target)"         = ari,
      "Cov"                            = s$covariance_spec$type,
      check.names      = FALSE,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)

  title_grob <- grid::textGrob(
    "genMPGMM -- Generated dataset summary",
    x = 0, hjust = 0,
    gp = grid::gpar(fontface = "bold", fontsize = 14)
  )

  noise_str <- if (s$add_noise)
    sprintf("Additive noise: TRUE (sd = %g)", s$noise_sd)
  else
    "Additive noise: FALSE"

  meta <- sprintf(
    "M = %d features  x  N = %d observations  |  P = %d profiles  |  Noise features: %d (%.0f %%)  |  %s",
    s$M, s$N, P,
    length(dat$noise_feature_indices),
    100 * s$noise_feature_fraction,
    noise_str
  )
  meta_grob <- grid::textGrob(
    meta, x = 0, hjust = 0,
    gp = grid::gpar(fontsize = 10, col = "grey35")
  )

  row_fills <- rep_len(c("white", "#f0f0f0"), P)
  tbl_grob <- gridExtra::tableGrob(
    df, rows = NULL,
    theme = gridExtra::ttheme_default(
      base_size = 10,
      core    = list(bg_params = list(fill = row_fills, col = "grey80"),
                     fg_params = list(hjust = 0, x = 0.05)),
      colhead = list(bg_params = list(fill = "white", col = "grey80"),
                     fg_params = list(col = "black", fontface = "bold",
                                      hjust = 0, x = 0.05))
    )
  )
  tbl_grob$widths <- grid::unit(
    c(0.07, 0.04, 0.20, 0.04, 0.20, 0.17, 0.17, 0.11), "npc"
  )

  gridExtra::arrangeGrob(
    title_grob,
    meta_grob,
    tbl_grob,
    ncol    = 1,
    heights = grid::unit(c(0.9, 0.6, 1, 1), c("cm", "cm", "null", "null"))
  )
}

# ---- landscape renderer ------------------------------------

.mp_render_landscape <- function(plots, dat, selected, page_w, page_h) {
  .mp_draw_margined(.mp_summary_table_grob(dat), page_w, page_h)

  list_plots_2pp <- c("pca", "umap", "silhouette", "distance_distributions")

  for (nm in selected) {
    if (nm %in% list_plots_2pp) {
      for (pg in .mp_paginate_plots(plots[[nm]], per_page = 2L, layout_ncol = 2L,
                                    title = .mp_list_title(nm)))
        .mp_draw_margined(pg, page_w, page_h)
    } else {
      .mp_draw_margined(.mp_entry_to_grob(plots[[nm]], nm), page_w, page_h)
    }
  }

  invisible(NULL)
}

# ============================================================
# Public: plot.MPObject
# ============================================================

#' Plot a multi-profile GMM object
#'
#' Displays diagnostic plots for an object of class \code{"MPObject"}.
#' Use \code{which} to select a subset.
#'
#' @section Plots available:
#' \describe{
#'   \item{1 / "feature_partitions"}{One coloured bar per profile showing which
#'     feature belongs to which group.}
#'   \item{2 / "observation_partitions"}{One coloured bar per profile showing
#'     the component assignment of each observation.}
#'   \item{3 / "pca"}{Global PCA scatter coloured by component assignment
#'     (one panel per profile).}
#'   \item{4 / "umap"}{Global UMAP scatter coloured by component assignment
#'     (one panel per profile).  Requires the \pkg{uwot} package.}
#'   \item{5 / "mahalanobis"}{Tile matrix of pairwise Mahalanobis distances
#'     between component means for each profile.}
#'   \item{6 / "silhouette"}{Silhouette width plot per profile, computed in the
#'     profile-specific mean-contrast subspace.  Requires \pkg{cluster}.}
#'   \item{7 / "distance_distributions"}{Density overlay of within- vs
#'     between-component pairwise distances, per profile, in the profile
#'     mean-contrast subspace.}
#'   \item{8 / "alluvial"}{Sankey/alluvial diagram of component assignment flows
#'     across profiles.  Requires \pkg{ggalluvial}.  Colour profile controlled
#'     by \code{alluvial_colour_profile}.}
#' }
#'
#' @param x An object of class \code{"MPObject"}.
#' @param which Plots to display.  \code{"all"} (default), an integer vector
#'   \code{1:8}, or a character vector of plot names (see section above).
#' @param ask If \code{TRUE} (default when interactive), pause between plots.
#' @param umap_n_neighbors Number of neighbours for UMAP (default 15).
#' @param umap_min_dist Minimum distance for UMAP (default 0.1).
#' @param umap_metric Distance metric for UMAP (default \code{"euclidean"}).
#' @param umap_n_epochs Number of UMAP epochs (default 300).
#' @param umap_seed Random seed for UMAP (default 42).
#' @param umap_scale_features Scale features before UMAP (default \code{TRUE}).
#' @param alluvial_colour_profile Integer 1..P: which profile's component
#'   assignments colour the alluvial flows (default 1).
#' @param ... Currently unused.
#'
#' @return Invisibly returns a named list of the built \code{ggplot} objects.
#'
#' @examples
#' \dontrun{
#' dat <- genMPGMM(
#'   n_feature_patterns = c(2, 2), n_components = c(2, 2),
#'   feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
#'   mixing_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
#'   dist_mahalanobis = 3, target_ari = c(1, 0.3),
#'   M = 60, N = 80, covariance_spec = list(type = "diagonal"),
#'   noise_feature_fraction = 0, seed = 1
#' )
#' plot(dat)                                    # all plots, interactive
#' plot(dat, which = c("pca", "silhouette"))    # subset by name
#' plot(dat, which = 4, umap_n_neighbors = 20) # UMAP with custom params
#' plot(dat, which = "alluvial", alluvial_colour_profile = 2)
#' }
#'
#' @export
#' @rawNamespace export(plot.MPObject)
plot.MPObject <- function(x,
                          which                   = "all",
                          ask                     = interactive(),
                          umap_n_neighbors        = 15,
                          umap_min_dist           = 0.1,
                          umap_metric             = "euclidean",
                          umap_n_epochs           = 300,
                          umap_seed               = 42,
                          umap_scale_features     = TRUE,
                          alluvial_colour_profile = 1L,
                          ...) {
  .mp_check_deps()
  selected  <- .mp_resolve_which(which)
  umap_args <- list(
    n_neighbors    = umap_n_neighbors,
    min_dist       = umap_min_dist,
    metric         = umap_metric,
    n_epochs       = umap_n_epochs,
    seed           = umap_seed,
    scale_features = umap_scale_features
  )
  plots <- .mp_build_all(x, selected, umap_args, alluvial_colour_profile)

  if (ask && length(selected) > 1L) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }

  for (nm in selected) {
    .mp_print_entry(plots[[nm]], nm)
  }

  invisible(plots)
}

# ============================================================
# Public: report_MPObject
# ============================================================

#' Save a diagnostic PDF report for a multi-profile GMM object
#'
#' Writes an A4 landscape PDF whose layout adapts automatically to any number
#' of profiles (P) and components (K).  A summary table is followed by one
#' page per plot section; multi-panel plots (PCA, UMAP, silhouette, distance
#' distributions) are paginated two per page.
#'
#' @param x An object of class \code{"MPObject"}.
#' @param file Path for the output PDF.  Defaults to \code{"mpgmm_report.pdf"}.
#' @param which Plots to include.  \code{"all"} (default), an integer vector
#'   \code{1:8}, or a character vector of plot names.
#'   See \code{\link{plot.MPObject}} for the full list.
#' @param open If \code{TRUE} (default when interactive), open the PDF after
#'   saving.
#' @param umap_n_neighbors,umap_min_dist,umap_metric,umap_n_epochs,umap_seed,umap_scale_features
#'   UMAP parameters; see \code{\link{plot.MPObject}}.
#' @param alluvial_colour_profile Integer 1..P; see \code{\link{plot.MPObject}}.
#' @param ... Currently unused.
#'
#' @return Invisibly returns the path to the saved PDF.
#'
#' @examples
#' \dontrun{
#' dat <- genMPGMM(
#'   n_feature_patterns = c(2, 2), n_components = c(2, 2),
#'   feature_group_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
#'   mixing_proportions = list(c(0.5, 0.5), c(0.5, 0.5)),
#'   dist_mahalanobis = 3, target_ari = c(1, 0.3),
#'   M = 60, N = 80, covariance_spec = list(type = "diagonal"),
#'   noise_feature_fraction = 0, seed = 1
#' )
#' report_MPObject(dat)
#' report_MPObject(dat, which = c("pca", "mahalanobis"), file = "partial.pdf")
#' }
#'
#' @export
report_MPObject <- function(x,
                            file                    = "mpgmm_report.pdf",
                            which                   = "all",
                            open                    = interactive(),
                            umap_n_neighbors        = 15,
                            umap_min_dist           = 0.1,
                            umap_metric             = "euclidean",
                            umap_n_epochs           = NULL,
                            umap_seed               = 42,
                            umap_scale_features     = TRUE,
                            alluvial_colour_profile = 1L,
                            ...) {
  .mp_check_deps()
  selected  <- .mp_resolve_which(which)
  umap_args <- list(
    n_neighbors    = umap_n_neighbors,
    min_dist       = umap_min_dist,
    metric         = umap_metric,
    n_epochs       = umap_n_epochs,
    seed           = umap_seed,
    scale_features = umap_scale_features
  )
  plots <- .mp_build_all(x, selected, umap_args, alluvial_colour_profile)
  file  <- normalizePath(file, mustWork = FALSE)

  grDevices::pdf(file, width = 11.69, height = 8.27, onefile = TRUE)
  tryCatch(
    .mp_render_landscape(plots, x, selected, 11.69, 8.27),
    finally = grDevices::dev.off()
  )

  message("Report saved to:\n  ", file)
  if (open)
    tryCatch(
      utils::browseURL(paste0("file:///", gsub("\\\\", "/", file))),
      error = function(e) NULL
    )

  invisible(file)
}
