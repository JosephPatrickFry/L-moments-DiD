# ================================================================
# Simultaneous inference for L-moments: Gaussian multiplier bootstrap
# (note, Section 3.3) applied to 2x2 DiD contrasts of independent cells.
#
# Single sample (as in the note): with jackknife influences l_{i,k},
#   Z_k^(b) = n^{-1/2} sum_i xi_i l_{i,k},  xi_i ~ N(0,1) iid,
#   T^(b)   = max_{k<=K} |Z_k^(b)| / sigma_hat_k,
#   c_hat   = empirical (1 - alpha) quantile of T^(1..B),
#   I_k     = lambda_hat_k +/- c_hat * sigma_hat_k / sqrt(n).
#
# Contrasts across independent cells (e.g. a 2x2 DiD with weights
# w = (+1, -1, -1, +1) on T_post, T_pre, C_post, C_pre): each cell gets
# its own multipliers and the bootstrap error of contrast r, order k, is
#   Z_{r,k}^(b) = sum_c w_{r,c} (1/n_c) sum_{i in c} xi_i l_{c,i,k},
# whose conditional covariance is exactly the jackknife covariance of
# the contrast (lmoment_contrast_vcov). Then
#   T^(b) = max_{r,k} |Z_{r,k}^(b)| / se_{r,k},
#   CI_{r,k} = est_{r,k} +/- c_hat * se_{r,k}  (simultaneous over r, k).
# With one contrast and one cell this is exactly the note's procedure
# (se = sigma_hat / sqrt(n)).
#
# Because Z^(b) is linear in Gaussian multipliers it is exactly
# N(0, V_hat) given the data; method = "gaussian" draws from that law
# directly (cost free of n), method = "multiplier" draws the xi_i.
# Requires lmoment_jackknife.R to be sourced first.
# ================================================================

# Critical value for simultaneous inference over all contrasts x orders.
#   cells:   named list of lmom_u_jackknife() results (independent samples)
#   weights: contrasts x cells matrix (column names = cell names)
# Returns c_hat, the bootstrap max-t draws, and the contrast estimates/SEs.
multiplier_critical_value <- function(cells, weights, alpha = 0.05, B = 10000,
                                      method = c("multiplier", "gaussian"),
                                      seed = NULL, chunk = 1000) {
  method <- match.arg(method)
  weights <- as.matrix(weights)
  if (!is.null(seed)) set.seed(seed)
  cv <- lmoment_contrast_vcov(cells, weights)
  se <- as.vector(t(cv$se))                    # contrast-major, matching cv$vcov
  stopifnot(all(se > 0))
  K <- length(cells[[1]]$lambda)
  R <- nrow(weights)

  Tb <- numeric(B)
  if (method == "gaussian") {
    L <- chol(cv$vcov + diag(1e-12 * max(diag(cv$vcov)), nrow(cv$vcov)))
    done <- 0
    while (done < B) {
      b <- min(chunk, B - done)
      Z <- matrix(rnorm(b * nrow(L)), b) %*% L      # b x (R*K), rows ~ N(0, V)
      Tb[done + seq_len(b)] <- apply(abs(sweep(Z, 2, se, `/`)), 1, max)
      done <- done + b
    }
  } else {
    done <- 0
    while (done < B) {
      b <- min(chunk, B - done)
      Z <- matrix(0, b, R * K)
      for (cn in colnames(weights)) {
        cl <- cells[[cn]]
        xi <- matrix(rnorm(cl$n * b), cl$n, b)
        S  <- crossprod(xi, cl$influence) / cl$n          # b x K: (1/n_c) sum xi_i l_i
        for (r in seq_len(R)) {
          if (weights[r, cn] == 0) next
          cols <- (r - 1) * K + seq_len(K)
          Z[, cols] <- Z[, cols] + weights[r, cn] * S
        }
      }
      Tb[done + seq_len(b)] <- apply(abs(sweep(Z, 2, se, `/`)), 1, max)
      done <- done + b
    }
  }
  # empirical conditional (1 - alpha) quantile: inf{c : F_B(c) >= 1 - alpha}
  c_hat <- unname(quantile(Tb, 1 - alpha, type = 1))
  list(crit = c_hat, T_draws = Tb, estimate = cv$estimate, se = cv$se, vcov = cv$vcov,
       alpha = alpha, B = B, method = method)
}

# 2x2 DiD of the first K U-statistic L-moments with jackknife SEs,
# pointwise and simultaneous (multiplier-bootstrap) confidence intervals.
#   y: outcome; treat, post: 0/1 vectors of the same length
# Cells (treat x post) are treated as independent iid samples.
did2x2_lmoment_ci <- function(y, treat, post, K = 4, alpha = 0.05, B = 10000,
                              method = c("multiplier", "gaussian"), seed = NULL) {
  method <- match.arg(method)
  ok <- !is.na(y) & !is.na(treat) & !is.na(post)
  y <- y[ok]; treat <- treat[ok]; post <- post[ok]
  cell <- function(tr, po) lmom_u_jackknife(y[treat == tr & post == po], K)
  cells <- list(T_post = cell(1, 1), T_pre = cell(1, 0), C_post = cell(0, 1), C_pre = cell(0, 0))
  w <- matrix(c(1, -1, -1, 1), 1, dimnames = list("DiD", names(cells)))
  mb <- multiplier_critical_value(cells, w, alpha = alpha, B = B, method = method, seed = seed)
  est <- mb$estimate[1, ]; se <- mb$se[1, ]; z <- qnorm(1 - alpha / 2)
  out <- data.frame(k = seq_len(K), estimate = est, se = se,
                    lo_pointwise = est - z * se, hi_pointwise = est + z * se,
                    lo_simultaneous = est - mb$crit * se, hi_simultaneous = est + mb$crit * se,
                    row.names = names(est))
  attr(out, "crit") <- mb$crit
  attr(out, "z") <- z
  attr(out, "cells") <- cells
  attr(out, "n") <- sapply(cells, `[[`, "n")
  out
}
