# ================================================================
# Jackknife covariance of the unbiased (U-statistic) sample L-moments
#
# Following the note (Section 3.3): with lambda_hat_k the U-statistic
# estimator and lambda_hat_k^(-i) the same estimator after deleting
# observation i,
#   l_{i,k} = (n - 1) * (lambda_hat_k - lambda_hat_k^(-i)),
#   Sigma_hat[j,k] = (1/n) sum_i l_{i,j} l_{i,k},
# where Sigma_hat estimates the asymptotic covariance of
# sqrt(n) * (lambda_hat - lambda), so Var(lambda_hat) ~ Sigma_hat / n.
# For a U-statistic the deleted estimates average back to the full
# estimate, so sum_i l_{i,k} = 0 exactly and Sigma_hat is PSD.
#
# Fast exact leave-one-out. lambda_hat_k = sum_r p*_{k-1,r} b_r with the
# unbiased PWMs
#   b_r = (1/n) sum_i a_r(i; n) x_(i),  a_r(i; m) = prod_{s=1}^r (i-s)/(m-s).
# Deleting the observation at sorted position j leaves an (n-1)-sample
# in which x_(i) keeps rank i for i < j and has rank i-1 for i > j, so
#   b_r^(-j) = 1/(n-1) * [ sum_{i<j} a_r(i; n-1) x_(i)
#                          + sum_{i>j} a_r(i-1; n-1) x_(i) ],
# a prefix sum plus a suffix sum. All n deletions cost O(n K) after one
# sort. Ties are harmless: deleting any of several tied values leaves
# the same multiset. The data are centered first (lambda_1's influence
# is x_i - mean(x) and lambda_k, k >= 2, are shift-invariant), which
# keeps the prefix sums small and the differences numerically accurate.
#
# Uses the probability-weighted-moment route, which is accurate for the
# small K used in the DiD tables (K <= ~15); see lmom_u_stable() in
# lmoment_choose_K.R for large K. Self-contained, except that the
# brute-force checker lmom_u_jackknife_bruteforce() needs lmom_u() from
# lmoment_did_utils.R (used only by the test scripts).
# ================================================================

# coefficients p*_{k-1,r}: lambda_k = sum_{r=0}^{k-1} P[k, r+1] * b_r
lmom_pwm_coef <- function(K) {
  P <- matrix(0, K, K)
  for (k in 1:K) {
    r <- 0:(k - 1)
    P[k, r + 1] <- (-1)^(k - 1 - r) * choose(k - 1, r) * choose(k - 1 + r, r)
  }
  P
}

# a_r(i; m) for i = idx, vectorized; zero whenever i - s = 0 for some s
pwm_weight <- function(idx, m, r) {
  w <- rep(1, length(idx))
  if (r > 0) for (s in 1:r) w <- w * (idx - s) / (m - s)
  w
}

# Main function.
#   x: numeric vector (NAs dropped); K: number of L-moments (n >= K + 1)
# Returns
#   lambda     K-vector of U-statistic L-moments
#   influence  n x K matrix l_{i,k}, rows in the ORIGINAL order of x
#   Sigma      K x K jackknife estimate of the asymptotic covariance
#   vcov       Sigma / n, the covariance of lambda_hat itself
#   se         sqrt(diag(vcov))
#   n
lmom_u_jackknife <- function(x, K = 4) {
  x <- x[!is.na(x)]
  n <- length(x)
  stopifnot(n >= K + 1)
  o  <- order(x)
  mu <- mean(x)
  xs <- x[o] - mu                       # sorted, centered
  i  <- seq_len(n)
  m  <- n - 1

  b_full <- numeric(K)
  b_del  <- matrix(0, n, K)             # row j: PWMs after deleting sorted position j
  for (r in 0:(K - 1)) {
    b_full[r + 1] <- mean(pwm_weight(i, n, r) * xs)
    lo <- pwm_weight(i, m, r) * xs              # contribution of x_(i) when i < j
    hi <- pwm_weight(i - 1, m, r) * xs          # contribution of x_(i) when i > j
    prefix <- c(0, cumsum(lo))[i]               # sum_{i < j}
    suffix <- c(rev(cumsum(rev(hi))), 0)[i + 1] # sum_{i > j}
    b_del[, r + 1] <- (prefix + suffix) / m
  }
  P <- lmom_pwm_coef(K)
  lam_full <- drop(P %*% b_full)
  lam_del  <- b_del %*% t(P)                    # n x K, centered scale

  infl_sorted <- m * sweep(-lam_del, 2, lam_full, `+`)   # (n-1)(lambda - lambda^(-j))
  infl <- matrix(0, n, K)
  infl[o, ] <- infl_sorted
  lam_full[1] <- lam_full[1] + mu
  nm <- paste0("l", seq_len(K))
  colnames(infl) <- nm
  Sigma <- crossprod(infl) / n
  dimnames(Sigma) <- list(nm, nm)
  vc <- Sigma / n
  list(lambda = setNames(lam_full, nm), influence = infl, Sigma = Sigma,
       vcov = vc, se = sqrt(diag(vc)), n = n)
}

# Brute-force jackknife (n calls to lmom_u); only for checking.
lmom_u_jackknife_bruteforce <- function(x, K = 4) {
  x <- x[!is.na(x)]
  n <- length(x)
  full <- lmom_u(x, K)
  del  <- t(sapply(seq_len(n), \(i) lmom_u(x[-i], K)))
  infl <- (n - 1) * sweep(-del, 2, full, `+`)
  list(lambda = full, influence = infl, Sigma = crossprod(infl) / n)
}

# ------------------------------------------------
# Covariance of linear contrasts of L-moments across independent cells
# ------------------------------------------------
# cells:   named list of lmom_u_jackknife() results (one per group x
#          period cell), treated as mutually independent samples
# weights: matrix with one row per contrast and one column per cell
#          (column names = cell names), e.g. a 2x2 DiD has weights
#          (+1, -1, -1, +1) on (T_post, T_pre, C_post, C_pre)
# Returns, for every L-moment order k, the contrast estimates, their
# covariance matrix (contrasts x contrasts) and standard errors, plus
# the full (contrast x order) joint covariance.
#   Var(sum_c w_c lambda_c) = sum_c w_c^2 Var(lambda_c)  (independence)
lmoment_contrast_vcov <- function(cells, weights) {
  weights <- as.matrix(weights)
  stopifnot(all(colnames(weights) %in% names(cells)))
  K  <- length(cells[[1]]$lambda)
  nm <- names(cells[[1]]$lambda)
  ct <- if (is.null(rownames(weights))) paste0("c", seq_len(nrow(weights))) else rownames(weights)
  est <- matrix(0, nrow(weights), K, dimnames = list(ct, nm))
  # joint covariance over (contrast, order), ordered contrast-major
  J <- matrix(0, nrow(weights) * K, nrow(weights) * K)
  for (cn in colnames(weights)) {
    w <- weights[, cn]
    est <- est + outer(w, cells[[cn]]$lambda)
    J <- J + kronecker(tcrossprod(w), cells[[cn]]$vcov)
  }
  lab <- as.vector(t(outer(ct, nm, paste, sep = ":")))
  dimnames(J) <- list(lab, lab)
  se <- matrix(sqrt(diag(J)), nrow(weights), K, byrow = TRUE, dimnames = list(ct, nm))
  list(estimate = est, se = se, vcov = J)
}
