# ================================================================
# Effect on the Gini coefficient from the L-moment DiD
#
# For a nonnegative outcome the Gini coefficient is the L-CV
#   G = lambda_2 / lambda_1,   lambda_2 = E|X - X'| / 2.
# With U-statistic L-moments, lambda_2_hat / lambda_1_hat is the Gini
# mean difference over distinct pairs divided by twice the mean, i.e.
# n/(n-1) times the usual plug-in sample Gini.
#
# 2x2 design with cells T1 (treated, post), T0, C1, C0. Under parallel
# trends in the L-moments, the treated group's counterfactual post
# L-moments are
#   lambda_k^cf = lambda_k(T0) + [lambda_k(C1) - lambda_k(C0)],  k = 1, 2,
# and the effect on the Gini coefficient is
#   Delta_G = lambda_2(T1)/lambda_1(T1) - lambda_2^cf/lambda_1^cf.
# Inference: delta method with the jackknife covariances of the four
# (independent) cells. With a = lambda(T1), c = lambda^cf,
#   dDelta/d lambda(T1) = ( -a2/a1^2, 1/a1 )
#   dDelta/d lambda(T0) = dDelta/d lambda(C1) = -( -c2/c1^2, 1/c1 )
#   dDelta/d lambda(C0) = +( -c2/c1^2, 1/c1 )
#   Var(Delta_G) = sum_cells g_cell' V_cell g_cell.
#
# did_lratio_cells() is the general version (any lambda_num/lambda_den,
# e.g. tau_3, tau_4) and is what eg_08 uses for all three ratio rows;
# did_gini_lmoment() is the Gini special case working from raw data.
# Requires lmoment_jackknife.R to be sourced first.
# ================================================================

# General L-moment ratio effect from jackknife cells (lmom_u_jackknife
# results named T1, T0, C1, C0, i.e. treated/control x post/pre):
#   rho = lambda_num / lambda_den, observed for T1 minus counterfactual
#   lambda^cf = lambda(T0) + lambda(C1) - lambda(C0).
# Examples: Gini = (2, 1); L-skewness tau_3 = (3, 2); L-kurtosis tau_4 = (4, 2).
# Delta method: d rho / d lambda_num = 1/lambda_den, d rho / d lambda_den
# = -lambda_num/lambda_den^2, applied to T1 (sign +) and to the
# counterfactual (T0, C1: -, C0: +); Var = sum_cells g' V_cell g.
did_lratio_cells <- function(cells, num, den, alpha = 0.05) {
  K <- length(cells$T1$lambda)
  stopifnot(max(num, den) <= K)
  grad <- function(l) { g <- numeric(K); g[num] <- 1 / l[den]; g[den] <- g[den] - l[num] / l[den]^2; g }
  a  <- cells$T1$lambda
  cf <- cells$T0$lambda + cells$C1$lambda - cells$C0$lambda
  r_obs <- a[num] / a[den]; r_cf <- cf[num] / cf[den]
  ga <- grad(a); gc <- grad(cf)
  G <- list(T1 = ga, T0 = -gc, C1 = -gc, C0 = gc)
  v <- sum(sapply(names(G), \(c) drop(t(G[[c]]) %*% cells[[c]]$vcov %*% G[[c]])))
  eff <- unname(r_obs - r_cf); se <- sqrt(v); z <- qnorm(1 - alpha / 2)
  list(effect = eff, se = se, ci = c(eff - z * se, eff + z * se),
       observed = unname(r_obs), counterfactual = unname(r_cf))
}

did_gini_lmoment <- function(y, treat, post, alpha = 0.05) {
  ok <- !is.na(y) & !is.na(treat) & !is.na(post)
  y <- y[ok]; treat <- treat[ok]; post <- post[ok]
  if (any(y < 0)) warning("negative outcomes: lambda_2/lambda_1 is then an L-CV, not a Gini coefficient")
  cell <- function(tr, po) lmom_u_jackknife(y[treat == tr & post == po], K = 2)
  cl <- list(T1 = cell(1, 1), T0 = cell(1, 0), C1 = cell(0, 1), C0 = cell(0, 0))
  a  <- cl$T1$lambda
  cf <- cl$T0$lambda + cl$C1$lambda - cl$C0$lambda
  stopifnot(cf[1] > 0, cf[2] > 0)
  g_obs <- a[2] / a[1]; g_cf <- cf[2] / cf[1]
  d_a  <- c(-a[2] / a[1]^2, 1 / a[1])
  d_cf <- c(-cf[2] / cf[1]^2, 1 / cf[1])
  grads <- list(T1 = d_a, T0 = -d_cf, C1 = -d_cf, C0 = d_cf)
  var_eff <- sum(sapply(names(cl), \(c) drop(t(grads[[c]]) %*% cl[[c]]$vcov %*% grads[[c]])))
  # SE of the observed and counterfactual Gini separately (for reporting)
  se_obs <- sqrt(drop(t(d_a) %*% cl$T1$vcov %*% d_a))
  se_cf  <- sqrt(sum(sapply(c("T0", "C1", "C0"), \(c) drop(t(d_cf) %*% cl[[c]]$vcov %*% d_cf))))
  eff <- unname(g_obs - g_cf); se <- sqrt(var_eff); z <- qnorm(1 - alpha / 2)
  list(effect = eff, se = se, ci = c(eff - z * se, eff + z * se),
       gini_observed = unname(g_obs), se_observed = se_obs,
       gini_counterfactual = unname(g_cf), se_counterfactual = se_cf,
       gini_cells = sapply(cl, \(c) unname(c$lambda[2] / c$lambda[1])),
       lambda_observed = a, lambda_counterfactual = cf,
       n = sapply(cl, `[[`, "n"))
}
