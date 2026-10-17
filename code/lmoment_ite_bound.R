# ================================================================
# Lower bound on the variance of individual treatment effects among
# the treated from the L-moment DiD of lambda_2:
#   V(Y_1(1) - Y_1(0) | G = 1) >= (3/4) {GMD(F_11^1) - GMD(F_11^0)}^2,
# with GMD = E|X - X'| = 2 lambda_2, so the bound is 3 * Delta_2^2,
#   Delta_2 = lambda_2(T1) - [lambda_2(T0) + lambda_2(C1) - lambda_2(C0)].
# Inference uses the jackknife variance of Delta_2 (independent cells):
#   bound_hat     = 3 Delta_hat^2            (plug-in; E = 3(Delta^2 + se^2))
#   bound_bc      = 3 max(0, Delta_hat^2 - se^2)   (bias-corrected)
#   se_bound      = 6 |Delta_hat| se          (delta method)
#   lower_conf    = 3 max(0, |Delta_hat| - z_{1-alpha/2} se)^2
# lower_conf is a valid one-sided (1 - alpha) lower confidence bound for
# the ITE variance: with prob. >= 1 - alpha, Delta lies in the two-sided
# CI, so V >= 3 Delta^2 >= 3 min_{d in CI} d^2. (Using z_{1-alpha} would
# under-cover near Delta = 0 because the sign of Delta_hat is random.)
# SD version, SD(ITE) >= sqrt(3) |Delta_2|:
#   sd_bound      = sqrt(3) |Delta_hat|  (bias E|Delta_hat| - |Delta| ~ se*phi(Delta/se),
#                                         negligible unless |Delta|/se is small)
#   se_sd_bound   = sqrt(3) se            (delta method, Delta != 0)
#   sd_lower_conf = sqrt(3) max(0, |Delta_hat| - z_{1-alpha/2} se) = sqrt(lower_conf)
# Non-coverage is alpha/2 + Phi(-2|Delta|/se - z_{1-alpha/2}) <= alpha, so the
# two-sided critical value gives coverage >= 1 - alpha uniformly in Delta;
# z_{1-alpha} would cover 1 - 2 alpha at Delta = 0.
# The OOP table (eg_08) reports the SD version; the variance quantities
# are returned for eg_09_ite_variance_bound.R.
# Requires lmoment_jackknife.R to be sourced first.
# ================================================================

ite_var_lower_bound <- function(y, treat, post, alpha = 0.05) {
  ok <- !is.na(y) & !is.na(treat) & !is.na(post)
  y <- y[ok]; treat <- treat[ok]; post <- post[ok]
  cell <- function(tr, po) lmom_u_jackknife(y[treat == tr & post == po], K = 2)
  cl <- list(T1 = cell(1, 1), T0 = cell(1, 0), C1 = cell(0, 1), C0 = cell(0, 0))
  l2 <- sapply(cl, \(c) unname(c$lambda[2]))
  cf <- l2[["T0"]] + l2[["C1"]] - l2[["C0"]]
  d  <- l2[["T1"]] - cf
  se <- sqrt(sum(sapply(cl, \(c) c$vcov[2, 2])))
  z  <- qnorm(1 - alpha / 2)
  lc <- 3 * max(0, abs(d) - z * se)^2
  list(delta_lambda2 = d, se_lambda2 = se,
       gmd_observed = 2 * l2[["T1"]], gmd_counterfactual = 2 * cf, gmd_diff = 2 * d, se_gmd_diff = 2 * se,
       bound = 3 * d^2, bound_bc = 3 * max(0, d^2 - se^2), se_bound = 6 * abs(d) * se,
       lower_conf = lc, sd_bound = sqrt(3) * abs(d), se_sd_bound = sqrt(3) * se, sd_lower_conf = sqrt(lc),
       alpha = alpha, n = sapply(cl, `[[`, "n"))
}
