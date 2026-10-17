# ================================================================
# Data-generating processes for the choose-K simulation (sim_K_config.R).
# Each DGP is a list with
#   sample(n)   draws n observations
#   theta       true orthonormal coefficients theta_j, j = 1..J
#               (theta_j = <Q, phi_j>, phi_j = sqrt(2j-1) P*_{j-1})
#   IQ2         int_0^1 Q(u)^2 du
#   dIQ(n)      int_{(i-1)/n}^{i/n} Q(u) du, i = 1..n (for the exact
#               L2 loss of the type-1 empirical quantile function)
#   par         calibrated parameters
# L2 loss of a K-term estimate (Parseval):
#   ||Q_hat_K - Q||^2 = sum_{j<=K} (theta_hat_j - theta_j)^2 + IQ2 - sum_{j<=K} theta_j^2
# L2 loss of the empirical quantile function Q_n(u) = x_(ceiling(n u)):
#   ||Q_n - Q||^2 = sum_i x_(i)^2 / n - 2 sum_i x_(i) dIQ_i + IQ2
# ================================================================

# Shifted Legendre polynomials P*_0..P*_{J-1} at u (length(u) x J)
shifted_legendre <- function(u, J) {
  P <- matrix(0, length(u), J); P[, 1] <- 1
  if (J >= 2) P[, 2] <- 2 * u - 1
  if (J >= 3) for (m in 1:(J - 2)) P[, m + 2] <- ((2 * m + 1) * (2 * u - 1) * P[, m + 1] - m * P[, m]) / (m + 1)
  P
}
phi_basis <- function(u, J) sweep(shifted_legendre(u, J), 2, sqrt(2 * seq_len(J) - 1), `*`)

# Antiderivatives of phi_j, j = 1..J (int_0^u P*_r = (P*_{r+1} - P*_{r-1}) / (2(2r+1)) + const)
phi_antideriv <- function(u, J) {
  P <- shifted_legendre(u, J + 1)
  F <- matrix(0, length(u), J); F[, 1] <- u
  if (J >= 2) for (j in 2:J) F[, j] <- (P[, j + 1] - P[, j - 1]) / (2 * (2 * j - 1))
  sweep(F, 2, sqrt(2 * seq_len(J) - 1), `*`)
}

# Gauss-Legendre nodes/weights on [0, 1] (Golub-Welsch)
gauss_legendre01 <- function(m) {
  b <- seq_len(m - 1) / sqrt(4 * seq_len(m - 1)^2 - 1)
  Tm <- matrix(0, m, m); Tm[cbind(1:(m - 1), 2:m)] <- b; Tm[cbind(2:m, 1:(m - 1))] <- b
  e <- eigen(Tm, symmetric = TRUE)
  list(x = (e$values + 1) / 2, w = e$vectors[1, ]^2)
}

# ---- A: Uniform[0, U] -------------------------------------------------
dgp_uniform <- function(U, J = 20) {
  th <- numeric(J); th[1] <- U / 2; th[2] <- U / (2 * sqrt(3))
  list(case = "A", label = "A: Uniform (exactly sparse)",
       sample = function(n) runif(n, 0, U), theta = th, IQ2 = U^2 / 3,
       dIQ = function(n) U * ((1:n)^2 - (0:(n - 1))^2) / (2 * n^2),
       par = c(U = U))
}

# ---- B: bounded convex, Q(u) = U (exp(c u) - 1) / (exp(c) - 1) ---------
dgp_expq <- function(U, c, J = 20, nodes = 200) {
  Q  <- function(u) U * expm1(c * u) / expm1(c)
  IQ <- function(u) U * (expm1(c * u) / c - u) / expm1(c)
  gl <- gauss_legendre01(nodes)
  th <- drop(crossprod(phi_basis(gl$x, J), gl$w * Q(gl$x)))
  IQ2 <- U^2 / expm1(c)^2 * (expm1(2 * c) / (2 * c) - 2 * expm1(c) / c + 1)
  list(case = "B", label = "B: Bounded, fast decay",
       sample = function(n) Q(runif(n)), theta = th, IQ2 = IQ2,
       dIQ = function(n) diff(IQ((0:n) / n)), par = c(U = U, c = c))
}

# ---- C: zero-inflated generalized Pareto --------------------------------
#   Q(u) = 0 for u <= p0, else (sigma/xi) ((1 - v)^-xi - 1), v = (u - p0)/(1 - p0); xi < 1/2
dgp_zigpd <- function(p0, sigma, xi, J = 20, nodes = 400) {
  stopifnot(xi > 0, xi < 0.5)
  Q  <- function(u) { v <- pmax(u - p0, 0) / (1 - p0); (sigma / xi) * ((1 - v)^(-xi) - 1) }
  IQ <- function(u) { v <- pmax(u - p0, 0) / (1 - p0)
    (1 - p0) * (sigma / xi) * ((1 - (1 - v)^(1 - xi)) / (1 - xi) - v) }
  # theta_j = (1-p0)(sigma/xi) int_0^1 [(1-v)^-xi - 1] phi_j(p0 + (1-p0) v) dv; the singular part
  # is made smooth by w = (1-v)^(1-xi): int (1-v)^-xi g(v) dv = int_0^1 g(1 - w^(1/(1-xi))) dw / (1-xi)
  gl <- gauss_legendre01(nodes)
  ph <- function(v) phi_basis(p0 + (1 - p0) * v, J)
  sing <- drop(crossprod(ph(1 - gl$x^(1 / (1 - xi))), gl$w)) / (1 - xi)
  reg  <- drop(crossprod(ph(gl$x), gl$w))
  th <- (1 - p0) * (sigma / xi) * (sing - reg)
  IQ2 <- (1 - p0) * (sigma / xi)^2 * (1 / (1 - 2 * xi) - 2 / (1 - xi) + 1)
  list(case = "C", label = "C: Heavy tail (zero-inflated GPD)",
       sample = function(n) Q(runif(n)), theta = th, IQ2 = IQ2,
       dIQ = function(n) diff(IQ((0:n) / n)), par = c(p0 = p0, sigma = sigma, xi = xi))
}

# ---- D: resampling a fixed data set --------------------------------------
dgp_resample <- function(y, J = 20) {
  y <- sort(y); m <- length(y)
  Fg <- phi_antideriv((0:m) / m, J)
  th <- drop(crossprod(diff(Fg), y))           # sum_i y_(i) int_{(i-1)/m}^{i/m} phi_j
  IQ <- c(0, cumsum(y)) / m                    # int_0^{k/m} Q
  list(case = "D", label = "D: Resampled data",
       sample = function(n) y[sample.int(m, n, replace = TRUE)], theta = th, IQ2 = mean(y^2),
       dIQ = function(n) { if (n == m) y / m else stop("dIQ for D implemented for n = m only") },
       par = c(m = m))
}

# ---- calibration to the 2007 treated cell (rules in sim_K_config.R) -----
calibrate_dgps <- function(x, cfg, J = max(cfg$K_max)) {
  x <- x[!is.na(x)]
  # A: match the mean
  A <- dgp_uniform(2 * mean(x), J)
  # B: bound at the B_q quantile, c matched to the mean below the bound
  U_B <- unname(quantile(x, cfg$B_q, type = 1)); mB <- mean(x[x <= U_B]) / U_B
  muc <- function(c) (expm1(c) / c - 1) / expm1(c)
  cB  <- uniroot(\(c) muc(c) - mB, c(1e-6, 200), tol = 1e-12)$root
  B <- dgp_expq(U_B, cB, J)
  # C: zero share, Hill tail index on top hill_frac of positives, median of positives
  xp <- sort(x[x > 0]); k <- floor(cfg$hill_frac * length(xp)); mp <- length(xp)
  xi <- mean(log(xp[(mp - k + 1):mp])) - log(xp[mp - k])
  sigma <- median(xp) * xi / (2^xi - 1)
  C <- dgp_zigpd(mean(x == 0), sigma, xi, J)
  # D: the cell itself
  D <- dgp_resample(x, J)
  list(A = A, B = B, C = C, D = D)
}

# L2 loss of the K-term estimates (K = 1..J) and of the empirical quantile function
sim_losses <- function(x, theta_hat, dgp, dIQ) {
  J <- length(theta_hat)
  lossK <- unname(cumsum((theta_hat - dgp$theta[1:J])^2) + dgp$IQ2 - cumsum(dgp$theta[1:J]^2))
  xs <- sort(x)
  c(lossK = lossK, loss_eqf = sum(xs^2) / length(xs) - 2 * sum(xs * dIQ) + dgp$IQ2)
}
