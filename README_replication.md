# Replication: Medicare Part D application and choose-K simulation

This reproduces two tables:

| Table | File |
|---|---|
| Effect of Medicare Part D on out-of-pocket prescription drug spending (L-moment DiD) | `output/eg/oop_did_table_K4_pre2005_age60plus.tex` |
| Choosing $K$ by the unbiased risk estimate (simulation) | `output/sim_K/sim_K_table.tex` |

## Data

Download the replication package for Engelhardt and Gruber (2011), "Medicare Part D and the Financial Protection of the Elderly," *AEJ: Economic Policy* 3(4): 77–102, from openICPSR, project 114772 (a free login is required, so the data are not redistributed here).

Unpack it into this folder so that these two files exist:

```
114772-V1/AEJ-Data-Appendix/meps_anal_indiv_60plus.dta   (used by both tables)
114772-V1/AEJ-Data-Appendix/meps_anal_indiv.dta          (read by the build script)
```

No other file from the package is needed.

## Software

- R 4.1 or later (tested with R 4.5.1 on Windows 11).
- Packages: `haven`, `dplyr`, `tidyr`, `foreach`, `doParallel`, `doRNG`, `ggplot2`, `scales`.

```r
install.packages(c("haven", "dplyr", "tidyr", "foreach", "doParallel", "doRNG", "ggplot2", "scales"))
```

`run_all.R` writes the versions used to `output/sessionInfo.txt`.

## Running

```
Rscript run_all.R
```

`run_all.R` treats the folder that contains it (the one with `code/` and `114772-V1/`) as the project root, whatever the working directory, so `Rscript path/to/run_all.R` also works. All outputs are written under that folder: `data/clean/` and `output/`. It can also be run with `source("path/to/run_all.R")` from an R session.

The scripts can also be run one at a time, in the order below, from the project folder (or from anywhere, with the environment variable `EA_ROOT` set to the project folder's path).

| Step | Script | Output |
|---|---|---|
| 1 | `code/eg_01_build.R` | `data/clean/eg_meps_2002_2007.rds`; prints a check against the paper's OOP decile table |
| 2 | `code/eg_08_oop_latex_table.R` | empirical table |
| 3 | `code/sim_K_run.R` | `output/sim_K/sim_K_raw_{A,B,C,D}.rds` |
| 4 | `code/sim_K_report.R` | simulation table, `sim_K_summary.csv`, `sim_K_khat.png/.pdf` |

Steps 3–4 need step 1, because the simulation designs are calibrated to the 2007 treated cell. Output folders are created as needed.

**Runtime.** Steps 1, 2 and 4 take under a minute each. Step 3 (4,000 replications) took about 4.5 minutes on 30 workers. It uses `parallel::detectCores() - 2` workers by default, so a 4-core machine needs roughly an hour; override with `Rscript code/sim_K_run.R cores=<n>`.

**Reproducibility.** Seeds are fixed in `code/eg_08_oop_latex_table.R` (multiplier bootstrap) and `code/sim_K_config.R` (one per design). The simulation uses `doRNG`, so results do not depend on the number of workers.

## Code files

| File | Contents |
|---|---|
| `code/setup.R` | Project root and path helpers, loaded by every script |
| `code/lmoment_jackknife.R` | U-statistic L-moments with exact leave-one-out jackknife covariance |
| `code/lmoment_inference.R` | Multiplier-bootstrap critical value; 2×2 DiD pointwise and simultaneous CIs |
| `code/lmoment_gini.R` | DiD effects on L-moment ratios (Gini, τ₃, τ₄), delta-method SEs |
| `code/lmoment_ite_bound.R` | Lower bound on the SD (and variance) of individual treatment effects |
| `code/lmoment_choose_K.R` | Numerically stable L-moments for large K; risk estimate R̂(K); choice of K |
| `code/sim_K_config.R` | Simulation settings, seeds and calibration rules |
| `code/sim_K_dgp.R` | Simulation designs A–D: samplers, true coefficients, exact losses |

## Using the estimation and inference functions in another application

The functions take plain vectors, so they apply to any outcome observed in independent group × period cells. The design is unweighted and has no covariates: each cell is treated as an independent iid sample, and the DiD is a contrast of cell-level U-statistic L-moments. Load them with:

```r
source("code/setup.R")   # or set EA_ROOT; defines ea_path()
for (f in c("lmoment_jackknife.R", "lmoment_inference.R", "lmoment_gini.R",
            "lmoment_ite_bound.R", "lmoment_choose_K.R"))
  source(ea_path("code", f))
```

The examples below use a data frame `d` with an outcome `y`, a 0/1 group indicator `treat`, a 0/1 period indicator `post`, and (for the event study) a `year` column:

```r
set.seed(1)
n <- 4000
d <- data.frame(treat = rbinom(n, 1, 0.5), year = sample(2004:2007, n, replace = TRUE))
d$post <- as.integer(d$year == 2007)
d$y <- rexp(n, rate = 1 / 500) * ifelse(d$treat == 1 & d$post == 1, 0.7, 1)
```

### `lmom_u_jackknife(x, K = 4)` — L-moments of one sample

Computes the unbiased (U-statistic) sample L-moments $\hat\lambda_1, \dots, \hat\lambda_K$ of `x` and their exact leave-one-out jackknife covariance, in $O(nK)$ time after one sort. `NA`s are dropped.

It returns a list with these components:
- `lambda`: the K estimates.
- `vcov`: the covariance of `lambda`, $\hat\Sigma/n$.
- `se`: the standard errors.
- `Sigma`: the asymptotic covariance of $\sqrt n(\hat\lambda - \lambda)$.
- `influence`: the $n \times K$ jackknife influence values, in the original order of `x`.
- `n`: the sample size.

It uses the probability-weighted-moment formulas, which are accurate up to about K = 15. For larger K, use `lmom_u_stable()`.

```r
jk <- lmom_u_jackknife(d$y[d$treat == 1 & d$post == 1], K = 4)
jk$lambda; jk$se
```

### `did2x2_lmoment_ci(y, treat, post, K = 4, alpha = 0.05, B = 10000, method = "multiplier", seed = NULL)` — 2×2 DiD with confidence intervals

The DiD of the first K L-moments, $\hat\Delta_k = [\hat\lambda_k(T_1) - \hat\lambda_k(T_0)] - [\hat\lambda_k(C_1) - \hat\lambda_k(C_0)]$, with jackknife standard errors. Each row of the returned data frame gives `k`, `estimate`, `se`, the pointwise CI (`lo_pointwise`, `hi_pointwise`) and the simultaneous CI over $k = 1, \dots, K$ (`lo_simultaneous`, `hi_simultaneous`).

The simultaneous intervals use a Gaussian multiplier-bootstrap critical value with B draws: `method = "multiplier"` draws one multiplier per observation, while `method = "gaussian"` draws directly from the implied normal distribution and gives the same result in distribution. Set `seed` for reproducibility.

The data frame also carries these attributes:
- `crit`: the critical value $\hat c$.
- `n`: the cell sizes.
- `cells`: the four `lmom_u_jackknife()` results, named `T_post`, `T_pre`, `C_post`, `C_pre`, for reuse below.

```r
ci <- did2x2_lmoment_ci(d$y, d$treat, d$post, K = 4, B = 10000, seed = 123)
ci                  # estimates, SEs, pointwise and simultaneous 95% CIs
attr(ci, "crit")    # shared critical value c-hat
```

### `did_lratio_cells(cells, num, den, alpha = 0.05)` — effects on L-moment ratios

Gives the effect on a ratio $\rho = \lambda_{num}/\lambda_{den}$ as the observed ratio for the treated post cell minus the counterfactual ratio. The counterfactual L-moments are $\lambda^{cf} = \lambda(T_0) + \lambda(C_1) - \lambda(C_0)$, i.e., parallel trends in the L-moments. The SE comes from the delta method using each cell's jackknife covariance.

`cells` must be a list of four `lmom_u_jackknife()` results named `T1`, `T0`, `C1`, `C0`, each with at least `max(num, den)` L-moments. It returns `effect`, `se`, `ci` (pointwise), `observed` and `counterfactual`.

Common choices:
- Gini coefficient: `num = 2, den = 1` (for a nonnegative outcome).
- L-skewness $\tau_3$: `num = 3, den = 2`.
- L-kurtosis $\tau_4$: `num = 4, den = 2`.

`did_gini_lmoment(y, treat, post)` computes the Gini case directly from the data.

```r
cells <- attr(ci, "cells")
names(cells) <- c("T1", "T0", "C1", "C0")   # T_post, T_pre, C_post, C_pre
gini <- did_lratio_cells(cells, num = 2, den = 1)
tau3 <- did_lratio_cells(cells, num = 3, den = 2)
c(gini$effect, gini$se); gini$ci
```

### `ite_var_lower_bound(y, treat, post, alpha = 0.05)` — heterogeneity of individual effects

Under parallel trends in $\lambda_2$, the Gini mean difference (GMD $= 2\lambda_2$) gives a lower bound on the dispersion of individual treatment effects among the treated:
$$\mathbb V(Y(1) - Y(0) \mid G = 1) \ge \tfrac34\{\mathrm{GMD}(F^1_{11}) - \mathrm{GMD}(F^0_{11})\}^2 = 3\Delta_{\lambda_2}^2,$$
so $\mathrm{SD} \ge \sqrt3\,|\Delta_{\lambda_2}|$.

The function returns these components:
- `sd_bound`: the estimate $\sqrt3|\hat\Delta_{\lambda_2}|$.
- `se_sd_bound`: its delta-method SE.
- `sd_lower_conf`: a one-sided $(1-\alpha)$ lower confidence bound, $\sqrt3\max\{0, |\hat\Delta_{\lambda_2}| - z_{1-\alpha/2}\,\mathrm{SE}\}$. It uses $z_{1-\alpha/2}$ so that coverage holds even when $\Delta_{\lambda_2}$ is near 0.
- The variance-scale versions: `bound`, `bound_bc` (bias-corrected), `se_bound` and `lower_conf`.

```r
ib <- ite_var_lower_bound(d$y, d$treat, d$post)
c(ib$sd_bound, ib$se_sd_bound, ib$sd_lower_conf)
```

### `lmoment_contrast_vcov(cells, weights)` and `multiplier_critical_value(cells, weights, alpha, B, method, seed)` — general contrasts

These generalize the 2×2 DiD to any set of linear contrasts of independent cells, such as an event study.

Arguments:
- `cells`: a named list of `lmom_u_jackknife()` results.
- `weights`: a matrix with one row per contrast and one column per cell; the column names must match the cell names.

What they return:
- `lmoment_contrast_vcov()` returns the contrast `estimate` and `se` (contrasts × orders) and the joint `vcov`.
- `multiplier_critical_value()` additionally returns `crit`, the critical value for confidence intervals that hold simultaneously over all contrasts and orders.

```r
cell <- function(tr, yr) lmom_u_jackknife(d$y[d$treat == tr & d$year == yr], K = 4)
cl <- list(T04 = cell(1, 2004), T05 = cell(1, 2005), T06 = cell(1, 2006), T07 = cell(1, 2007),
           C04 = cell(0, 2004), C05 = cell(0, 2005), C06 = cell(0, 2006), C07 = cell(0, 2007))
# event study relative to 2005: (T_t - T_2005) - (C_t - C_2005)
w <- rbind(`2004` = c(1, -1, 0, 0, -1, 1,  0,  0),
           `2006` = c(0, -1, 1, 0,  0, 1, -1,  0),
           `2007` = c(0, -1, 0, 1,  0, 1,  0, -1))
colnames(w) <- names(cl)
cv <- lmoment_contrast_vcov(cl, w)
mb <- multiplier_critical_value(cl, w, B = 10000, seed = 123)
cv$estimate                                   # years x L-moment orders
lo <- cv$estimate - mb$crit * cv$se           # simultaneous over years and orders
hi <- cv$estimate + mb$crit * cv$se
```

### `choose_K_lmoment(x, B, K_max, N = 200, seed = NULL)` — choosing the number of L-moments

Selects $\hat K$ by minimizing the unbiased risk estimate $\hat R(K) = \sum_{j \le K}(\hat\theta_j^2 - 2M_j)$ for the $L^2$ distance between the K-term quantile approximation and the true quantile function. Here $\theta_j = \sqrt{2j-1}\,\lambda_j$, and $M_j$ averages $\hat\theta_j(A)\hat\theta_j(B)$ over N random half-partitions $(A, B)$ of the sample.

Arguments:
- `B`: a bound on $|X|$. It only enters the constant penalty $B^2 K_{\max}/\sqrt n$, so it does not affect $\hat K$.
- `K_max`: the largest K considered. It warns if `K_max` exceeds 10% of the half-sample size, where high-order estimates become unstable.

It returns `K` and a `table` giving $\hat\lambda_j$, $\hat\theta_j$, $M_j$ and $\hat R(K)$ for each K.

```r
ck <- choose_K_lmoment(d$y[d$treat == 1 & d$post == 1], B = max(d$y), K_max = 20, N = 200, seed = 1)
ck$K
```

The other scripts in `code/` belong to other applications or exploratory analyses and are not needed for these tables.
