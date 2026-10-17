# ================================================================
# Simulation design for choosing K (fixed before running).
# Single-distribution designs calibrated to the 2007 treated cell of the
# Engelhardt & Gruber (2011) application: OOP Rx spending, age 65+,
# 2007, age 60+ sample, complete cases (n = 3,470).
#
#   A  exactly sparse      Uniform[0, U_A], U_A = 2 * cell mean
#                          (theta_j = 0 for j >= 3)
#   B  bounded, fast decay Q(u) = U_B (exp(c u) - 1) / (exp(c) - 1),
#                          U_B = cell 90th percentile (type 1), c matched
#                          to the cell mean of spending <= U_B
#   C  heavy tail          zero-inflated generalized Pareto: zero share =
#                          cell share of $0, shape xi = Hill estimate on
#                          the top 5% of positive spending, scale matched
#                          to the median of positive spending
#   D  the data            resample the cell with replacement
#                          (Q = the cell's type-1 empirical quantile function)
#
# Each replication draws n = 3,470 observations, computes R_hat(K) for
# K = 1..50 once (choose_K_lmoment, orthonormal basis, N = 200 random
# half-partitions), and reads off K_hat for K_max = 10, 20 and 50.
# Risk = average over replications of the realized L2 loss; the oracle
# K* is sample-specific (minimizes the realized loss in each replication).
#
# Pipeline (run from the EA folder, after code/eg_01_build.R):
#   Rscript code/sim_K_run.R      -> output/sim_K/sim_K_raw_<case>.rds
#   Rscript code/sim_K_report.R   -> output/sim_K/sim_K_table.tex, .csv, figure
# Sourced after code/setup.R (uses ea_path()).
# ================================================================

sim_cfg <- list(
  cases     = c("A", "B", "C", "D"),
  R         = 1000,                # replications per case
  K_max     = c(10, 20, 50),
  K_fixed   = 4,
  N_part    = 200,                 # random half-partitions for M_j
  basis     = "orthonormal",
  B_q       = 0.90,                # case B support bound: this quantile of the cell
  hill_frac = 0.05,                # case C: Hill estimate on the top 5% of positives
  seed      = c(A = 20261011, B = 20261012, C = 20261013, D = 20261014),
  n_cores   = max(1, parallel::detectCores() - 2),
  out_dir   = ea_path("output", "sim_K")
)

# The calibration cell (OOP Rx, age 65+, 2007, age 60+ sample, complete cases)
sim_cell <- function() {
  d <- readRDS(ea_path("data", "clean", "eg_meps_2002_2007.rds"))
  d$rrxoop[d$sample == "age60plus" & d$complete == 1 & d$year == 2007 & d$treat == 1]
}
