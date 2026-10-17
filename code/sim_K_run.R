# ================================================================
# Choose-K simulation: run the replications in parallel
# (foreach + doParallel; doRNG gives reproducible streams that do not
# depend on the number of workers).
# Run from the EA folder (needs data/clean/eg_meps_2002_2007.rds):
#   Rscript code/sim_K_run.R [cases=A,B,C,D] [R=1000] [cores=<n>]
# Writes output/sim_K/sim_K_raw_<case>.rds with one row per replication:
#   Khat_<K_max>   K_hat for each K_max (argmin of R_hat(K), K <= K_max)
#   Rhat1..        R_hat(K) path, K = 1..max(K_max) (saved for diagnostics;
#                  not used by sim_K_report.R)
#   lossK1..       L2 loss of the K-term estimate, K = 1..max(K_max)
#   loss_eqf       L2 loss of the type-1 empirical quantile function
# Results do not depend on the number of workers (doRNG streams). The
# default uses detectCores() - 2 workers; about 4.5 minutes on 30 workers.
# ================================================================

suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
source(file.path(Sys.getenv("EA_ROOT", "."), "code", "setup.R"))
source(ea_path("code", "sim_K_config.R"))
code_files <- ea_path("code", c("lmoment_choose_K.R", "sim_K_dgp.R"))
for (f in code_files) source(f)

args <- commandArgs(trailingOnly = TRUE)
arg  <- function(key, default) { a <- sub(paste0("^", key, "="), "", grep(paste0("^", key, "="), args, value = TRUE)); if (length(a)) a else default }
cases   <- strsplit(arg("cases", paste(sim_cfg$cases, collapse = ",")), ",")[[1]]
R       <- as.integer(arg("R", sim_cfg$R))
n_cores <- as.integer(arg("cores", sim_cfg$n_cores))
dir.create(sim_cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

x_cell <- sim_cell()
n      <- length(x_cell)
J      <- max(sim_cfg$K_max)
dgps   <- calibrate_dgps(x_cell, sim_cfg, J)

cl <- makeCluster(n_cores)
registerDoParallel(cl)
invisible(clusterCall(cl, function(files) { for (f in files) source(f); NULL }, code_files))
cat(sprintf("n = %d, R = %d, K_max = {%s}, N = %d, workers = %d\n",
            n, R, paste(sim_cfg$K_max, collapse = ", "), sim_cfg$N_part, n_cores))

for (cs in cases) {
  dgp <- dgps[[cs]]
  dIQ <- dgp$dIQ(n)
  K_max <- sim_cfg$K_max; N_part <- sim_cfg$N_part; basis <- sim_cfg$basis
  t0 <- Sys.time()
  res <- foreach(r = seq_len(R), .combine = rbind, .inorder = TRUE,
                 .options.RNG = sim_cfg$seed[[cs]]) %dorng% {
    x  <- dgp$sample(n)
    # B only enters the K-invariant penalty B^2 K_max / sqrt(n), so it does not affect K_hat
    ck <- choose_K_lmoment(x, B = max(abs(x)), K_max = J, N = N_part, basis = basis)
    Rhat <- ck$table$Rhat
    Khat <- sapply(K_max, \(km) which.min(Rhat[1:km]))
    c(setNames(Khat, paste0("Khat_", K_max)), setNames(Rhat, paste0("Rhat", 1:J)),
      sim_losses(x, ck$table$theta, dgp, dIQ))
  }
  attr(res, "rng") <- NULL
  rownames(res) <- NULL
  saveRDS(list(case = cs, label = dgp$label, par = dgp$par, theta = dgp$theta, IQ2 = dgp$IQ2,
               n = n, R = R, cfg = sim_cfg, res = res),
          file.path(sim_cfg$out_dir, sprintf("sim_K_raw_%s.rds", cs)))
  cat(sprintf("case %s done: %.1f s\n", cs, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
stopCluster(cl)
