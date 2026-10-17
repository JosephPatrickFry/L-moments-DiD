# ================================================================
# Reproduce the two tables from the Engelhardt & Gruber (2011)
# replication package (see README_replication.md):
#   output/eg/oop_did_table_K4_pre2005_age60plus.tex   (empirical application)
#   output/sim_K/sim_K_table.tex                       (choose-K simulation)
# All paths (inputs, data/clean/, output/) are relative to the folder
# containing this file, wherever it is and whatever the working directory:
#   Rscript path/to/run_all.R            all steps
#   Rscript path/to/run_all.R skip_sim   skip the simulation replications
#                                        (the slow step); sim_K_report.R then
#                                        uses existing output/sim_K/sim_K_raw_*.rds
# (or source("path/to/run_all.R") from an R session).
# Each step runs in its own R process; the run stops at the first failure.
# ================================================================

# project root = the folder containing run_all.R
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (!length(this_file)) this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)   # source()
if (!length(this_file)) stop("Run this file with Rscript or source() so its location can be found.")
Sys.setenv(EA_ROOT = dirname(normalizePath(this_file[1], winslash = "/", mustWork = TRUE)))
source(file.path(Sys.getenv("EA_ROOT"), "code", "setup.R"))

pkgs <- c("haven", "dplyr", "tidyr", "foreach", "doParallel", "doRNG", "ggplot2", "scales")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  stop("Install the missing packages first: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")
if (getRversion() < "4.1.0") stop("R >= 4.1 is required (native pipe and \\(x) lambdas).")

skip_sim <- "skip_sim" %in% commandArgs(trailingOnly = TRUE)
steps <- c("eg_01_build.R",              # data/clean/eg_meps_2002_2007.rds
           "eg_08_oop_latex_table.R",    # empirical table
           if (!skip_sim) "sim_K_run.R", # simulation replications
           "sim_K_report.R")             # simulation table

rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type != "windows") rscript <- shQuote(rscript)   # system2 quotes the command itself on Windows
for (s in steps) {   # child processes inherit EA_ROOT
  cat(sprintf("\n==== %s ====\n", s))
  t0 <- Sys.time()
  status <- system2(rscript, shQuote(ea_path("code", s)))
  if (status != 0) stop(sprintf("%s failed (exit status %d)", s, status))
  cat(sprintf("---- %s done in %.0f s\n", s, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

writeLines(capture.output(sessionInfo()), file.path(ea_dir("output"), "sessionInfo.txt"))
cat("\nTables:\n ", ea_path("output", "eg", "oop_did_table_K4_pre2005_age60plus.tex"),
    "\n ", ea_path("output", "sim_K", "sim_K_table.tex"), "\n")
