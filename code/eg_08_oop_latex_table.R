# ================================================================
# LaTeX table: effect of Medicare Part D on out-of-pocket Rx spending,
# Engelhardt & Gruber (2011) data. 2x2 DiD of 2007 vs the pre-period
# (default 2005), age 65+ vs 60-64, age 60+ sample, complete cases,
# unweighted. Rows:
#   lambda_1..lambda_4  U-statistic DiD estimate, jackknife SE, pointwise
#                       95% CI, simultaneous 95% CI (Gaussian multiplier
#                       bootstrap, 10,000 draws, shared critical value)
#   tau_3, tau_4, Gini  L-skewness (lambda_3/lambda_2), L-kurtosis
#                       (lambda_4/lambda_2), Gini (lambda_2/lambda_1):
#                       observed minus counterfactual ratio, delta-method
#                       SE, pointwise 95% CI
#   SD lower bound      SD of individual effects among the treated is at
#                       least sqrt(3)|Delta lambda_2|; delta-method SE and
#                       one-sided 95% lower confidence bound with z_0.975
#                       (lmoment_ite_bound.R)
# Input:  data/clean/eg_meps_2002_2007.rds (eg_01_build.R)
# Output: output/eg/oop_did_table_K4_pre<pre>_<sample>.tex
# Run from the EA folder:
#   Rscript code/eg_08_oop_latex_table.R [pre=2005] [sample=age60plus]
# (pre may be a range, e.g. pre=2002-2005)
# ================================================================

library(dplyr)
source(file.path(Sys.getenv("EA_ROOT", "."), "code", "setup.R"))
for (f in c("lmoment_jackknife.R", "lmoment_inference.R", "lmoment_gini.R", "lmoment_ite_bound.R"))
  source(ea_path("code", f))

args <- commandArgs(trailingOnly = TRUE)
arg  <- function(key, default) { a <- sub(paste0("^", key, "="), "", grep(paste0("^", key, "="), args, value = TRUE)); if (length(a)) a else default }
pre_spec <- arg("pre", "2005"); samp <- arg("sample", "age60plus")
pre_years <- if (grepl("-", pre_spec)) do.call(seq, as.list(as.integer(strsplit(pre_spec, "-")[[1]]))) else as.integer(pre_spec)
K <- 4; alpha <- 0.05; B <- 10000

d <- readRDS(ea_path("data", "clean", "eg_meps_2002_2007.rds")) |>
  filter(sample == samp, complete == 1, year %in% c(pre_years, 2007)) |>
  mutate(post = as.integer(year == 2007))

ci <- did2x2_lmoment_ci(d$rrxoop, d$treat, d$post, K = K, alpha = alpha, B = B,
                        method = "multiplier", seed = 20261008)
g  <- did_gini_lmoment(d$rrxoop, d$treat, d$post, alpha = alpha)
crit <- attr(ci, "crit")

neg <- function(s) sub("^-", "$-$", s)
f1 <- function(x) neg(formatC(x, format = "f", digits = 1, big.mark = "{,}"))
f3 <- function(x) neg(formatC(x, format = "f", digits = 3))
br <- function(lo, hi, f) sprintf("[%s, %s]", f(lo), f(hi))
rows <- sprintf("$\\lambda_%d$ & %s & %s & %s & %s \\\\", ci$k, f1(ci$estimate), f1(ci$se),
                br(ci$lo_pointwise, ci$hi_pointwise, f1), br(ci$lo_simultaneous, ci$hi_simultaneous, f1))
# ratio effects (observed minus counterfactual), delta-method SEs, from the K = 4 jackknife cells
cells <- attr(ci, "cells"); names(cells) <- c("T1", "T0", "C1", "C0")   # T_post, T_pre, C_post, C_pre
ratio_row <- function(label, num, den) {
  r <- did_lratio_cells(cells, num, den, alpha = alpha)
  sprintf("%s & %s & %s & %s & --- \\\\", label, f3(r$effect), f3(r$se), br(r$ci[1], r$ci[2], f3))
}
tau3_row <- ratio_row("L-skewness ($\\tau_3$)", 3, 2)
tau4_row <- ratio_row("L-kurtosis ($\\tau_4$)", 4, 2)
gini_row <- ratio_row("Gini coefficient", 2, 1)
stopifnot(abs(did_lratio_cells(cells, 2, 1)$se - g$se) < 1e-12)   # same as did_gini_lmoment
# lower bound on the SD of individual effects among the treated: sqrt(3) |Delta lambda_2|,
# one-sided 95% lower confidence bound with z_{0.975} (uniformly valid; see lmoment_ite_bound.R)
ib <- ite_var_lower_bound(d$rrxoop, d$treat, d$post, alpha = alpha)
stopifnot(abs(ib$delta_lambda2 - ci$estimate[2]) < 1e-8)
sd_row <- sprintf("SD of individual effects, lower bound & %s & %s & $[%s, \\infty)$ & --- \\\\",
                  f1(ib$sd_bound), f1(ib$se_sd_bound), f1(ib$sd_lower_conf))
pre_lab <- if (length(pre_years) > 1) sprintf("%d--%d", min(pre_years), max(pre_years)) else pre_years
n <- attr(ci, "n")

tex <- c(
  "% Requires \\usepackage{booktabs} and \\usepackage{threeparttable}",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\begin{threeparttable}",
  "\\caption{Effect of Medicare Part D on Out-of-Pocket Prescription Drug Spending}",
  "\\label{tab:oop_lmoment_did}",
  "\\begin{tabular}{lrrcc}",
  "\\toprule",
  " & Estimate & SE & Pointwise 95\\% CI & Simultaneous 95\\% CI \\\\",
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{L-moments (2007 dollars)}} \\\\",
  rows,
  "\\addlinespace",
  "\\multicolumn{5}{l}{\\textit{L-moment ratios}} \\\\",
  tau3_row,
  tau4_row,
  "\\addlinespace",
  "\\multicolumn{5}{l}{\\textit{Inequality}} \\\\",
  gini_row,
  "\\addlinespace",
  "\\multicolumn{5}{l}{\\textit{Effect heterogeneity (2007 dollars)}} \\\\",
  sd_row,
  "\\midrule",
  sprintf("\\multicolumn{5}{l}{Shared multiplier-bootstrap critical value: $\\hat c = %.3f$} \\\\", crit),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]\\footnotesize",
  sprintf("\\item \\textit{Notes:} %s Gaussian multiplier draws used in calculating $\\hat{c}$. The SD lower bound is $\\sqrt{3}\\,|\\hat\\Delta_{\\lambda_2}|$; its confidence interval is one-sided, $\\sqrt{3}\\max\\{0, |\\hat\\Delta_{\\lambda_2}| - 1.96\\,\\mathrm{SE}\\}$. MEPS data of Engelhardt and Gruber (2011), %s sample, 2007 vs. %s comparison. Cell sizes: treated %s / %s, control %s / %s (pre / post).",
          format(B, big.mark = "{,}"), ifelse(samp == "age60plus", "age 60+", "age 60--70"), pre_lab,
          format(n[["T_pre"]], big.mark = "{,}"), format(n[["T_post"]], big.mark = "{,}"),
          format(n[["C_pre"]], big.mark = "{,}"), format(n[["C_post"]], big.mark = "{,}")),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}")

out <- file.path(ea_dir("output", "eg"), sprintf("oop_did_table_K4_pre%s_%s.tex", gsub("--", "-", pre_lab), samp))
writeLines(tex, out)
cat(tex, sep = "\n")
