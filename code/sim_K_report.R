# ================================================================
# Choose-K simulation: summary table (CSV + LaTeX) and K_hat / K* figure
# from output/sim_K/sim_K_raw_<case>.rds (sim_K_run.R).
# Risk = mean L2 loss over replications, in thousands of dollars^2.
#   K*        sample-specific oracle: in each replication,
#             argmin_{K <= K_max} ||Q_hat_K - Q||^2 (the realized loss)
#   K_hat     argmin_{K <= K_max} R_hat(K)
#   K_fixed*  (CSV only) best fixed K: argmin_{K <= K_max} of the simulated risk
# Risk columns in the LaTeX table are rounded to 3 significant figures.
# Run from the EA folder after sim_K_run.R:
#   Rscript code/sim_K_report.R
# Outputs (output/sim_K/): sim_K_table.tex, sim_K_summary.csv,
#   sim_K_khat.png/.pdf (distributions of K_hat and K*)
# ================================================================

library(ggplot2)
source(file.path(Sys.getenv("EA_ROOT", "."), "code", "setup.R"))
source(ea_path("code", "sim_K_config.R"))
out_dir <- sim_cfg$out_dir
runs <- lapply(sim_cfg$cases, \(cs) readRDS(file.path(out_dir, sprintf("sim_K_raw_%s.rds", cs))))
names(runs) <- sim_cfg$cases
q1 <- function(x, p) unname(quantile(x, p, type = 1))

tab <- list(); kdist <- list()
for (cs in names(runs)) {
  r <- runs[[cs]]; res <- r$res
  L <- res[, grep("^lossK", colnames(res)), drop = FALSE]
  risk <- colMeans(L); idx <- seq_len(nrow(L))
  for (km in sim_cfg$K_max) {
    kh <- res[, paste0("Khat_", km)]
    ks <- apply(L[, 1:km, drop = FALSE], 1, which.min)        # sample-specific oracle
    kf <- which.min(risk[1:km])                               # best fixed K
    tab[[paste(cs, km)]] <- data.frame(
      case = cs, label = r$label, K_max = km,
      Kstar_median = median(ks), Kstar_q25 = q1(ks, 0.25), Kstar_q75 = q1(ks, 0.75),
      Khat_median = median(kh), Khat_q25 = q1(kh, 0.25), Khat_q75 = q1(kh, 0.75),
      pct_Khat_eq_Kstar = 100 * mean(kh == ks), mean_abs_Khat_minus_Kstar = mean(abs(kh - ks)),
      pct_Kstar_eq_Kmax = 100 * mean(ks == km), pct_Khat_eq_Kmax = 100 * mean(kh == km),
      risk_Khat = mean(L[cbind(idx, kh)]) / 1e3, risk_Kstar = mean(L[cbind(idx, ks)]) / 1e3,
      risk_K4 = risk[sim_cfg$K_fixed] / 1e3, risk_eqf = mean(res[, "loss_eqf"]) / 1e3,
      K_fixed_star = kf, risk_K_fixed_star = risk[kf] / 1e3,
      R = nrow(res), n = r$n)
    kdist[[paste(cs, km)]] <- rbind(data.frame(case = r$label, K_max = km, which = "hat(K)", K = kh),
                                    data.frame(case = r$label, K_max = km, which = "K^'*'", K = ks))
  }
}
tab <- do.call(rbind, tab); rownames(tab) <- NULL
write.csv(tab, file.path(out_dir, "sim_K_summary.csv"), row.names = FALSE)
print(tab[, setdiff(names(tab), c("label", "R", "n"))], digits = 4, row.names = FALSE)

# ---- LaTeX table --------------------------------------------------------
fr <- function(x) formatC(signif(x, 3), format = "fg", digits = 3, flag = "#", big.mark = "{,}")   # 3 significant figures
mi <- function(m, a, b) sprintf("%g [%g, %g]", m, a, b)
rows <- character(0)
for (cs in names(runs)) {
  tb <- tab[tab$case == cs, ]
  rows <- c(rows, sprintf("\\multicolumn{9}{l}{\\textit{%s}} \\\\", sub("^[A-D]: ", paste0(cs, ". "), tb$label[1])))
  rows <- c(rows, sprintf("& %d & %s & %s & %.1f & %s & %s & %s & %s \\\\",
                          tb$K_max, mi(tb$Kstar_median, tb$Kstar_q25, tb$Kstar_q75),
                          mi(tb$Khat_median, tb$Khat_q25, tb$Khat_q75), tb$pct_Khat_eq_Kstar,
                          fr(tb$risk_Khat), fr(tb$risk_Kstar), fr(tb$risk_K4), fr(tb$risk_eqf)))
  if (cs != tail(names(runs), 1)) rows <- c(rows, "\\addlinespace")
}
tex <- c(
  "% Requires \\usepackage{booktabs} and \\usepackage{threeparttable}",
  "\\begin{table}[htbp]", "\\centering", "\\begin{threeparttable}",
  "\\caption{Choosing $K$ by the Unbiased Risk Estimate}",
  "\\label{tab:sim_choose_K}",
  "\\begin{tabular}{lrccrrrrr}", "\\toprule",
  " & & & & & \\multicolumn{4}{c}{Risk ($\\times 10^{3}$)} \\\\", "\\cmidrule(lr){6-9}",
  " & $K_{\\max}$ & $K^{*}$: median [IQR] & $\\hat K$: median [IQR] & \\% $\\hat K = K^{*}$ & $\\hat K$ & $K^{*}$ & $K = 4$ & $\\hat Q_n$ \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]\\footnotesize",
  sprintf(paste0("\\item \\textit{Notes:} %s replications of samples of size $n = %s$. Risk is the mean over replications of ",
                 "$\\|\\hat Q - Q\\|_{2}^{2}$ (2007 dollars squared) for the $K$-term L-moment estimate at $\\hat K$, at the ",
                 "sample-specific oracle $K^{*}$ (the $K \\le K_{\\max}$ minimizing the realized loss in each replication), and at $K = 4$, ",
                 "and for the empirical quantile function $\\hat Q_n$. $\\hat K$ minimizes $\\hat R(K)$ over ",
                 "$K \\le K_{\\max}$ with $M_j$ from %d random half-partitions. Designs calibrated to out-of-pocket prescription drug spending ",
                 "of the age 65+ group in 2007 (Engelhardt and Gruber 2011 data, age 60+ sample): A uniform on $[0, 2\\bar{x}]$; ",
                 "B $Q(u) = U(e^{cu}-1)/(e^{c}-1)$ with $U$ the 90th percentile; C zero-inflated generalized Pareto with Hill tail index; ",
                 "D resampling the observed data."),
          format(sim_cfg$R, big.mark = "{,}"), format(runs[[1]]$n, big.mark = "{,}"), sim_cfg$N_part),
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
writeLines(tex, file.path(out_dir, "sim_K_table.tex"))

# ---- figure: distributions of K_hat and the sample-specific K* ------------
kdist <- do.call(rbind, kdist)
kdist$K_max <- factor(kdist$K_max, levels = sim_cfg$K_max, labels = paste0("K[max]==", sim_cfg$K_max))
kdist$case  <- factor(kdist$case, levels = sapply(runs, `[[`, "label"))
kdist$which <- factor(kdist$which, levels = c("hat(K)", "K^'*'"))
p <- ggplot(kdist, aes(K, fill = which)) +
  geom_bar(aes(y = after_stat(prop)), position = position_dodge(preserve = "single"), width = 0.8) +
  facet_grid(case ~ K_max, scales = "free", labeller = labeller(K_max = label_parsed, case = label_wrap_gen(14))) +
  scale_x_continuous(breaks = \(l) { b <- scales::breaks_extended(6)(l); b[b == round(b) & b >= 1] }) +
  scale_fill_manual(values = c("#2a78d6", "#eb6834"), labels = scales::label_parse()) +
  labs(x = "K", y = "Share of replications", fill = NULL,
       title = expression("Distributions of " * hat(K) * " and the sample-specific oracle " * K^"*"),
       subtitle = sprintf("n = %s, %s replications", format(runs[[1]]$n, big.mark = ","), format(sim_cfg$R, big.mark = ","))) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", legend.justification = "left", panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0),
        plot.subtitle = element_text(colour = "#52514e", size = 9), plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(out_dir, "sim_K_khat.png"), p, width = 10, height = 8, dpi = 180)
ggsave(file.path(out_dir, "sim_K_khat.pdf"), p, width = 10, height = 8)
