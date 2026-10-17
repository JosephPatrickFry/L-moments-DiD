# ================================================================
# Engelhardt & Gruber (2011), "Medicare Part D and the Financial
# Protection of the Elderly", AEJ: Economic Policy 3(4):77-102.
# Replication package: openICPSR 114772 (folder 114772-V1/ in EA).
#
# The package ships the authors' analysis files, built from MEPS by
# primary_data_assembly_program.do:
#   meps_anal_indiv.dta         ages 60-70 (Tables 2-4, 6)
#   meps_anal_indiv_60plus.dta  ages 60+   (Tables 7-8, IVQR)
# Already applied there: MEPS 2002-2007; drop persons without end-of-
# year Rx payment measures, under-65s on Medicare (DI), and missing
# marital status; dollars in 2007 $ (CPI-U).
#
# Paper's analysis conventions (table*.do), applied here:
#   treated = age 65+ (rdage65), control = age 60-64
#   pre = 2002-2005, post = 2007; 2006 (transition) dropped in all
#     estimates -- kept here with a flag for the event study
#   complete cases on the table's variable list
#   unweighted; clusters = household x (age 65+)
# Outcomes: rrxoop (OOP Rx), rrxtot (total Rx), rrxgov2 (Medicare +
#   Medicaid Rx), rrxpriv2 (all other payers, i.e. private + other).
# Validation: unweighted deciles of OOP Rx in the 60+ sample against the
#   working-paper Table 7 (= published Table 8).
#
# Input:  114772-V1/AEJ-Data-Appendix/meps_anal_indiv.dta and
#         meps_anal_indiv_60plus.dta (download from openICPSR 114772)
# Output: data/clean/eg_meps_2002_2007.rds, both samples stacked
#         (column `sample`); the OOP table and the choose-K simulation
#         use only sample == "age60plus", complete == 1.
# Run from the EA folder:  Rscript code/eg_01_build.R
# ================================================================

library(haven)
library(dplyr)
library(tidyr)

source(file.path(Sys.getenv("EA_ROOT", "."), "code", "setup.R"))
pkg     <- ea_path("114772-V1", "AEJ-Data-Appendix")
out_dir <- ea_dir("data", "clean")
for (f in c("meps_anal_indiv.dta", "meps_anal_indiv_60plus.dta"))
  if (!file.exists(file.path(pkg, f)))
    stop("Missing ", file.path(pkg, f), ": unpack the openICPSR 114772 package into the EA folder.")

# variables kept by table8.do / table7.do, plus outcomes and weights
vars <- c("rdupersid", "rwho_ref", "year", "rage", "rdage65", "dpost", "rperwt",
          "rrxoop", "rrxtot", "rrxgov2", "rrxpriv2", "rrxnum",
          "rdrxacov", "rdrxgcov2", "rdrxprivcov2", "rdrxdoublecov",
          "rrace", "reduc", "rdfemale", "rmarstat", "rcenreg", "rshlt", "hincquint")
cc_vars <- c("rrxoop", "rrxtot", "rrxgov2", "rrxpriv2", "rage", "year", "rdrxacov", "rdrxgcov2",
             "rdrxprivcov2", "rdrxdoublecov", "rrace", "reduc", "rdfemale", "rmarstat",
             "rcenreg", "rshlt", "hincquint")

load_sample <- function(file, label) {
  read_dta(file.path(pkg, file), col_select = all_of(vars)) |>
    zap_labels() |>
    mutate(sample = label,
           treat = as.integer(rdage65),
           period = case_when(year <= 2005 ~ "pre", year == 2006 ~ "transition", TRUE ~ "post"),
           # paper's clusters: household x age-65 group
           clustid = as.integer(factor(paste(rwho_ref, rdage65))),
           complete = as.integer(complete.cases(across(all_of(cc_vars)))),
           # analysis sample used by the paper's tables: complete, 2006 dropped
           in_paper_sample = as.integer(complete == 1 & year != 2006))
}

eg <- bind_rows(load_sample("meps_anal_indiv.dta", "age60_70"),
                load_sample("meps_anal_indiv_60plus.dta", "age60plus"))
saveRDS(eg, file.path(out_dir, "eg_meps_2002_2007.rds"))

# ------------------------------------------------
# validation
# ------------------------------------------------
cat("Person-years (complete cases; 2006 shown separately):\n")
print(as.data.frame(eg |> filter(complete == 1) |>
  count(sample, treat, period) |> pivot_wider(names_from = period, values_from = n)))

cat("\nUnweighted OOP Rx deciles, 60+ sample, 2006 dropped (2007 $)\n",
    "vs working-paper Table 7 (published Table 8):\n", sep = "")
paper <- rbind(`65+ pre`  = c(1, 73, 170, 290, 453, 679, 990, 1492, 2494),
               `65+ post` = c(0, 45, 105, 174, 267, 395, 568, 846, 1415),
               `60-64 pre`  = c(0, 5, 46, 115, 210, 340, 505, 765, 1350),
               `60-64 post` = c(0, 5, 40, 102, 188, 290, 446, 648, 1178))
s <- filter(eg, sample == "age60plus", in_paper_sample == 1)
q <- function(tr, po) quantile(s$rrxoop[s$treat == tr & s$dpost == po], seq(.1, .9, .1), type = 2)
mine <- rbind(`65+ pre` = q(1, 0), `65+ post` = q(1, 1), `60-64 pre` = q(0, 0), `60-64 post` = q(0, 1))
cmp <- rbind(mine, paper); rownames(cmp) <- c(paste(rownames(mine), "(this)"), paste(rownames(paper), "(paper)"))
colnames(cmp) <- paste0("q", 1:9 * 10)
print(round(cmp[c(1, 5, 2, 6, 3, 7, 4, 8), ]))
dd <- (mine[2, ] - mine[1, ]) - (mine[4, ] - mine[3, ])
cat("\nQuantile DiD (this):  ", round(dd), "\n")
cat("Quantile DiD (paper): ", c(-1, -28, -59, -103, -164, -235, -363, -529, -907), "\n")
