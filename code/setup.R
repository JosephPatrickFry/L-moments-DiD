# ================================================================
# Project paths. The project root is the EA folder (the folder that
# contains code/ and 114772-V1/). Run scripts from it, e.g.
#   Rscript code/eg_08_oop_latex_table.R
# or set the environment variable EA_ROOT to its path (run_all.R sets
# EA_ROOT to the folder containing run_all.R). Scripts load this
# file with
#   source(file.path(Sys.getenv("EA_ROOT", "."), "code", "setup.R"))
# ================================================================

ea_root <- normalizePath(Sys.getenv("EA_ROOT", "."), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(ea_root, "code", "setup.R")))
  stop("Could not find code/setup.R under '", ea_root,
       "'. Run from the EA project folder or set EA_ROOT to it.")

# path under the project root
ea_path <- function(...) file.path(ea_root, ...)
# path under the project root, creating the directory if needed
ea_dir <- function(...) {
  p <- ea_path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}
