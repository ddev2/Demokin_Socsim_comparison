# Step 1. Run the Socsim simulations of the paper.
#
# Each scenario of scenarios/scenarios.csv is run with its seed. A scenario whose
# results already exist for this seed is skipped (runSocsim(..., overwrite = TRUE) to
# run it again). Each run takes from a few minutes (high mortality) to several hours
# (low mortality and high fertility: populations of several million persons).
# Requires the R package rsocsim (https://github.com/MPIDR/rsocsim).

# ==== 0. Setup ====

# Folder of the repository: the folder that contains config.R, found from the script
# open in RStudio, from Rscript, or from the working directory
local({
  d <- tryCatch(dirname(rstudioapi::getActiveDocumentContext()$path), error = function(e) "")
  if (!nzchar(d)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    d <- if (length(f) == 1) dirname(normalizePath(f)) else getwd()
  }
  d <- normalizePath(d)
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("cannot find config.R: setwd() to the repository")
  REPO_DIR <<- d
})
source(file.path(REPO_DIR, "config.R"))


# ==== 1. Scenarios ====

# All the scenarios of the paper; replace by a subset to run only some of them, e.g.
# todo <- c("TFR=4_HM", "TFR=4_HM_bint12")
todo  <- SCENARIOS$scenario
seeds <- SCENARIOS$seed[match(todo, SCENARIOS$scenario)]


# ==== 2. Run ====

status <- runSocsimScenarios(todo, seeds)
