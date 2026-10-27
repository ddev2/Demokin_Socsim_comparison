# Step 3. DemoKin with the rates of each Socsim run.
#
# DemoKin (one-sex, time-invariant model) is run with the female rates estimated from
# each Socsim run over its last RATE_YEARS years and the proportion of female births of
# the run (section 4.3 of the paper). The rates are read from results/empirical_rates
# (written by step 2) or, if step 2 has not been run, from data/empirical_rates (the
# rates of the runs of the paper). With RATES = "input", DemoKin is run instead with
# the rates written in the Socsim rate files of each scenario.
#
# Output:
#   results/demokin_counts.csv        kin ever born at age EGO_AGE, formula [1] of GKP,
#                                     TFR, GRR, NRR and e0 of the rates, by scenario
#   results/demokin_kin_summary.rds   DemoKin kin_summary of each scenario (all ages)
# >>> rates saved for the later steps, start
#   results/demokin_inputs.rds        rates given to DemoKin, by scenario (used by step 4,
#                                     so that tau2 is computed with the same rates)
# <<< end

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

RATES <- "empirical"   # or "input"


# ==== 1. Run DemoKin ====

summaries <- list()
# >>> rates saved for the later steps, start
inputs <- list()
# <<< end
rows <- lapply(SCENARIOS$scenario, function(s) {
  if (RATES == "empirical") {
    f <- file.path(RESULTS_DIR, "empirical_rates", paste0(s, ".csv"))
    if (!file.exists(f)) f <- file.path(RATES_DIR, paste0(s, ".csv"))
    if (!file.exists(f)) { message(s, ": no rates"); return(NULL) }
    inp <- demokinInputs(utils::read.csv(f), MAX_AGE_DK)
  } else {
    inp <- inputsFromScenario(file.path(SCENARIO_DIR, s), MAX_AGE_DK)
  }
# >>> rates saved for the later steps, start
  inputs[[s]] <<- inp
# <<< end
  sdk <- runDemoKin(inp)
  summaries[[s]] <<- sdk
  g <- gkpSisters(inp)
  cbind(scenario = s, rates = RATES, ratesSummary(inp), demokinAtAge(sdk, EGO_AGE),
        GKP_sisters = unname(g["GKP_sisters"]), r = unname(g["r"]))
})
dk <- do.call(rbind, rows)


# ==== 2. Save ====

utils::write.csv(dk, file.path(RESULTS_DIR, "demokin_counts.csv"), row.names = FALSE)
saveRDS(summaries, file.path(RESULTS_DIR, "demokin_kin_summary.rds"))
# >>> rates saved for the later steps, start
saveRDS(inputs, file.path(RESULTS_DIR, "demokin_inputs.rds"))
# <<< end
print(dk[, c("scenario", "TFR", "e0", "dk_sisters", "dk_cousins", "GKP_sisters")], digits = 3)
