# Step 2. Kin of the egos in the Socsim runs, and female rates of the runs.
#
# For each scenario (see socsimKinScenario in R/kin_count_lib.R):
#   results/egos/<scenario>.rds              egos (rule of section 4.2 of the paper)
#   results/socsim_counts/<scenario>.csv     kin at age EGO_AGE, DemoKin definition and
#                                            all kin ever born; mean and variance of the
#                                            numbers of children and daughters
#   results/empirical_rates/<scenario>.csv   female rates over the last RATE_YEARS years
# and, if GENEALOGIES is TRUE, the genealogies of a sample of egos and their kin by age
# of ego in DemoKin format (results/genealogies/). Then all the counts are gathered in
# results/socsim_counts.csv.

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

# Full genealogies and kin by age of ego (slow: about an hour per run of several
# million persons); not needed for the tables and figures of the paper
GENEALOGIES <- FALSE
# TRUE to process again the scenarios already processed
OVERWRITE <- FALSE


# ==== 1. Scenarios ====

todo  <- SCENARIOS$scenario
seeds <- SCENARIOS$seed[match(todo, SCENARIOS$scenario)]


# ==== 2. Kin counts and rates ====

status <- forEachScenario(todo, seeds, "socsimKinScenario",
                          genealogies = GENEALOGIES, overwrite = OVERWRITE)


# ==== 3. All scenarios in one file ====

counts <- collectCsv(file.path(RESULTS_DIR, "socsim_counts"))
# no file is written before at least one scenario has been processed (step 4 then uses
# data/socsim_counts.csv)
if (!is.null(counts) && nrow(counts) > 0) {
  utils::write.csv(counts, file.path(RESULTS_DIR, "socsim_counts.csv"), row.names = FALSE)
  message("results/socsim_counts.csv: ", nrow(counts), " scenarios")
} else {
  message("No scenario processed yet: results/socsim_counts.csv not written")
}
