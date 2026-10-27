# Settings shared by all the scripts in scripts/. Each script finds the folder of the
# repository (REPO_DIR), then sources this file, which sets the paths and the
# parameters of the comparison and loads the four libraries in R/.
#
# To use Socsim results stored elsewhere (for instance the original project folder),
# set the environment variable before running a script:
#   Sys.setenv(KIN_SCENARIO_DIR = "/path/to/data/ socsim")
# or edit SCENARIO_DIR below.


# ==== 1. Paths ====

if (!exists("REPO_DIR")) REPO_DIR <- getwd()

# Scenario folders: one folder per Socsim scenario, with its .sup file, its rate files
# and its initial population. The Socsim results are written into these folders.
SCENARIO_DIR <- Sys.getenv("KIN_SCENARIO_DIR", file.path(REPO_DIR, "scenarios"))

# Description of the scenarios (population, mortality, bint, hetfert, seed)
SCENARIO_FILE <- file.path(REPO_DIR, "scenarios", "scenarios.csv")

# Results of the scripts (not under version control, see .gitignore)
RESULTS_DIR <- Sys.getenv("KIN_RESULTS_DIR", file.path(REPO_DIR, "results"))

# Female rates estimated from the Socsim runs used in the paper, so that DemoKin can
# be run without running Socsim first (scripts/03_run_demokin.R)
RATES_DIR <- file.path(REPO_DIR, "data", "empirical_rates")


# ==== 2. Parameters of the comparison ====

EGO_AGE        <- 65      # age of ego at which kin are counted
MIN_AGE_AT_END <- 70      # egos are at least this old at the end of the run
WINDOW_YEARS   <- 60      # egos are born in the last WINDOW_YEARS eligible years
RATE_YEARS     <- 100     # rates given to DemoKin: last RATE_YEARS years of each run
MAX_AGE_DK     <- 100     # last age of the rates given to DemoKin
MAX_EGOS_GENEA <- 20000   # egos sampled for the full genealogies (optional step)


# ==== 3. Libraries ====

for (f in c("socsim_lib.R", "genealogy_lib.R", "kin_count_lib.R", "demokin_lib.R")) {
  source(file.path(REPO_DIR, "R", f))
}

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

SCENARIOS <- utils::read.csv(SCENARIO_FILE, stringsAsFactors = FALSE)
