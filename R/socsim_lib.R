# Socsim library: read a scenario folder, prepare the initial population, run Socsim
# through rsocsim, read the results and select the egos.
#
# A scenario folder <SCENARIO_DIR>/<name> contains:
#   <name>.sup         supervisory file (bint, hetfert, duration, rate files)
#   rates/...          fertility, mortality and marriage rates named in the .sup
#   presim.opop(.gz)   initial population (built from the rates if absent)
# Socsim writes its results to <name>/sim_results_<seed>_ (rsocsim up to 1.7) or
# <name>/sim_results_<name>.sup_<seed>_ (rsocsim 1.9 and later).
#
# Main functions:
#   runSocsim("TFR=4_HM", seed = 600)            run one scenario (skipped if done)
#   runSocsimScenarios(c("TFR=4_HM", ...), seeds) several, each in its own R process
#   readSocsimResults(findResultDir(...))        result.opop and result.omar
#   selectEgos(opop)                             egos of the paper (section 4.2)
#   readSocsimRates(folder)                      rates as applied by Socsim

library(data.table)


# ==== 1. Reading a scenario folder ====

# Settings of a .sup file
readSup <- function(folder) {
  sup <- file.path(folder, paste0(basename(folder), ".sup"))
  if (!file.exists(sup)) stop(sup, " does not exist")
  txt <- readLines(sup, warn = FALSE)
  txt <- txt[!grepl("^\\s*\\*", txt)]                      # drop comment lines
  val <- function(key) {
    x <- sub(paste0("^\\s*", key, "\\s+"), "",
             grep(paste0("^\\s*", key, "\\s"), txt, value = TRUE))
    trimws(sub("\\s.*$", "", x))
  }
  list(supFile  = basename(sup),
       bint     = as.numeric(tail(c("0", val("bint")), 1)),
       hetfert  = as.numeric(tail(c("0", val("hetfert")), 1)),
       duration = as.numeric(val("duration")),              # one value per segment
       includes = unique(val("include")))
}

# One block of a Socsim rate file as a vector of monthly rates by month of age, read
# as Socsim reads it: each line "y m r" gives the rate r up to the upper bound of y
# years and m months, from the previous bound. Months after the last bound get NA.
readRateBlock <- function(file, event = "birth", sex = "F", state = "single",
                          nMonths = 1332) {
  txt  <- readLines(file, warn = FALSE)
  txt  <- trimws(txt[!grepl("^\\s*\\*", txt) & nzchar(trimws(txt))])
  head <- grep(paste0("^", event, "\\s+\\d+\\s+", sex, "\\s+", state), txt)
  if (length(head) == 0) stop("no '", event, " ", sex, " ", state, "' block in ", file)
  out  <- rep(NA_real_, nMonths)
  prev <- 0
  for (ln in txt[-(1:head[1])]) {
    if (!grepl("^\\d", ln)) break                           # next block
    v <- as.numeric(strsplit(ln, "\\s+")[[1]][1:3])
    b <- min(12 * v[1] + v[2], nMonths)
    if (b > prev) out[(prev + 1):b] <- v[3]
    prev <- max(prev, b)
  }
  out
}

# Annual fertility (births per woman at each age 0..maxAge) and survival
# probabilities by sex, as applied by Socsim, from the rate files named in the .sup
readSocsimRates <- function(folder, maxAge = 110) {
  s      <- readSup(folder)
  fFile  <- file.path(folder, grep("fert", s$includes, ignore.case = TRUE, value = TRUE)[1])
  mFile  <- file.path(folder, grep("mort", s$includes, ignore.case = TRUE, value = TRUE)[1])
  n      <- 12 * (maxAge + 1)
  byYear <- function(v, f) sapply(0:maxAge, function(a) f(v[(12 * a + 1):(12 * a + 12)]))
  fm <- readRateBlock(fFile, "birth", "F", "single", n)
  fm[is.na(fm)] <- 0
  fx <- byYear(fm, sum)
  px <- function(sex) {
    q <- readRateBlock(mFile, "death", sex, "single", n)
    q[is.na(q)] <- 1          # no rate after the last line: close the life table
    byYear(q, function(x) prod(1 - pmin(x, 1)))
  }
  pxF <- px("F")
  pxM <- px("M")
  e0  <- function(p) sum(c(1, cumprod(p))[1:(maxAge + 1)] * (1 + c(p, 0)[1:(maxAge + 1)]) / 2)
  list(fx = fx, px_F = pxF, px_M = pxM, TFR = sum(fx), e0_F = e0(pxF), e0_M = e0(pxM),
       fertFile = fFile, mortFile = mFile)
}


# ==== 2. Initial population ====

# Stable population by single year of age for both sexes, by projection of a
# stationary population with the rates (fx: births per woman; px: survival)
stablePopulation <- function(fx, px_w, px_m, years = 200, initialPop = 20000,
                             pFemale = 0.488) {
  nAges <- length(px_w)
  w <- pFemale * initialPop * c(1, cumprod(px_w[-nAges]))
  m <- (1 - pFemale) * initialPop * c(1, cumprod(px_m[-nAges]))
  for (y in seq_len(years)) {
    b <- sum(w * fx, na.rm = TRUE)
    w <- c(b * pFemale, w[-nAges] * px_w[-nAges])
    m <- c(b * (1 - pFemale), m[-nAges] * px_m[-nAges])
  }
  data.frame(age = 0:(nAges - 1), women = w, men = m)
}

# Socsim initial population (opop format) of about initPop persons with the age and
# sex structure of a stable population; dates of birth in months, ego of age a born
# in the year yBirth - a
socsimStablePop <- function(pop_women, pop_men, initPop = 50000, yBirth = 100) {
  tot <- sum(pop_women) + sum(pop_men)
  nW  <- round(pop_women * initPop / tot)
  nM  <- round(pop_men * initPop / tot)
  ages <- seq_along(nW) - 1
  one <- function(n, fem) {
    a <- rep(ages, n)
    data.frame(fem = fem, dob = (yBirth - a + 1) * 12 + sample(1:12, length(a), replace = TRUE))
  }
  x <- rbind(one(nW, 1), one(nM, 0))
  x <- x[order(-x$dob), ]
  opop <- data.frame(pid = seq_len(nrow(x)), fem = x$fem, group = 1, nev = 0, dob = x$dob,
                     mom = 0, pop = 0, nesibm = 0, nesibp = 0, lborn = 0, marid = 0,
                     mstat = 0, dod = 0, fmult = 0)
  opop
}

# Size of the initial population such that the population alive at the end of the
# run is about targetEnd persons (growth rate from the net reproduction rate)
autoInitialSize <- function(rates, years, targetEnd = 3e5, minSize = 1000, maxSize = 1e5) {
  lx  <- c(1, cumprod(rates$px_F))[seq_along(rates$fx)]
  nrr <- 0.488 * sum(lx * rates$fx)
  Tg  <- sum((seq_along(rates$fx) - 0.5) * lx * rates$fx) / sum(lx * rates$fx)
  r   <- log(nrr) / Tg
  round(min(max(targetEnd * exp(-r * years), minSize), maxSize))
}

# Make sure that presim.opop exists: use the existing file, decompress presim.opop.gz,
# copy presim_from_<scenario>.opop, or build a stable population from the rates.
# presim.omar (no marriages at the start) is always written.
prepareInitialPop <- function(folder, rates = readSocsimRates(folder), years = NULL,
                              initSize = NULL, seed = 1) {
  opopFile <- file.path(folder, "presim.opop")
  gzFile   <- paste0(opopFile, ".gz")
  sib      <- list.files(folder, pattern = "^presim_from_.*\\.opop$", full.names = TRUE)
  if (file.exists(opopFile) && file.size(opopFile) > 0) {
    message("Initial population: existing presim.opop")
  } else if (file.exists(gzFile)) {
    con <- gzfile(gzFile, "rt")
    writeLines(readLines(con), opopFile)
    close(con)
    message("Initial population: presim.opop.gz decompressed")
  } else if (length(sib) == 1) {
    file.copy(sib, opopFile)
    message("Initial population copied from ", basename(sib))
  } else {
    if (is.null(years)) years <- sum(readSup(folder)$duration) / 12
    if (is.null(initSize)) initSize <- autoInitialSize(rates, years)
    set.seed(seed)
    st   <- stablePopulation(rates$fx, rates$px_F, rates$px_M, initialPop = 5e6)
    opop <- socsimStablePop(st$women, st$men, initPop = initSize, yBirth = 100)
    utils::write.table(opop, opopFile, row.names = FALSE, col.names = FALSE)
    message("Initial population: stable population of ", nrow(opop), " persons written")
  }
  writeLines(character(0), file.path(folder, "presim.omar"))
  invisible(opopFile)
}


# ==== 3. Results ====

# Folder of the results of a run with a given seed, or NA. An empty result.opop (a run
# that did not write its results) counts as no results.
findResultDir <- function(folder, supFile = paste0(basename(folder), ".sup"), seed) {
  cand <- file.path(folder, c(paste0("sim_results_", supFile, "_", seed, "_"),
                              paste0("sim_results_", seed, "_")))
  f  <- file.path(cand, "result.opop")
  ok <- cand[file.exists(f) & file.size(f) > 0]
  if (length(ok) == 0) NA_character_ else ok[1]
}

# result.opop (one line per person) and result.omar (one line per marriage)
readSocsimResults <- function(resDir, omar = TRUE) {
  opop <- data.table::fread(file.path(resDir, "result.opop"), header = FALSE)
  data.table::setnames(opop, c("pid", "fem", "group", "nev", "dob", "mom", "pop",
                               "nesibm", "nesibp", "lborn", "marid", "mstat", "dod",
                               "fmult"))
  omarFile <- file.path(resDir, "result.omar")
  om <- data.table::data.table()
  if (omar && file.exists(omarFile) && file.size(omarFile) > 0) {
    om <- data.table::fread(omarFile, header = FALSE)
    data.table::setnames(om, c("mid", "wpid", "hpid", "dstart", "dend", "rend",
                               "wprior", "hprior"))
  }
  list(opop = as.data.frame(opop), omar = as.data.frame(om))
}

# Results of a scenario for a given seed (error if the scenario has not been run)
readScenario <- function(scenario, seed, scenarioDir = SCENARIO_DIR, omar = TRUE) {
  folder <- file.path(scenarioDir, scenario)
  resDir <- findResultDir(folder, seed = seed)
  if (is.na(resDir)) stop(scenario, ": no results for seed ", seed, " (run Socsim first)")
  c(readSocsimResults(resDir, omar = omar), list(resultDir = resDir))
}


# ==== 4. Running Socsim ====

# Run one scenario. If results for this seed already exist the run is skipped, unless
# overwrite = TRUE. Socsim runs in a local working folder (workDir): Socsim creates an
# empty result.opop at the start and writes it at the end, which cloud-synchronised
# folders (Google Drive, Dropbox) do not always handle. The results are copied back
# to the scenario folder at the end, and their size is checked.
# processMethod: "callr" (default when callr is installed: Socsim runs in a separate
# R process) or "inprocess" (safe for one run per R session only).
runSocsim <- function(scenario, seed = 600, scenarioDir = SCENARIO_DIR, overwrite = FALSE,
                      workDir = file.path(tempdir(), "socsim_work"), processMethod = NULL,
                      initSize = NULL) {
  if (is.null(processMethod)) {
    processMethod <- if (requireNamespace("callr", quietly = TRUE)) "callr" else "inprocess"
  }
  folder <- normalizePath(file.path(path.expand(scenarioDir), scenario), mustWork = TRUE)
  sup    <- readSup(folder)
  rates  <- readSocsimRates(folder)
  message(sprintf("%s: bint %g, hetfert %g, %g months; TFR %.2f, e0 women %.1f",
                  scenario, sup$bint, sup$hetfert, sum(sup$duration), rates$TFR, rates$e0_F))
  resDir <- findResultDir(folder, sup$supFile, seed)
  if (!is.na(resDir) && !overwrite) {
    message(scenario, ": results for seed ", seed, " already exist in ", basename(resDir))
    return(invisible(resDir))
  }
  prepareInitialPop(folder, rates, sum(sup$duration) / 12, initSize)

  # --- copy the inputs to the local working folder ---
  runDir <- file.path(workDir, scenario)
  unlink(runDir, recursive = TRUE)
  dir.create(runDir, recursive = TRUE)
  for (f in c(sup$supFile, sup$includes, "presim.opop", "presim.omar")) {
    dir.create(dirname(file.path(runDir, f)), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(folder, f), file.path(runDir, f)))
      stop("cannot copy ", f, " to ", runDir)
  }
  runDir <- normalizePath(runDir)

  # --- run ---
  message(scenario, ": Socsim started ", format(Sys.time()))
  if (processMethod == "callr") {
    callr::r(function(folder, supFile, seed) {
      rsocsim::socsim(folder, supFile, seed, process_method = "inprocess")
    }, args = list(runDir, sup$supFile, as.character(seed)), show = TRUE)
  } else {
    rsocsim::socsim(runDir, sup$supFile, as.character(seed), process_method = processMethod)
  }
  message(scenario, ": Socsim ended ", format(Sys.time()))
  runRes <- findResultDir(runDir, sup$supFile, seed)
  if (is.na(runRes)) stop("Socsim produced no result.opop in ", runDir)

  # --- copy the results back to the scenario folder ---
  dest <- file.path(folder, basename(runRes))
  unlink(dest, recursive = TRUE)
  if (!file.copy(runRes, folder, recursive = TRUE)) stop("cannot copy ", runRes, " to ", folder)
  if (file.size(file.path(dest, "result.opop")) != file.size(file.path(runRes, "result.opop")))
    stop("result.opop was not copied completely to ", dest)
  unlink(runDir, recursive = TRUE)
  message(scenario, ": results in ", dest)
  invisible(dest)
}

# Apply a function of this repository to several scenarios, each in its own R process
# when callr is installed (separate = TRUE): the memory used by a large population is
# freed after each scenario, and a failure (for instance out of memory) stops only
# that scenario, not the R session. fun: name of a function whose first two arguments
# are the scenario and the seed; seeds: one seed, or one per scenario.
forEachScenario <- function(scenarios, seeds, fun, ...,
                            separate = requireNamespace("callr", quietly = TRUE)) {
  seeds  <- rep_len(seeds, length(scenarios))
  status <- character(0)
  for (i in seq_along(scenarios)) {
    s <- scenarios[i]
    r <- tryCatch({
      if (separate) {
        callr::r(function(repo, fun, s, seed, args) {
          REPO_DIR <<- repo
          source(file.path(repo, "config.R"))
          do.call(fun, c(list(s, seed), args))
          invisible(NULL)
        }, args = list(REPO_DIR, fun, s, seeds[i], list(...)), show = TRUE)
      } else {
        do.call(fun, c(list(s, seeds[i]), list(...)))
      }
      "ok"
    }, error = function(e) sub("(?s).*\n! ", "", conditionMessage(e), perl = TRUE))
    message("==== ", s, ": ", r)
    status[s] <- r
  }
  message("\n==== Summary")
  for (s in names(status)) message(sprintf("  %-32s %s", s, substr(status[s], 1, 100)))
  invisible(status)
}

# Run Socsim for several scenarios (see forEachScenario)
runSocsimScenarios <- function(scenarios, seeds = 600, ...) {
  forEachScenario(scenarios, seeds, "runSocsim", ...)
}

# New scenario folder that copies the rates and the initial population of an existing
# one, with a .sup that differs only by bint, hetfert and (optionally) duration
newScenarioFolder <- function(scenario, from, bint = 0, hetfert = 0, duration = NULL,
                              scenarioDir = SCENARIO_DIR) {
  src <- file.path(path.expand(scenarioDir), from)
  dst <- file.path(path.expand(scenarioDir), scenario)
  if (dir.exists(dst)) stop(dst, " already exists")
  s <- readSup(src)
  if (is.null(duration)) duration <- sum(s$duration)
  dir.create(file.path(dst, "rates"), recursive = TRUE)
  for (f in s$includes) file.copy(file.path(src, f), file.path(dst, f))
  for (f in c("presim.opop", "presim.opop.gz")) {
    if (file.exists(file.path(src, f))) file.copy(file.path(src, f), file.path(dst, f))
  }
  writeLines(c(
    "*************************************************************************",
    paste0("* Supervisory file for ", scenario),
    paste0("* Rates and initial population copied from ", from, " on ", Sys.Date()),
    "*************************************************************************",
    "segments  1", "input_file presim", "",
    "random_father 1", "marriage_after_childbirth 0", "marriage_eval distribution",
    "marriage_queues 1", paste("hetfert", hetfert), paste("bint", bint), "",
    paste("duration", duration),
    paste("include", s$includes),
    "run"), file.path(dst, paste0(scenario, ".sup")))
  message("Created ", dst)
  invisible(dst)
}


# ==== 5. Choice of egos ====

# Egos (section 4.2 of the paper): women (fem = 1) born in the simulation whose mother
# and maternal grandmother were also born in it (mom > 0 at each step: the initial
# population has mom = 0), so that the fertility histories of both are complete, who
# are at least minAgeAtEnd years old at the end of the run, and born in the last
# windowYears years for which these conditions are met. If there are more than maxEgos, a
# random sample is drawn (seed sampleSeed). Returns the pids, with the selection in
# attribute "info".
selectEgos <- function(opop, minAgeAtEnd = 70, windowYears = 60, egoFem = 1,
                       maxEgos = Inf, sampleSeed = 1, minEgos = 5000) {
  endMonth <- max(opop$dob)
  momOf    <- function(p) { m <- opop$mom[match(p, opop$pid)]; m[is.na(m)] <- 0; m }
  e  <- opop[opop$fem == egoFem & opop$mom > 0 & opop$dob <= endMonth - 12 * minAgeAtEnd, ]
  gm <- momOf(e$mom)                         # maternal grandmother (0 if mother not born in run)
  e  <- e[gm > 0 & momOf(gm) > 0, ]          # mother and grandmother born in the run
  if (nrow(e) == 0) stop("no ego meets the conditions: the run is too short")
  last  <- max(e$dob)
  first <- last - 12 * windowYears
  e <- e[e$dob >= first, ]
  info <- list(endMonth = endMonth, firstDob = min(e$dob), lastDob = last,
               nEligible = nrow(e), minAgeAtEnd = minAgeAtEnd, windowYears = windowYears)
  pid <- e$pid
  if (length(pid) > maxEgos) {
    set.seed(sampleSeed)
    pid <- sort(sample(pid, maxEgos))
  }
  info$nEgos <- length(pid)
  message(sprintf("Egos: %d eligible, %d kept, born in months %d to %d (end of run: month %d)",
                  info$nEligible, info$nEgos, info$firstDob, info$lastDob, endMonth))
  if (info$nEligible < minEgos)
    warning("only ", info$nEligible,
            " eligible egos: consider a longer run or a larger initial population")
  attr(pid, "info") <- info
  pid
}
