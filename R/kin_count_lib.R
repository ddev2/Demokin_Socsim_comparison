# Counts of kin in Socsim populations, with the definition of kin ever born of
# DemoKin, and female rates estimated from the runs.
#
# Definition (section 4.3 and Appendix of the paper). In DemoKin, kin ever born at
# age x of ego are the kin alive at that exact age plus the kin who died between the
# birth of ego and that age: the kin alive at the birth of ego or born after it. Kin
# who died before ego was born (older sisters who died in childhood, for instance) are
# left out. Kin are counted through the female line, as in the one-sex model:
#   sisters           other daughters of ego's mother
#   maternal aunts    other daughters of the maternal grandmother
#   maternal cousins  daughters of these aunts
#   nieces            daughters of ego's sisters
#
# 1. countKinSocsim(opop, egos)    kin of all egos at age EGO_AGE, under the DemoKin
#                                  definition and in all, with the mean and variance of
#                                  the numbers of children and daughters of the women of
#                                  the mothers' generation (used in Tables 1 and A1)
# 2. empiricalRates(opop)          female exposure, deaths, births and births of
#                                  daughters by single year of age over the last years
#                                  of the run (the rates given to DemoKin)
# 3. socsimKinScenario(scenario)   1 and 2 for one scenario, written to the results
# 4. genealogiesToDemoKin(dfKin)   numbers of kin by age of ego in the format of
#                                  DemoKin (kin_full, kin_summary), from the genealogies
#                                  of genealogy_lib.R (all ages of ego, slower)

library(data.table)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


# ==== 1. Kin of all egos at a given age ====

# opop:  result.opop of a run (data frame)
# egos:  pids of the egos (default: selectEgos with the rule of the paper)
# Returns one row: the mean numbers of kin per ego, "_dk" with the DemoKin definition
# and "_all" for all kin ever born, born before ego reaches egoAge.
countKinSocsim <- function(opop, egos = NULL, egoAge = 65, minAgeAtEnd = 70,
                           windowYears = 60) {
  if (is.null(egos)) egos <- selectEgos(opop, minAgeAtEnd, windowYears)
  d <- data.table::as.data.table(opop)[, list(pid, fem, dob, mom, dod)]
  data.table::setkey(d, pid)

  # numbers of children and daughters of each person
  kids <- d[mom > 0]
  nchT  <- kids[, .N, by = mom]
  ndauT <- kids[fem == 1, .N, by = mom]
  d[, nch := 0L]
  d[nchT, on = list(pid = mom), nch := i.N]
  d[, ndau := 0L]
  d[ndauT, on = list(pid = mom), ndau := i.N]
  pMale <- mean(kids$fem == 0)

  e  <- d[list(egos)]                         # egos, mothers, grandmothers
  m  <- d[list(e$mom)]
  E  <- data.table::data.table(ego = e$pid, edob = e$dob, mother = e$mom, gm = m$mom)
  n  <- nrow(E)
  fd <- d[fem == 1 & mom > 0, list(pid, mom, dob, dod)]   # women born in the run
  limit   <- 12 * egoAge
  aliveAt <- function(dod, edob) dod == 0 | dod > edob    # alive at ego's birth or born after
  dk      <- function(x) sum(aliveAt(x$dod, x$edob)) / n

  # sisters: daughters of ego's mother other than ego
  S <- merge(E, fd, by.x = "mother", by.y = "mom")
  S <- S[pid != ego & dob < edob + limit]
  # maternal aunts: daughters of the grandmother other than ego's mother
  A <- merge(E, fd, by.x = "gm", by.y = "mom")
  A <- A[pid != mother & dob < edob + limit]
  # maternal cousins: daughters of the aunts, born before ego reaches egoAge
  C <- merge(A[, list(ego, edob, aunt = pid)], fd, by.x = "aunt", by.y = "mom")
  C <- C[dob < edob + limit]
  # nieces: daughters of the sisters, born before ego reaches egoAge
  N <- merge(S[, list(ego, edob, sister = pid)], fd, by.x = "sister", by.y = "mom")
  N <- N[dob < edob + limit]

  # women of the mothers' generation: born in the period of birth of the mothers
  w   <- d[fem == 1 & mom > 0 & dob >= min(m$dob) & dob <= max(m$dob)]
  w50 <- w[dod == 0 | dod - dob >= 600]                   # survivors to age 50
  mv  <- function(x) c(mean(x), mean((x - mean(x))^2))    # mean and variance

  data.frame(
    n_ego = n, ego_dob_first = min(e$dob), ego_dob_last = max(e$dob),
    end_month = max(d$dob), p_male = pMale,
    sisters_all = nrow(S) / n,           sisters_dk = dk(S),
    older_sisters_all = sum(S$dob < S$edob) / n,
    older_sisters_dk  = sum(S$dob < S$edob & aliveAt(S$dod, S$edob)) / n,
    aunts_all   = nrow(A) / n,           aunts_dk   = dk(A),
    cousins_all = nrow(C) / n,           cousins_dk = dk(C),
    nieces_all  = nrow(N) / n,           nieces_dk  = dk(N),
    daughters_ego = mean(e$ndau),
    F_ch    = mv(w$nch)[1],    var_ch    = mv(w$nch)[2],
    F_dau   = mv(w$ndau)[1],   var_dau   = mv(w$ndau)[2],
    F50_ch  = mv(w50$nch)[1],  var50_ch  = mv(w50$nch)[2],
    F50_dau = mv(w50$ndau)[1], var50_dau = mv(w50$ndau)[2],
    alive50 = nrow(w50) / nrow(w))
}


# ==== 2. Female rates estimated from a run ====

# Female exposure (person-years), deaths, births (both sexes) and births of daughters
# by single year of age of the woman, over the last `years` years of the run. Women
# alive at the end of the run (dod = 0) are exposed until the end. Births are counted
# by the age of the mother at the birth.
empiricalRates <- function(opop, years = 100, maxAge = 110) {
  d   <- data.table::as.data.table(opop)[, list(pid, fem, dob, mom, dod)]
  P1  <- max(d$dob)
  P0  <- P1 - 12 * years
  w   <- d[fem == 1]
  dob <- as.numeric(w$dob)
  dod <- as.numeric(ifelse(w$dod == 0, P1, pmin(w$dod, P1)))
  A   <- maxAge + 1
  expo <- vapply(0:maxAge, function(a) {
    s <- pmax(dob + 12 * a, P0)
    e <- pmin(dob + 12 * (a + 1), dod, P1)
    sum(pmax(e - s, 0)) / 12
  }, numeric(1))
  dd <- w[dod > 0 & dod >= P0 & dod < P1]
  ad <- pmin((dd$dod - dd$dob) %/% 12, maxAge)
  deaths <- tabulate(ad + 1, A)
  k    <- d[mom > 0 & dob >= P0 & dob < P1]
  mdob <- d$dob[match(k$mom, d$pid)]
  am   <- (k$dob - mdob) %/% 12
  ok   <- am >= 0 & am <= maxAge
  births <- tabulate(am[ok] + 1, A)
  bdau   <- tabulate(am[ok & k$fem == 1] + 1, A)
  data.frame(age = 0:maxAge, expo = expo, deaths = deaths, births = births, bdau = bdau)
}


# ==== 3. Processing of one scenario ====

# Everything that is computed from the results of one Socsim run, written to
# resultsDir:
#   egos/<scenario>.rds              pids of the egos and the selection rule
#   socsim_counts/<scenario>.csv     countKinSocsim() at age egoAge
#   empirical_rates/<scenario>.csv   empiricalRates() over the last rateYears years
# and, with genealogies = TRUE (slow with large populations; sample of maxEgosGenea egos):
#   genealogies/<scenario>_genealogies.rds   socsimToGenealogies()
#   genealogies/<scenario>_kin_by_age.rds    genealogiesToDemoKin()
# A scenario already processed is skipped unless overwrite = TRUE.
socsimKinScenario <- function(scenario, seed, scenarioDir = SCENARIO_DIR,
                              resultsDir = RESULTS_DIR, egoAge = EGO_AGE,
                              minAgeAtEnd = MIN_AGE_AT_END, windowYears = WINDOW_YEARS,
                              rateYears = RATE_YEARS, genealogies = FALSE,
                              maxEgosGenea = MAX_EGOS_GENEA, overwrite = FALSE) {
  out <- function(sub, ext, create = TRUE) {
    if (create) dir.create(file.path(resultsDir, sub), showWarnings = FALSE, recursive = TRUE)
    file.path(resultsDir, sub, paste0(scenario, ext))
  }
  fCount <- out("socsim_counts", ".csv", create = FALSE)
  fGenea <- out("genealogies", "_genealogies.rds", create = FALSE)
# >>> a scenario is reprocessed if its rates are missing or its seed has changed, start
  fRates  <- out("empirical_rates", ".csv", create = FALSE)
  sameSeed <- file.exists(fCount) &&
    isTRUE(utils::read.csv(fCount, nrows = 1)$seed == seed)
  doCount <- overwrite || !file.exists(fCount) || !file.exists(fRates) || !sameSeed
# <<< end
  doGenea <- genealogies && (overwrite || !file.exists(fGenea))
  if (!doCount && !doGenea) {
    message(scenario, ": already processed")
    return(invisible(NULL))
  }
  res  <- readScenario(scenario, seed, scenarioDir, omar = doGenea)
  message(scenario, ": ", nrow(res$opop), " persons in ", basename(res$resultDir))
  egos <- selectEgos(res$opop, minAgeAtEnd, windowYears)
  if (doCount) {
    saveRDS(egos, out("egos", ".rds"))
    k <- countKinSocsim(res$opop, egos, egoAge)
    utils::write.csv(cbind(scenario = scenario, seed = seed, k), out("socsim_counts", ".csv"),
                     row.names = FALSE)
    er <- empiricalRates(res$opop, rateYears)
    utils::write.csv(er, out("empirical_rates", ".csv"), row.names = FALSE)
    message(scenario, sprintf(": sisters %.3f, cousins %.3f (DemoKin definition)",
                              k$sisters_dk, k$cousins_dk))
  }
  if (doGenea) {
    eg    <- selectEgos(res$opop, minAgeAtEnd, windowYears, maxEgos = maxEgosGenea)
    kin   <- retrieveKin(res$opop, res$omar, as.vector(eg))
    dfKin <- socsimToGenealogies(res$opop, as.vector(eg), kin)
    saveRDS(dfKin, out("genealogies", "_genealogies.rds"))
    saveRDS(genealogiesToDemoKin(dfKin), out("genealogies", "_kin_by_age.rds"))
    message(scenario, ": genealogies of ", length(eg), " egos saved")
  }
  invisible(NULL)
}

# All the csv files of a folder in one data frame
collectCsv <- function(dir) {
  f <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  do.call(rbind, lapply(f, utils::read.csv, stringsAsFactors = FALSE))
}


# ==== 4. Kin by age of ego, in DemoKin format, from genealogies ====

# Numbers of kin of female egos (female line) by age of ego, from the genealogies of
# socsimToGenealogies(), in the format of DemoKin::kin(): kin_full (by age of kin and
# of ego) and kin_summary (count_living, count_dead, count_cum_dead, count_everborn
# by age of ego). Kin codes are those of DemoKin: m, gm, ggm, d, gd, ggd, os, ys, oa,
# ya, coa, cya, nos, nys. Counts are divided by the number of egos alive at each age.
genealogiesToDemoKin <- function(dfKin, aSex = "F", nCols = 100) {
  df <- createAndAddDemoKinFull(dfKin, aSex = aSex, nCols = nCols)
  kin_full <- createDemoKinFull(df, divideByEgo = TRUE, changeKinNames = TRUE, nCols = nCols)
  createDemoKinSummary(kin_full)
}

# --- Functions of the KINFERT pipeline (KinFertLib.R and demokin_functions.R) ---

# first step for creating a DemoKin dataframe from a KINFERT genealogy file
processGenea <- function(df, aSex) {
  selectDownKin <- function (df_tmp, kinTypeSel, aSex) {
# >>> female line through mothers, male line through fathers, start
    if (aSex == "F") {
# <<< end
      df_tmp$select <- ifelse((df_tmp$kinType==kinTypeSel)&df_tmp[match(df_tmp$idMother, df_tmp$id),]$sex==aSex,TRUE,df_tmp$select) 
    } else {
      df_tmp$select <- ifelse((df_tmp$kinType==kinTypeSel)&df_tmp[match(df_tmp$idFather, df_tmp$id),]$sex==aSex,TRUE,df_tmp$select) 
    }
    df_tmp <- subset(df_tmp, !((df_tmp$kinType==kinTypeSel)&(is.na(df_tmp$select))))
    return (df_tmp)
  }
  
  df <- df[,intersect(c("Family", "id", "kinOf", "kinType", "sex", "idEgo", "idFather", "idMother", "yBirthFloat", "yDeathFloat", "nChildren", "ageRelEgo", "censored"), names(df))]
  #select the correct kin set
  kinTypesBoth <- c("aunt-uncle", "child", "ego",
                    "father", "first cousin", "grand child",
                    "grand father", "grand mother",
                    "great grand child", "great grand father", "great grand mother", "mother",
                    "niece-nephew", "sibling")
  kinTypesFemale <- c("aunt-uncle", "child", "ego",
                      "first cousin", "grand child",
                      "grand mother",
                      "great grand child", "great grand mother", "mother",
                      "niece-nephew", "sibling")
  kinTypesMale <- c("aunt-uncle", "child", "ego",
                    "father", "first cousin", "grand child",
                    "grand father",
                    "great grand child", "great grand father",
                    "niece-nephew", "sibling")
  
  if (length (aSex)>1) {
    df <- subset(df, df$kinType %in% kinTypesBoth)
  } else {
    if (aSex == "F") {
      df <- subset(df, df$kinType %in% kinTypesFemale)
    } else {
      df <- subset(df, df$kinType %in% kinTypesMale)
    }
    df$kinType <- factor(df$kinType)
    #select families with ego of the selected sex
    df <- df %>%
      dplyr::group_by(Family) %>%
      dplyr::filter((any((kinType == "ego") & (sex %in% aSex)))) %>%
      dplyr::ungroup()
    #select kin of the correct sex
    df <- subset(df, sex %in% aSex)
    df$select <- NA
    df$select <- ifelse(df$kinType=="ego",TRUE,df$select)
    df$select <- ifelse(df$kinType=="mother",TRUE,df$select)
    df$select <- ifelse(df$kinType=="child",TRUE,df$select)
    
    df <- selectDownKin (df, "grand child", aSex) #only grandchildren who are children of ego's children of the correct sex
    df <- selectDownKin (df, "great grand child", aSex) #only great-grandchildren who are children of ego's grand-children of the correct sex
    #siblings should have the same mother or the same father than ego (here we exclude half-siblings)
    if (aSex == "F") {
      egoMothers <- subset(df, (df$kinType=="ego"))$idMother
      df$select <- ifelse((df$kinType=="sibling")&(df$idMother %in% egoMothers),TRUE,df$select) 
      rm(egoMothers)
    } else {
      egoFathers <- subset(df, (df$kinType=="ego"))$idFather
      df$select <- ifelse((df$kinType=="sibling")&(df$idFather %in% egoFathers),TRUE,df$select)
      rm(egoFathers)
    }
    df <- subset(df, !((df$kinType=="sibling")&(is.na(df$select))))
    df <- selectDownKin (df, "niece-nephew", aSex) #only niece-nephew who are children of ego's sibling of the correct sex
    #select grandparents who are parents of either mother or father
    if (aSex == "F") {
      egoGrandMothers <- subset(df, (df$kinType=="mother"))$idMother
      df$select <- ifelse((df$kinType=="grand mother")&(df$id %in% egoGrandMothers),TRUE,df$select) 
      df$select <- ifelse((df$kinType=="aunt-uncle")&(df$idMother %in% egoGrandMothers),TRUE,df$select) 
      rm(egoGrandMothers)
      df <- subset(df, !((df$kinType=="grand mother")&(is.na(df$select))))
      df <- subset(df, !((df$kinType=="aunt-uncle")&(is.na(df$select))))
    } else {
      egoGrandFathers <- subset(df, (df$kinType=="father"))$idFather
      df$select <- ifelse((df$kinType=="grand father")&(df$id %in% egoGrandFathers),TRUE,df$select) 
      df$select <- ifelse((df$kinType=="aunt-uncle")&(df$idFather %in% egoGrandFathers),TRUE,df$select) 
      rm(egoGrandFathers)
      df <- subset(df, !((df$kinType=="grand father")&(is.na(df$select))))
      df <- subset(df, !((df$kinType=="aunt-uncle")&(is.na(df$select))))
    }
    df <- selectDownKin (df, "first cousin", aSex) #only first cousin who are children of ego's sibling of the correct sex
    if (aSex == "F") {
      egoGreatGrandMothers <- subset(df, (df$kinType=="grand mother"))$idMother
      df$select <- ifelse((df$kinType=="great grand mother")&(df$id %in% egoGreatGrandMothers),TRUE,df$select) 
      rm(egoGreatGrandMothers)
      df <- subset(df, !((df$kinType=="great grand mother")&(is.na(df$select))))
    } else {
      egoGreatGrandFathers <- subset(df, (df$kinType=="grand father"))$idFather
      df$select <- ifelse((df$kinType=="great grand father")&(df$id %in% egoGreatGrandFathers),TRUE,df$select) 
      rm(egoGreatGrandFathers)
      df <- subset(df, !((df$kinType=="great grand father")&(is.na(df$select))))
    }
    df$select <- NULL
  } # (length (aSex)>1)
  
  ## older/younger siblings, aunts-uncles and their respective children
  df$younger <- NA
  #1) siblings
  df$yBirthEgo <- df[match(df$idEgo, df$id),]$yBirthFloat
  df$sexEgo <- df[match(df$idEgo, df$id),]$sex
  df$younger <- ifelse((df$kinType=="sibling")&(df$sex==df$sexEgo)&(df$yBirthFloat>df$yBirthEgo),TRUE,df$younger)
# >>> a sibling born in the same month as ego is counted as older, so that os + ys = s, start
  df$younger <- ifelse((df$kinType=="sibling")&(df$sex==df$sexEgo)&(df$yBirthFloat<=df$yBirthEgo),FALSE,df$younger)
# <<< end
  df$yBirthEgo <- NULL
  df$sexEgo <- NULL
  #2) aunts-uncles
  if ("F" %in% aSex) {
    df$idMotherEgo <- df[match(df$idEgo, df$id),]$idMother
    df$yBirthMotherEgo <- df[match(df$idMotherEgo, df$id),]$yBirthFloat
    df$younger <- ifelse((df$kinType=="aunt-uncle")&(df$sex=="F")&(df$yBirthFloat>df$yBirthMotherEgo),TRUE,df$younger)
# >>> same-month aunts counted as older, start
    df$younger <- ifelse((df$kinType=="aunt-uncle")&(df$sex=="F")&(df$yBirthFloat<=df$yBirthMotherEgo),FALSE,df$younger)
# <<< end
    df$idMotherEgo <- NULL
    df$yBirthMotherEgo <- NULL
  }
  if ("M" %in% aSex) {
    df$idFatherEgo <- df[match(df$idEgo, df$id),]$idFather
    df$yBirthFatherEgo <- df[match(df$idFatherEgo, df$id),]$yBirthFloat
    df$younger <- ifelse((df$kinType=="aunt-uncle")&(df$sex=="M")&(df$yBirthFloat>df$yBirthFatherEgo),TRUE,df$younger)
# >>> same-month uncles counted as older, start
    df$younger <- ifelse((df$kinType=="aunt-uncle")&(df$sex=="M")&(df$yBirthFloat<=df$yBirthFatherEgo),FALSE,df$younger)
# <<< end
    df$idFatherEgo <- NULL
    df$yBirthFatherEgo <- NULL
  }
  #3) nieces-nephews and first cousins
  if ("F" %in% aSex) {
    df$motherYounger <- df[match(df$idMother, df$id),]$younger
    df$younger <- ifelse(df$kinType %in% c("niece-nephew","first cousin"),df$motherYounger,df$younger)
    df$motherYounger <- NULL
  }
  if ("M" %in% aSex) {
    df$fatherYounger <- df[match(df$idFather, df$id),]$younger
    df$younger <- ifelse(df$kinType %in% c("niece-nephew","first cousin"),df$fatherYounger,df$younger)
    df$fatherYounger <- NULL
  }
  
  return (df)
}

# >>> status of kin computed from the dates at each exact age of ego, start
# Status of each kin at each exact age a = 0..nCols of ego, computed from the dates
# (in years, with month precision) rather than from ages truncated to whole years.
# At exact age a of ego (date t = birth of ego + a), if ego is alive at t:
#   living    the kin is born and not yet dead at t
#   everBorn  the kin is born by t and was alive at the birth of ego or born after it
#             (the definition of kin ever born of DemoKin)
#   age       completed age of the kin at t, if born by t
#   death     the kin dies between exact ages a and a + 1 of ego (not if censored),
#             as DemoKin records at age a the deaths of the following year
# Returns a list of four matrices with one row per kin and nCols + 1 columns.
kinStatusByEgoAge <- function(df, nCols = 100) {
  a   <- 0:nCols
  t   <- outer(df$yBirthEgo, a, "+")
  egoAlive <- t < df$yDeathEgo
  born     <- df$yBirthFloat <= t
  living   <- born & (t < df$yDeathFloat) & egoAlive
  everBorn <- born & (df$yDeathFloat > df$yBirthEgo) & egoAlive
  age      <- floor(t - df$yBirthFloat)
  age[!(everBorn | living)] <- NA
  cens     <- if ("censored" %in% names(df)) as.logical(df$censored) else rep(FALSE, nrow(df))
  cens[is.na(cens)] <- FALSE
  aDeath   <- floor(df$yDeathFloat - df$yBirthEgo)
  death    <- (matrix(aDeath, nrow(df), length(a)) == matrix(a, nrow(df), length(a), byrow = TRUE)) &
    (df$yDeathFloat > df$yBirthEgo) & !cens & egoAlive
  list(living = living * 1, everBorn = everBorn * 1, age = age, death = death * 1)
}

# The four row-wise functions below are no longer used by createAndAddDemoKinFull(),
# which calls kinStatusByEgoAge(); they are kept for reference.
# <<< end
#whether kin is alive during ego lifetime, for how many years and from what age of ego
row_wise_kin_alive <- function(row, names_egoAgeRange, nCols=100) {
  # Current values of the nCols+ columns
  current_vals1 <- as.numeric(row[names_egoAgeRange])
  if (!as.logical(row["aliveDuringEgoLifeTime"])) return (current_vals1)
  ageRelEgo <- as.numeric(row["ageRelEgo"])
  lengthLife <- as.numeric(row["lengthLife"])
  egoLengthLife <- min(nCols, trunc (as.numeric(row["egoLengthLife"])))
  
  # ageRelEgo is negative when a kin is younger than ego
  if (ageRelEgo < 0) {
    start_pos <- - ageRelEgo # ... therefore we convert to positive
    count <- min(egoLengthLife + ageRelEgo, lengthLife) + 1
  } else {
    start_pos <- 0
    count <- min(egoLengthLife, lengthLife - ageRelEgo) + 1
  }
  
  # Set values to 1 for the specified range
  if (!is.na(start_pos) && !is.na(count) && count > 0) {
    end_pos <- min(start_pos + count - 1, nCols)  # Don't exceed column nCols
    if (start_pos >= 0 && start_pos <= nCols) {
      cols_to_change <- (start_pos:end_pos) + 1  # +1 because R is 1-indexed
      current_vals1[cols_to_change] <- 1
    }
  }
  
  return(as.numeric(current_vals1))
}

#whether kin is ever-born during ego lifetime, and from what age of ego
row_wise_kin_everBorn <- function(row, names_egoAgeRange_EB, nCols=100) {
  # Current values of the nCols+ columns
  current_vals1 <- as.numeric(row[names_egoAgeRange_EB])
  if (!as.logical(row["aliveDuringEgoLifeTime"])) return (current_vals1)
  ageRelEgo <- as.numeric(row["ageRelEgo"])
  lengthLife <- as.numeric(row["lengthLife"])
  egoLengthLife <- min(nCols, trunc (as.numeric(row["egoLengthLife"])))
  
  # ageRelEgo is negative when a kin is younger than ego
  if (ageRelEgo < 0) {
    start_pos <- - ageRelEgo # ... therefore we convert to positive
    count <- egoLengthLife + ageRelEgo + 1
  } else {
    start_pos <- 0
    count <- egoLengthLife + 1
  }
  
  # Set values to 1 for the specified range
  if (!is.na(start_pos) && !is.na(count) && count > 0) {
    end_pos <- min(start_pos + count - 1, nCols)  # Don't exceed column nCols
    if (start_pos >= 0 && start_pos <= nCols) {
      cols_to_change <- (start_pos:end_pos) + 1  # +1 because R is 1-indexed
      current_vals1[cols_to_change] <- 1
    }
  }
  
  return (as.numeric(current_vals1))
}

#ages of kin during ego lifetime
row_wise_kin_age <- function(row, names_kinAges, nCols=100, kinAlive=FALSE) {
  if (!as.logical(row["aliveDuringEgoLifeTime"])) return (as.numeric(row[names_kinAges]))
  ageRelEgo <- as.numeric(row["ageRelEgo"])
  lengthLife <- as.numeric(row["lengthLife"])
  egoLengthLife <- min(nCols, trunc (as.numeric(row["egoLengthLife"])))
  
  start_pos <- ageRelEgo
  if (start_pos < 0) {
    first_age <- 0
    if (kinAlive) {
      count <- min (egoLengthLife + start_pos, lengthLife) + 1
    } else {
      count <- egoLengthLife + ageRelEgo + 1
    }
    start_pos <- - start_pos
  } else {
    first_age <- start_pos
    if (kinAlive) {
      count <- min (egoLengthLife, lengthLife - start_pos) + 1
    } else {
      count <- egoLengthLife + 1
    }
    start_pos <- 0
  }
  
  # Current values of the nCols+1 columns
  current_vals2 <- as.numeric(row[names_kinAges])
  
  # Set values to ages of kin for the specified range
  if (!is.na(start_pos) && !is.na(count) && count > 0) {
    end_pos <- min(start_pos + count - 1, nCols)  # Don't exceed column nCols
    if (start_pos >= 0 && start_pos <= nCols) {
      cols_to_change <- (start_pos:end_pos) + 1  # +1 because R is 1-indexed
      current_vals2[cols_to_change] <- (0:(end_pos-start_pos)) + first_age
    }
  }
  
  return(as.numeric(current_vals2))
}

#death of kin during ego lifetime
row_wise_kin_death <- function(row, names_kinDeath_egoAgeRange, nCols=100) {
  if (!as.logical(row["aliveDuringEgoLifeTime"])) return (as.numeric(row[names_kinDeath_egoAgeRange]))
  if (("censored" %in% names(row)) && isTRUE(as.logical(trimws(row["censored"])))) return (as.numeric(row[names_kinDeath_egoAgeRange]))
  ageRelEgo <- as.numeric(row["ageRelEgo"])
  lengthLife <- as.numeric(row["lengthLife"])
  egoLengthLife <- min(nCols, trunc (as.numeric(row["egoLengthLife"])))
  
  age_death_rel <- min(nCols, min (nCols, lengthLife) - ageRelEgo)
  
  # Current values of the nCols+1 columns
  current_vals3 <- as.numeric(row[names_kinDeath_egoAgeRange])
  
  # Set values to 1 for the specified range
  if ((age_death_rel >= 0) & (age_death_rel <= egoLengthLife)) {
    current_vals3[age_death_rel + 1] <- 1
  }
  
  return(as.numeric(current_vals3))
}

# Function to duplicate rows with younger or older kinType
younger_older_rows <- function(df) {
  
  # Create the modified dataframe
  df_modified <- df %>%
    # Add row numbers to maintain original order
    dplyr::mutate(original_row = dplyr::row_number())
  # Create duplicated rows for non-NA 'younger' values
  # Create the duplicated rows separately
  duplicated_rows <- df_modified %>%
    dplyr::filter(!is.na(younger)) %>%
    dplyr::mutate(kinType = paste0(ifelse(younger, "y_", "o_"), kinType))
  
  # Combine them
  df_modified <- dplyr::bind_rows(df_modified, duplicated_rows) %>%
    dplyr::arrange(original_row, is.na(younger)) %>%
    dplyr::select(-original_row)
  
  return (df_modified)
}

createDemoKinLevels <- function(kinTypes) {
  levelsKin <- levels(kinTypes)
  named_vector <- c(
    "ego" = "ego",
    "o_first cousin" = "coa",
    "y_first cousin" = "cya",
    "first cousin" = "c",
    "child" = "d",
    "grand child" = "gd",
    "great grand child" = "ggd",
    "great grand father" = "ggm",
    "great grand mother" = "ggm",
    "grand father" = "gm",
    "grand mother" = "gm",
    "father" = "m",
    "mother" = "m",
    "o_niece-nephew" = "nos",
    "y_niece-nephew" = "nys",
    "niece-nephew" = "n",
    "o_aunt-uncle" = "oa",
    "y_aunt-uncle" = "ya",
    "aunt-uncle" = "a",
    "o_sibling" = "os",
    "y_sibling" = "ys",
    "sibling" = "s"
  )
  
  return (named_vector[levelsKin])
}

createDemoKinFull <- function(df=df_DemoKin, aCohort=NULL, divideByEgo=TRUE, changeKinNames=TRUE, nCols=100) {
  
  # Column names for living during ego lifetime (0-nCols)
  names_egoAgeRange <- as.character(0:nCols)
  # Column names for ever-born during ego lifetime (0-nCols)
  names_egoAgeRange_EB <- paste0("EB_", 0:nCols)
  # Column names for kin ages
  names_kinAges <- paste0("age_", 0:nCols)
  # Column names for death of kin by ego ages
  names_kinDeath_egoAgeRange <- paste0("death_", 0:nCols)
  
  # 0) sum egos by age alive and divide all cols by this number
  if (divideByEgo) {
    df_egos <- subset(df, kinType=="ego")
    sumEgos <- colSums(df_egos[names_egoAgeRange])
    rm(df_egos)
    # remove 0 from sumEgos
    sumEgos[sumEgos == 0] <- 1e-30
    # divide kin counts by egos numbers
    df[names_egoAgeRange] <- sweep(df[names_egoAgeRange], 2, sumEgos, FUN = "/")
    df[names_egoAgeRange_EB] <- sweep(df[names_egoAgeRange_EB], 2, sumEgos, FUN = "/")
    df[names_kinDeath_egoAgeRange] <- sweep(df[names_kinDeath_egoAgeRange], 2, sumEgos, FUN = "/")
    df <- subset(df, kinType!="ego")
  }
  
  # 1) ego age and alive status to long format
  df_long <- df [,c("kinType", "sex", names_egoAgeRange)] %>%
    # Pivot ego age columns
    tidyr::pivot_longer(
      cols = all_of(names_egoAgeRange),
      names_to = "ego_age",
      values_to = "isAlive"
    )
  # Convert ego_age to numeric
  df_long <- dplyr::mutate(df_long, ego_age = as.numeric(ego_age))
  
  # 2) ego age and ever born status to long format
  df_long_everBorn <- df [,c("kinType", "sex", names_egoAgeRange_EB)] %>%
    # Pivot ego age columns
    tidyr::pivot_longer(
      cols = all_of(names_egoAgeRange_EB),
      names_to = "ego_age",
      values_to = "everBorn"
    )
  df_long$everBorn <- df_long_everBorn$everBorn
  rm(df_long_everBorn)
  
  # 3) ego age and death status to long format
  df_long_death <- df [,c("kinType", "sex", names_kinDeath_egoAgeRange)] %>%
    # Pivot ego age columns
    tidyr::pivot_longer(
      cols = all_of(names_kinDeath_egoAgeRange),
      names_to = "kin_death",
      values_to = "death"
    )
  df_long$death <- df_long_death$death
  rm(df_long_death)
  
  # 4) kin age to long format
  df_age_kin <- df [,c("kinType", "sex", names_kinAges)] %>%
    # Pivot ego age columns
    tidyr::pivot_longer(
      cols = all_of(names_kinAges),
      names_to = "kin_ages",
      values_to = "age_kin"
    )
  df_long$age_kin <- df_age_kin$age_kin
  rm(df_age_kin)
  
  df_long <- df_long[,c("kinType", "sex", "age_kin", "ego_age", "isAlive", "everBorn", "death")]
  
  df_full <- df_long %>%
    dplyr::group_by(kinType, sex, age_kin, ego_age) %>%
    dplyr::summarise(
      count_alive = sum(isAlive, na.rm = TRUE),
      count_everBorn = sum(everBorn, na.rm = TRUE),
      count_death = sum(death, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Complete all combinations and fill missing with 0
    tidyr::complete(kinType, sex, age_kin, ego_age, 
                    fill = list(count_alive = 0, count_everBorn = 0, count_death = 0))
  
  if (changeKinNames) {
    df_full$kinType <- factor(df_full$kinType)
    demoKinLevels <- createDemoKinLevels (df_full$kinType)
    levels(df_full$kinType) <- demoKinLevels
    df_full$kinType <- factor(df_full$kinType)
    df_full$kinType <- as.character(df_full$kinType)
  }  
  df_full <- subset(df_full,(!(is.na(age_kin))))
  changeSex <- c("F"="f", "M"="m")
  df_full$sex <- changeSex[df_full$sex]
  names(df_full) <- c("kin", "sex_kin", "age_kin", "age_focal", "living", "everBorn", "dead")
  
  if (is.null(aCohort)) {
    df_full$year <- NA
    df_full$cohort <- NA
  } else {
    df_full$year <- aCohort + df_full$age_focal
    df_full$cohort <- aCohort
  }
  
  return (df_full)
}

# take a KinFert dataframe and first step of creating a DemoKin one:
# columns of ages to keep tab of whether ego and/or kin are alive, dead or ever-born
createAndAddDemoKinFull <- function (df, aSex=c("F"), nCols=100) {
  #aSex=c("F"), c("M"), c"F","M"): one of those three, the last option meaning both sexes
  
  df <- processGenea(df, aSex)
  
  ## length of life during ego life
  df$lengthLife <- trunc (df$yDeathFloat-df$yBirthFloat)
  df$egoLengthLife <- df[match(df$idEgo, df$id),]$lengthLife
  df$yBirthEgo <- df[match(df$idEgo, df$id),]$yBirthFloat
  df$yDeathEgo <- df[match(df$idEgo, df$id),]$yDeathFloat
  df$ageRelEgo <- trunc (df$yBirthEgo - df$yBirthFloat)
  # kin alive during ego's lifetime (ever-born): kin contemporary of ego
  df$aliveDuringEgoLifeTime <- ((df$yBirthFloat >= df$yBirthEgo) & (df$yBirthFloat < df$yDeathEgo)) | # kin born during ego's lifetime
    ((df$yDeathFloat > df$yBirthEgo) & (df$yDeathFloat <= df$yDeathEgo)) | # kin deceased during ego's lifetime
    ((df$yBirthFloat >= df$yBirthEgo) & (df$yDeathFloat <= df$yDeathEgo)) | # kin born and died during ego's lifetime
    ((df$yBirthFloat < df$yBirthEgo) & (df$yDeathFloat > df$yDeathEgo)) # kin born before and deceased after ego's lifetime
  
  # Create nCols columns filled with 0 for alive/dead_notBorn status during ego lifetime
  names_egoAgeRange <- as.character(0:nCols) # ego age range
  new_cols <- setNames(data.frame(matrix(0, nrow = nrow(df), ncol = nCols+1)), 
                       names_egoAgeRange)
  # Add to existing dataframe
  df <- cbind(df, new_cols)
  
  # Create nCols columns filled with 0 for ever-born status during ego lifetime
  names_egoAgeRange_EB <- paste0("EB_", 0:nCols)
  new_cols <- setNames(data.frame(matrix(0, nrow = nrow(df), ncol = nCols+1)), 
                       names_egoAgeRange_EB)
  # Add to existing dataframe
  df <- cbind(df, new_cols)
  
  # Create nCols columns with ages of kin during ego lifetime
  names_kinAges <- paste0("age_", 0:nCols)
  new_cols <- setNames(data.frame(matrix(NA, nrow = nrow(df), ncol = nCols+1)), 
                       names_kinAges)
  # Add to existing dataframe
  df <- cbind(df, new_cols)
  
  # Create nCols columns for kin death during ego lifetime
  names_kinDeath_egoAgeRange <- paste0("death_", 0:nCols)
  new_cols <- setNames(data.frame(matrix(0, nrow = nrow(df), ncol = nCols+1)), 
                       names_kinDeath_egoAgeRange)
  # Add to existing dataframe
  df <- cbind(df, new_cols)
  
# >>> status of kin from the dates, with the requested number of columns, start
  st <- kinStatusByEgoAge(df, nCols)
  df[names_egoAgeRange]          <- st$living
  df[names_egoAgeRange_EB]       <- st$everBorn
  df[names_kinAges]              <- st$age
  df[names_kinDeath_egoAgeRange] <- st$death
# <<< end
  
  # add younger and older kin
  df <- younger_older_rows (df)
  
  return (df)
}

createDemoKinSummary <- function (kin_full, aCohort=NULL) {
  # copied from https://github.com/IvanWilli/DemoKin/blob/main/R/kin.R
  # with some modifications
  agrupar <- c()
  if(!is.null(aCohort)){
    agrupar <- c("year", "cohort")
  }
  agrupar_no_age_focal <- c("kin", agrupar)
  agrupar <- c("age_focal", "kin", agrupar)
  
  kin_summary <- dplyr::bind_rows(
    kin_full %>%
      dplyr::rename(total=living) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(agrupar))) %>%
      dplyr::summarise(count_living = sum(total),
                       mean_age = sum(total*age_kin)/sum(total),
                       sd_age  = (sum(total*age_kin^2)/sum(total)-mean_age^2)^.5) %>%
      tidyr::pivot_longer(count_living:sd_age, names_to = "indicator", values_to = "value"),
    kin_full %>%
      dplyr::rename(total=dead) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(agrupar))) %>%
      dplyr::summarise(count_dead = sum(total)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(agrupar_no_age_focal))) %>%
      dplyr::mutate(count_cum_dead = cumsum(count_dead),
                    mean_age_lost = cumsum(count_dead * age_focal)/cumsum(count_dead)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_longer(count_dead:mean_age_lost, names_to = "indicator", values_to = "value")) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = indicator, values_from = value)
  
  if ("everBorn" %in% names(kin_full)) {
    kin_summary <- kin_summary %>%
      dplyr::left_join(
        kin_full %>%
          dplyr::rename(total=everBorn) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(agrupar))) %>%
          dplyr::summarise(count_ever_born = sum(total)) %>%
          tidyr::pivot_longer(count_ever_born, names_to = "indicator", values_to = "value")
      )
    kin_summary$indicator <- NULL
    names(kin_summary)[ncol(kin_summary)] <- "count_everborn"
  }
  
  # return
  kin_out <- list(kin_full = kin_full, kin_summary = kin_summary)
  
  return (kin_out)
}
