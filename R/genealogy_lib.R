# Genealogies of egos from a Socsim population.
#
#   kin   <- retrieveKin(opop, omar, egos)      pids of each type of kin of each ego
#   dfKin <- socsimToGenealogies(opop, egos, kin)
#
# retrieveKin follows retrieve_kin of the rsocsim documentation (as adapted in the
# project, with the kin types "unclesaunts", "firstcousins" and "niblings"), and gives
# the same result; the person-by-person searches of the population are replaced by
# vectorised lookups, which makes it much faster with large populations.
#
# socsimToGenealogies returns one data frame with one row per ego and per kin of ego
# (the KinFert format used by kin_count_lib.R):
#   Family       number of the genealogy (one per ego)
#   id           identifier of the person within the file
#   kinType      "ego", "mother", "sibling", "aunt-uncle", "first cousin", ...
#   sex          "F" or "M"
#   idEgo, idFather, idMother   identifiers within the file (0 if not in the genealogy)
#   yBirthFloat, yDeathFloat    dates in years (Socsim months / 12, shifted by yearBase)
#   ageRelEgo    age difference with ego, in completed years
#   censored     TRUE if alive at the end of the run (yDeathFloat = end of the run)


# ==== 1. Kin of each ego ====

retrieveKin <- function(opop, omar = NULL, pid,
                        extra_kintypes = c("unclesaunts", "firstcousins", "niblings")) {
  `%ni%` <- Negate(`%in%`)
  nGenea <- length(pid)
  message("Building ", nGenea, " genealogies")
  zna <- function(x) ifelse(x == 0, NA, x)

  # --- spouses (last marriage) ---
  if (!is.null(omar) && nrow(omar) > 0) {
    opop$spouse <- ifelse(opop$fem == 1, omar[zna(opop$marid), "hpid"],
                          omar[zna(opop$marid), "wpid"])
  } else {
    opop$spouse <- NA
  }

  # --- ancestors ---
  mt <- function(p) match(p, opop$pid)
  opop$FM  <- opop$mom[mt(opop$pop)]
  opop$MM  <- opop$mom[mt(opop$mom)]
  opop$MF  <- opop$pop[mt(opop$mom)]
  opop$FF  <- opop$pop[mt(opop$pop)]
  opop$FFM <- opop$FM[mt(opop$pop)]
  opop$FMM <- opop$MM[mt(opop$pop)]
  opop$FMF <- opop$MF[mt(opop$pop)]
  opop$FFF <- opop$FF[mt(opop$pop)]
  opop$MFM <- opop$FM[mt(opop$mom)]
  opop$MMM <- opop$MM[mt(opop$mom)]
  opop$MMF <- opop$MF[mt(opop$mom)]
  opop$MFF <- opop$FF[mt(opop$mom)]

  # --- children of each person, indexed by pid ---
  allChildren <- with(opop, c(tapply(pid, mom, c), tapply(pid, pop, c)))
  allChildren["0"] <- NULL; allChildren["0"] <- NULL
  KidsOf <- list()
  KidsOf[as.numeric(names(allChildren))] <- allChildren
  rm(allChildren)
  ko <- function(p) lapply(p, function(x) unique(as.vector(unlist(KidsOf[x]))))

  # columns of the egos' rows, as a list with one vector per ego
  rowsOf <- function(cols) {
    m <- as.matrix(opop[mt(pid), cols, drop = FALSE])
    lapply(seq_len(nrow(m)), function(i) as.vector(m[i, ]))
  }

  res <- list()
  res$ggparents <- rowsOf(c("MMM", "MMF", "MFM", "MFF", "FMM", "FMF", "FFM", "FFF"))
  res$gparents  <- rowsOf(c("MM", "MF", "FM", "FF"))
  res$parents   <- rowsOf(c("mom", "pop"))
  if ("unclesaunts" %in% extra_kintypes) {
    g2 <- ko(res$gparents)
    res$unclesaunts <- lapply(seq_along(g2), function(i) g2[[i]][g2[[i]] %ni% res$parents[[i]]])
  }
  res$siblings <- ko(res$parents)
  res$siblings <- lapply(seq_along(res$siblings),
                         function(i) res$siblings[[i]][res$siblings[[i]] %ni% pid[[i]]])
  if ("firstcousins" %in% extra_kintypes) res$firstcousins <- ko(res$unclesaunts)
  res$children   <- ko(pid)
  res$gchildren  <- ko(res$children)
  res$ggchildren <- ko(res$gchildren)
  if ("niblings" %in% extra_kintypes) res$niblings <- ko(res$siblings)
  res$spouse <- as.list(opop$spouse[mt(pid)])
  message(sprintf("Kin found: %s",
                  paste(names(res), sapply(res, function(l) sum(unlist(l) != 0, na.rm = TRUE)),
                        sep = " ", collapse = ", ")))

  # empty sets as NA (as in retrieve_kin); great-grandparents are kept as they are
  for (i in seq_along(res)) {
    res[[i]][sapply(res[[i]], function(x) length(x) == 0)] <- NA
  }
  res$ggparents <- rowsOf(c("MMM", "MMF", "MFM", "MFF", "FMM", "FMF", "FFM", "FFF"))
  res
}


# ==== 2. Genealogies in KinFert format ====

kinRows <- function(kinInfo, tmpDf, kinTypeBySex) {
  n <- nrow(kinInfo)
  result <- tmpDf[rep(1, n), ]
  result$id          <- tmpDf$id + seq_len(n)
  result$kinType     <- kinTypeBySex[1 + kinInfo$fem]
  result$sex         <- ifelse(kinInfo$fem == 1, "F", "M")
  result$yBirth      <- trunc(kinInfo$yBirthFloat)
  result$yBirthFloat <- kinInfo$yBirthFloat
  result$yDeathFloat <- kinInfo$yDeathFloat
  result$ageRelEgo   <- trunc(tmpDf$yBirthFloat - kinInfo$yBirthFloat)
  result$idFather    <- kinInfo$pop
  result$idMother    <- kinInfo$mom
  result$SocSimId    <- kinInfo$pid
  result$censored    <- kinInfo$censored
  result
}

# Socsim pids of parents replaced by identifiers within the genealogy (0 if absent)
updateIdParents <- function(dfFam) {
  lookup <- stats::setNames(dfFam$id, dfFam$SocSimId)
  dfFam$idFather <- lookup[as.character(dfFam$idFather)]
  dfFam$idFather[is.na(dfFam$idFather)] <- 0
  dfFam$idMother <- lookup[as.character(dfFam$idMother)]
  dfFam$idMother[is.na(dfFam$idMother)] <- 0
  dfFam$SocSimId <- NULL
  dfFam
}

socsimToGenealogies <- function(opop, egos_pid, kin, yearBase = 1933,
                                lastDobInitialPop = 1212) {
  # persons alive at the end of the run (dod = 0) are censored at the end of the run
  opop$censored <- (opop$dod == 0)
  opop$dod[opop$censored] <- max(opop$dob) + 1
  opop$yBirthFloat <- yearBase + (opop$dob - lastDobInitialPop - 1) / 12
  opop$yDeathFloat <- yearBase + (opop$dod - lastDobInitialPop - 1) / 12

  # only the persons who appear in the genealogies, one data frame per pid
  needed <- unique(c(egos_pid, unlist(kin, use.names = FALSE)))
  sub    <- opop[opop$pid %in% needed, ]
  lookup <- split(sub, sub$pid)
  rm(sub)

  kinTypeMapping <- list(
    ggparents    = c("M" = "great grand father", "F" = "great grand mother"),
    gparents     = c("M" = "grand father", "F" = "grand mother"),
    parents      = c("M" = "father", "F" = "mother"),
    unclesaunts  = c("M" = "aunt-uncle", "F" = "aunt-uncle"),
    siblings     = c("M" = "sibling", "F" = "sibling"),
    firstcousins = c("M" = "first cousin", "F" = "first cousin"),
    children     = c("M" = "child", "F" = "child"),
    gchildren    = c("M" = "grand child", "F" = "grand child"),
    ggchildren   = c("M" = "great grand child", "F" = "great grand child"),
    niblings     = c("M" = "niece-nephew", "F" = "niece-nephew"),
    spouse       = c("M" = "partner", "F" = "partner"))

  nGenea   <- length(egos_pid)
  kinTypes <- names(kin)
  out      <- vector("list", nGenea)
  idKin    <- 0
  for (fam in seq_len(nGenea)) {
    if (fam %% 500 == 0) message(sprintf("Genealogy %d of %d", fam, nGenea))
    idKin <- idKin + 1
    idEgo <- idKin
    egoInfo <- lookup[[as.character(egos_pid[fam])]]
    if (is.null(egoInfo) || nrow(egoInfo) == 0) {
      warning(sprintf("Ego pid %s not found", egos_pid[fam]))
      next
    }
    egoInfo <- egoInfo[1, ]
    famList <- vector("list", length(kin) + 1)
    k <- 1
    famList[[k]] <- data.frame(
      Family = fam, id = idKin, kinOf = "ego", kinType = "ego",
      sex = ifelse(egoInfo$fem == 1, "F", "M"), idEgo = idEgo,
      idFather = egoInfo$pop, idMother = egoInfo$mom, nUnion = 0, idPartners = "",
      yBirth = trunc(egoInfo$yBirthFloat), yBirthFloat = egoInfo$yBirthFloat,
      yDeathFloat = egoInfo$yDeathFloat, nChildren = 0, ageRelEgo = 0,
      SocSimId = egoInfo$pid, censored = egoInfo$censored, stringsAsFactors = FALSE)
    k <- k + 1
    tmpDf <- famList[[1]]
    nChildren <- 0
    nUnion    <- 0
    for (indKin in seq_along(kin)) {
      kinPid <- kin[[indKin]][fam]
      if (is.na(kinPid)) next
      v <- unlist(kinPid)
      v <- v[!is.na(v) & v != 0]
      if (length(v) == 0) next
      kinInfo <- do.call(rbind, lookup[as.character(v)])
      if (is.null(kinInfo) || nrow(kinInfo) == 0) next
      tmp <- kinRows(kinInfo, tmpDf, kinTypeMapping[[kinTypes[indKin]]])
      famList[[k]] <- tmp
      k <- k + 1
      if (tmp$kinType[1] == "child")   nChildren <- nrow(tmp)
      if (tmp$kinType[1] == "partner") nUnion    <- nrow(tmp)
      tmpDf$id <- tmpDf$id + nrow(tmp)
      idKin    <- tmpDf$id
    }
    dfFam <- do.call(rbind, famList[1:(k - 1)])
    dfFam[1, "nChildren"] <- nChildren
    dfFam[1, "nUnion"]    <- nUnion
    out[[fam]] <- updateIdParents(dfFam)
  }
  dfKin <- do.call(rbind, out)
  rownames(dfKin) <- NULL
  dfKin
}
