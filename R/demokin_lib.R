# DemoKin library: rates given to DemoKin, runs of the one-sex time-invariant model,
# kin ever born at a given age of ego, and the two formulas of the paper.
#
#   inp <- demokinInputs(read.csv("data/empirical_rates/TFR=4_HM.csv"))
#   s   <- runDemoKin(inp)                  kin_summary with count_everborn
#   demokinAtAge(s, 65)                     sisters, aunts, cousins, ... at age 65
#   gkpSisters(inp)                         formula [1] of GKP (all sisters ever born)
#   pearsonSisters(F_dau, var_dau)          formula [2]: G + var / G - 1
#
# The rates of the paper are the female rates estimated from each Socsim run over its
# last 100 years (empiricalRates() in kin_count_lib.R; files in data/empirical_rates),
# with the proportion of female births observed in the run. inputsFromScenario() gives
# instead the rates written in the Socsim rate files of a scenario.

suppressPackageStartupMessages(library(DemoKin))
suppressPackageStartupMessages(library(dplyr))


# ==== 1. Rates ====

# From a table of exposure, deaths, births and births of daughters by age (columns
# age, expo, deaths, births, bdau): survival probabilities px = exp(-mx), closed at the
# last age (px = 0), fertility fx (births of both sexes per woman), and the proportion
# of female births pf
demokinInputs <- function(r, maxAge = 100) {
  r  <- r[r$age <= maxAge, ]
  mx <- ifelse(r$expo > 0, r$deaths / r$expo, 1)
  px <- exp(-mx)
  px[length(px)] <- 0
  fx <- ifelse(r$expo > 0, r$births / r$expo, 0)
  pf <- sum(r$bdau) / sum(r$births)
  list(age = r$age, mx = mx, px = px, fx = fx, pf = pf)
}

# From the rate files of a Socsim scenario (rates as applied by Socsim, see
# readSocsimRates in socsim_lib.R), with a given proportion of female births
inputsFromScenario <- function(folder, maxAge = 100, pf = 0.488) {
  r  <- readSocsimRates(folder, maxAge = maxAge)
  px <- r$px_F
  px[length(px)] <- 0
  list(age = 0:maxAge, mx = -log(pmax(r$px_F, 1e-300)), px = px, fx = r$fx, pf = pf)
}

# Summary measures of the rates: TFR, GRR, NRR and female life expectancy at birth
ratesSummary <- function(inp) {
  lx <- c(1, cumprod(inp$px))[seq_along(inp$px)]
# >>> NRR with the person-years lived at each age, as in Table A1 and lambdaMoments(), start
  Lx <- (lx + c(lx[-1], 0)) / 2
  data.frame(TFR = sum(inp$fx), pf = inp$pf, GRR = sum(inp$fx) * inp$pf,
             NRR = sum(inp$fx * inp$pf * Lx), e0 = sum(lx) - 0.5)
# <<< end
}


# ==== 2. DemoKin ====

# One-sex, time-invariant model (DemoKin::kin). Kin ever born at exact age x of ego:
# count_living(x) + count_cum_dead(x - 1), because DemoKin records at age x the
# deaths that occur during the following year. This counts the kin alive at the birth
# of ego or born after it (the definition of kin_count_lib.R).
runDemoKin <- function(inp, output_kin = c("m", "gm", "d", "os", "ys", "oa", "ya",
                                           "coa", "cya", "nos", "nys")) {
  k <- DemoKin::kin(p = inp$px, f = inp$fx, time_invariant = TRUE,
                    birth_female = inp$pf, output_kin = output_kin)
  k$kin_summary %>%
    dplyr::arrange(kin, age_focal) %>%
    dplyr::group_by(kin) %>%
    dplyr::mutate(count_everborn = count_living + dplyr::lag(count_cum_dead, default = 0)) %>%
    dplyr::ungroup()
}

# Kin ever born (and living) at a given age of ego, older and younger kin together
demokinAtAge <- function(s, age = 65) {
  g <- function(kins, v) sum(s[[v]][s$kin %in% kins & s$age_focal == age])
  data.frame(dk_sisters  = g(c("os", "ys"), "count_everborn"),
             dk_aunts    = g(c("oa", "ya"), "count_everborn"),
             dk_cousins  = g(c("coa", "cya"), "count_everborn"),
             dk_nieces   = g(c("nos", "nys"), "count_everborn"),
             dk_daughters = g("d", "count_everborn"),
             dk_sisters_living = g(c("os", "ys"), "count_living"),
             dk_cousins_living = g(c("coa", "cya"), "count_living"))
}


# ==== 3. Formulas [1] and [2] ====

# Formula [1] of Goodman, Keyfitz and Pullum (1974): expected number of sisters ever
# born of a woman, averaged over the ages of her mother at her birth in the stable
# population. Older sisters: all daughters born to the mother before age x; younger
# sisters: daughters born after x, weighted by the survival of the mother from x.
# Ages are single years, with half of the daughters of the year of ego's birth
# counted as older and half as younger.
gkpSisters <- function(inp) {
  fx   <- inp$fx * inp$pf                      # daughters per woman
  mx   <- inp$mx
  A    <- length(fx)
  lmid <- exp(-(cumsum(mx) - mx / 2))          # survival to the middle of each age
  ages <- 0:(A - 1) + 0.5
  r    <- stats::uniroot(function(r) sum(exp(-r * ages) * lmid * fx) - 1, c(-0.2, 0.2))$root
  W    <- exp(-r * ages) * lmid * fx           # ages of mothers at the birth of ego
  W    <- W / sum(W)
  S1   <- 0
  for (x in seq_len(A)) {
    older   <- sum(fx[seq_len(x - 1)]) + fx[x] / 2
    younger <- if (x < A) sum((lmid[(x + 1):A] / lmid[x]) * fx[(x + 1):A]) else 0
    S1 <- S1 + W[x] * (older + younger + fx[x] / 2)
  }
  c(r = r, GKP_sisters = S1)
}

# Formula [2] (Pearson): mean number of sisters of a woman, from the mean G and the
# variance of the number of daughters of the women of the mothers' generation
pearsonSisters <- function(F_dau, var_dau) F_dau + var_dau / F_dau - 1

# >>> moments of the expected number of daughters by age at death, start
# Appendix B of the paper. Lambda(t) is the expected number of daughters born to a woman
# by age t (zero before the first age of fertility). For a woman of the mothers'
# generation who dies at age T, with the life table of inp:
#   mu   = E[Lambda(T)], the net reproduction rate,
#   tau2 = Var[Lambda(T)], the part of the variance of the number of daughters due to
#          the differences in age at death,
#   S0   = mu + tau2 / mu, the number of sisters of formula [1] with r = 0.
# Ages are single years: a woman who dies during age x is given half of the daughters
# of that age.
lambdaMoments <- function(inp) {
  fx  <- inp$fx * inp$pf                       # daughters per woman
  A   <- length(fx)
  l   <- cumprod(c(1, inp$px))                 # survival at exact ages 0..A (l = 0 at A)
  Lam <- c(0, cumsum(fx))                      # expected daughters by exact ages 0..A
  d   <- -diff(l)                              # deaths during each age
  LT  <- (Lam[-1] + Lam[-(A + 1)]) / 2         # expected daughters of women dying in it
  mu  <- sum(d * LT)
  tau2 <- sum(d * LT^2) - mu^2
  c(mu = mu, tau2 = tau2, S0 = mu + tau2 / mu)
}
# <<< end
