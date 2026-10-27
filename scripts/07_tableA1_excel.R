# >>> new script, Table A1 of the paper in an Excel file, start
# Step 7. Table A1 of the paper (all 32 runs) in an Excel file.
#
# Gathers, for each scenario:
#   - the description of the run: scenarios/scenarios.csv;
#   - the Socsim results of step 2: results/socsim_counts.csv, or the counts of the runs
#     of the paper in data/socsim_counts.csv if step 2 has not been run;
#   - DemoKin run here with the female rates of each Socsim run (data/empirical_rates or
#     results/empirical_rates), as in step 3 and in the scripts of demokin/ (step 6);
# and writes the table in the layout of the paper.
#
# Columns of Table A1:
#   e0, TFR        female life expectancy at birth and total fertility rate of the rates
#                  given to DemoKin (estimated from the last 100 years of the run)
# >>> NRR column, start
#   NRR            net reproduction rate of the same rates
# <<< end
#   bint, hetfert  minimum birth interval (months) and heterogeneity in fecundability
#   V/F ch.        variance to mean ratio of the number of children ever born of the
#                  women of the mothers' generation who survive to age 50
#   s2/mu dau.     variance to mean ratio of the number of daughters of all the women of
#                  the mothers' generation
# >>> kappa column, start
#   kappa          (s2 - tau2) / mu: the same ratio once the part due to the mortality of
#                  mothers (tau2, Appendix B) is removed; it gives the error of formula [4]
# <<< end
#   Sisters, Maternal cousins: kin ever born (alive at the birth of ego or born after
#                  it) of a woman aged 65, female line, in DemoKin and Socsim, and the
#                  ratio Socsim / DemoKin
#
# >>> file name distinct from results/sisters_formulas.csv of step 4, start
# Output: results/paper_TableA1.xlsx (and results/paper_TableA1.csv). Needs DemoKin and
# openxlsx.
# <<< end


# ==== 0. Setup ====

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

suppressPackageStartupMessages(library(DemoKin))
library(openxlsx)

MAX_AGE     <- 100
EGO_AGE     <- 65
# >>> results folder as in config.R (environment variable KIN_RESULTS_DIR), start
RESULTS_DIR <- Sys.getenv("KIN_RESULTS_DIR", file.path(REPO_DIR, "results"))
# <<< end
dir.create(RESULTS_DIR, showWarnings = FALSE)

# Label of the populations in the paper. The scenario folders and rate objects keep
# their original names (KOREA_2022...): the Korean schedule is labelled 2018 in the paper.
POP_LABEL <- c("TFR 3" = "TFR 3", "TFR 4" = "TFR 4", "TFR 5" = "TFR 5",
               "Korea 2022" = "Korea 2018", "Korea 2018" = "Korea 2018", "USA 1933" = "USA 1933")


# ==== 1. Socsim results ====

SCENARIOS <- utils::read.csv(file.path(REPO_DIR, "scenarios", "scenarios.csv"),
                             stringsAsFactors = FALSE)
fss <- file.path(RESULTS_DIR, "socsim_counts.csv")
if (!file.exists(fss) || file.size(fss) < 10) fss <- file.path(REPO_DIR, "data", "socsim_counts.csv")
message("Socsim counts read from ", fss)
ss <- utils::read.csv(fss, stringsAsFactors = FALSE)


# ==== 2. DemoKin with the rates of each run ====

ratesFile <- function(s) {
  f <- file.path(RESULTS_DIR, "empirical_rates", paste0(s, ".csv"))
  if (!file.exists(f)) f <- file.path(REPO_DIR, "data", "empirical_rates", paste0(s, ".csv"))
  f
}

demokinScenario <- function(s) {
  r  <- utils::read.csv(ratesFile(s))
  r  <- r[r$age <= MAX_AGE, ]
  mx <- ifelse(r$expo > 0, r$deaths / r$expo, 1)
  px <- exp(-mx)
  px[length(px)] <- 0
  fx <- ifelse(r$expo > 0, r$births / r$expo, 0)
  pf <- sum(r$bdau) / sum(r$births)
  k  <- DemoKin::kin(p = px, f = fx, time_invariant = TRUE, birth_female = pf,
                     output_kin = c("os", "ys", "coa", "cya"))
  s2 <- k$kin_summary
  s2 <- s2[order(s2$kin, s2$age_focal), ]
  # kin ever born at exact age x: count_living(x) + count_cum_dead(x - 1)
  s2$count_everborn <- s2$count_living +
    stats::ave(s2$count_cum_dead, s2$kin, FUN = function(z) c(0, z[-length(z)]))
  eb <- function(kins) sum(s2$count_everborn[s2$kin %in% kins & s2$age_focal == EGO_AGE])
  lx <- cumprod(c(1, px))[seq_along(px)]
# >>> NRR column, start
  Lx <- (lx + c(lx[-1], 0)) / 2            # person-years lived at each age
# >>> kappa column, start
  # tau2: variance of the expected number of daughters by age at death (lambdaMoments in R/demokin_lib.R)
  lE  <- c(lx, 0)                            # survival at exact ages 0..A, closed at A
  Lam <- c(0, cumsum(fx * pf))
  dth <- -diff(lE)
  LT  <- (Lam[-1] + Lam[-length(Lam)]) / 2
  tau2 <- sum(dth * LT^2) - sum(dth * LT)^2
# <<< end
  data.frame(scenario = s, TFR = sum(fx), e0 = sum(lx) - 0.5, NRR = sum(Lx * fx) * pf,
             tau2 = tau2,
             dk_sisters = eb(c("os", "ys")), dk_cousins = eb(c("coa", "cya")))
}

dk <- do.call(rbind, lapply(SCENARIOS$scenario, demokinScenario))


# ==== 3. Table A1 ====

t <- merge(SCENARIOS, ss[, setdiff(names(ss), "seed")], by = "scenario")
t <- merge(t, dk, by = "scenario")
# >>> every scenario must be present once, start
missing <- setdiff(SCENARIOS$scenario, t$scenario)
if (length(missing) > 0) stop("scenarios missing: ", paste(missing, collapse = ", "))
stopifnot(!anyDuplicated(t$scenario))
# <<< end
t$mortality <- factor(t$mortality, c("low", "intermediate", "high"))
t$population <- factor(POP_LABEL[t$population], c("TFR 3", "TFR 4", "TFR 5", "Korea 2018", "USA 1933"))
t <- t[order(t$mortality, t$population, t$hetfert, t$bint), ]

tab <- data.frame(
  Mortality       = as.character(t$mortality),
  Scenario        = as.character(t$population),
  e0              = round(t$e0),
  TFR             = round(t$TFR, 2),
# >>> NRR column, start
  NRR             = round(t$NRR, 2),
# <<< end
  bint            = t$bint,
  hetfert         = ifelse(t$hetfert == 1, "yes", "no"),
  VF_ch           = round(t$var50_ch / t$F50_ch, 2),
  s2mu_dau        = round(t$var_dau / t$F_dau, 2),
# >>> kappa column, start
  kappa           = round((t$var_dau - t$tau2) / t$F_dau, 2),
# <<< end
  Sisters_DemoKin = round(t$dk_sisters, 2),
  Sisters_Socsim  = round(t$sisters_dk, 2),
  Sisters_ratio   = round(t$sisters_dk / t$dk_sisters, 2),
  Cousins_DemoKin = round(t$dk_cousins, 2),
  Cousins_Socsim  = round(t$cousins_dk, 2),
  Cousins_ratio   = round(t$cousins_dk / t$dk_cousins, 2),
  stringsAsFactors = FALSE)
utils::write.csv(tab, file.path(RESULTS_DIR, "paper_TableA1.csv"), row.names = FALSE)


# ==== 4. Excel file in the layout of the paper ====

wb <- createWorkbook()
addWorksheet(wb, "Table A1")
title <- paste("Table A1. Expected number of sisters and maternal cousins ever born of a",
               "woman aged 65, female line, in DemoKin and Socsim, by mortality level,",
               "minimum birth interval (bint) and heterogeneity (hetfert); all 32 runs")
writeData(wb, 1, title, startRow = 1, startCol = 1)

# two header rows
# >>> NRR column, start
# >>> kappa column, start
h1 <- c("Scenario", "e0", "TFR", "NRR", "bint", "hetfert", "σ²/mean", "", "", "Sisters", "", "",
        "Maternal cousins", "", "")
h2 <- c("", "", "", "", "", "", "ch.", "dau.", "κ", "DemoKin", "Socsim", "Ratio",
        "DemoKin", "Socsim", "Ratio")
# <<< end
# <<< end
writeData(wb, 1, t(h1), startRow = 3, colNames = FALSE)
writeData(wb, 1, t(h2), startRow = 4, colNames = FALSE)
# >>> NRR column, start
# >>> kappa column, start
mergeCells(wb, 1, cols = 7:9, rows = 3)
mergeCells(wb, 1, cols = 10:12, rows = 3)
mergeCells(wb, 1, cols = 13:15, rows = 3)
# <<< end
# <<< end

# body, with a row for each mortality level
row <- 5
groupRows <- integer(0)
for (m in c("low", "intermediate", "high")) {
  x <- tab[tab$Mortality == m, -1]
  if (nrow(x) == 0) next
  writeData(wb, 1, paste0(toupper(substring(m, 1, 1)), substring(m, 2), " mortality"),
            startRow = row, startCol = 1)
  groupRows <- c(groupRows, row)
  writeData(wb, 1, x, startRow = row + 1, colNames = FALSE)
  row <- row + 1 + nrow(x)
}
note <- paste("Note: e0: life expectancy at birth of women; TFR: total fertility rate;",
              "both from the rates given to DemoKin, estimated from the last 100 years of",
              "each run. NRR: net reproduction rate of the same rates; at low mortality it",
              "is close to the number of sisters in DemoKin. σ²/mean ch.: variance to mean ratio of the number of children ever",
              "born of women of the mothers' generation who survive to age 50. dau.: the",
              "same ratio for the number of daughters of all women of the mothers'",
              "generation. κ: (σ² − τ²)/μ, the same ratio once the part due to the mortality",
              "of mothers is removed (formula [4], Appendix B). Kin ever born: kin alive at the birth of ego or born after it.",
              "Cousins are counted up to age 65 of ego. Ratio: Socsim divided by DemoKin.")
writeData(wb, 1, note, startRow = row + 1, startCol = 1)

# styles
bold   <- createStyle(textDecoration = "bold")
head   <- createStyle(textDecoration = "bold", halign = "center", border = "bottom")
dec2   <- createStyle(numFmt = "0.00", halign = "right")
dec0   <- createStyle(numFmt = "0", halign = "right")
centre <- createStyle(halign = "center")
body   <- 5:(row - 1)
addStyle(wb, 1, bold, rows = 1, cols = 1)
# >>> NRR column, start
addStyle(wb, 1, head, rows = 3:4, cols = 1:15, gridExpand = TRUE)   # 15 columns
addStyle(wb, 1, dec0, rows = body, cols = 2, gridExpand = TRUE)
addStyle(wb, 1, dec2, rows = body, cols = c(3:4, 7:15), gridExpand = TRUE)
addStyle(wb, 1, centre, rows = body, cols = 5:6, gridExpand = TRUE)
# <<< end
addStyle(wb, 1, bold, rows = groupRows, cols = 1, gridExpand = TRUE)
setColWidths(wb, 1, cols = 1, widths = 20)
setColWidths(wb, 1, cols = 2:15, widths = 9)   # 15 columns with NRR and kappa

f <- file.path(RESULTS_DIR, "paper_TableA1.xlsx")
saveWorkbook(wb, f, overwrite = TRUE)
message("Table A1 written to ", f)
print(tab[, -1], row.names = FALSE)
# <<< end
