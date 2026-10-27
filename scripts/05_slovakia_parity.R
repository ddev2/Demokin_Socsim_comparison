# Step 5. Slovakia 1980: DemoKin without and with parity.
#
# >>> rates derived from the parity data, start
# The same population is run in DemoKin twice:
#   - the age model of Caswell (2019), kin(), with the age-specific fertility rates
#     implied by the rates by age and parity, and the survival probabilities of Slovakia
#     in 1980;
#   - the age x parity model of Caswell (2020), kin_multi_stage(..., parity = TRUE),
#     with the rates by age and parity of the same population (data svk_Uxs, svk_fxs,
#     svk_pxs and svk_Hxs of the DemoKin package).
# The numbers of living kin by age of ego are compared (section 3.4 of the paper). With
# the same fertility by age, both models give the same numbers of kin in the direct line
# and of older sisters and aunts; the parity model gives about 4 per cent fewer younger
# sisters and aunts (4 per cent more with parity = FALSE).
#
# Output:
#   results/slovakia1980_kin.csv            living kin by age of ego, both models
#   results/figures/slovakia1980.pdf (.png) the two models by type of kin
#   results/figures/fig_slovakia.pdf (.png) the figure of the paper (eight kin types)
# <<< end
# >>> benchmark without mortality of section 3.3, start
#   results/slovakia1980_no_mortality.csv   distribution of births, mean and variance of
#                                           daughters, Pearson's formula and the sisters
#                                           of both models, without mortality
# <<< end
# No Socsim run is needed.

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
suppressPackageStartupMessages(library(ggplot2))


# >>> rates derived from the parity data, start
# ==== 1. Rates ====

# >>> source of the Slovak data of DemoKin, start
# Source of the data (checked against the HFD file SVKmi.txt, version of 9 May 2018):
# the rates by age and parity distributed with DemoKin are those of Slovakia in 1980,
# although the DemoKin documentation says 1990. They are the conditional fertility
# rates m_k(x) of the Human Fertility Database turned into probabilities,
# u_k = m_k / (1 + 0.5 m_k) (Caswell 2020, equation 37): svk_Uxs holds the u_k and
# svk_fxs = 0.5 u_k (daughters). Row r of these data holds the HFD rates of age r, so
# that, with row 1 read as age 0, DemoKin places fertility one year younger than the HFD.
# The survival probabilities svk_pxs are 1 - qx of the HMD female life table of 1980
# (e0 = 74.3), with no shift in age (row 1 = age 0).
# <<< end


# Survival probabilities: the same for all parities in svk_pxs
surv <- svk_pxs[, 1]

# Age-specific fertility (daughters per woman) implied by the rates by age and parity:
# the parity distribution of a cohort of women at each age, obtained by applying the
# parity transitions svk_Uxs from parity 0 at birth (mortality does not depend on parity,
# so that it does not change this distribution), weights the fertility of each parity.
# Both models thus have exactly the same fertility by age. (An earlier version used the
# rates of the Human Fertility Database, whose level and timing differ slightly from
# those of the parity data.)
nAges <- nrow(svk_fxs)
parityDist <- matrix(0, nAges, ncol(svk_fxs))
parityDist[1, 1] <- 1
for (a in 1:(nAges - 1)) parityDist[a + 1, ] <- svk_Uxs[[a]] %*% parityDist[a, ]
fx <- rowSums(parityDist * svk_fxs)
message(sprintf("Daughters per woman (GRR): %.3f; mean age at birth: %.1f",
                sum(fx), sum((0:(nAges - 1) + 0.5) * fx) / sum(fx)))
# <<< end


# ==== 2. Age model (Caswell 2019) ====

# >>> rates derived from the parity data, start
# fx are daughters per woman, as in the parity data: birth_female = 1
ageModel <- kin(p = surv, f = fx, birth_female = 1, time_invariant = TRUE)
# <<< end


# ==== 3. Age x parity model (Caswell 2020) ====

outputKin <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya",
               "os", "ys")
parityModel <- kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs,
                               output_kin = outputKin,
                               birth_female = 1,    # one-sex model: all births are female
                               parity = TRUE)       # stages are parities


# ==== 4. Comparison ====

parityTotal <- parityModel %>%
  dplyr::group_by(age_focal, kin) %>%
  dplyr::summarise(count_living = sum(living), .groups = "drop") %>%
  dplyr::mutate(model = "Age x parity (Caswell 2020)")
ageTotal <- ageModel$kin_summary %>%
  dplyr::filter(kin %in% outputKin) %>%
  dplyr::select(age_focal, kin, count_living) %>%
  dplyr::mutate(model = "Age (Caswell 2019)")
both <- dplyr::bind_rows(ageTotal, parityTotal)
utils::write.csv(both, file.path(RESULTS_DIR, "slovakia1980_kin.csv"), row.names = FALSE)

# >>> rates derived from the parity data, start
# relative difference at ages 30 and EGO_AGE of ego, by type of kin
cmpAge <- both %>%
  dplyr::filter(age_focal %in% c(30, EGO_AGE)) %>%
# <<< end
  tidyr::pivot_wider(names_from = model, values_from = count_living)
cmpAge$ratio_parity_age <- cmpAge[["Age x parity (Caswell 2020)"]] / cmpAge[["Age (Caswell 2019)"]]
print(cmpAge, n = Inf)


# ==== 5. Figure ====

figDir <- file.path(RESULTS_DIR, "figures")
dir.create(figDir, showWarnings = FALSE)
pSlovakia <- suppressMessages(DemoKin::rename_kin(both)) %>%
  ggplot(aes(age_focal, count_living, colour = model, linetype = model)) +
  geom_line() +
  facet_wrap(~kin_label) +
  labs(x = "Age of ego", y = "Number of living female kin",
       title = "Slovakia 1980: age model and age x parity model of DemoKin") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figDir, "slovakia1980.pdf"), pSlovakia, width = 29.7, height = 21, units = "cm")
ggsave(file.path(figDir, "slovakia1980.png"), pSlovakia, width = 29.7, height = 21, units = "cm",
       dpi = 200)

# >>> rates derived from the parity data, start
# --- figure of the paper: kin in the direct line and lateral kin ---
paperKin <- c(m = "Mother", d = "Daughters", os = "Older sisters", ys = "Younger sisters",
              oa = "Older aunts", ya = "Younger aunts", coa = "Cousins (older aunts)",
              cya = "Cousins (younger aunts)")
fp <- both[both$kin %in% names(paperKin) & both$age_focal <= 90, ]
fp$kin <- factor(paperKin[fp$kin], paperKin)
fp$model <- factor(fp$model, c("Age (Caswell 2019)", "Age x parity (Caswell 2020)"))
pPaper <- ggplot(fp, aes(age_focal, count_living, colour = model, linetype = model)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~kin, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = c("#2a6fb0", "#555555")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  expand_limits(y = 0) +
  labs(x = "Age of ego", y = "Living kin") +
  theme_classic(base_size = 8, base_family = "serif") +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(size = 7),
        legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figDir, "fig_slovakia.pdf"), pPaper, width = 16, height = 8.5, units = "cm")
ggsave(file.path(figDir, "fig_slovakia.png"), pPaper, width = 16, height = 8.5, units = "cm",
       dpi = 300)
# <<< end


# >>> benchmark without mortality of section 3.3, start
# ==== 6. Without mortality: completed parity and Pearson's formula ====

# Distribution of the number of births of a cohort of women who all survive to the end
# of the reproductive ages, with the rates by age and parity. The number of births is
# followed exactly, including the births of women already in the open class 5+. Each
# birth is a daughter with probability one half, so that the number of daughters has
# mean G = E(B)/2 and variance E(B)/4 + Var(B)/4, where B is the number of births.
pBirth <- 2 * svk_fxs                        # probability of a birth, by age and stage
maxB   <- ncol(pBirth) + nrow(pBirth)
prob   <- c(1, numeric(maxB))                  # prob[n + 1]: probability of n births
for (a in seq_len(nrow(pBirth))) {
  new <- numeric(length(prob))
  for (n in 0:(maxB - 1)) {
    p <- pBirth[a, min(n + 1, ncol(pBirth))]
    new[n + 1] <- new[n + 1] + prob[n + 1] * (1 - p)
    new[n + 2] <- new[n + 2] + prob[n + 1] * p
  }
  prob <- new
}
nb      <- 0:maxB
meanB   <- sum(nb * prob)
varB    <- sum((nb - meanB)^2 * prob)
G       <- meanB / 2
varDau  <- meanB / 4 + varB / 4
pearson <- G + varDau / G - 1

# DemoKin without mortality, age model and age x parity model: sisters ever born at
# age EGO_AGE are the living sisters, since nobody dies
noDeath <- c(rep(1, length(fx) - 1), 0)
ka <- kin(p = noDeath, f = fx, birth_female = 1, time_invariant = TRUE,
          output_kin = c("os", "ys"))$kin_summary
kp <- kin_multi_stage(U = svk_Uxs, f = svk_fxs,
                      D = matrix(noDeath, nrow = length(fx), ncol = ncol(svk_fxs)),
                      H = svk_Hxs, birth_female = 1, parity = TRUE, output_kin = c("os", "ys"))
bench <- data.frame(
  births_0 = prob[1], births_1 = prob[2], births_2 = prob[3], births_3 = prob[4],
  births_4plus = sum(prob[-(1:4)]), mean_daughters = G, var_daughters = varDau,
  ratio = varDau / G, pearson_sisters = pearson,
  demokin_age_sisters = sum(ka$count_living[ka$age_focal == EGO_AGE]),
  demokin_parity_sisters = sum(kp$living[kp$age_focal == EGO_AGE]))
utils::write.csv(bench, file.path(RESULTS_DIR, "slovakia1980_no_mortality.csv"), row.names = FALSE)
print(round(t(bench), 4))
# <<< end
