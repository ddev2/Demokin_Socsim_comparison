# Step 4. Tables and figures of the paper, from the results of steps 2 and 3.
#
#   results/comparison.csv     all the measures, one row per scenario
# >>> file names distinct from the tables of the paper (written by step 7), start
#   results/runs_detail.csv (.md)       all the measures of each run: rates, variance
#                                       ratios, sisters and maternal cousins
#   results/sisters_formulas.csv (.md)  the two definitions of kin ever born, formula
#                                       [1] of GKP and Pearson's formula [2]
# <<< end
#   results/figures/fig1.pdf   Figure 1: Socsim against DemoKin, runs without minimum
#                              birth interval or heterogeneity
#   results/figures/fig2.pdf   Figure 2: (a) sisters, Socsim - DemoKin against
#                              sigma2/G - 1; (b) cousins against sisters, ratios

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


# ==== 1. Data ====

# Socsim counts of step 2, or the counts of the paper shipped in data/ if step 2 has
# not been run
fss <- file.path(RESULTS_DIR, "socsim_counts.csv")
if (!file.exists(fss) || file.size(fss) < 10) fss <- file.path(REPO_DIR, "data", "socsim_counts.csv")
message("Socsim counts read from ", fss)
ss <- utils::read.csv(fss, stringsAsFactors = FALSE)
dk <- utils::read.csv(file.path(RESULTS_DIR, "demokin_counts.csv"), stringsAsFactors = FALSE)
cmp <- merge(SCENARIOS, ss[, setdiff(names(ss), "seed")], by = "scenario")
cmp <- merge(cmp, dk, by = "scenario")
# >>> every scenario must be present once, start
missing <- setdiff(SCENARIOS$scenario, cmp$scenario)
if (length(missing) > 0) stop("scenarios missing from the Socsim or DemoKin results: ",
                              paste(missing, collapse = ", "))
stopifnot(!anyDuplicated(cmp$scenario))
# <<< end

cmp$VF_ch       <- cmp$var50_ch / cmp$F50_ch      # children, women surviving to 50
cmp$s2G_dau     <- cmp$var_dau / cmp$F_dau        # daughters, all women of mothers' generation
cmp$pearson     <- pearsonSisters(cmp$F_dau, cmp$var_dau)
# >>> tau2 and kappa = (sigma2 - tau2) / mu (formula [4], Appendix B), start
# >>> tau2 from the rates given to DemoKin in step 3 (empirical or input rates), start
fInp <- file.path(RESULTS_DIR, "demokin_inputs.rds")
if (!file.exists(fInp)) stop("run scripts/03_run_demokin.R first (", fInp, " is missing)")
inputs <- readRDS(fInp)
mom <- t(vapply(cmp$scenario, function(s) lambdaMoments(inputs[[s]]), numeric(3)))
# <<< end
cmp$tau2      <- mom[, "tau2"]
cmp$S0        <- mom[, "S0"]                     # formula [1] with r = 0
cmp$kappa_dau <- (cmp$var_dau - cmp$tau2) / cmp$F_dau
# <<< end
cmp$ratio_sis   <- cmp$sisters_dk / cmp$dk_sisters
cmp$ratio_cous  <- cmp$cousins_dk / cmp$dk_cousins
cmp$ratio_aunts <- cmp$aunts_dk / cmp$dk_aunts
cmp$mortality   <- factor(cmp$mortality, c("low", "intermediate", "high"))
cmp$design      <- factor(paste0("bint ", cmp$bint, ", ", ifelse(cmp$hetfert == 1, "hetfert", "no hetfert")),
                        c("bint 0, no hetfert", "bint 12, no hetfert",
                          "bint 0, hetfert", "bint 12, hetfert"))
cmp <- cmp[order(cmp$mortality, cmp$population, cmp$bint, cmp$hetfert), ]
utils::write.csv(cmp, file.path(RESULTS_DIR, "comparison.csv"), row.names = FALSE)


# ==== 2. Tables ====

r2 <- function(x) sprintf("%.2f", x)
mdTable <- function(df, file) {
  df[] <- lapply(df, as.character)
  txt <- c(paste("|", paste(names(df), collapse = " | "), "|"),
           paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
           apply(df, 1, function(x) paste("|", paste(x, collapse = " | "), "|")))
  writeLines(txt, file)
}

tab1 <- data.frame(
  Mortality = cmp$mortality, Scenario = cmp$population, e0 = sprintf("%.0f", cmp$e0),
  TFR = r2(cmp$TFR), bint = cmp$bint, hetfert = ifelse(cmp$hetfert == 1, "yes", "no"),
  "V/F ch." = r2(cmp$VF_ch), "s2/mu dau." = r2(cmp$s2G_dau),   # label s2/mu
  "kappa dau." = r2(cmp$kappa_dau),                             # (sigma2 - tau2) / mu
  "Sisters DemoKin" = r2(cmp$dk_sisters), "Sisters Socsim" = r2(cmp$sisters_dk),
  "Sisters ratio" = sprintf("%.3f", cmp$ratio_sis),
  "Cousins DemoKin" = r2(cmp$dk_cousins), "Cousins Socsim" = r2(cmp$cousins_dk),
  "Cousins ratio" = sprintf("%.3f", cmp$ratio_cous), check.names = FALSE)
utils::write.csv(tab1, file.path(RESULTS_DIR, "runs_detail.csv"), row.names = FALSE)
mdTable(tab1, file.path(RESULTS_DIR, "runs_detail.md"))

tabA1 <- data.frame(
  Mortality = cmp$mortality, Scenario = cmp$population, bint = cmp$bint,
  hetfert = ifelse(cmp$hetfert == 1, "yes", "no"),
  "Sisters GKP" = r2(cmp$GKP_sisters), "Sisters DemoKin" = r2(cmp$dk_sisters),
  "Sisters Socsim (DK)" = r2(cmp$sisters_dk), "Sisters Socsim (all)" = r2(cmp$sisters_all),
  "Sisters Pearson" = r2(cmp$pearson),
  "Cousins DemoKin" = r2(cmp$dk_cousins), "Cousins Socsim (DK)" = r2(cmp$cousins_dk),
  "Cousins Socsim (all)" = r2(cmp$cousins_all), check.names = FALSE)
utils::write.csv(tabA1, file.path(RESULTS_DIR, "sisters_formulas.csv"), row.names = FALSE)
mdTable(tabA1, file.path(RESULTS_DIR, "sisters_formulas.md"))
print(tab1, row.names = FALSE)


# ==== 3. Figures ====

# >>> figures of the paper drawn with ggplot2 (Figures 2 and 3), start
# Figure 1 of the paper (Slovakia 1980) is drawn by scripts/05_slovakia_parity.R.
figDir <- file.path(RESULTS_DIR, "figures")
dir.create(figDir, showWarnings = FALSE)
blue <- "#2a6fb0"; grey <- "#555555"; orange <- "#c8782a"
theme_paper <- theme_classic(base_size = 9, base_family = "serif") +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(size = 9),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = unit(0.4, "cm"))
saveFig <- function(p, name, width, height) {
  ggsave(file.path(figDir, paste0(name, ".pdf")), p, width = width, height = height, units = "cm")
  ggsave(file.path(figDir, paste0(name, ".png")), p, width = width, height = height, units = "cm",
         dpi = 300)
}

# --- Figure 2: Socsim against DemoKin, runs without minimum interval or heterogeneity ---
A <- cmp[cmp$bint == 0 & cmp$hetfert == 0, ]
F2 <- rbind(data.frame(kin = "Sisters", DemoKin = A$dk_sisters, Socsim = A$sisters_dk,
                       mortality = A$mortality),
            data.frame(kin = "Maternal cousins", DemoKin = A$dk_cousins, Socsim = A$cousins_dk,
                       mortality = A$mortality))
F2$kin <- factor(F2$kin, c("Sisters", "Maternal cousins"))
F2$mortality <- factor(F2$mortality, c("low", "intermediate", "high"),
                       c("Low mortality", "United States 1933", "High mortality"))
# square panels starting at zero, with the line of equality
lim2 <- do.call(rbind, lapply(split(F2, F2$kin), function(d) {
  m <- 1.08 * max(d$DemoKin, d$Socsim)
  data.frame(kin = d$kin[1], DemoKin = c(0, m), Socsim = c(0, m))
}))
fig2 <- ggplot(F2, aes(DemoKin, Socsim)) +
  geom_blank(data = lim2) +
  geom_abline(colour = "grey60", linewidth = 0.3) +
  geom_point(aes(shape = mortality, colour = mortality), size = 1.8) +
  facet_wrap(~kin, scales = "free") +
  scale_colour_manual(values = c(blue, orange, grey)) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "DemoKin", y = "Socsim") +
  theme_paper + theme(aspect.ratio = 1)
saveFig(fig2, "fig2_sisters_cousins", 14, 7.5)

# --- Figure 3: (a) sisters, Socsim - DemoKin against sigma2/G - 1, low mortality;
#     (b) ratios Socsim / DemoKin, cousins against sisters, all runs ---
designPoints <- list(
  geom_point(aes(shape = design, colour = design), size = 1.8, stroke = 0.6),
  scale_shape_manual(values = c(16, 17, 5, 0), drop = FALSE),
  scale_colour_manual(values = c(blue, grey, blue, grey), drop = FALSE),
  scale_x_continuous(expand = c(0, 0)), scale_y_continuous(expand = c(0, 0)),
  theme_paper, theme(aspect.ratio = 1, plot.title = element_text(size = 9, hjust = 0.5)))
lowM <- cmp[cmp$mortality == "low", ]
fig3a <- ggplot(lowM, aes(s2G_dau - 1, sisters_dk - dk_sisters)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(colour = "grey60", linewidth = 0.3) +
  coord_cartesian(xlim = c(-0.3, 1.1), ylim = c(-0.3, 1.1)) +
  # plotmath for sigma squared: the pdf() device does not draw these characters as text
  # >>> mean of daughters written mu (it estimates the NRR), start
  labs(x = expression(sigma^2 / mu - 1 ~ "(daughters, mothers' generation)"),
  # <<< end
       y = "Sisters: Socsim - DemoKin", title = "(a) Sisters, low mortality") +
  designPoints +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_blank())
fig3b <- ggplot(cmp, aes(ratio_sis, ratio_cous)) +
  geom_hline(yintercept = 1, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 1, colour = "grey80", linewidth = 0.3) +
  geom_abline(colour = "grey60", linewidth = 0.3) +
  coord_cartesian(xlim = c(0.8, 1.55), ylim = c(0.8, 1.55)) +
  labs(x = "Sisters: Socsim / DemoKin", y = "Maternal cousins: Socsim / DemoKin",
       title = "(b) Cousins against sisters, all runs") +
  designPoints + theme(legend.position = "none")
# the two panels side by side, with the grid package only (no other package needed)
sideBySide <- function(p1, p2) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}
# >>> new Figure 3a, formula [4] tested on all 32 runs, start
# All sisters ever born (not DemoKin's narrower count): Socsim minus formula [1] of GKP,
# against kappa - 1 = (sigma2 - tau2) / mu - 1. Formula [4] predicts that the points lie
# on the diagonal of slope one, at all levels of mortality.
fig3aNew <- ggplot(cmp, aes(kappa_dau - 1, sisters_all - GKP_sisters)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(colour = "grey60", linewidth = 0.3) +
  coord_cartesian(xlim = c(-0.3, 1.1), ylim = c(-0.3, 1.1)) +
  labs(x = expression((sigma^2 - tau^2) / mu - 1 ~ "(daughters, mothers' generation)"),
       y = "All sisters ever born: Socsim - formula [1]",
       title = "(a) Sisters, all runs") +
  designPoints +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_blank())
# <<< end
# >>> Figure 3 of the paper with the new panel (a); the former panel (a),
# low mortality runs only, is kept in fig3_differences_lowmortality, start
for (ext in c("pdf", "png")) {
  for (v in list(list("fig3_differences", fig3aNew), list("fig3_differences_lowmortality", fig3a))) {
    f <- file.path(figDir, paste0(v[[1]], ".", ext))
    if (ext == "pdf") pdf(f, width = 18 / 2.54, height = 9 / 2.54) else
      png(f, width = 18, height = 9, units = "cm", res = 300)
    sideBySide(v[[2]], fig3b)
    invisible(dev.off())
  }
}
# <<< end
# <<< end
message("Tables and figures written to ", RESULTS_DIR)
