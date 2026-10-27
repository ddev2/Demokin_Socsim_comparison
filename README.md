# Lateral kin in matrix kinship models: DemoKin and Socsim

Code and data for the comparison of the numbers of lateral kin (sisters, maternal
aunts, maternal cousins) given by the one-sex matrix kinship model of Caswell (2019),
as implemented in the R package DemoKin, with those of the Socsim microsimulation
model, in the paper <!-- >>> new title, start -->*How many sisters? Bias in the
counts of lateral kin from integral and matrix kinship models*<!-- <<< end -->
(Devolder, in preparation).

<!-- >>> paper withdrawn until acceptance, start -->
**Paper:** under review; it will be added to the repository once accepted.
<!-- <<< end -->

The work has five steps:

1. Socsim simulations with constant rates, one per scenario (32 scenarios).
2. In each simulation, selection of the egos and count of their kin at age 65 under
   the definition of kin ever born used by DemoKin; estimation of the female rates of
   the simulation.
3. DemoKin with the rates of each simulation, and formula [1] of Goodman, Keyfitz and
   Pullum (1974).
4. Tables and figures of the paper.
5. A separate check of DemoKin with Slovakia 1980: the age model of Caswell (2019) and
   the age x parity model of Caswell (2020) run on the same population.

Steps 3, 4 and 5 can be run without Socsim, with the rates and kin counts of the
simulations of the paper stored in `data/`.


## Contents

```
config.R                    paths and parameters, loads the libraries
R/
  socsim_lib.R              scenario folders, initial population, running Socsim,
                            reading results, selection of egos
  genealogy_lib.R           genealogies of egos (kin of each type, KinFert format)
  kin_count_lib.R           counts of kin with the DemoKin definition, rates of a run,
                            kin by age of ego in DemoKin format from genealogies
  demokin_lib.R             rates for DemoKin, DemoKin runs, formulas [1] and [2]
scripts/
  01_run_socsim.R           step 1
  02_socsim_kin.R           step 2
  03_run_demokin.R          step 3
  04_tables_figures.R       step 4
  05_slovakia_parity.R      step 5
<!-- >>> steps 6 and 7, start -->
  06_demokin_scripts.R      step 6: one self-contained DemoKin script per run
  07_tableA1_excel.R        step 7: Table A1 of the paper in an Excel file
demokin/                    the 32 DemoKin scripts written by step 6, and run_all.R
<!-- <<< end -->
scenarios/
  scenarios.csv             the 32 scenarios: population, mortality, bint, hetfert, seed
  <scenario>/               .sup file, rates/ and initial population (presim.opop.gz)
data/
  empirical_rates/          female rates of the 32 simulations of the paper (step 2)
  socsim_counts.csv         kin counts of the 32 simulations of the paper (step 2)
results/                    created by the scripts (not under version control)
```


## Requirements

- R 4.1 or later.
- [rsocsim](https://github.com/MPIDR/rsocsim) for step 1 (tested with versions 1.7.3
  and 1.9.18):
  `remotes::install_github("MPIDR/rsocsim")`.
- [DemoKin](https://cran.r-project.org/package=DemoKin) 1.0.3 for step 3.
- data.table, dplyr, tidyr, ggplot2, and callr (strongly recommended: each simulation
  then runs in its own R process):
  `install.packages(c("DemoKin", "data.table", "dplyr", "tidyr", "ggplot2", "callr"))`.

Open `lateral-kin.Rproj` in RStudio, or set the working directory to this folder, and
run the scripts in order. Each script finds the folder of the repository by itself
(from RStudio, from `Rscript scripts/0x_....R`, or from the working directory).


## Step 1. Socsim simulations: `scripts/01_run_socsim.R`

Runs every scenario of `scenarios/scenarios.csv` with its seed. Socsim writes the
results to `scenarios/<scenario>/sim_results_<seed>_/` (rsocsim 1.7) or
`scenarios/<scenario>/sim_results_<scenario>.sup_<seed>_/` (rsocsim 1.9 and later); the
code reads both. A scenario whose results already exist is skipped.

- Each run starts from the initial population stored in `presim.opop.gz`, so that the
  runs of the paper can be reproduced with the same seeds. For a new scenario without
  an initial population, a stable population is built from its rates.
- Socsim runs in a local temporary folder and the results are then copied to the
  scenario folder. Socsim creates `result.opop` empty at the start of a run and fills it
  at the end; in a folder synchronised by Google Drive or Dropbox the empty file can
  replace the results.
- Run times range from a few minutes (high mortality) to several hours (low mortality
  and high fertility, with populations of several million persons and result files of
  up to 650 MB). Several R sessions can run different scenarios at the same time, but
  not the same scenario.
- `newScenarioFolder("TFR=4_LM_bint12", from = "TFR=4_LM", bint = 12)` creates a new
  scenario with the rates and initial population of an existing one.

To use simulations that are stored elsewhere (for instance the original project folder),
set the folder of the scenarios before running the scripts:

```r
Sys.setenv(KIN_SCENARIO_DIR = "/path/to/data/ socsim")
```


## Step 2. Kin in Socsim: `scripts/02_socsim_kin.R`

For each simulation (function `socsimKinScenario()`):

- **Egos** (`selectEgos()`): women born during the simulation whose mother and maternal
  grandmother were also born during it, so that the fertility histories of both are
  complete, who are at least 70 years old at the end of the run, and born in the last
  60 years for which these conditions are met (about 9,000 to 445,000 egos per run).
- **Kin** (`countKinSocsim()`), through the female line, as in the one-sex model:
  sisters (other daughters of ego's mother), maternal aunts (other daughters of the
  maternal grandmother), maternal cousins (daughters of these aunts) and nieces
  (daughters of the sisters), born before ego reaches age 65. Two counts are given:
  - `_dk`: the definition of kin ever born of DemoKin, that is, the kin alive at the
    birth of ego or born after it (kin who died before ego was born are left out);
  - `_all`: all kin ever born.
  Counts are means over all egos, whether or not ego survives to 65 (the kin of ego do
  not depend on her survival).
- **Variance ratios**, for the women of the mothers' generation (born during the period
  of birth of the egos' mothers): mean and variance of the numbers of children
  (`F50_ch`, `var50_ch`, women surviving to 50) and of daughters (`F_dau`, `var_dau`, all
  women). V/F for children describes the birth process; sigma2/G for daughters enters
  formula [2].
- **Rates** (`empiricalRates()`): female exposure, deaths, births and births of daughters
  by single year of age over the last 100 years of the run.

Output: `results/socsim_counts.csv` (one row per scenario), `results/socsim_counts/`,
`results/empirical_rates/` and `results/egos/`. Most of the time is spent reading
`result.opop`: from a few seconds to a few minutes per run, with up to a few GB of memory
for the largest runs (about 5 million persons). The files `data/empirical_rates/` and
`data/socsim_counts.csv` are copies of these results for the runs of the paper.

With `GENEALOGIES <- TRUE`, the script also builds the full genealogies of a sample of
20,000 egos (`retrieveKin()`, `socsimToGenealogies()`) and their kin by age of ego in
the format of DemoKin (`genealogiesToDemoKin()`: `kin_full` and `kin_summary` with
`count_living`, `count_dead`, `count_cum_dead` and `count_everborn`). These counts are
divided by the number of egos alive at each age, so that, with high mortality, they
are conditional on the survival of ego and differ somewhat from the counts above.
They are not needed for the tables and figures.


## Step 3. DemoKin: `scripts/03_run_demokin.R`

DemoKin (`kin()`, one-sex, time-invariant) is run for each scenario with the female
rates of the simulation (`results/empirical_rates/`, or `data/empirical_rates/` if step 2
has not been run) and the proportion of female births observed in it. Kin ever born
at exact age x of ego are `count_living(x) + count_cum_dead(x - 1)`, since DemoKin
records at age x the deaths that occur during the following year; this is the same
definition as the `_dk` counts of step 2. The script also computes formula [1] of
Goodman, Keyfitz and Pullum (1974) for sisters with the same rates (`gkpSisters()`),
and the TFR, GRR, NRR and female life expectancy of the rates.

With `RATES <- "input"`, DemoKin is run instead with the rates written in the Socsim
rate files of each scenario (`inputsFromScenario()`).

<!-- >>> rates saved for step 4, start -->
Output: `results/demokin_counts.csv`, `results/demokin_kin_summary.rds` and
`results/demokin_inputs.rds` (the rates given to DemoKin, used by step 4 to compute tau2
with the same rates). The NRR is computed with the person-years lived at each age, as
in Table A1.
<!-- <<< end -->


## Step 4. Tables and figures: `scripts/04_tables_figures.R`

The Socsim counts are read from `results/socsim_counts.csv`, or from
`data/socsim_counts.csv` if step 2 has not been run.

- `results/comparison.csv`: all measures, one row per scenario.
<!-- >>> file names distinct from the tables of the paper, start -->
- `results/runs_detail.csv` and `.md`: one row per run, sisters and maternal cousins at
  age 65 in DemoKin and Socsim, with V/F (children), sigma2/mu and kappa (daughters).
- `results/sisters_formulas.csv` and `.md`: the two definitions of kin ever born, formula
  [1] and Pearson's formula [2] (`pearsonSisters()`: mu + sigma2/mu - 1).
- The script stops if a scenario is missing from the Socsim or DemoKin results.
<!-- <<< end -->
<!-- >>> figures of the paper drawn with ggplot2, start -->
- `results/figures/fig2_sisters_cousins.pdf` (and `.png`): Figure 2 of the paper, Socsim
  against DemoKin, runs without minimum birth interval or heterogeneity.
<!-- >>> new panel (a) of Figure 3, start -->
- `results/figures/fig3_differences.pdf` (and `.png`): Figure 3 of the paper, (a) all
  sisters ever born, Socsim minus formula [1], against kappa - 1 = (sigma2 - tau2)/mu - 1,
  all runs (the check of formula [4]); (b) ratios Socsim / DemoKin, cousins against
  sisters, all runs. The former panel (a), low mortality runs only, is kept in
  `fig3_differences_lowmortality.pdf`.
<!-- <<< end -->

All the figures of the paper are drawn with ggplot2: Figure 1 by step 5, Figures 2 and 3
by step 4.
<!-- <<< end -->


## Step 5. Slovakia 1980, DemoKin with and without parity: `scripts/05_slovakia_parity.R`

The same population, Slovakia in 1980, is run in DemoKin with the age model of Caswell
(2019), `kin()`, and with the age x parity model of Caswell (2020),
`kin_multi_stage(..., parity = TRUE)`, using the rates by age and parity distributed with
<!-- >>> rates derived from the parity data, start -->
DemoKin (`svk_Uxs`, `svk_fxs`, `svk_pxs`, `svk_Hxs`). The age model uses the
age-specific fertility rates implied by these data (the parity distribution of a cohort at
each age, applied to the fertility of each parity), so that both models have exactly the
same fertility by age. The script compares the numbers of living kin by age of ego: kin in
the direct line and older sisters and aunts are the same in both models; younger sisters
and aunts are about 4 per cent fewer in the parity model (4 per cent more with
`parity = FALSE`).

<!-- >>> source of the Slovak data of DemoKin, start -->
The rates by age and parity distributed with DemoKin are those of Slovakia in **1980**
(the DemoKin documentation says 1990): they reproduce exactly the conditional fertility
rates of the Human Fertility Database for 1980 (file `SVKmi.txt`, version of 2018),
converted to probabilities u = m / (1 + 0.5 m), with fertility placed one year younger
than in the HFD. The survival probabilities are those of the HMD female life table of 1980
(e0 = 74.3 years), with ages not shifted.
<!-- <<< end -->

Output: `results/slovakia1980_kin.csv`, `results/figures/slovakia1980.pdf` (all kin types)
and `results/figures/fig_slovakia.pdf` (Figure 1 of the paper), with `.png` versions.
<!-- >>> benchmark without mortality, start -->
`results/slovakia1980_no_mortality.csv` gives the numbers of section 3.3 without
mortality: the distribution of the number of births (followed exactly, including births
in the open class 5+), the mean and variance of the number of daughters, Pearson's
formula, and the sisters ever born at age 65 in the age model and the age x parity model.
<!-- <<< end -->
<!-- <<< end -->
It takes less than a minute.


<!-- >>> steps 6 and 7, start -->
## Step 6. One DemoKin script per run: `scripts/06_demokin_scripts.R`

Writes `demokin/<scenario>.R` for each of the 32 runs: a short script that contains the
rates given to DemoKin and runs it, so that the DemoKin side of the comparison can be
read and reproduced with the DemoKin package alone, without Socsim or the rest of the
repository. `demokin/run_all.R` runs the 32 scripts and writes `demokin/demokin_results.csv`.

The rates are those of step 3. They are estimated from each Socsim run over its last 100
years (`empiricalRates()` in `R/kin_count_lib.R`), by single year of age:

- exposure: years lived by women at each age during these years, from their dates of birth
  and death in months;
- deaths of women by age at death, and births (of both sexes and of daughters) by the age
  of the mother at the birth, in completed years;
- death rate m(x) = deaths / exposure, and survival probability p(x) = exp(-m(x)), set to 0
  at age 100 to close the life table;
- fertility rate f(x) = births of both sexes / exposure, and proportion of female births
  of the same years, which DemoKin applies to f(x) (`birth_female`).

Kin ever born at age 65 are `count_living(65) + count_cum_dead(64)`, as in step 3.


## Step 7. Table A1 in Excel: `scripts/07_tableA1_excel.R`

Gathers the description of the runs (`scenarios/scenarios.csv`), the Socsim counts of step
2 (`results/socsim_counts.csv`, or `data/socsim_counts.csv` for the runs of the paper) and
DemoKin run with the rates of each run, and writes Table A1 of the paper in the layout of
the paper: `results/paper_TableA1.xlsx` (and `results/paper_TableA1.csv`; the names differ
from `sisters_formulas.csv` of step 4, since macOS does not distinguish upper and lower
case in file names). It needs DemoKin and
openxlsx (`install.packages("openxlsx")`) and takes less than a minute. The Korean
scenarios keep their original folder names (`KOREA_2022...`) but are labelled Korea 2018
in `scenarios/scenarios.csv` and in all the tables, as in the paper. Scripts 06 and 07
use the results folder of `config.R` (environment variable `KIN_RESULTS_DIR`).
<!-- <<< end -->

## Scenarios

`scenarios/scenarios.csv` gives, for each scenario, the population (synthetic schedules
with TFR 3, 4 and 5; Korea 2018; United States 1933), the mortality level, the minimum
interval between births (`bint`, in months: 0 or 12), heterogeneity in fecundability
(`hetfert`: 0 or 1), the seed and the length of the run. The female life expectancy at
birth is about 95 years with low mortality (85 for Korea), 63 for the United States in
1933, and 31 to 51 years with high mortality.

The names of the folders are those of the original runs and do not follow one rule: for
instance `TFR=3_LM_HetFert` has `bint` 12 and `hetfert` 1, while `TFR=3_HM_HetFert` has
`bint` 0 and `hetfert` 1. The values in `scenarios.csv`, read from the `.sup` files, are
the reference.

The rate files follow the Socsim format: each line `y m r` gives the monthly rate `r` up
to the upper age bound of `y` years and `m` months. `readSocsimRates()` reads them as
Socsim applies them.


## Git

<!-- >>> paper withdrawn until acceptance, start -->
The paper is not part of the repository until it is accepted for publication
(`.gitignore`).
<!-- <<< end -->
The results of the scripts (`results/`) and the Socsim outputs are not
under version control (`.gitignore`).


## Reference

Caswell, H. (2019). The formal demography of kinship: A matrix formulation.
*Demographic Research* 41(24): 679-712.

Caswell, H. (2020). The formal demography of kinship II: Multistate models, parity, and
sibship. *Demographic Research* 42(38): 1097-1146.

Goodman, L.A., Keyfitz, N., and Pullum, T.W. (1974). Family formation and the frequency
of various kinship relationships. *Theoretical Population Biology* 5(1): 1-27.
