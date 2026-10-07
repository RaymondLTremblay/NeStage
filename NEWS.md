# NeStage (development version)

## Bug fix: internally computed generation time L

- When `L = NULL`, `Ne_clonal_Y2000()`, `Ne_sexual_Y2000()`, `Ne_mixed_Y2000()`
  and `Ne_sensitivity_L()` (reference line) computed L incorrectly. The term
  `sum(F * T^x)` already equals l_x * m_x, but it was multiplied by l_x a second
  time, and the cohort was started from the stable stage distribution instead
  of from newborns. The returned L was about one third of the correct value
  (4.728 instead of 13.399 for *Fritillaria camtschatcensis* Miz), so Ne/N was
  overestimated and minimum census size underestimated.
- L is now the mean age of reproduction of a newborn cohort entering
  `recruit_row`, and reproduces Table 4 of Yonezawa et al. (2000) exactly
  (13.399 Miz, 8.353 Nan). Results obtained by supplying `L` directly are
  unchanged.
- `Ne_sensitivity_L()` now calls the same internal function rather than a
  separate copy of the calculation.
- New tests in `tests/testthat/test-generation-time.R`; the Yonezawa vignette
  replication check now compares the computed L with Table 4.

# NeStage 0.8.0 (2026-03-07) — First CRAN submission

- First submission to CRAN
- Package passes R CMD check with 0 errors | 0 warnings | 0 notes
- Tested on macOS (R 4.5.2) and Windows (win-builder devel)

# NeStage 0.7.3 (2026-03-07)

- Removed `Ne_compadre_sexual_filter()` and `Ne_compadre_sexual()` — moved to
  the dedicated research repository `NeStage_COMPADRE`
- Removed `vignettes/NeStage_COMPADRE.Rmd` — moved to `NeStage_COMPADRE`
- Removed `inst/scripts/Ne_compadre_sexual_analysis.R` — moved to `NeStage_COMPADRE`
- Removed `Rcompadre` from Suggests
- Fixed `License` field to reference LICENSE file
- Fixed `ggplot2` declared in Imports
- Removed duplicate `expm` from Suggests
- devtools::check(): 0 errors | 0 warnings | 0 notes

# NeStage 0.7.2 (2026-03-06)

- Added `Ne_compadre_sexual_filter()` and `Ne_compadre_sexual()`
- Added vignette `NeStage_COMPADRE.Rmd`
- Fixed `source()` blocks in all three existing vignettes
- Added `Rcompadre` to Suggests
- devtools::check(): 0 errors | 0 warnings | 0 notes

# NeStage 0.7.1 (2026-03-06)

- Added Quick Start vignette (`NeStage_quickstart.Rmd`)

# NeStage 0.7.0 (2026-03-06)

- Complete roxygen2 documentation for all 9 exported functions
- Added `Ne_sexual_Y2000()`, `Ne_clonal_Y2000()`, `Ne_mixed_Y2000()`
- Added sensitivity/elasticity functions
- devtools::check(): 0 errors | 0 warnings | 0 notes
