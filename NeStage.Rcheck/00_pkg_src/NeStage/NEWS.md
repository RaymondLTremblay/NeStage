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
