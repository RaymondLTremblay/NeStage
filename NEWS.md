# NeStage 0.7.1 (2026-03-06)
- Added Quick Start vignette (`NeStage_quickstart.Rmd`): minimum working
  examples for all three core functions with data format checklist.

# NeStage 0.7.0 (2026-03-06)
- Complete roxygen2 documentation for all 9 exported functions.
- Fixed non-ASCII characters (em-dashes) in R/Ne_clonal_Y2000.R,
  R/Ne_mixed_Y2000.R, and R/Ne_sexual_Y2000.R.
- Fixed DOI in all references: Yonezawa et al. (2000) tb01243.x.
- Fixed unused argument `Ne_at_census` in `Ne_clonal_Y2000_both()`
  and `Ne_mixed_Y2000_both()`.
- Added `utils::globalVariables()` to suppress R CMD check notes
  for ggplot2 aesthetic variable names in Ne_sensitivity.R.
- Fixed `[j, i]` roxygen cross-reference warning in three files.
- Fixed LICENSE stub format; added LICENSE.md.
- Corrected ORCID for Raymond L. Tremblay (0000-0002-8588-4372).
- Bumped RoxygenNote to 7.3.3.
- devtools::check() result: 0 errors | 0 warnings | 0 notes.

# NeStage 0.6.1 (2026-03-03)
- Added estimator dispatch with `estimator="paper_v1"` to reproduce
  Yonezawa (2000) Table 4 exactly (Miz/Nan).
- Kept proxy estimators as default to avoid breaking changes.
- Fixed lambda assembly (row 1 = fecundity only).
- Added examples for Table 4 replication and sex-only scenarios.
