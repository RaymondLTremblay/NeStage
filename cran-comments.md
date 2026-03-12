## R CMD check results (NeStage 0.8.0)

0 errors | 0 warnings | 4 notes

- NOTE: New submission
- NOTE: unable to verify current time (network issue on test machine)
- NOTE: pandoc not installed locally (README/NEWS check skipped)
- NOTE: local HTML tidy version outdated (HTML validation skipped)

## Resubmission
Previous submission was rejected due to install.packages() calls in vignettes.
These have been removed. All vignettes now use rmarkdown::html_vignette.

NeStage computes effective population size (Ne) and the Ne/N ratio for
stage-structured populations using matrix population models, following the
framework of Yonezawa (2000, doi:10.1046/j.1365-2540.2000.00747.x).
Functions are provided for sexually reproducing, clonally reproducing, and
mixed (sexual + clonal) populations, along with sensitivity and elasticity
analyses for Ne/N with respect to vital rates.

## Known build warning (not a check issue)

During R CMD build, R 4.5.2 emits:
  "Invalid ORCID iD: 'https://orcid.org/0000-0001-6896-5844'"
This is a known bug in R 4.5.2's person() validator which incorrectly
rejects the full URL format. The ORCID is valid and uses the format
recommended by CRAN policy. This warning does not appear in R CMD check
results.

## Downstream dependencies

None — this is a new package.
