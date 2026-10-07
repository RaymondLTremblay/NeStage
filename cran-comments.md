## Update (NeStage 0.8.1)

This is a patch release fixing a bug in the internally computed generation
time L (used when the user does not supply `L`). Survival was counted twice
and the cohort was started from the stable stage distribution instead of
from newborns, so L was underestimated and Ne/N overestimated. L now
reproduces Table 4 of Yonezawa et al. (2000) exactly. Results obtained with
a user-supplied `L` are unchanged. New regression tests cover the fix.

The DESCRIPTION License field is GPL-3, as in the published 0.8.0.

## R CMD check results

Local, macOS Tahoe 26.6.2, R 4.6.1 (`devtools::check()`, --as-cran):
0 errors | 0 warnings | 0 notes

win-builder, R-devel (2026-10-05 r90641 ucrt):
0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are currently no downstream dependencies for this package.
