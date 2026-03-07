# NeStage CRAN Submission Checklist

## Step 1 — Spell check
Run in RStudio console:
```r
devtools::spell_check()
```
Fix any real typos in documentation or DESCRIPTION. Add false positives
(e.g. "Ne", "Yonezawa", "elasticity") to the `inst/WORDLIST` file.

## Step 2 — Check on Windows (win-builder)
```r
devtools::check_win_devel()
```
This emails results to raymond.tremblay@upr.edu within ~30 minutes.
Wait for the email and confirm 0 errors | 0 warnings | 0 notes before proceeding.

## Step 3 — Check on R-hub (multiple platforms)
```r
# Install rhub if needed
install.packages("rhub")

# Check on multiple platforms
rhub::check_for_cran()
```
This tests on Linux and Windows simultaneously.

## Step 4 — Bump version to 0.8.0 for CRAN
CRAN prefers x.y.z versions where x.y.z >= 0.1.0.
Update DESCRIPTION: Version: 0.8.0
Update NEWS.md with a # NeStage 0.8.0 (CRAN submission) entry.

## Step 5 — Final local check
```r
devtools::check(remote = TRUE, manual = TRUE)
```

## Step 6 — Submit
```r
devtools::submit_cran()
```
This builds the package, runs final checks, and opens the CRAN web form.
On the form:
- Confirm maintainer email
- Paste contents of cran-comments.md into the "Comments" box
- Submit

## Step 7 — Confirm email
CRAN sends a confirmation email — click the link to complete submission.

## Step 8 — Wait and respond
- Initial automated checks: within hours
- Human review: 2–10 days
- If reviewers request changes: make them promptly and resubmit with
  updated cran-comments.md noting what was changed

## Files added to NeStage root for submission
- cran-comments.md   ← this checklist's companion file
