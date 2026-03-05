# GitHub Repository Setup Checklist for NeStage

## Files to commit (replace existing)

1. **README.md** — replace with the new README.md provided
2. **.gitignore** — replace with the new .gitignore provided

## Files to DELETE from the repository

These should not be tracked in a public R package repo.
Remove them with:



---

## Manual steps on GitHub.com (cannot be done via file commits)

### 1. Add repository description and website

Go to: https://github.com/RaymondLTremblay/NeStage
Click the gear icon ⚙️ next to "About" (top right of the repo page).

- **Description**: `R package for computing variance effective population size (Ne) in stage-structured populations — implements Yonezawa et al. (2000)`
- **Website**: leave blank for now, or add pkgdown site URL once built
- **Topics**: add all of the following tags:
  - `r`
  - `r-package`
  - `population-genetics`
  - `effective-population-size`
  - `conservation-biology`
  - `stage-structured`
  - `population-ecology`
  - `fritillaria`
  - `demographic-matrix`

### 2. Tag a release (v0.6.1)

Go to: https://github.com/RaymondLTremblay/NeStage/releases/new

- **Tag**: `v0.6.1`
- **Release title**: `NeStage v0.6.1 — Initial public release`
- **Description**:

```
## NeStage v0.6.1

Initial public release of the NeStage R package.

### What's included
- Full implementation of Yonezawa et al. (2000) variance effective population
  size framework for stage-structured populations
- Validated replication of all Table 4 values from the paper
- `paper_v1` exact estimator using the correct overline{u²} formula
- General model (Eq. 6) for mixed sexual/clonal reproduction
- Vignette: Ne_Yonezawa2000 — step-by-step replication with formula derivations
- ggplot2 plotting utilities

### Reference
Yonezawa K., Kinoshita E., Watano Y., and Zentoh H. (2000).
Evolution 54(6): 2007–2013. https://doi.org/10.1111/j.0014-3820.2000.tb01243.x
```

### 3. Enable GitHub Pages (optional but recommended)

Once you run `pkgdown::build_site()` locally and push the `docs/` folder:

Go to: https://github.com/RaymondLTremblay/NeStage/settings/pages
- Source: Deploy from branch
- Branch: `main` / folder: `/docs`

This will give you a package website at:
https://raymondltremblay.github.io/NeStage/

---

## Suggested DESCRIPTION updates

Open `DESCRIPTION` in the repo and ensure these fields are set:

```
Package: NeStage
Title: Variance Effective Population Size for Stage-Structured Populations
Version: 0.6.1
Authors@R: person("Raymond L.", "Tremblay",
    email = "your@email.com",
    role = c("aut", "cre"),
    comment = c(ORCID = "YOUR-ORCID-HERE"))
Description: Computes the variance effective population size (Ne) for
    stage-structured populations following the framework of Yonezawa et al.
    (2000) <doi:10.1111/j.0014-3820.2000.tb01243.x>. Provides exact
    replication of the paper's Table 4 results and a general model for
    populations with mixed sexual and clonal reproduction.
License: MIT + file LICENSE
URL: https://github.com/RaymondLTremblay/NeStage
BugReports: https://github.com/RaymondLTremblay/NeStage/issues
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.x.x
Suggests:
    knitr,
    rmarkdown,
    ggplot2,
    gt,
    popbio,
    popdemo,
    testthat (>= 3.0.0)
VignetteBuilder: knitr
```
