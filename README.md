# NeStage

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R package](https://img.shields.io/badge/R-package-blue)](https://github.com/RaymondLTremblay/NeStage)
[![lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

## Overview

**NeStage** computes the **variance effective population size** (*N*e) for
stage-structured plant and animal populations, implementing the framework of
Yonezawa et al. (2000). It is designed for conservation biologists and
population geneticists who need to translate demographic transition matrices
into genetic drift metrics.

In stage-structured populations, individuals move through distinct life-history
stages — for example, seedling → juvenile → reproductive adult — each with
different survival and fecundity rates. Standard *N*e formulas that assume
simple age structure are not appropriate for these systems. NeStage provides the
correct formulation derived directly from the variance of allele-frequency change
across stages.

> **Key biological insight**: *N*e is often only 20–30% of the census size *N*
> in stage-structured populations. To maintain *N*e ≥ 5,000 — the threshold
> recommended for long-term gene diversity conservation (Lande 1995) — census
> populations must typically exceed 17,000–23,000 individuals, far larger than
> naive counts suggest.

---

## The Yonezawa (2000) model

For a population with *s* stages, let *u*ⱼᵢ be the annual transition rate
from stage *i* to stage *j*, and *D*ᵢ the fraction of the population in stage *i*.
The model operates in two tiers:

**General model (Eq. 6)** — for populations with mixed sexual and clonal
reproduction, non-Poisson variance in reproductive output, and arbitrary
deviations from Hardy–Weinberg proportions:

$$N_e = \frac{2N}{V \cdot L}$$

where *V* captures the full variance structure of survival and reproduction
across stages, and *L* is the mean generation time.

**Poisson clonal-dominant case (Eqs. 10–11)** — the simplification that applies
when reproduction is primarily clonal and reproductive contributions follow a
Poisson distribution (as in *Fritillaria camtschatcensis*). Here the key
intermediate quantity is the **stage-weighted mean of squared survival rates**:

$$\overline{u^2} = \sum_{i=1}^{s} D_i \, u_{\cdot i}^2, \qquad
u_{\cdot i} = \sum_j u_{ji}$$

which yields the effective size ratios:

$$\frac{N_y}{N} = \frac{1}{1 - \overline{u^2}}, \qquad
\frac{N_e}{N} = \frac{1}{(1 - \overline{u^2})\,L}$$

> Equations render on GitHub (LaTeX support added 2022). See the vignette for
> a fully typeset and annotated derivation.

### Key quantities computed

| Symbol | Name | Meaning |
|--------|------|---------|
| *ū²* | Stage-weighted mean squared survival | Central intermediate; drives all *N*e calculations |
| *N*y/*N* | Annual effective size ratio | Genetic drift per year relative to census size |
| *N*e/*N* | Generation-time effective size ratio | Drift per generation relative to census size |
| *L* | Generation time | Mean age of reproduction under stage structure (years) |
| Min *N* | Minimum viable census size | Census size for *N*e ≥ 5,000 (Lande 1995) |

---

## Installation

NeStage is not yet on CRAN. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("RaymondLTremblay/NeStage")
```

---

## Quick start

### Replicate Yonezawa et al. (2000) Table 4

The package ships with the *Fritillaria camtschatcensis* data from the
foundational paper. One function call reproduces all values in Table 4:

```r
library(NeStage)

nestage_example_table4()
#> === Table 4 replication ===
#> Miz: L=13.399  Ny/N=2.932 (2.977)  Ne/N=0.219 (0.222)
#> Nan: L=8.353   Ny/N=2.428 (2.444)  Ne/N=0.291 (0.293)
#> # values in () use the stable-stage (expected) distribution instead of observed counts
```

### Single population — observed stage fractions

```r
miz <- get_table2_inputs("Miz")

out <- nestage_exact(
  A_obs = NULL,
  meta  = list(
    estimator  = "paper_v1",
    table2_pop = "Miz",
    Ds         = miz$D_obs,
    population = "Miz"
  )
)

out$L    # 13.399  — generation time (years)
out$NyN  #  2.932  — annual effective size ratio
out$NeN  #  0.219  — generation-time effective size ratio
5000 / out$NeN  # 22,831 — minimum census size for Ne >= 5,000 (Lande 1995)
```

### Using your own transition matrix

```r
# U = survival/transition matrix (rows = destination stage, cols = source stage)
# F = fecundity vector (newborns produced per individual per stage)

U <- matrix(c(
  0.50, 0.10, 0.00,
  0.20, 0.60, 0.05,
  0.00, 0.15, 0.70
), nrow = 3, byrow = TRUE)

F_vec <- c(0.0, 0.5, 2.0)

# Derive stable stage distribution as D
A <- U; A[1, ] <- A[1, ] + F_vec
D <- Re(eigen(A)$vectors[, 1])
D <- abs(D) / sum(abs(D))

# Note: L must be supplied externally. It is computed by iterating the
# matrix to age 500 following the formula in Yonezawa (2000, p. 2008):
#   L = sum(x * sum_j(Fj * uj1x)) / sum(sum_j(Fj * uj1x))
# This iteration is computationally intensive; use a published L from
# your study or derive it from the transition matrix before calling nestage_exact().

out <- nestage_exact(
  A_obs = NULL,
  meta  = list(
    estimator  = "paper_v1",
    U          = U,
    F          = F_vec,
    Ds         = D,
    L          = 10.0,       # generation time in years — supply from your study
    population = "my_species"
  )
)

cat("Ne/N =", round(out$NeN, 3), "\n")
cat("Min N for Ne >= 5,000:", ceiling(5000 / out$NeN), "\n")
```

---

## Functions

| Function | Description |
|----------|-------------|
| `nestage_exact()` | Core *N*e calculation from transition matrix and stage fractions |
| `nestage_example_table4()` | Replicates Table 4 of Yonezawa et al. (2000) |
| `get_table2_inputs()` | Returns *Fritillaria* demographic data for `"Miz"` or `"Nan"` |

---

## Biological system: *Fritillaria camtschatcensis*

The package validation is based on *Fritillaria camtschatcensis* (Liliaceae),
a perennial alpine herb studied at two populations on Mount Hakusan, central
Honshu, Japan (Yonezawa et al. 2000). Plants persist through three
life-history stages:

| Stage | Description |
|-------|-------------|
| Stage 1 | One-leaf, non-flowering |
| Stage 2 | Multi-leaf, non-flowering |
| Stage 3 | Multi-leaf, flowering |

All three stages reproduce clonally via bulblets; no seedling recruitment was
observed in the field, so the populations are maintained almost entirely by
clonal reproduction. Two populations were studied:

| Population | Elevation | *L* (years) | *N*e/*N* | Min *N* for *N*e ≥ 5,000 |
|------------|-----------|-------------|----------|--------------------------|
| **Miz** (Mizuyajiri) | 2,450 m | 13.40 | 0.219 | 22,831 |
| **Nan** (Nanryu)     | 2,050 m |  8.35 | 0.291 | 17,183 |

Miz, at higher elevation with a drier and colder habitat, has a longer
generation time and lower *N*e/*N*, but is subject to less annual drift
(*N*y/*N* = 2.932 vs. 2.428). Both populations require census sizes
well above 17,000 individuals to conserve normal levels of gene diversity
under the Ne ≥ 5,000 criterion of Lande (1995).

---

## Vignettes

```r
vignette("Ne_Yonezawa2000", package = "NeStage")
```

| Vignette | Description |
|----------|-------------|
| `Ne_Yonezawa2000` | Step-by-step replication of Table 4, with full formula derivations, notation reference, and validation tables |

---

## Citation

If you use NeStage in published research, please cite both the package and the
foundational paper:

**Package:**

> Tremblay, R.L. (2026). *NeStage: Effective population size for
> stage-structured populations*. R package version 0.6.1.
> https://github.com/RaymondLTremblay/NeStage

**Foundational paper:**

> Yonezawa, K., Kinoshita, E., Watano, Y., and Zentoh, H. (2000).
> Formulation and estimation of the effective size of stage-structured
> populations in *Fritillaria camtschatcensis*, a perennial herb with a
> complex life history. *Evolution* **54**(6): 2007–2013.
> https://doi.org/10.1111/j.0014-3820.2000.tb01243.x

**Conservation threshold:**

> Lande, R. (1995). Mutation and conservation. *Conservation Biology*
> **9**: 728–791.

---

## Contributing

Bug reports and feature requests are welcome via
[GitHub Issues](https://github.com/RaymondLTremblay/NeStage/issues).
Pull requests are also welcome — please open an issue first to discuss
proposed changes.

---

## License

MIT © Raymond L. Tremblay
