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
stages — seedling → juvenile → reproductive adult — each with different
survival and fecundity rates. Standard *N*e formulas that assume simple age
structure are not appropriate for these systems. NeStage derives *N*e directly
from the variance of allele-frequency change across stages.

> **Key biological insight**: *N*e is often only 10–30% of the census size *N*
> in stage-structured populations (Frankham 1995). A population of 40
> *Lepanthes* orchids — already large for this endangered genus — typically
> achieves Ne of only 7–37, well below the short-term inbreeding threshold
> of Ne ≥ 50 (Franklin 1980).

---

## The Yonezawa (2000) model

For a population with *s* stages, the **generation-time effective size** is:

$$N_e = \frac{2N}{V \cdot L}$$

where *V* captures the full variance structure of survival and reproduction
across stages and *L* is the mean generation time. Three special cases are
implemented:

| Model | When to use |
|-------|------------|
| `Ne_clonal_Y2000()` | All reproduction is clonal (e.g. *Fritillaria*) |
| `Ne_sexual_Y2000()` | All reproduction is sexual (e.g. *Lepanthes* orchids) |
| `Ne_mixed_Y2000()` | Mixed sexual + clonal reproduction; full general model (Eq. 6) |

---

## Installation

NeStage is not yet on CRAN. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("RaymondLTremblay/NeStage")
```

---

## Quick start

### Sexual reproduction — *Lepanthes rupestris* population 1

```r
library(NeStage)

# Survival-transition matrix (MatU) — monthly time step
T_mat <- matrix(c(
  0,     0,     0,     0,     0,     0,
  0.738, 0.738, 0,     0,     0,     0,
  0,     0,     0.515, 0,     0.076, 0.013,
  0,     0.038, 0,     0.777, 0,     0,
  0,     0.002, 0.368, 0.011, 0.730, 0.171,
  0,     0,     0.037, 0,     0.169, 0.790
), nrow = 6, byrow = TRUE)

# Fecundity vector (row 1 of MatF) — monthly newborns per individual
F_vec <- c(0, 0, 0, 0, 0.120, 0.414)

# Stage frequency vector D (from observed counts)
D <- c(22, 22, 36, 36, 48.5, 48.5)
D <- D / sum(D)

out <- Ne_sexual_Y2000(
  T_mat      = T_mat,
  F_vec      = F_vec,
  D          = D,
  Ne_target  = 50,     # Franklin (1980) short-term inbreeding threshold
  census_N   = 40,     # actual or expected census population size
  population = "L. rupestris pop 1"
)

print(out)
#>
#>  ─────────────────────────────────────────────────────
#>   Ne_sexual_Y2000 results  ·  L. rupestris pop 1
#>  ─────────────────────────────────────────────────────
#>   Stages s                           = 6
#>   Stage-weighted survival (u_bar)    = 0.741
#>   Generation time L                  = 14.6 months
#>   V (total variance)                 = 0.420
#>   Ne/N                               = 0.572
#>   --- Conservation threshold ---
#>   Ne target                          = 50
#>   Minimum census size N              = 88
#>   Ne at your census size (N = 40)    = 22.9
#>   WARNING: Ne (22.9) < Ne target (50) at N = 40
#>  ─────────────────────────────────────────────────────
```

### Clonal reproduction — *Fritillaria camtschatcensis* (replicates Yonezawa 2000)

```r
T_Miz <- matrix(c(
  0.789, 0.121, 0.054,
  0.007, 0.621, 0.335,
  0.001, 0.258, 0.611
), nrow = 3, byrow = TRUE)

out_Miz <- Ne_clonal_Y2000(
  T_mat      = T_Miz,
  F_vec      = c(0.055, 1.328, 2.398),
  D          = c(0.935, 0.038, 0.027),
  L          = 13.399,   # from Yonezawa et al. (2000) Table 4
  Ne_target  = 5000,     # Lande (1995) long-term evolutionary threshold
  census_N   = 200,
  population = "Fritillaria Miz"
)
# Ne/N = 0.219  — matches Table 4 exactly
```

---

## Sensitivity analysis

Four functions sweep one parameter at a time and return a data frame,
a ggplot2 figure, and elasticity statistics — telling you which parameter
drives Ne/N most and where management effort has the greatest genetic return.

```r
sens <- Ne_sensitivity_Vk(
  model_fn    = Ne_sexual_Y2000,
  T_mat       = T_mat,
  F_vec       = F_vec,
  D           = D,
  stage_index = 5,                       # reproductive adults
  Vk_range    = seq(0.5, 8, by = 0.5),
  Ne_target   = 50,
  population  = "L. rupestris pop 1"
)

print(sens)    # elasticity table
sens$plot      # ggplot2 figure
```

| Function | Parameter swept | Management question |
|---|---|---|
| `Ne_sensitivity_Vk()` | Sexual reproductive variance Vk/k̄ | How unequal is seed parentage? |
| `Ne_sensitivity_Vc()` | Clonal reproductive variance Vc/c̄ | Do dominant clonal genets reduce Ne? |
| `Ne_sensitivity_d()` | Clonal fraction *d* | Does shifting from sexual to clonal matter? |
| `Ne_sensitivity_L()` | Generation time *L* | How much does L uncertainty affect Ne/N? |

---

## Choosing your Ne target

The `Ne_target` parameter controls what minimum Ne is considered viable.

| Ne target | Criterion | Source |
|-----------|-----------|--------|
| **50** | Avoid short-term inbreeding depression | Franklin (1980) |
| **500** | Maintain long-term quantitative genetic variation | Franklin (1980) |
| **5000** | Preserve long-term evolutionary potential | Lande (1995) |

For small, range-restricted species such as *Lepanthes* (census N typically
< 100), the Ne ≥ 50 threshold is the relevant near-term goal. The Ne ≥ 5,000
criterion is more appropriate for large-population conservation planning.

---

## Functions

| Function | Model | Description |
|----------|-------|-------------|
| `Ne_clonal_Y2000()` | Clonal | Eqs. 10–11 of Yonezawa et al. (2000) |
| `Ne_clonal_Y2000_both()` | Clonal | Observed + expected D in one call |
| `Ne_sexual_Y2000()` | Sexual | Full sexual model |
| `Ne_mixed_Y2000()` | Mixed | General model (Eq. 6); any d, Vk, Vc, a |
| `Ne_mixed_Y2000_both()` | Mixed | Observed + expected D in one call |
| `Ne_sensitivity_Vk()` | Sexual, mixed | Sensitivity to sexual repro variance |
| `Ne_sensitivity_Vc()` | Mixed | Sensitivity to clonal repro variance |
| `Ne_sensitivity_d()` | Mixed | Sensitivity to clonal fraction |
| `Ne_sensitivity_L()` | All | Sensitivity to generation time |

---

## Vignettes

| Vignette | Description |
|----------|-------------|
| `Ne_Yonezawa2000` | Step-by-step replication of Table 4 of Yonezawa et al. (2000); formula derivations and notation reference |
| `NeStage_functions` | Applied user guide: Ne/N for 17 populations of three Puerto Rican *Lepanthes* orchid species (Tremblay & Ackerman 2001) |
| `NeStage_sensitivity` | Conservation decision-making: elasticity analysis and a management prioritisation framework |

```r
vignette("Ne_Yonezawa2000",     package = "NeStage")
vignette("NeStage_functions",   package = "NeStage")
vignette("NeStage_sensitivity", package = "NeStage")
```

---

## Citation

**Package:**

> Tremblay, R.L. (2026). *NeStage: Effective population size for
> stage-structured populations*. R package version 0.1.0.
> https://github.com/RaymondLTremblay/NeStage

**Foundational paper:**

> Yonezawa, K., Kinoshita, E., Watano, Y., and Zentoh, H. (2000).
> Formulation and estimation of the effective size of stage-structured
> populations in *Fritillaria camtschatcensis*. *Evolution* **54**(6):
> 2007–2013. https://doi.org/10.1111/j.0014-3820.2000.tb01244.x

**Field data:**

> Tremblay, R.L. and Ackerman, J.D. (2001). Gene flow and effective
> population size in *Lepanthes* (Orchidaceae): a case for genetic drift.
> *Biological Journal of the Linnean Society* **72**: 47–62.

---

## References

- Franklin I.R. (1980). Evolutionary change in small populations. In: Soulé
  & Wilcox, eds. *Conservation Biology*. Sinauer. pp. 135–149.
- Frankham R. (1995). Effective population size/adult population size ratios
  in wildlife. *Genetical Research* **66**: 95–107.
- Lande R. (1995). Mutation and conservation. *Conservation Biology* **9**:
  782–791.
- Orive M.E. (1993). Effective population size in organisms with complex
  life histories. *Theoretical Population Biology* **44**: 316–340.

---

## License

MIT © Raymond L. Tremblay
