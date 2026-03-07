---
title: 'NeStage: An R Package for Computing Effective Population Size from Stage-Structured Matrix Population Models'
tags:
  - R
  - population genetics
  - conservation biology
  - effective population size
  - matrix population models
  - stage-structured populations
authors:
  - name: Raymond L. Tremblay
    orcid: 0000-0001-6896-5844
    affiliation: 1
affiliations:
  - name: Department of Biology, University of Puerto Rico at Humacao, Puerto Rico, USA
    index: 1

date: 2026-03-07
bibliography: paper.bib
---

# Summary

`NeStage` is an R package that computes the variance effective population size
($N_e$) and the $N_e/N$ ratio for stage-structured populations using matrix
population models (MPMs). The package implements the analytical framework of
@Yonezawa2000, which derives $N_e$ from the demographic parameters encoded in
the survival/transition matrix ($U$) and the fecundity matrix ($F$) of an MPM.
Functions are provided for three reproductive systems: sexually reproducing
populations, clonally reproducing populations, and mixed (sexual + clonal)
populations. Sensitivity and elasticity analyses allow users to identify which
vital rates most strongly influence $N_e/N$, directly informing conservation
management decisions.

# Statement of Need

Effective population size ($N_e$) is one of the most important parameters in
conservation biology, genetics and evolutionary processes. It determines the rate of inbreeding, the loss of
genetic diversity through genetic drift, and the efficacy of natural selection
[@Frankham1995; @Lande1987]. In practice, conservation biologists must estimate
the minimum census population size ($N_{min}$) needed to maintain a target $N_e$
(commonly $N_e \geq 500$ to preserve long-term evolutionary potential;
@Frankham2014).

For most plant and animal species, populations are stage-structured: individuals
differ in their survival, growth, and reproductive rates depending on their life
stage (seedling, juvenile, adult). MPMs are the standard tool for capturing
this demographic heterogeneity [@Caswell2001]. However, most existing software
for $N_e$ estimation either assumes a simple age-structured demography or
requires individual-level genotypic data. No general-purpose R package existed
for computing $N_e$ directly from the MPM matrices that ecologists already
routinely construct.

`NeStage` fills this gap. Given the $U$ and $F$ matrices from any MPM — whether
constructed from field data or obtained from databases such as COMPADRE
[@Salguero2015] — `NeStage` computes $N_e/N$, generation time $L$, the
variance in lifetime reproductive success $V_k$, and the minimum census size
$N_{min}$ needed to achieve any user-specified $N_e$ target. The package is
designed for conservation practitioners, population ecologists, and evolutionary
biologists who work with stage-structured organisms.

# State of the Field

Several R packages address related problems. `popbio` [@Stubben2007] and
`popdemo` [@Stott2012] provide general tools for MPM analysis — computing
lambda, stable stage distributions, and sensitivity analyses for population
growth — but do not compute $N_e$. The `demogR` package computes $N_e$ for
age-classified populations but does not handle stage-classified models or
clonal reproduction. Individual-based approaches to $N_e$ estimation (e.g.,
`NeEstimator`, @Do2014) require genetic data and cannot leverage the large
existing databases of published MPMs.

`NeStage` is unique in combining three features: (1) direct computation of $N_e$
from stage-classified MPMs without genetic data; (2) support for sexual,
clonal, and mixed reproductive systems; and (3) sensitivity/elasticity analyses
that identify the vital rates most influential for $N_e/N$, enabling targeted
management interventions.

# Software Design

The package exports nine functions organized around three reproductive models:

**Sexual reproduction** (`Ne_sexual_Y2000`): Implements equations 3–5 of
@Yonezawa2000. Requires the survival/transition matrix `T_mat`, the fecundity
vector `F_vec` (row 1 of the $F$ matrix), the stage frequency distribution `D`,
and optionally the inbreeding coefficient $F_{IS}$ and variance in offspring
number $V_k$.

**Clonal reproduction** (`Ne_clonal_Y2000`): Implements the clonal model of
@Yonezawa2000, requiring `T_mat`, a clonal propagule vector `C_vec`, and the
stage distribution `D`.

**Mixed reproduction** (`Ne_mixed_Y2000`): Handles populations with both sexual
and clonal pathways, combining the sexual and clonal variance components.

**Sensitivity and elasticity** (`Ne_sensitivity_Vk`, `Ne_sensitivity_L`,
`Ne_sensitivity_Vc`, `Ne_sensitivity_d`): Compute the partial derivatives and
proportional sensitivities of $N_e/N$ with respect to variance in reproductive
success ($V_k$), generation time ($L$), clonal variance ($V_c$), and clonal
fraction ($d$). Results are visualized as `ggplot2` figures.

All functions return tidy data frames and `ggplot2` graphics objects, making
results easy to tabulate or incorporate into downstream analyses. The package
includes four vignettes demonstrating the full workflow from raw matrix data to
publication-ready figures, including a reproduction of Table 4 from @Yonezawa2000
using empirical data for *Fritillaria camtschatcensis*.

# Research Impact

`NeStage` enables two complementary research applications. First, it allows
conservation biologists to compute $N_{min}$ — the minimum census population
size required to maintain a genetically viable population — from demographic
data alone, without requiring genetic sampling. This is particularly valuable
for rare or endangered species where genetic data are unavailable. Second, the
sensitivity analyses identify the specific life-history transitions (e.g.,
adult survival, seedling recruitment) where management effort will have the
greatest positive effect on $N_e/N$, enabling evidence-based prioritization of
conservation actions.

The package has been validated against the analytical results of @Yonezawa2000
for *Fritillaria camtschatcensis* (a clonally reproducing lily) and against
empirical data for *Lepanthes eltoroensis* (a lithophytic orchid endemic to
Puerto Rico), studied by the author [@Tremblay2002]. Integration with the COMPADRE Plant Matrix Database
[@Salguero2015] via the `Rcompadre` package [@Jones2022] allows batch
computation of $N_e/N$ across hundreds of plant species, enabling macroecological
analyses of how life history strategy relates to genetic vulnerability.

# AI Usage Disclosure

The author used AI assistance (Claude, Anthropic) for code review, 
roxygen2 documentation formatting, and vignette drafting. All scientific 
content, mathematical derivations, biological interpretations, and design 
decisions were made by the author.

# Acknowledgements

The author thanks the maintainers of the COMPADRE Plant Matrix Database for
making demographic data freely available.

# References
