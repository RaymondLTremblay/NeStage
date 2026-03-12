pkgname <- "NeStage"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "NeStage-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('NeStage')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Ne_clonal_Y2000")
### * Ne_clonal_Y2000

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Ne_clonal_Y2000
### Title: Variance effective population size for 100% clonal
###   stage-structured populations
### Aliases: Ne_clonal_Y2000

### ** Examples

# ---------------------------------------------------------
# Example 1: Replicate Yonezawa Table 4 for population Miz
# ---------------------------------------------------------
# Transition matrix from Table 2 of Yonezawa et al. (2000)
T_Miz <- matrix(c(
  0.789, 0.121, 0.054,
  0.007, 0.621, 0.335,
  0.001, 0.258, 0.611
), nrow = 3, byrow = TRUE)

# Observed stage fractions (Table 2)
D_Miz <- c(0.935, 0.038, 0.027)

# Fecundity vector (Table 2) -- clonal propagules per plant per year
F_Miz <- c(0.055, 1.328, 2.398)

# Using L from Table 4 directly (recommended for exact replication)
result <- Ne_clonal_Y2000(
  T_mat      = T_Miz,
  F_vec      = F_Miz,
  D          = D_Miz,
  L          = 13.399,
  population = "Miz"
)
print(result)
# Expected: Ny/N = 2.932, Ne/N = 0.219, Min_N = 22,831

# ---------------------------------------------------------
# Example 2: Let the function compute L internally
# ---------------------------------------------------------
result2 <- Ne_clonal_Y2000(
  T_mat      = T_Miz,
  F_vec      = F_Miz,
  D          = D_Miz,
  population = "Miz (computed L)"
)
print(result2)

# ---------------------------------------------------------
# Example 3: Replicate Table 4 for population Nan
# ---------------------------------------------------------
T_Nan <- matrix(c(
  0.748, 0.137, 0.138,
  0.006, 0.669, 0.374,
  0.001, 0.194, 0.488
), nrow = 3, byrow = TRUE)

D_Nan <- c(0.958, 0.027, 0.015)
F_Nan <- c(0.138, 2.773, 5.016)

result3 <- Ne_clonal_Y2000(
  T_mat      = T_Nan,
  F_vec      = F_Nan,
  D          = D_Nan,
  L          = 8.353,
  population = "Nan"
)
print(result3)
# Expected: Ny/N = 2.428, Ne/N = 0.291, Min_N = 17,183




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Ne_clonal_Y2000", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Ne_clonal_Y2000_both")
### * Ne_clonal_Y2000_both

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Ne_clonal_Y2000_both
### Title: Run Ne_clonal_Y2000 with both observed and expected stage
###   fractions
### Aliases: Ne_clonal_Y2000_both

### ** Examples

T_Miz <- matrix(c(
  0.789, 0.121, 0.054,
  0.007, 0.621, 0.335,
  0.001, 0.258, 0.611
), nrow = 3, byrow = TRUE)

both <- Ne_clonal_Y2000_both(
  T_mat      = T_Miz,
  F_vec      = c(0.055, 1.328, 2.398),
  D_obs      = c(0.935, 0.038, 0.027),
  D_exp      = c(0.921, 0.046, 0.033),
  L          = 13.399,
  population = "Miz"
)
print(both$observed)   # main values in Table 4
print(both$expected)   # parenthetical values in Table 4




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Ne_clonal_Y2000_both", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Ne_mixed_Y2000")
### * Ne_mixed_Y2000

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Ne_mixed_Y2000
### Title: Variance effective population size for mixed sexual/clonal
###   stage-structured populations
### Aliases: Ne_mixed_Y2000

### ** Examples

# ---------------------------------------------------------
# Example 1: 3-stage perennial herb, stage 3 reproduces both
# sexually (seeds) and clonally (rhizomes). Stages 1 and 2
# reproduce sexually only.
# ---------------------------------------------------------
T_herb <- matrix(c(
  0.30, 0.05, 0.00,
  0.40, 0.65, 0.10,
  0.00, 0.20, 0.80
), nrow = 3, byrow = TRUE)

result <- Ne_mixed_Y2000(
  T_mat      = T_herb,
  F_vec      = c(0.0, 0.5, 3.0),
  D          = c(0.60, 0.25, 0.15),
  d          = c(0.0, 0.0, 0.7),   # stage 3 is 70% clonal
  population = "mixed herb"
)
print(result)

# ---------------------------------------------------------
# Example 2: Show Ny/N for a predominantly clonal population
# ---------------------------------------------------------
result_Ny <- Ne_mixed_Y2000(
  T_mat      = T_herb,
  F_vec      = c(0.0, 0.5, 3.0),
  D          = c(0.60, 0.25, 0.15),
  d          = c(0.5, 0.8, 0.9),
  show_Ny    = TRUE,
  population = "predominantly clonal herb"
)
print(result_Ny)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Ne_mixed_Y2000", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Ne_sensitivity_L")
### * Ne_sensitivity_L

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Ne_sensitivity_L
### Title: Sensitivity of Ne/N to generation time L
### Aliases: Ne_sensitivity_L

### ** Examples

# Clonal model -- how sensitive is Ne/N to L for Fritillaria Miz?
T_Miz <- matrix(c(
  0.789, 0.121, 0.054,
  0.007, 0.621, 0.335,
  0.001, 0.258, 0.611
), nrow = 3, byrow = TRUE)

sens <- Ne_sensitivity_L(
  model_fn   = Ne_clonal_Y2000,
  T_mat      = T_Miz,
  F_vec      = c(0.055, 1.328, 2.398),
  D          = c(0.935, 0.038, 0.027),
  L_range    = seq(5, 25, by = 1),
  L_ref      = 13.399,
  population = "Fritillaria Miz"
)
print(sens$data)
sens$plot




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Ne_sensitivity_L", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Ne_sensitivity_Vc")
### * Ne_sensitivity_Vc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Ne_sensitivity_Vc
### Title: Sensitivity of Ne/N to clonal reproductive variance in one stage
### Aliases: Ne_sensitivity_Vc

### ** Examples

T_herb <- matrix(c(
  0.30, 0.05, 0.00,
  0.40, 0.65, 0.10,
  0.00, 0.20, 0.80
), nrow = 3, byrow = TRUE)

sens <- Ne_sensitivity_Vc(
  T_mat       = T_herb,
  F_vec       = c(0.0, 0.5, 3.0),
  D           = c(0.60, 0.25, 0.15),
  d           = c(0.0, 0.0, 0.7),
  stage_index = 3,
  Vc_range    = seq(0.5, 6, by = 0.5),
  population  = "mixed herb"
)
print(sens$data)
sens$plot




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Ne_sensitivity_Vc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Ne_sensitivity_Vk")
### * Ne_sensitivity_Vk

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Ne_sensitivity_Vk
### Title: Sensitivity of Ne/N to sexual reproductive variance in one stage
### Aliases: Ne_sensitivity_Vk

### ** Examples

T_herb <- matrix(c(
  0.30, 0.05, 0.00,
  0.40, 0.65, 0.10,
  0.00, 0.20, 0.80
), nrow = 3, byrow = TRUE)

sens <- Ne_sensitivity_Vk(
  model_fn    = Ne_sexual_Y2000,
  T_mat       = T_herb,
  F_vec       = c(0.0, 0.5, 3.0),
  D           = c(0.60, 0.25, 0.15),
  stage_index = 3,
  Vk_range    = seq(0.5, 6, by = 0.5),
  population  = "hypothetical herb"
)
print(sens$data)
sens$plot




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Ne_sensitivity_Vk", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Ne_sensitivity_d")
### * Ne_sensitivity_d

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Ne_sensitivity_d
### Title: Sensitivity of Ne/N to clonal fraction d_i in one stage
### Aliases: Ne_sensitivity_d

### ** Examples

T_herb <- matrix(c(
  0.30, 0.05, 0.00,
  0.40, 0.65, 0.10,
  0.00, 0.20, 0.80
), nrow = 3, byrow = TRUE)

# Under Poisson defaults: flat curve (d_i has no effect)
sens_poisson <- Ne_sensitivity_d(
  T_mat       = T_herb,
  F_vec       = c(0.0, 0.5, 3.0),
  D           = c(0.60, 0.25, 0.15),
  d_fixed     = c(0.0, 0.0, 0.5),
  stage_index = 3,
  population  = "mixed herb (Poisson)"
)
sens_poisson$plot

# With non-Poisson clonal variance: d_i matters
sens_nonpoisson <- Ne_sensitivity_d(
  T_mat       = T_herb,
  F_vec       = c(0.0, 0.5, 3.0),
  D           = c(0.60, 0.25, 0.15),
  d_fixed     = c(0.0, 0.0, 0.5),
  stage_index = 3,
  Vc_over_c   = c(1, 1, 3),
  population  = "mixed herb (Vc/c_bar = 3 in stage 3)"
)
sens_nonpoisson$plot




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Ne_sensitivity_d", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Ne_sexual_Y2000")
### * Ne_sexual_Y2000

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Ne_sexual_Y2000
### Title: Variance effective population size for 100% sexual
###   stage-structured populations
### Aliases: Ne_sexual_Y2000

### ** Examples

# ---------------------------------------------------------
# Example 1: Simple 3-stage plant, Poisson defaults
# A hypothetical outcrossing perennial herb with three stages:
#   Stage 1: seedling/juvenile (high mortality, no reproduction)
#   Stage 2: vegetative adult  (moderate survival, low reproduction)
#   Stage 3: reproductive adult (high survival, high reproduction)
# ---------------------------------------------------------
T_plant <- matrix(c(
  0.30, 0.05, 0.00,
  0.40, 0.65, 0.10,
  0.00, 0.20, 0.80
), nrow = 3, byrow = TRUE)

F_plant <- c(0.0, 0.5, 3.0)   # seeds per individual per stage per year
D_plant <- c(0.60, 0.25, 0.15) # observed stage fractions

result <- Ne_sexual_Y2000(
  T_mat      = T_plant,
  F_vec      = F_plant,
  D          = D_plant,
  population = "hypothetical herb"
)
print(result)

# ---------------------------------------------------------
# Example 2: Same population, high reproductive variance in stage 3
# Vk/k_bar = 3 for reproductive adults -- pollinator-limited, so only
# a few adults contribute most of the seeds in any given year.
# This should reduce Ne relative to Example 1.
# ---------------------------------------------------------
result_highvar <- Ne_sexual_Y2000(
  T_mat       = T_plant,
  F_vec       = F_plant,
  D           = D_plant,
  Vk_over_k   = c(1, 1, 3),
  population  = "hypothetical herb (high repro variance stage 3)"
)
print(result_highvar)
# Ne/N should be lower than Example 1.

# ---------------------------------------------------------
# Example 3: Supply L directly from a published source
# ---------------------------------------------------------
result_Luser <- Ne_sexual_Y2000(
  T_mat      = T_plant,
  F_vec      = F_plant,
  D          = D_plant,
  L          = 8.5,
  population = "hypothetical herb (published L)"
)
print(result_Luser)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Ne_sexual_Y2000", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
