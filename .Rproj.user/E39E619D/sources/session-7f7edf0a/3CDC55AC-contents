## scripts/00_bootstrap.R

## 1) Make sure no installed package version is attached/loaded
if ("package:NeStage" %in% search()) {
  detach("package:NeStage", unload = TRUE, character.only = TRUE)
}
if ("NeStage" %in% loadedNamespaces()) {
  try(unloadNamespace("NeStage"), silent = TRUE)
}

## 2) Minimal deps
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

## 3) Utility to source all scripts in this folder (except self)
this_file <- "00_bootstrap.R"
src <- function(file) source(file, chdir = TRUE, encoding = "UTF-8")
all_files <- list.files("scripts", pattern = "\\.R$", full.names = TRUE)
for (f in setdiff(all_files, file.path("scripts", this_file))) src(f)

## 4) Create outputs folder
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

## 5) Run the Yonezawa demo end-to-end (produces tables + ggplots)
run_yonezawa_demo(save_plots = TRUE)
