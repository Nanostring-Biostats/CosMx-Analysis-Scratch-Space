# InSituTree

## Motivation

InSituTree is an R package that builds on the InSituType package to perform cell typing on 
single cell spatial trasncriptomics data. It performs supervised cell typing in a hierarchical
manner, only considering genes likely to be good at differentiating the cell types or subtypes
being assessed.

## Installation

```
remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space", 
                        subdir = "_code/InSituTree", ref = "Main")
```

## Usage:

See the vignette.