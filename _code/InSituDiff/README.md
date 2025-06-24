# InSituDiff

## Motivation

This is an R package for exploring spatial transcriptomics datasets containing control and disease samples. 
It evaluate all disease cellular neighborhoods for how they are perturbed from the 
most comparable control neighborhoods. Applications include:

- Single-gene analyses, e.g. plotting gene perturbations across space or comparing perturbations across samples
- Obtaining lists of highly-perturbed genes
- Finding clusters of co-perturbed genes
- Finding tissue regions with distinct perturbation patterns (i.e. spatial clustering of perturbations)

## Installation

```
remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space", 
                        subdir = "_code/InSituDiff", ref = "Main")
```

## Usage:

See the vignette.   
