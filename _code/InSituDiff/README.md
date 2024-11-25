![image](https://github.com/Nanostring-Biostats/InSituDiff/assets/4357938/609a70d1-8bb2-40bc-a0e1-d04672e90ee6)
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
devtools::install_github("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/InSituDiff")
```

## Usage:

See the vignette.   