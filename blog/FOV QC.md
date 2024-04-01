# FOV QC

## Background

In most experiments, all FOVs will perform comparably, and data analyses need not consider FOV.
However, FOVs can suffer a variety of technical effects, sometimes causing obvious bias in the 
data (e.g. all the cells in an FOV will be clustered as the same cell type), and sometimes more subtle. 
We recommend that FOV QC be performed early in analyses. Should misbehaving FOVs be detected,
we almost always recommend they should be excluded. 

Here we'll describe known FOV-level artifacts, and we'll show use of R code for detecting 
impacted FOVs.

## FOV artifacts



## Approach to FOV QC

Because all the above artifacts impact *bits*, i.e. reporter cycle/color pairs, 
we will look for artifacts at the level of bits, not genes. 
Specifically, for each bit in our barcodes, we'll
