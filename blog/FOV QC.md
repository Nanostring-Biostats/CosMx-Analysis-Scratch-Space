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

We have observed the below artifacts:

#### Registration failure: 
the images from each reporter cycle must be "registered", i.e. 
aligned to the images from the other cycles, in both horizontal and vertical position. 
This process can go wrong in various ways, but all with the same impact: the barcode bits 
from that reporter cycle are assigned to the wrong positions, and they no longer be used
to identify the RNA transcript they came from. This phenomenon drives down expression
for all genes with a barcode bit in the impacted reporter cycle. 
The CosMx instrument performs 8 "cycles" (as opposed to "reporter cycles") of data acquisition;
registration failure can impact a reporter cycle across one or all of these cycles, causing
either a slight decrease or a total loss of signal for the impacted genes. 

#### Autofluorescence: 
if the tissue in an FOV is autofluorescent, it can make fluorescent 
signal from CosMx reporter probes harder to detect. When this happens, all genes with barcode
bits in the impacted color will have lower expression.


## Approach to FOV QC

Because all the above artifacts impact *bits*, i.e. reporter cycle/color pairs, 
we will look for artifacts at the level of bits, not genes. 
Specifically, for each barcode bit, we'll look for FOVs where genes using the bit are
underexpressed compared to comparable regions elsewhere. 

#### Technical details: 

we place a 7x7 grid across each FOV. For each grid square, we find the 10 most similar squares
in other FOVs, with "similar" being based on the square's expression profile. (We also only accept one 
neighbor per other FOV.)
Then, for each barcode bit, we take the genes using the bit, and We contrast their 
expression in the square vs. in the average of the 10 most similar squares elsewhere. 
For a given FOV and barcode bit, this gives us 49 constrasts. 
WHen and FOV's grid squares consistently underexpress the relevant gene set, we flag the FOV.

## Code

Functions for FOV QC can be found [here](code/FOVQC/FOV%20QC/FOV%20QC%20utils.md). 

This approach is new as of April 2024, and as-yet lightly tested. 
Use thoughtfully. 


