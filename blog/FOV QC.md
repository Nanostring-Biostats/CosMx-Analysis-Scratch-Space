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

All known FOV artifacts act by modulating our ability to detect reporter probes.
Most commonly, we see a single reporter cycle in which all 4 colors of probes lose efficiency;
that is, 4 "bits" of our color barcode are impacted, and in turn, so are all the genes sharing those barcode bits. 

Thus we see phenomena like the below, where genes with impacted bits are muted, while other genes behave normally:

![image](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/assets/4357938/3e3dbed6-5469-4bab-885e-ad1534dba420)


We have observed the below root causes of FOV artifacts:

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
bits in the impacted color will be harder to detect. And at the same time, they will
suffer higher rates of FalseCode style background events - i.e., their barcode will 
more often be spuriously observed in the absense of hyb probes for the gene. 


## Approach to FOV QC

Because all the above artifacts impact *bits*, i.e. reporter cycle/color pairs, 
we will look for artifacts at the level of bits, not genes. 
Specifically, for each barcode bit, we'll look for FOVs where genes using the bit are
underexpressed compared to comparable regions elsewhere. 

#### Technical details: 

We place a 7x7 grid across each FOV. For each grid square, we find the 10 most similar squares
in other FOVs, with "similar" being based on the square's expression profile. (We also only accept one 
neighbor per other FOV.)
Then, for each barcode bit, we take the genes using the bit, and We contrast their 
expression in the square vs. in the average of the 10 most similar squares elsewhere. 
For a given FOV and barcode bit, this gives us 49 contrasts. 
When an FOV's grid squares consistently underexpress the relevant gene set, we flag the FOV.

Below we demonstrate this approach, looking at a tissue with particularly dramatic FOV effects. 

On the left, we plot expression of a single barcode bit impacted by FOV effects. 
FOV 19 has almost entirely lost expression of the genes from this barcode bit,
and FOV 16 looks as though it could be losing some expression. 

On the right, we show the results of our FOV QC approach: for a 7x7 grid within 
each FOV, we see estimated change in barcode bit expression compared to similar 
grid squares in other FOVs. FOV 19 still stands out as an obvious failure. 
In contrast, the low expression in FOV 16 is shown to be similar - sometimes higher, sometimes lower - 
than biologically similar regions elsewhere in the tissue. 

## Code

Functions for FOV QC can be found [here](../code/FOV%20QC). 

We advise this approach be applied separately to each slide or tissue in a study. 

This approach is new as of April 2024, and as-yet lightly tested. 
Use thoughtfully. 



