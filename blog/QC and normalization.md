# QC and normalization of CosMx RNA data

We've tried a lot of options here, and we've settled on very simple procedures for most cases. 

## QC

QC in CosMx is motivated by known error modes. Here's a list of major things that can go wrong:
- A cell might be undersampled, leading to excessively low counts (Either only a tip of it is in the slide, or detection efficiency is poor within it.) Solution: remove the cell. 


QC logic would then proceed as follows:


 

