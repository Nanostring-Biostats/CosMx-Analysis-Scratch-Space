# QC and normalization of CosMx RNA data

We've tried a lot of options here, and we've settled on very simple procedures for most cases. 

## QC

QC in CosMx is motivated by known error modes. Here's a list of major things that can go wrong:
- A cell might be undersampled, leading to excessively low counts (Either only a tip of it is in the slide, or detection efficiency is poor within it.) Solution: remove the cell. 
- A cell might suffer extremely high background, either due to intrinsic tissue stickiness (e.g. associated with necrosis) or due to optical artifacts. Solution: remove the cell.
- Errors in cell segmentation might assign multiple cells to the same "cell". Solution: remove these multiplets. 
- A FOV's expression profile can be distorted by image registration errors or by imaging artifacts, e.g. fluorescence hiding spots of one color. These FOVs can be analyzable if you're careful, but given the uncertainty they pose it's usually best to remove them. 



QC logic would then proceed as follows:

1. Remove cells with too few counts. For our 1000plex assay, we use a pretty generous threshold of 20 counts. A higher threshold would be reasonable. 

\code{
# counts is the matrix of raw expression profiles, cells in rows, genes in columns
totalcounts <- Matrix::rowSums(counts)  
drop <- totalcounts < 20
}

2. Remove cells with high outlier areas. You can use Grubb's test to detect outliers, 
or you can draw a histogram of cell areas and choose a cutoff on your own. 





 

