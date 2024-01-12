# Updating cell typing results using alternative data

### Background

The single cell world has produced a wealth of algorithms for cell typing, many of 
which can work in spatial data. 
However, spatial datasets contain additional information useful for cell typing, 
and algorithms that only use gene expression miss out on this data. 
For example, a cell's spatial context can provide important insight into its cell type,
as can its image, especially in platforms like CosMx where several channels of protein imagery are collected. 

Here we describe a simple algorithm for using space and/or cell images to update
cell typing results from any algorithm. This approach harnesses our Insitutype 
cell typing algorithm.

### Approach

Inputs:
- Cell type assignments from any algorithm
- CosMx expression data (raw counts + mean negprobes per cell)
- Altenative data types - a matrix of cells' neighborhood characteristics, immunofluorescence staining,
 or morphology features. Alternatively, we can take simple xy positions and automatically
 derive a matrix of neighborhood attributes. 

The algorithm proceeds as follows:
\enumerate{
 \item Step 1: Get the mean background-subtracted expression profile of each cell type.
 \item Step 2: Cluster the alternative data to assign cells to pre-clusters, or "cohorts" in Insitutype's phrasing.
 \item Step 3: Call Insitutype to reclassify the cells, assigning cells to whichever profile
  best fits their data, while using the pre-clusters as an informative prior. 
  ([the Insitutype paper](https://www.biorxiv.org/content/10.1101/2022.10.19.512902v1) describes this logic in detail.)
}





