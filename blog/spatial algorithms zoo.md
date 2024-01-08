# Recommended algorithms for spatial biology

Spatial statistics is a well-developed field, with deep statistical methodology and highly efficient open-source tools. 
In CosMx data, where a single study can contain millions of cells, computational efficiency is vital. 
Here we recommend some toolkits we've found useful:


### Fast nearest-neighbors search:

#### Returning the indexes and distances to a cell's K-nearest neighbors:
```
neighbors <- FNN::get.knnx(data = xy, # 2-column matrix of xy locations
                           query = xy, 
                           k = 50)
# returns 2 outputs: a matrix of each cell's nearest neighbor indices (including itself),
#  and a matrix of distances to these neighbors.
```

#### Returning a sparse matrix of cells' K-nearest neighbors

See the function [nearestNeighborGraph in the Insitucor package](https://github.com/Nanostring-Biostats/InSituCor/blob/main/R/NeighborhoodCalculations.R)

```
# xy is a 2-column matrix of xy locations
neighbors <- InSituCor:::nearestNeighborGraph(x = xy[, 1], y = xy[, 2], n=50)
```

#### Returning a sparse matrix of cells' neighbors within a radius

See the function [radiusBasedGraph in the Insitucor package](https://github.com/Nanostring-Biostats/InSituCor/blob/main/R/NeighborhoodCalculations.R)

```
# xy is a 2-column matrix of xy locations
neighbors <- InSituCor:::radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = 0.1)
```


### Measuring a gene's degree of spatial correlation

Our goal here is to measure how much a gene's expression depends on spatial location. 
Genes with strong spatial dependence are presumably more interesting, deserving human attention. 
A much-less-than-comprehensive list of methods is below.

Methods:
- Moran's I statistic: This is a time-honored method in spatial statistics, published in 1950. Using the analytical rather than the permutation p-value speeds it up greatly, and we find their performance to be similar. 
- [SpatialDE](https://github.com/Teichlab/SpatialDE): the first attempt to measure spatial autocorrelation in spatial transcriptomics. Can be slow. 
- [Maxspin](https://github.com/dcjones/maxspin): A more recent method using machine learning and information theory to get performance improvements. Can be slow. 
-[SPARK-X](https://github.com/xzhoulab/SPARK): Runs at speed similar to Moran's I.


### Measuring spatial correlation between two genes

When two or more genes are spatially correlated it can be of high biological interest. 
These genes might regulate each other, or they could be jointly regulated by some latent variable. 

Methods for measuring spatial correlation between genes include:
- Lee's L: another spatial statistics classic. 
- [SpatialDE](https://github.com/Teichlab/SpatialDE)

However, we have found methods like the above to be unsatisfying, since genes with cell-type-specific expression
end up sharing strong spatial correlations. E.g. CD19 and MS4A1 are expressed mainly by B-cells, so if B-cells are 
spatially clustered, then these genes will be spatially correlated, but for biologically trivial reasons. 
To isolate more interesting spatial correlations, we developed [InSituCor](https://github.com/Nanostring-Biostats/insitucor). This is our recommended approach. 
It can analyze hundreds of thousands of cells and thousands of genes in minutes. 


### Counting occurrences within cell neighborhoods

Analysts will often want to score cells for how often something occurs in their neighborhoods.
For example, you might want to know how many T-cell neighbors each cell has, or 
how many transcripts of a gene surround it. 

The below code demonstrates how to use the spatstat::marktable function to do this. 

```
# "xy"" is a 2-column matrix of cell locations
# "clust"" is a vector of cell type assignments
# create a point process object:
pp <- spatstat.geom::ppp(xy[, 1], xy[, 2], xrange = range(xy[, 1]), yrange = range(xy[, 2]))
marks(pp) <- clust
marks(pp) <- as.factor(marks(pp))
# count neighbors of each db cluster:
mt05 <- spatstat::marktable(X = pp, R = 0.05, N = NULL, exclude=TRUE, collapse=FALSE)
rownames(mt05) <- names(which(use))
```


