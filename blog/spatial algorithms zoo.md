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
