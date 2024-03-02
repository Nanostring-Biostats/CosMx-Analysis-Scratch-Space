# Strategies for big datasets

CosMx data can be truly huge: millions of cells and thousands of genes. 
This prevents many typical analysis strategies. 
Here we'll discuss ways to work with big datasets.

### Data types:

CosMx data comes several varieties:

#### Large sparse matrices:
These are matrices of cells * genes or cells * cells, populated mainly by 0 values. 
Sparse matrix formats allow us to only store information for non-zero values, greatly reducing memory demands.
When working with sparse matrices, try to use methods that can act on this data type. 

Examples of sparse matrices:
- raw counts (sparse matrix, integers)
- norm counts (sparse, but now decimals. can round to 3 or 4 decimal places to control size a bit)
- cells' neighbor relationships (e.g. 50 entries per cell for 50 nearest neighbors)

#### Large dense datasets:
Some data is inevitably dense. Ideally, only pull this data into memory when you need it.

Examples of dense data:
- Cell metadata. Storing as a data table is most efficient. Since this usually has dozens of variables
  that are unnecessary for most analyses, you can also keep in memory only the columns you need for a given analysis. 
- Principal components. Unavoidably large. To save memory, store only the top 20-50 PCs, throwing out the information-light remaining PCs.

#### Small enough to not be a problem:
- umap
- xy locations

Other data:
- transcript locations (huge)
- polygons (very big)
