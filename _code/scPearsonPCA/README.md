# scPearsonPCA

An R Package for computing PCA from pearson residual normalization, without needing to compute dense pearson residuals for a full cells x genes counts matrix.

Uses summary statistics and algebra shortcuts to derive the gene x gene covariance matrix of pearson residuals 
(with default options for centering, scaling, and clipping extreme values) rather than the large and dense cells x genes matrix.

### Installation
```
remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space", 
                        subdir = "_code/scPearsonPCA", ref = "Main")

```
