### On cell typing with marker genes

Our basic recommendation is this: relying on marker genes do not produce successful cell typing.

Spatial transcriptomics data has two features that make marker genes challenging to use.
1. Background: single cell profiles include two kinds of false counts: these platforms 
sometimes see transcripts that aren't present (false detections), and errors in cell
segmentation lead transcripts from one cell to be assigned to its neighbor. 
Both these phenomena lead to marker genes being counted in cells where they aren't truly present. 
2. Variable signal strength: tissues and cells vary widely in how efficiently existing RNA 
molecules are read. Thus genes with low expression are easily missed in many cells.

Applying the above phenomena to FOXP3, the canonical marker for Treg cells, we can envision 
non-Treg cells with spurious FOXP3 coming from false detections or contamination from a neighboring Treg,
and we can imagine Treg cells where FOXP3 isn't detected. A cell typing regime that applied an expression
threshold to FOXP3 would be unacceptably error-prone.

Instead of using marker genes, we recommend cell typing using most or all of cells' 
expression profiles. The data for a single gene in a single cell is noisy, but the
evidence from a complete expression profile is much more stable. Marker genes are useful 
for determining what cell type a cluster has captured. E.g., if a cluster is enriched in 
FOXP3, you can safely label it Tregs.

As an advanced approach, we have had success cell typing using *smoothed* expression 
of marker genes. We replace each cell's observed profile with the average profile of
the 20 or cells that have the most similar expression profiles to it. This essentially
performs a variance-bias tradeoff: we bias a cell to look like its neighbors in expression space,
but we greatly cut down the noise in the expression level. Cell typing based on marker genes in
this smoothed data can be successful. 
