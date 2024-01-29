# Basic CosMx analysis workflow

## Introduction

This directory contains a complete demo analysis of a CosMx dataset. 
It's intended to be used as a template for other analyses to follow. 

We'll be analyzing a 6000-plex dataset from a _______.
We begin with the files written by the AtoMx flat file export module. 
The full dataset we've used is too big to be saved on Github. 
To follow along, we recommend using your own data, also as output by the AtoMx flat file export module.

## File structure

We'll organize the files for this analysis as follows:

```insert snapshot of the folders```

- The "raw_data" folder holds the exports from AtoMx
- "processed_data" holds data objects generated during analysis, meant to be used by later analyses. 
- "results" holds results intended for human consumption. 
- analysis scripts are in the top-level directory
- "utils" holds R scripts containing functions used by analyses
- scripts are numbered by the order in which they should be run. 
  Most scripts require output created by earlier scripts. 
  
## Data structure

Our flat file exports contain the following:

```insert cartoon of initial data types```

These contain:
- Raw counts
- Cell metadata: other attributes of cells, e.g. size, immunofluorescence values, tissue and FOV IDs,...
- Spatial locations: xy locations given in mm. **Warning**: studies containing multiple slides may initially have overlapping xy locations.
- Transcript data: for all RNA transcripts detected, location, gene ID, and cell ID. Most analyses use cell-level data, not this transcript-level data, but it can make compelling plots.
- Tissue images. Not used by most analyses, but useful for Figures. 

Our analyses will append lots of new information to this starting point, ending here:

```insert cartoon of all data types```

New data types include:
- UMAP coordinates
- Data acting as new metadata columns, e.g. cell type assignment and spatial cluster
- Special results objects from analyses: cell typing, differential expression, InSituCor, ...



## Workflow:

Our workflow performs the below steps:

First, the fundamentals:

- [Parse and format data export by AtoMx]()
- [Custom tiling of tissues and FOVs in space]()
- [QC and normalization]()
- [Use transcript locations to assess/refine cell segmentation]()
- [Dimension reduction (PCA and UMAP)]()
- [Cell typing]()

Then, we go after biology:

- [Defining cells' spatial context]()
- [Hypothesis-driven analyses, i.e. differential expression: how do cells change behavior based on their spatial context?]()
- [Hypothesis-generating analyses: identifying spatially correlated genes with InSituCor]()


## General analysis advice

```list of topics here, some statements, some links to CASS articles```






