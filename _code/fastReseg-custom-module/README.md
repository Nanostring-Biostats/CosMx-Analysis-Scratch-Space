# Custom Module: FastReseg RNA Custom Module

Evaluate cell segmentation error and trim off contaminating transcripts based on the spatial profile of RNA transcripts using FastReseg package. 

This module would first check the existence of FastReseg outcomes and run the evaluation if no existing data or different input hash. The FastReseg evaluation outcomes are saved as new soma called `RNA_trimmed` in current study. 

User can specify whether to overwrite the `RNA` soma with the new soma for rerunning AtoMx pipeline from the step of QC module. It’s recommended to export your current `RNA` soma before overwriting it with new data. 

This module also adds a column called `lrtest_nlog10P` to the `RNA` obs and can be used in DE analysis.

## Setup (in the AtoMx custom module creation interface)

Please refer to the `CosMxDAFastResegSetup.docx` document on how to set up this custom module in AtoMx platform. 


### Outputs:

- A histogram plot for the segmentation contamination score calculated by FastReseg module.

- One new soma `soma_RNA_trimmed` attached to current study if `overwrite_RNA_soma` is FALSE.

- Existing soma `soma_RNA` is replaced with post-resegmented data in current study if `overwrite_RNA_soma` is TRUE.

### Notes:

It's recommended to run the fundation RNA pipeline on CosMx run to get cell types first before proceed to FastReseg custom module. 

If the module ran into memory issue, try to decrease the `maxCoreNum` to process fewer FOVs in parallel and thus reduce the amount of memory needed.

### Resources

- [Manuscript preprint](https://www.biorxiv.org/content/10.1101/2024.12.05.627051v1.full)

- [The FastReseg R package](https://github.com/Nanostring-Biostats/FastReseg)

- [Package documents and tutorial](https://nanostring-biostats.github.io/FastReseg/articles/tutorial.html)

