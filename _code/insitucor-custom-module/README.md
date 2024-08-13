# InSituCor-custom-module


## Setup (in the AtoMx custom module creation interface)

- Give it 128Gb of RAM and 16 cores. 
- Specify a max run time of 4 hours
- Only one argument to specify: "celltype", display name = "celltype", range = 1-50, Default value = "Insitutype" (unless you want to use Leiden clustering), argument type = "String", and check "required"
- Put it in the pipeline after either leiden clustering or Insitutype

![image](https://github.com/user-attachments/assets/a8311d34-e0a4-416a-ba30-ecd224546ed8)

## Interpreting results

See the [Insitucor preprint](https://www.biorxiv.org/content/10.1101/2023.09.19.558514v1) 
and [Github page](https://github.com/Nanostring-Biostats/InSituCor)
for everything there is to know about the method. 

Briefly, we recommend Insitucor as a tool for hypothesis generation. 
It discovers modules of spatially correlated genes, and these modules 
often represent interesting biology.  

Once you've got results, we recommend starting with the cell type attribution heatmaps (output as .pdf files.)
They give a quick view of which genes and which cell types drive most modules. 
When a module strikes your interest, go explore its spatial plot (output as a .png file). 
Finally, to more thoroughly understand what a module is doing, use the AtoMx viewer to plot 
 its genes in space. For a good starting view, show individual transcripts from module genes, and color cells by cell type. 

 ### Attribution heatmap:
![image](https://github.com/user-attachments/assets/0fa10221-7bc6-4eca-98bb-72423a7aa203)


### Advanced users

For greater control over InSituCor, you can modify the below lines of the plugin .R script.
Run ?insitucor to see details of how to use these arguments. 

```
  k = 50,
  radius = NULL,
  tissue = tissuevec,
  min_module_size = 3,
  max_module_size = 30,
  resolution = 0.02,
  corthresh = 0.1,
  min_module_cor = 0.1,
  gene_weighting_rule = "inverse_sqrt",
  roundcortozero = 0.01,
  max_cells = 5000,
  attribution_subset_size = 5000
```
