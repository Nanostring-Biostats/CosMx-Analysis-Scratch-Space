
![image](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/assets/4357938/24ab79ab-e9c5-486e-9595-68f48797d757)

# CosMx Analysis Scratch Space

This repo is here to accelerate analysis of CosMx datasets. It contains code and
 writeups addressing challenges we've faced in analyzing CosMx data. 
After dozens of analyses, we have developed a sense of what works and what doesn't, 
 and we have developed (sometimes partial) solutions for many frequent problems.

*Code in this repo is provided without warranty or promise of support.*

### How to navigate

See the [table of contents](#table-of-contents) below.

The repo contains:
1. Blog posts of commentary & advice 
2. Code, sometimes usable, sometimes just for use as a template
3. CASSutils, a slowly-growing R package of useful functions (tbd)

### Contributing

- For questions or suggestions, write an Issue.
- Contributions are welcome from the community - just reach out to get write permission. 


## Table of contents

Note:
- [ ] Items with a checkbox like this are awaiting material. 
- Items with links will take you to a .md writeup, which may in turn link to code.
 
### Analysis strategies
- [What is spatial data for?](blog/what%20is%20high%20plex%20spatial%20data%20for.md)
- [ ] The impact of segmentation error on differential expression analyses
- [ ] Spatial statistics resources - fast nearest neighbors search etc
- [ ] A generally satisfying set of UMAP parameters for CosMx data
- [ ] Quick & comprehensive searches for interesting trends with "Everything vs. everything DE"
- [ ] Smoothing single cell gene expression for enhanced plotting


### Visualization
- [Functions for condensing FOVs and tissues to minimize whitespace](blog/condensing%20FOVs%20and%20tissues%20in%20XY%20space.md)
- [ ] Inferring cell polygons from transcript locations
- [ ] (For fun) Spatial transcriptomics plots in stained glass 

### Cell typing
- [Cell typing: what we've found to work](blog/cell%20typing%20basics.md)
- [On the use of marker genes](blog/on%20cell%20typing%20with%20marker%20genes.md)
- [ ] Hierarchical ("plinko") cell typing
- [ ] Cell typing with smoothed marker genes

### Tissue-specific solutions
- [ ] A workflow for kidney samples: cell typing and glomerulus definitions
- [ ] Scoring brain cells for distance to plaques



 
